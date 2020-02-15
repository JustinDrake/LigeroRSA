#include <thread>
#include <future>
#include <boost/serialization/utility.hpp>

#include "LatticeEncryption.hpp"
#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"
#include "TestHelpers.hpp"
#include "gtest/gtest.h"

using namespace ligero;
using namespace std::chrono_literals;


constexpr auto degree = 1 << 16;
constexpr auto p = 9;
constexpr auto q = 21;
constexpr auto parties = 2;
constexpr auto sigma = 8;
constexpr auto lambda = 128;
constexpr auto tau_limit_bit = 175; // 128
constexpr auto tau = 1000;
constexpr auto pbs = 1000;
constexpr auto throughput_cutoff = 100 * 1024; /* 200 KB/s */

constexpr auto ip = "tcp://127.0.0.1:5555";

void coordinate() {
    constexpr auto duration = 2s;
    constexpr auto wait_time = 12s;
    constexpr auto cleanup_time = 5s;
    constexpr auto data_size = 1 * 1024; /* 32 KB */
    constexpr auto throughput_cutoff = 100 * 1024; /* 200 KB/s */


    ProtocolConfig<uint64_t> config(parties, p, q, degree, sigma, lambda, tau_limit_bit, pbs);
    auto coordinator = EncryptedCoordinator<uint64_t, degree, p, q>();
    ZeroMQCoordinatorTransport transport(ip, parties);
    auto e = lattice::LatticeEncryption<uint64_t, degree, p, q>(config);


    protocolBeginningHelper(transport, e, config);
    
    auto presieve_success = hostPreSieving(transport, e, config); 
    if (!presieve_success) {
        LOG(FATAL) << "Presieving timed out!";
    }
    auto [eit_pruned, alphasPS_tick, c_can, c_gcd] = *presieve_success;

    
    auto candidates_success = hostModulusCandidates (transport, eit_pruned, alphasPS_tick, c_can, config);
    if (!candidates_success) {
        LOG(FATAL) << "Modulus candidate timed out!";
    }
    auto candidates = *candidates_success;

    auto [Bs, Ms] = ligero::math::compute_m_b_vec(4096, 104729, 2);
    auto nb_threads = 90;
    auto postSieve = std::bind(ligero::math::test_factorizable_threaded, std::placeholders::_1, Bs, Ms, nb_threads);
    {
        DBG("Start sanity check");
        for (int i = 0; i < 100; ++i) {
            auto x = boost::math::prime(i);
            for (int j = 0; j < candidates.size(); ++j) {
                if (candidates[j] % x == 0) {
                    DBG("candidates["<<j<<"] = " << candidates[j] << " x = " << x << "candidates[j] mod x = " << (candidates[j] % x));
                    assert(false);
                }
            }
        }
        DBG("Passed sanity check");
    }

    boost::dynamic_bitset<> discardFlags; 
    {
        DBG("Running postSieve");
        discardFlags = postSieve(std::vector<mpz_class>(candidates.data(), candidates.data() + candidates.size()));
        transport.broadcast(MessageType::POST_SIEVE, ZeroMQCoordinatorTransport::ids, discardFlags);
        candidates = discardCandidates(candidates, discardFlags);
        DBG("Completed postSieve");
    }

    bool foundModuli = true; 

    using PairOfVec = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>;
    if (foundModuli) {
        auto maybe_pq = transport.awaitAggregateVectorInput<PairOfVec>(
                MessageType::P_CLEAR_DEBUG,
                ZeroMQCoordinatorTransport::ids,
                [](const PairOfVec& a, const PairOfVec& b) {
                assert(a.first.size() == b.first.size());
                assert(a.second.size() == b.second.size());

                std::vector<mpz_class> fsum(a.first.size());
                std::vector<mpz_class> ssum(b.second.size());

                for (size_t i = 0; i < a.first.size(); ++i) {
                fsum[i] = a.first[i] + b.first[i];
                ssum[i] = a.second[i] + b.second[i];
                }
                return std::pair{fsum, ssum};
                },
                std::pair{std::vector<mpz_class>(candidates.size(), mpz_class(0)),
                std::vector<mpz_class>(candidates.size(), mpz_class(0))}
                );

        if (!maybe_pq) {
            LOG(FATAL) << "PQ timed out!";
        }
        auto pq = *maybe_pq;

        // TODO: Check why miller rabin fails
        /*
        for(size_t i = 0; i < pq.first.size(); ++i)
        {
            assert(boost::multiprecision::miller_rabin_test(MPInt(pq.first[i].get_mpz_t()), 4));
            assert(boost::multiprecision::miller_rabin_test(MPInt(pq.second[i].get_mpz_t()), 4));
        }
        */

        LOG(INFO) << "Prime factorization of the candidates are: ";
        for (size_t i = 0; i < pq.first.size(); ++i) {
            LOG(INFO) << pq.first[i];
            LOG(INFO) << pq.second[i];
        }

        assert(pq.first.size() == candidates.size());

        for(size_t i = 0; i< pq.first.size(); ++i)
        {
            assert(pq.first[i] * pq.second[i] == candidates[i]);
        }
        LOG(INFO) << "All candidate moduli checked.";
        DBG("Done.");
    }

}


void participate() {
    constexpr auto wait_timeout = 5 * 1000; /* wait 5 sec when start */
    constexpr auto test_timeout = 25 * 1000; /* a little longer than test */
    constexpr auto data_size = 32 * 1024; /* 32 KB */
    constexpr auto nb_max_send = 1024;


    auto id = boost::uuids::random_generator()();
    ZeroMQClientTransport transport(ip, id);

    auto proceed = transport.joinThroughputTest(wait_timeout, test_timeout, data_size, nb_max_send);

    if (!proceed) {
        LOG(FATAL) << "Kicked out due to poor network connection";
    }

    EncryptedClient<uint64_t, uint64_t, degree, p, q> ec(id);
    auto config = ec.registerAndAwaitConfiguration(transport, ip);
    auto e = lattice::LatticeEncryption<uint64_t, degree, p, q>(config);


    auto maybe_keys = ec.generateKeyPair(transport, e.chi());
    if (!maybe_keys) {
        LOG(FATAL) << "Kill/Restart received when keygen";
    }
    auto [publicKey, secretKey] = *maybe_keys;


    // Await designation
    DBG("Awaiting designation.");

    auto maybe_assignment = transport.awaitReply();
    if (!maybe_assignment) {
        LOG(FATAL) << "Kill/Restart when wait assignment";
    }
    MessageType assignment = *maybe_assignment;
    bool special = assignment == MessageType::ASSIGNMENT_P1;
    ec.speciald() = special;
    DBG("Participant (" << ec.socketId() << ") got my designation.");

    LOG(INFO) << "Connected to coordinator.";

    auto maybe_shares = ec.generatePreSievingShares(e, transport, publicKey, secretKey, special, config);
    if (!maybe_shares) {
        LOG(FATAL) << "Kill/Restart received when generate presieve shares";
    }
    auto [xcan, ycan, zcan, xgcd, ygcd, zgcd, both_shares] = *maybe_shares;

    auto maybe_modulus = ec.performModulusGeneration(transport, xcan, ycan, zcan, both_shares, special, config);
    if (!maybe_modulus) {
        LOG(FATAL) << "Kill/Restart received when perform modulus generation";
    }
    auto [candidates, p_shares, q_shares] = *maybe_modulus;

    {
        for (int i = 0; i < p_shares.size(); ++i) {
            if (special) {
                if (p_shares[i] % 4 != 3 || q_shares[i] % 4 != 3) {
                    DBG("special");
                    DBG("i = " << i);
                    DBG("p_shares[i] " << p_shares[i]);
                    DBG("p_shares[i] " << q_shares[i]);
                    assert(false);
                }
            } else {
                if (p_shares[i] % 4 != 0 || q_shares[i] % 4 != 0) {
                    DBG("ordinary");
                    DBG("i = " << i);
                    DBG("p_shares[i] " << p_shares[i]);
                    DBG("p_shares[i] " << q_shares[i]);
                    assert(false);
                }
            }
        }
    }

    auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>(MessageType::POST_SIEVE);
    if (!maybe_discardFlags) {
        LOG(FATAL) << "Kill/Restart when wait discard flags";
    }
    boost::dynamic_bitset<> discardFlags = *maybe_discardFlags;
    candidates = discardCandidates(candidates, discardFlags);
    p_shares = discardCandidates(p_shares, discardFlags);
    q_shares = discardCandidates(q_shares, discardFlags);

    bool foundModuli = true;
    if (foundModuli) {
        LOG(INFO) << "Found " << candidates.size() << " valid moduli:";

        for (int i = 0; i < candidates.size(); ++i) {
            LOG(INFO) << candidates[i];
        }

        // The following code transmits the p shares and q shares in the clear to the coordinator.
        // This is ONLY for checking the result. We need to remove this part finally.
        DBG("Participant (" << ec.socketId() << ") Transmitting the p and q values in the clear for debugging purposes only.");
        transport.send(MessageType::P_CLEAR_DEBUG, std::pair{p_shares, q_shares});
        DBG("Participant (" << ec.socketId() << ") done.");
    }

}

TEST(EncryptedProtocol, CandidatesGeneration) {

    std::vector<std::thread> threads;

    for (auto i = 0; i < parties; i++) {
        threads.emplace_back(participate);
    }

    threads.emplace_back(coordinate);

    for (auto i = 0; i < parties + 1; i++) {
        threads[i].join();
    }
}


int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

