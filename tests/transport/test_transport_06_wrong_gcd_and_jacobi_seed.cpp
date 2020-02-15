#include <thread>
#include <future>
#include <boost/serialization/utility.hpp>
#include <boost/multiprecision/miller_rabin.hpp>

#include "LatticeEncryption.hpp"
#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"
#include "TestHelpers.hpp"
#include "gtest/gtest.h"

using namespace ligero;
using namespace std::chrono_literals;


constexpr auto degree = 1 << 16;
constexpr auto p = 9; //8; //6; 
constexpr auto q = 21; // 19; // 14;
constexpr auto parties = 3;
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
    ZeroMQCoordinatorTransport transport(ip, config.numParties());

    auto coordinator = EncryptedCoordinator<uint64_t, degree, p, q>();
    coordinator.host(transport, config, {duration, wait_time, cleanup_time, data_size}, throughput_cutoff);
}


void participate() {
    constexpr auto wait_timeout = 5 * 1000; /* wait 5 sec when start */
    constexpr auto test_timeout = 25 * 1000; /* a little longer than test */
    constexpr auto data_size = 32 * 1024; /* 32 KB */
    constexpr auto nb_max_send = 1024;

    SocketId socketId = boost::uuids::random_generator()();
    ZeroMQClientTransport transport(ip, socketId);
    auto client = EncryptedClient<uint64_t, uint64_t, degree, p, q>(socketId);
    client.start(transport, wait_timeout, test_timeout, data_size, nb_max_send, ip);
}

void participate_failing() {
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


    // corrupt jacobiAndGCDSeedShares
    ec.jacobiAndGCDSeedShares() = ec.jacobiAndGCDSeedShares() + mpz_class(1);

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



    bool foundModuli = false;
    int numberOfTrials = 0;
    while (!foundModuli && (numberOfTrials < kTrialsThreshold)) {
        numberOfTrials++;

        LOG(INFO) << "Generating shares for pre-sieving.";
        auto maybe_shares = ec.generatePreSievingShares(e, transport, publicKey, secretKey, special, config);
        if (!maybe_shares) {
            LOG(ERROR) << "Kill/Restart received when generate presieve shares";
            LOG(FATAL) << "Unexpected failure";
        }
        auto [xcan, ycan, zcan, xgcd, ygcd, zgcd, both_shares] = *maybe_shares;

        LOG(INFO) << "Generated shares for " << both_shares.size() / 2 << " candidates.";

        LOG(INFO) << "Using pre-sieved shares for candidate generation.";
        auto maybe_modulus = ec.performModulusGeneration(transport, xcan, ycan, zcan,
                both_shares, special, config);
        if (!maybe_modulus) {
            LOG(ERROR) << "Kill/Restart received when perform modulus generation";
            LOG(FATAL) << "Unexpected failure";
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

        auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>
            (MessageType::POST_SIEVE);
        if (!maybe_discardFlags) {
            LOG(ERROR) << "Kill/Restart when wait discard flags";
            LOG(FATAL) << "Unexpected failure";
        }
        boost::dynamic_bitset<> discardFlags = *maybe_discardFlags;
        candidates = discardCandidates(candidates, discardFlags);
        p_shares = discardCandidates(p_shares, discardFlags);
        q_shares = discardCandidates(q_shares, discardFlags);

        // NOTE: here we need to adjust
        //indicesArray = discardCandidates(indicesArray, discardFlags);

        LOG(INFO) << "Running distributed primality test on the surviving " <<
            candidates.size() << " candidates.";

        //Finally, we move onto the biprimality test, which we implement here
        //in a loop because the coordinator may decide to run each test a variable
        //number of times.
        bool stillRunningTests = true;
        bool ranJacobi = false;

        DBG("Failing party started");
        while (stillRunningTests) {
            auto success = transport.awaitReply();
            if (!success) {
                LOG(ERROR) << "Kill/Restart received";
                LOG(FATAL) << "Unexpected failure";
            }
            switch (*success) {
                case MessageType::GAMMA_SHARES: {
                                                    if (!ranJacobi) {
                                                        LOG(INFO) << "Running Jacobi test on candidates.";
                                                        ranJacobi = true;
                                                    }
                                                    auto maybe_JacobiResult = ec.performJacobiTest(transport, candidates, p_shares, q_shares, special);
                                                    if (!maybe_JacobiResult) {
                                                        LOG(FATAL) << "Unexpected failure";
                                                    }

                                                    auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>(MessageType::DISCARD_FLAGS);
                                                    if (!maybe_discardFlags) {
                                                        LOG(ERROR) << "Kill/Restart received during wait discard flags";
                                                        LOG(FATAL) << "Unexpected failure";
                                                    }
                                                    boost::dynamic_bitset<> discardFlags = *maybe_discardFlags;

                                                    candidates = discardCandidates(candidates, discardFlags);
                                                    p_shares = discardCandidates(p_shares, discardFlags);
                                                    q_shares = discardCandidates(q_shares, discardFlags);

                                                    break;
                                                }
                case MessageType::GCD_RAND_SHARES: {
                                                       LOG(INFO) << "Running GCD test on candidates.";

                                                       if (candidates.size() > 1) {
                                                           candidates.resize(1);
                                                           p_shares.resize(1);
                                                           q_shares.resize(1);
                                                       }


                                                       std::vector<mpz_class> ax_shares(primesGCD);
                                                       std::vector<mpz_class> by_shares(primesGCD);
                                                       for (size_t axi = 0; axi < ax_shares.size(); ++axi) {
                                                           ax_shares[axi] = mpz_class(1);
                                                           by_shares[axi] = mpz_class(1);
                                                       }
                                                       transport.send(MessageType::GCD_AX_BY_SHARES, std::pair{std::pair{ax_shares, by_shares}, ec.jacobiAndGCDSeedShares()});
                                                       LOG(INFO) << "Corrupted party finished nicely";
                                                       return;
                                                       /*
                                                       auto maybe_JacobiGCDResult = ec.performGCDandJacobiTest(
                                                               transport,
                                                               candidates,
                                                               p_shares,
                                                               q_shares,
                                                               xgcd,
                                                               ygcd,
                                                               zgcd,
                                                               special,
                                                               config
                                                               );

                                                       if (!maybe_JacobiGCDResult) {
                                                           LOG(FATAL) << "Unexpected failure";
                                                       }

                                                       auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>(MessageType::DISCARD_FLAGS);
                                                       if (!maybe_discardFlags) {
                                                           LOG(ERROR) << "Kill/Restart received during wait discard flags";
                                                           LOG(FATAL) << "Unexpected failure";
                                                       }
                                                       boost::dynamic_bitset<> discardFlags = *maybe_discardFlags;

                                                       candidates = discardCandidates(candidates, discardFlags);
                                                       p_shares = discardCandidates(p_shares, discardFlags);
                                                       q_shares = discardCandidates(q_shares, discardFlags);

                                                       break;
                                                       */
                                                   }
                case MessageType::FOUND_MODULI:
                                                   foundModuli = true;
                                                   stillRunningTests = false;
                                                   break;
                case MessageType::NO_MODULI: {
                                                 LOG(INFO) << "Did not find a modulus. Restarting";
                                                 stillRunningTests = false;
                                                 break;
                                             }
                default:
                                             throw std::runtime_error("Received a message type out of order.");
            }
        }

        if (foundModuli) {
            LOG(INFO) << "Found " << candidates.size() << " valid moduli:";

            for (int i = 0; i < candidates.size(); ++i) {
                LOG(INFO) << candidates[i];
            }

            /*
            size_t idx = 0;
            bool rec = false;
            for (; idx < candidatesCAN().size(); idx++) {
                if (candidatesCAN()[idx] == candidates[0]) {
                    rec = true;
                    break;
                }
            }

            if (!rec) {
                // handle error
            }

            candidateIndices() = getCandidateIndices(idx);
            final_index() = idx;
            //std::cout<< speciald() << ": p_shares[0] = " << p_shares[0] << std::endl;
            */

            // The following code transmits the p shares and q shares in the clear to the coordinator.
            // This is ONLY for checking the result. We need to remove this part finally.
            DBG("Participant (" << ec.socketId() <<
                    ") Transmitting the p and q values in the clear for debugging purposes only.");
            transport.send(MessageType::P_CLEAR_DEBUG, std::pair{p_shares, q_shares});
            DBG("Participant (" << ec.socketId() << ") done.");

        }
    }
    
}

TEST(EncryptedProtocol, JacobiTest) {

    std::vector<std::thread> threads;

    for (auto i = 0; i < parties - 1; i++) {
        threads.emplace_back(participate);
    }

    threads.emplace_back(participate_failing);

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
