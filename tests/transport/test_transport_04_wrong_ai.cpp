#include <thread>
#include <future>
#include <boost/serialization/utility.hpp>

#include "LatticeEncryption.hpp"
#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"
#include "gtest/gtest.h"
#include "TestHelpers.hpp"

using namespace ligero;
using namespace std::chrono_literals;

constexpr auto degree = 1 << 16;
constexpr auto p = 9;
constexpr auto q = 21;
constexpr auto parties = 3;
constexpr auto sigma = 8;
constexpr auto lambda = 128;
constexpr auto tau_limit_bit = 175;
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


    DBG("Host Throughput Test");
    auto survivor = transport.hostThroughputTest({duration, wait_time, cleanup_time, data_size}, [&](size_t t) { return t > throughput_cutoff; });
    if (survivor.size() < 2) {
        LOG(FATAL) << "too few survivor, quit";
    }
    DBG("Host Registration");
    hostRegistration(transport, config);

    DBG("Generate Keys");
    auto keygen_success = hostGenerateKeyPair<uint64_t, degree, q>(transport);

    if (!keygen_success) {
        LOG(INFO) << "Received restart (GOOD), starting coordinator second time";
        keygen_success = hostGenerateKeyPair<uint64_t, degree, q>(transport);
        if (!keygen_success) {
            LOG(FATAL) << "Keygen failed";
        } else {
            LOG(INFO) << "Succesfully finished coordinate";
        }
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


    try {
        auto maybe_keys = ec.generateKeyPair(transport, e.chi());
        if (!maybe_keys) {
            LOG(INFO) << "Received restart (Good), starting second run";
            maybe_keys = ec.generateKeyPair(transport, e.chi());
            if (!maybe_keys) {
                LOG(FATAL) << "Second run failed, aborting...";
            } else {
                LOG(INFO) << "Second run finished succesfully for Participant(" << id << ")";

            }
        } else {
            LOG(FATAL) << "Should have detected wrong ai and restart";
        }
    } catch (...) {
        LOG(FATAL) << "Caught exception in client...";
    }

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


    // change ai
    using Q = nfl::poly_p<uint64_t, degree, q>;
    Q ai(nfl::uniform{});
    ai.ntt_pow_phi();

    transport.send(MessageType::PUBLIC_KEY_A_SHARES, ai);
}

TEST(EncryptedProtocol, CandidatesGeneration) {

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

