#include <thread>
#include <future>
#include <thread>
#include <chrono>
#include <boost/serialization/utility.hpp>

#include "Ligero.hpp"
#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"
#include "gtest/gtest.h"

using namespace ligero;
using namespace std::chrono_literals;

constexpr auto degree = 1 << 16;
constexpr auto p = 9; //8; //6; 
constexpr auto q = 21; // 19; // 14;
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
    ZeroMQCoordinatorTransport transport(ip, parties);
    auto coordinator = EncryptedCoordinator<uint64_t, degree, p, q>();
    coordinator.host(transport, config, {duration, wait_time, cleanup_time, data_size}, throughput_cutoff);
}

void participate(int index) {
    constexpr auto wait_timeout = 5 * 1000; /* wait 5 sec when start */
    constexpr auto timeout = 25 * 1000; /* a little longer than test */
    constexpr auto data_size = 32 * 1024; /* 32 KB */
    constexpr auto nb_max_send = 1024;

    auto id = boost::uuids::random_generator()();
    auto client = EncryptedClient<uint64_t, uint64_t, degree, p, q>(id);
    ZeroMQClientTransport transport(ip, id);
    client.start(transport, wait_timeout, timeout, data_size, nb_max_send, ip);
}

TEST(DistributedEncryption, 2PartiesProduct) {

    std::vector<std::thread> threads;

    for (auto i = 0; i < parties; i++) {
        threads.emplace_back(participate, i);
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
