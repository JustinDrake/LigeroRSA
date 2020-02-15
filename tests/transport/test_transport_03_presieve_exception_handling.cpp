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
        LOG(FATAL) << "Keygen timed out!";
    }

    DBG("Hosting");

    // Assign P1
    transport.broadcast(MessageType::ASSIGNMENT_P1, std::vector<SocketId>{ZeroMQCoordinatorTransport::ids[0]});
    DBG("Sending ASSIGNMENT_P1");

    // Assign Pn
    transport.broadcast(MessageType::ASSIGNMENT_PN, std::vector<SocketId>(ZeroMQCoordinatorTransport::ids.begin() + 1, ZeroMQCoordinatorTransport::ids.end()));
    DBG("Sending ASSIGNMENT_PN");

    
    auto presieve_success = hostPreSieving(transport, e, config); 
    if (!presieve_success) {
        LOG(INFO) << "Presieving restarted (Good), starting second run";
        presieve_success = hostPreSieving(transport, e, config); 
        if (!presieve_success) {
            LOG(FATAL) << "Presieving failed";
        }
    }
    auto [eit_pruned, alphasPS_tick, c_can, c_gcd] = *presieve_success;

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
        LOG(INFO) << "Kill/Restart received when generate presieve shares (Good), starting second run";
        maybe_shares = ec.generatePreSievingShares(e, transport, publicKey, secretKey, special, config);
        if (!maybe_shares) {
            LOG(FATAL) << "Kill/restart received on second run, aborting...";
        }
    }
    auto [xcan, ycan, zcan, xgcd, ygcd, zgcd, both_shares] = *maybe_shares;

}

void participate_failing() {
    constexpr auto wait_timeout = 5 * 1000; /* wait 5 sec when start */
    constexpr auto test_timeout = 25 * 1000; /* a little longer than test */
    constexpr auto data_size = 32 * 1024; /* 32 KB */
    constexpr auto nb_max_send = 1024;


    auto id = boost::uuids::random_generator()();
    CorruptedClient transport(ip, id);

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

    using Q = nfl::poly_p<uint64_t, degree, q>;
    // first uniformly choose a
    Q enc_x_shares(nfl::uniform{});
    enc_x_shares.ntt_pow_phi();

    transport.corruptedDataSend(MessageType::ENCRYPTED_X_SHARES, enc_x_shares);

    /*
    auto maybe_shares = ec.generatePreSievingShares(e, transport, publicKey, secretKey, special, config);
    if (!maybe_shares) {
        LOG(FATAL) << "Kill/Restart received when generate presieve shares";
    }
    auto [xcan, ycan, zcan, xgcd, ygcd, zgcd, both_shares] = *maybe_shares;
    */

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

