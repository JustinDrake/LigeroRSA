#pragma once

#include "LatticeEncryption.hpp"
#include "EncryptedClient.hpp"

using namespace ligero;
using namespace std::chrono_literals;

template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
void protocolBeginningHelper(
        ZeroMQCoordinatorTransport& transport, 
        lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& e,
        ProtocolConfig<T>& config) {
    constexpr auto duration = 2s;
    constexpr auto wait_time = 12s;
    constexpr auto cleanup_time = 5s;
    constexpr auto data_size = 1 * 1024; /* 32 KB */
    constexpr auto throughput_cutoff = 100 * 1024; /* 200 KB/s */


    DBG("Host Throughput Test");
    auto survivor = transport.hostThroughputTest({duration, wait_time, cleanup_time, data_size}, [&](size_t t) { return t > throughput_cutoff; });
    if (survivor.size() < 2) {
        LOG(FATAL) << "too few survivor, quit";
    }
    DBG("Host Registration");
    hostRegistration(transport, config);

    DBG("Generate Keys");
    auto keygen_success = hostGenerateKeyPair<T, Degree, NbPrimesQ>(transport);

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
}

