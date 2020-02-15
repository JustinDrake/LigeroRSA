#pragma once

#include "LatticeEncryption.hpp"
#include "EncryptedClient.hpp"

using namespace ligero;
using namespace std::chrono_literals;

class CorruptedClient: public ZeroMQClientTransport {
public:
    using ZeroMQClientTransport::ZeroMQClientTransport;

    template <typename T>
    void corruptedDataSend (MessageType t, const T& message) {
        std::stringstream ss;

        // Serilize column vector to string
        boost::archive::binary_oarchive oa(ss);
        oa << t;
        oa << message;

        // Push it through the socket
        std::string msg = ss.str();

        // corrupt
        for (size_t i = 0; i < msg.length(); ++i) {
            msg[i] = 1;
        }

        bool success = socket.send(msg.c_str(), msg.size());
        if (!success) throw std::runtime_error("Failing to send message.");
    }


    template <typename T>
    void corruptedSizeSend (MessageType t, const T& message) {
        std::stringstream ss;

        // Serilize column vector to string
        boost::archive::binary_oarchive oa(ss);
        oa << t;
        oa << message;

        // Push it through the socket
        std::string msg = ss.str();

        // corrupt by taking only half of a message
        size_t halfSize = size_t(msg.size() / 2);
        std::string corruptedMsg = msg.substr(0, halfSize);
        DBG("origin msg size = " << msg.size());
        DBG("halfSize = " << halfSize);
        

        bool success = socket.send(corruptedMsg.c_str(), corruptedMsg.size());
        if (!success) throw std::runtime_error("Failing to send message.");
    }

};

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

