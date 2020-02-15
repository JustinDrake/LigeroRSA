#include <thread>
#include <thread>
#include <future>

#include <iostream>
#include <tuple>
#include <array>

#include <gmpxx.h>
#include <nfl.hpp>

#include "Transport.hpp"
#include "Ligero.hpp"
#include "RegistrationClient.hpp"
#include "RegistrationCoordinator.hpp"

using C = nfl::poly_from_modulus<uint16_t, 16, 28>;
static const C referencePoly(nfl::uniform{});
const std::string addr("tcp://127.0.0.1:5555");

using namespace ligero;

//==============================================================================
int main(int argc, char** argv)
{
    if (argc != 2) {
        LOG(FATAL) << "Usage: benchmark_party_01_serialization <ipaddress>";
    }

    std::string ipAddr(argv[1]);
    ligero::SocketId socketId = boost::uuids::random_generator()();
    ligero::ZeroMQClientTransport transport (std::string("tcp://")+ipAddr+std::string(":5555"), socketId);
    
    auto config = ligero::client::registerAndAwaitConfiguration<uint64_t>(transport, "n/a");
    transport.send(MessageType::MUTHU_ACK, referencePoly);

}