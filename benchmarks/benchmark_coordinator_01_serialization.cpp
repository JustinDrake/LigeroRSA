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

#include "RegistrationCoordinator.hpp"

using C = nfl::poly_from_modulus<uint16_t, 16, 28>;
static const C referencePoly(nfl::uniform{});

using namespace ligero;

//==============================================================================
int main(int argc, char** argv)
{
    if (argc != 4) {
        LOG(FATAL) << "Usage: benchmark_coordinator_01_serialization <ipaddress> <numParties> <logfile>";
    }

    std::string ipAddr(argv[1]);
    int numParties = atoi(argv[2]);
    std::string logfile(argv[3]);

    // Configuring Easylogging
    el::Configurations conf("logging");
    conf.set(el::Level::Global,el::ConfigurationType::Filename, logfile);
    el::Loggers::reconfigureAllLoggers(conf);

    ligero::ZeroMQCoordinatorTransport transport (std::string("tcp://")+ipAddr+std::string(":5555"), numParties);
    ligero::ProtocolConfig<C> protocolConfig;

    std::vector<SocketId> ids = coordinator::hostRegistration(transport, protocolConfig);

    timers myTimers;
    myTimers.initialize_timer();

    myTimers.begin(1,"serialization",0);

    C r = transport.awaitAggregateVectorInput<C>(
            MessageType::MUTHU_ACK, 
            ids,
            [](const C& a, const C& b) { std::cout << "a: " << a << std::endl << " b:" << b << std::endl; return a + b; },
            C());

    // assert(r == C(numParties));
    myTimers.end(1,"serialization",0);

    return 1;
}