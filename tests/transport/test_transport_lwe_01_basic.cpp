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

#include "gtest/gtest.h"

using C = nfl::poly_from_modulus<uint16_t, 16, 28>;
static const C referencePoly(nfl::uniform{});

using namespace ligero;

void participate (
    ligero::SocketId socketId
) {
    ligero::ZeroMQClientTransport transport ("tcp://127.0.0.1:5555", socketId);

    transport.send(MessageType::MUTHU_ACK, referencePoly);
}

void coordinate (
    int numParties,
    std::vector<ligero::SocketId> ids
) {
    ligero::ZeroMQCoordinatorTransport transport ("tcp://127.0.0.1:5555", numParties);
    
    C r = transport.awaitAggregateVectorInput<C>(
            MessageType::MUTHU_ACK, 
            ids,
            [](const C& a, const C& b) { std::cout << "a: " << a << std::endl << " b:" << b << std::endl; return a + b; },
            C());

    ASSERT_EQ(referencePoly, r);
}

TEST(BasicTransportTest, BroadcastsPoly) {
    constexpr int kNumParties = 1;

    // Spin up the test protocol
    std::vector<ligero::SocketId> socketIds;
    std::vector<std::thread> threads;

    // Spin up numParties parties
    for (int i = 0; i < kNumParties; ++i) {
        ligero::SocketId socketId = boost::uuids::random_generator()();

        threads.emplace_back(
            participate,
            socketId
        );

        socketIds.push_back(socketId);
    }

    // Spin up a coordinator
    threads.emplace_back(
        coordinate,
        kNumParties,
        socketIds
    );


    // Wait for them each to finish
    while (threads.size() > 0) {
        threads.back().join();
        threads.pop_back();
    }

}

//==============================================================================
int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
