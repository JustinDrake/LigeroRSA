/* Ligero */
#include "gtest/gtest.h"

// Boost
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include  <boost/serialization/unordered_set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

// STL
#include <vector>
#include <thread>
#include <future>

// Ligero
#include "Common.hpp"

#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"

using namespace std::chrono_literals;
using namespace ligero;

// Parameters RSA Ceremony
// ===========================

constexpr char *ip = "127.0.0.1";
constexpr char *full_addr = "tcp://127.0.0.1:5555";

constexpr auto degree = 1 << 16;
constexpr auto p = 9;
constexpr auto q = 21;
constexpr auto sigma = 8;
constexpr auto lambda = 128;
constexpr auto tau_limit_bit = 175;
constexpr auto tau = 1000;
constexpr auto numCandidates = 2048;

// throughput parties
constexpr auto wait_timeout = 300 * 1000; /* wait 5 min when start */
constexpr auto timeout = 25 * 1000; /* a little longer than test */
constexpr auto data_size = 32 * 1024; /* 32 KB */
constexpr auto nb_max_send = 1024;
constexpr auto parties = 2;

// throughtput coordinator
constexpr auto pbs = 1000;
constexpr auto duration = 20s;
constexpr auto wait_time = 120s;
constexpr auto cleanup_time = 5s;
constexpr auto throughput_cutoff = 100 * 1024; /* 200 KB/s */

size_t l = 1024;                                // Bloc size
size_t t = 256;                                 // Columns to Open Up (By default half the columns)
constexpr auto ptNumberEvaluation = 3;

//==============================================================================
void coordinate() {

    ProtocolConfig<uint64_t> config(parties, p, q, degree, sigma, lambda, tau_limit_bit, pbs);
    auto coordinator = EncryptedCoordinator<uint64_t, degree, p, q>();
    ZeroMQCoordinatorTransport transport(full_addr, config.numParties());
    coordinator.host(transport, config, {duration, wait_time, cleanup_time, data_size}, throughput_cutoff);
}

// For testing purposes, this test will generate data with a live run of the RSA Ceremony

void participate(int index) {

        // RSA Ceremony
        SocketId socketId = boost::uuids::random_generator()();
        ZeroMQClientTransport transport(full_addr, socketId);

        // RSA Ceremony
        auto client = EncryptedClient<uint64_t, uint64_t, degree, p, q>(socketId);

        client.start(transport, wait_timeout, timeout, data_size, nb_max_send, ip);

        // Zero-knowledge
        PublicData pdata;
        SecretData sdata;

        client.gatherData(pdata, sdata);

        // Set all parameters for active and passive protocols
            // core
        pdata.degree = degree;
        pdata.p = p;
        pdata.q = q; 
        pdata.sigma = sigma;
        pdata.lambda = lambda;
        pdata.tau_limit_bit = tau_limit_bit; 
        pdata.tau = tau;
        
            // active security
        pdata.ptNumberEvaluation = ptNumberEvaluation;

        // Export
        if (index == 1) {
            std::stringstream ssi;
            ssi << "RSACeremonyData";
            std::ofstream dataFile(ssi.str());

            std::stringstream ss;
            boost::archive::binary_oarchive oa(ss);
            oa << std::pair<PublicData, SecretData>(pdata, sdata);

            // Push it through the socket
            std::string msg = ss.str();

            dataFile << msg;
            dataFile.close();
            }
}

bool testZKceremony() {
    std::vector<std::thread> threads;

    for (auto i = 0; i < parties; i++) {
        threads.emplace_back(participate, i);
    }

    threads.emplace_back(coordinate);

    for (auto i = 0; i < parties + 1; i++) {
        threads[i].join();
    }

    return true;
}

//==============================================================================
/** Testing specific RSA Ceremony components for any ZK argument of knowledge */

// ==============================================================================
TEST(TestRSACeremonyComponentsZK, dataGeneration)               {
    ASSERT_EQ(testZKceremony(), true);
}

int main(int argc, char** argv)
{

    el::Configurations c;
    c.setToDefault();
    std::string config = std::string("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s.%g} %msg\n LOG_FLUSH_THRESHOLD = 1\n ENABLED = true\n FILENAME = \"test") + std::string(".log\" \n TO_FILE = true \n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n LOG_FLUSH_THRESHOLD = 1\n *INFO:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false");
    c.parseFromText(config.c_str());
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
    el::Loggers::setDefaultConfigurations(c, true);
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
