/* Ligero */
#include "gtest/gtest.h"

// Boost
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
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

#include "ZkArgument.hpp"
#include "ZkTests.hpp"
#include "protocol/ExpressNPStatement.hpp"

#include "ZkArgument.hpp"
#include "SecretSharingNTT.hpp"

// Settings Zero-Knowledge
RegularInteger p_zk(4611686018326724609ULL);
typedef Fpp_Fixed<p_zk> FieldT;

void g(int a) {

    std::cout << "exec" << std::endl;
}

int main(int argc, char** argv)
{
    FieldT::modulusIdx_ = 0;

    std::vector<std::thread> l;

    l.emplace_back(g, 1);
    l.emplace_back(g, 1);

    l[0].join();
    l[1].join();

    std::cout << l.size();

    l.clear();

    l.emplace_back(g, 1);
    l.emplace_back(g, 1);

    l[0].join();
    l[1].join();

    return 1;
    std::vector<ligero::zksnark::FullFSTranscript<FieldT>> data;

    std::cout << argv[1] << std::endl;
    fetchFromFile(data, argv[1]);

    std::cout << "done" << std::endl;

    return 0;
}
