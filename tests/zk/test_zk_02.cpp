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

#include "ZkArgument.hpp"
#include "ZkTests.hpp"
#include "protocol/ExpressNPStatement.hpp"

#include "ZkArgument.hpp"
#include "SecretSharingNTT.hpp"

using namespace std::chrono_literals;
using namespace ligero;

// Parameters RSA Ceremony
// ===========================

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
constexpr auto parties = 1;

// throughtput coordinator
constexpr auto pbs = 1000;
constexpr auto duration = 20s;
constexpr auto wait_time = 120s;
constexpr auto cleanup_time = 5s;
constexpr auto throughput_cutoff = 100 * 1024; /* 200 KB/s */

// debug
constexpr auto actQ = 21;

// com
SocketId _socketId;

// Parameters Zero-Knowledge
// ============================

// Settings Zero-Knowledge
RegularInteger p_zk(10514644122014433281ULL);
typedef Fpp_Fixed<p_zk> FieldT;

zksnark::InterleavedRSCode irs;                 // Reed-Solomon Interleaved Code
ligero::MTParameters<FieldT> mtParams;          // Merkle Tree Parameters
constexpr size_t l = 4096;                                // Bloc size
constexpr size_t k = 8192;                                // Degree
constexpr size_t n = 16384;                               // # of Shares

constexpr size_t t = 256;                                 // Columns to Open Up (By default half the columns)
constexpr auto ptNumberEvaluation = 3;

// Random Generator
const hash::HashIntegerSeedGeneratorF<FieldT>& fieldGenerator = hash::blake2bIntegerSeedGenerator<FieldT>;
const hash::HashStringSeedGeneratorF& stringGenerator = hash::blake2bStringHash;

//==============================================================================
void coordinate(Statements statement) {
    testGatherer<FieldT>(parties, actQ); // we only run one party in these tests for simplification
}

// For testing purposes, we will generate data with a live run first 
// And then test various sections of the statement against it

template<typename FieldT>
void participate(int index, Statements statement) {

        LOG(INFO) << "BROAD,KICKOFF"; 

        // Ceremony Data
        std::ifstream dataStream("RSACeremonyData");
        boost::archive::binary_iarchive ia(dataStream);
        std::pair<PublicData, SecretData> deserializer;
        ia >> deserializer;
        dataStream.close();

        PublicData& pdata = deserializer.first;
        SecretData& sdata = deserializer.second;

        _socketId = boost::uuids::random_generator()();
        ligero::ZeroMQClientTransport transport(full_addr, _socketId);

        // Primitive Roots of Unity for All Moduli
        nfl::poly_p<uint64_t, degree, q> roots({0});
        for (size_t idx = 0; idx < q; idx++) {roots(idx,1) = 1ull;}
        roots.ntt_pow_phi();

        for (pdata.modulusIdx = 0; pdata.modulusIdx < actQ; pdata.modulusIdx++) {
            Statements localStatement = statement;
            LOG(INFO) << "Modulus idx:"; 

            // Roots of Unity
            pdata.roots.assign(roots.poly_obj().data()+pdata.modulusIdx*(size_t)roots.degree,roots.poly_obj().data()+(pdata.modulusIdx+1)*(size_t)roots.degree);
            pdata.ptNumberEvaluation = ptNumberEvaluation;

            // Zero-Knowledge
            p_zk = nfl::params<uint64_t>::P[pdata.modulusIdx];
            FieldT::modulusIdx_ = pdata.modulusIdx;
            FieldT::preprocessing(log2(n) + 1); 

            // For portions 7 to 12, only run the first 9 moduli
            if ((statement == Statements::RSACeremony)&&(pdata.modulusIdx>8)) localStatement = Statements::RSACeremony_Rounds3to6;

            zksnark::Prover<FieldT> prover(localStatement, irs, t, l, fieldGenerator, stringGenerator, mtParams, pdata, sdata, permut<k>::compute, permut<n>::compute, transport);
            prover.produceArgumentOfKnowledge(1);
        }

        LOG(INFO) << "BROAD,DONE_PROVING"; 
}

template<typename FieldT>
bool testZKceremony(Statements statement) {

    // Producing the proof
    // ==============================================================================
    std::vector<std::thread> threads;

    for (auto i = 0; i < parties; i++) {
        threads.emplace_back(participate<FieldT>, i, statement);
    }

    threads.emplace_back(coordinate, statement);

    for (auto i = 0; i < parties + 1; i++) {
        threads[i].join();
    }

    // Verifying
    // ==============================================================================
    size_t aggregateSuccess = 0;
    Statements localStatement = statement;

    // Primitive Roots of Unity for All Moduli
    nfl::poly_p<uint64_t, degree, q> roots({0});
    for (size_t idx = 0; idx < q; idx++) {roots(idx,1) = 1ull;}
    roots.ntt_pow_phi();

    SecretData sdata;
    PublicData pdata;
    
    std::string fileRoot("Test");
    std::string stagPublicData("_PublicData_");
    std::string stagProof("_Proof_");
    
    std::stringstream ss;
    ss << _socketId;
    
    std::string ssocketId(ss.str());

    // Import Public Data
    std::string publicDataFile = fileRoot + stagPublicData + std::string("test");
    std::string argumentFile = fileRoot + stagProof + std::string("test");

    fetchFromFile<PublicData>(pdata, publicDataFile);

    for (size_t modulusIdx = 0; modulusIdx < actQ; modulusIdx++) {
        bool success = true;
        ligero::zksnark::FullFSTranscript<FieldT> transcript;
        
        fetchFromFile<ligero::zksnark::FullFSTranscript<FieldT>>(transcript, argumentFile + std::string("_") + std::to_string(modulusIdx+1));

        LOG(INFO) << "Modulus idx:"; 

        // Zero-Knowledge
        p_zk = nfl::params<uint64_t>::P[modulusIdx];
        pdata.modulusIdx = modulusIdx;
        FieldT::modulusIdx_ = modulusIdx;
        FieldT::preprocessing(log2(n) + 1); 

        // For portions 7 to 12, only run the first 9 moduli
        if ((localStatement == Statements::RSACeremony)&&(pdata.modulusIdx>8)) localStatement = Statements::RSACeremony_Rounds3to6;

        // Roots of Unity
        pdata.roots.assign(roots.poly_obj().data()+pdata.modulusIdx*(size_t)roots.degree,roots.poly_obj().data()+(pdata.modulusIdx+1)*(size_t)roots.degree);
        pdata.ptNumberEvaluation = ptNumberEvaluation;

        zksnark::Verifier<FieldT> verifier(localStatement, irs, t, l, fieldGenerator, stringGenerator, mtParams, pdata, sdata, permut<k>::compute, permut<n>::compute);

        try {
            verifier.verifyConsolidatedZKSnarkFile(transcript);
        }
        catch(failedTest& e) {
            DBG("Failed Test: " << e.what());
            success = false;
        }

        if (success) {
            LOG(INFO) << "Verified Argument Of Knowledge for Party " << "test, modulus " << modulusIdx;
            aggregateSuccess++;
        } else LOG(INFO) << "!!Argument Of Knowledge Failed for Party " << "test, modulus " << modulusIdx;
    }

    LOG(INFO) << "Verified All Arguments Of Knowledge for Party " << "test " << aggregateSuccess << "/" << actQ;
    LOG(INFO) << "BROAD,DONE_VERIFYING"; 

    return (aggregateSuccess == actQ);
}

//==============================================================================
/** Testing specific RSA Ceremony components for any ZK argument of knowledge */
 
// ==============================================================================
TEST(TestRSACeremonyComponentsZK, fullRSACeremony)              {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSACeremony),                     true);}
// TEST(TestRSACeremonyComponentsZK, RSACeremonyP)                 {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSACeremony_Rounds3to6),          true);}
// TEST(TestRSACeremonyComponentsZK, round3_testKeygen)            {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound3_keyGen),                true);}
// TEST(TestRSACeremonyComponentsZK, round4_preSieveTripleGen)     {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound4_preSieving),            true);}
// TEST(TestRSACeremonyComponentsZK, round5_aggregation)           {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound5_preSieving),            true);}
// TEST(TestRSACeremonyComponentsZK, round6_partialDecrypt)        {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound6_partialDecrypt),        true);}
// TEST(TestRSACeremonyComponentsZK, round7_candidateGeneration)   {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound7_candidateGeneration),   true);}
// TEST(TestRSACeremonyComponentsZK, round8_beaversTriples)        {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound8_beaversTriples),        true);}
// TEST(TestRSACeremonyComponentsZK, round10_jacobiTest)           {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound10_jacobiTest),           true);}
// TEST(TestRSACeremonyComponentsZK, round11_jacobiGCD)            {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSARound11_12_jacobiGCD),         true);}
// TEST(TestRSACeremonyComponentsZK, connecting_Proofs)           {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSACeremony_Connecting_Proofs),           true);}
// TEST(TestRSACeremonyComponentsZK, bounding_Variables)           {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSACeremony_Bounding_Variables),           true);}
// TEST(TestRSACeremonyComponentsZK, equating_Variables)           {ASSERT_EQ(testZKceremony<FieldT>(Statements::RSACeremony_Equating_Variables),           true);}

int main(int argc, char** argv)
{
    // Set Zero-Knowledge Parameters 
    irs.k                       = k;       // Degree of the Polynomial
    irs.n                       = n;       // Number of Shares

    mtParams.digestLength       = 32;
    mtParams.zkByteSize         = 8;
    mtParams.hashInnerNodeF     = hash::blake2bTwoToOneHash;
    mtParams.hashLeafContentF   = hash::blake2bFieldElementHash<FieldT>;
    mtParams.hashZKRandomnessF  = hash::blake2bZKElementHash;
    mtParams.leavesNumber       = irs.n;

    // Configure easylogging
    configureEasyLogging();

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
