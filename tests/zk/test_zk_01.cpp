/* Ligero */
#include "gtest/gtest.h"
#include <vector>

#include "Common.hpp"
#include "ZkTests.hpp"
#include "protocol/ExpressNPStatement.hpp"
#include "SecretSharingNTT.hpp"

// Parameters RSA Ceremony
// ===========================

// Core Protocol
constexpr auto l = 8;//1024;
constexpr auto k = 16;//2048;
constexpr auto n = 32;//4096;
constexpr auto t = 4;//256;

// Parameters Zero-Knowledge
// ============================

// Settings Zero-Knowledge
RegularInteger p_zk(4611686018326724609ULL);
typedef Fpp_Fixed<p_zk> FieldT;

// Modulus Index
size_t modulusIdx = 0;

template<typename FieldT>
bool test(size_t numberProofs, Statements st, size_t l, size_t k, size_t n, size_t t, void (*computeSmall) (uint64_t*,uint64_t const*), void (*computeLarge) (uint64_t*,uint64_t const*)) {

    // Reed-Solomon Interleaved Code
    zksnark::InterleavedRSCode irs;

    irs.k       = k;       // Degree of the Polynomial
    irs.n       = n;       // Number of Shares

    // Merkle Tree Parameters
    ligero::MTParameters<FieldT> mtParams;

    mtParams.digestLength       = 32;
    mtParams.zkByteSize         = 8;
    mtParams.hashInnerNodeF     = hash::blake2bTwoToOneHash;
    mtParams.hashLeafContentF   = hash::blake2bFieldElementHash<FieldT>;
    mtParams.hashZKRandomnessF  = hash::blake2bZKElementHash;
    mtParams.leavesNumber       = irs.n;

    // Random Generator
    const hash::HashIntegerSeedGeneratorF<FieldT>& fieldGenerator = hash::blake2bIntegerSeedGenerator<FieldT>;
    const hash::HashStringSeedGeneratorF& stringGenerator = hash::blake2bStringHash;

    // Empty Public Data
    PublicData pdata;

    bool success = true;

    // Running the prover
    // ==============================
        DBG("Running component " << componentsDescription[static_cast<uint64_t>(st)]);

        // Generate the proof
        std::thread g(testGatherer<FieldT>, numberProofs);

        pdata.modulusIdx = modulusIdx;
        p_zk = nfl::params<uint64_t>::P[pdata.modulusIdx]; // Zero-Knowledge
        FieldT::modulusIdx_ = pdata.modulusIdx;
        FieldT::preprocessing(log2(n) + 1); 

        testProver(st, irs, l, t, fieldGenerator, stringGenerator, mtParams, pdata, computeSmall, computeLarge);
        g.join();

    // Running the verifier
    // ==============================

    SecretData sdata;

    std::ifstream argumentStream("Proof");
    boost::archive::binary_iarchive ia(argumentStream);


    try {
            std::pair<PublicData, ligero::zksnark::FullFSTranscript<FieldT>> deserializer;
            ia >> deserializer;
            PublicData& pdata = deserializer.first;

            pdata.modulusIdx = modulusIdx;
            p_zk = nfl::params<uint64_t>::P[pdata.modulusIdx]; // Zero-Knowledge
            FieldT::modulusIdx_ = pdata.modulusIdx;
            FieldT::preprocessing(log2(n) + 1); 

            zksnark::Verifier<FieldT> verifier(st, irs, t, l, fieldGenerator, mtParams, pdata, sdata, computeSmall, computeLarge);
            verifier.verifyConsolidatedZKSnarkFile(deserializer.second);
    }
    catch(internalError& e) {
        DBG("Internal Error: " << e.what());
        success = false;
    }
    catch(failedTest& e) {
        DBG("Failed Test: " << e.what());
        success = false;
    }

    DBG("Verified Argument Of Knowledge.");

    return success;
}

//==============================================================================
/** Testing generic building blocs for any ZK argument of knowledge */

TEST(TestBuildingBlocsZK, testLinearConstraint)         {ASSERT_EQ(test<FieldT>(1, Statements::Generic_Addition,                    l, k, n, t, permut<k>::compute, permut<n>::compute),       true);}
TEST(TestBuildingBlocsZK, testQuadraticConstraint)      {ASSERT_EQ(test<FieldT>(1, Statements::Generic_Multiplication,              l, k, n, t, permut<k>::compute, permut<n>::compute),       true);}
TEST(TestBuildingBlocsZK, testIntraBlocConstraint)      {ASSERT_EQ(test<FieldT>(1, Statements::Generic_IntraBloc,                   l, k, n, t, permut<k>::compute, permut<n>::compute),       true);}
TEST(TestBuildingBlocsZK, testInterblocBoundedVariable) {ASSERT_EQ(test<FieldT>(1, Statements::Generic_InterblocBoundingVariables,  l, k, n, t, permut<k>::compute, permut<n>::compute),       true);}
TEST(TestBuildingBlocsZK, testIntrablocBoundedVariable) {ASSERT_EQ(test<FieldT>(1, Statements::Generic_IntrablocBoundingVariables,  l, k, n, t, permut<k>::compute, permut<n>::compute),       true);}
TEST(TestBuildingBlocsZK, testEvalPoly)                 {ASSERT_EQ(test<FieldT>(1, Statements::Generic_PolyEval,                    l, k, n, t, permut<k>::compute, permut<n>::compute),       true);}

int main(int argc, char** argv)
{

    // Configure easylogging
    configureEasyLogging();

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}