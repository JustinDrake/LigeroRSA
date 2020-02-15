/* Ligero */
#include "gtest/gtest.h"

#include <vector>

#include "Common.hpp"
#include "ZkArgument.hpp"
#include "FiniteFields.hpp"
#include "SecretSharing.hpp"
#include "Hash.hpp"
#include "Transport.hpp"
#include "GatheringCoordinator.hpp"

#include "protocol/Constraint.hpp"
#include "protocol/Builder.hpp"
#include "protocol/ExpressNPStatement.hpp"

using std::vector;
using namespace ligero;

constexpr char *ip = "127.0.0.1";
constexpr char *full_addr = "tcp://127.0.0.1:5555";

//==============================================================================
/** Generic procedure for testing ZK components */

//==============================================================================

std::vector<std::string> componentsDescription = {
                "Generic_Addition",
                "Generic_Multiplication",
                "Generic_IntraBloc",
                "Generic_InterblocBoundingVariables",
                "Generic_Convolution",
                "RSAKeyGen",
                "RSApreSieveTripleGen",
                "RSAaggregation",
                "RSApartialDecrypt",
                "RSAcandidateGeneration",
                "RSAbeaversTriples",
                "RSAjacobiTest",
                "GenerateData"
            };

template<typename FieldT>
bool testVerifier(Statements statement, zksnark::InterleavedRSCode& irs, size_t l, size_t t, const hash::HashIntegerSeedGeneratorF<FieldT>& generator, const hash::HashStringSeedGeneratorF& stringGenerator, ligero::MTParameters<FieldT>& mtParams, void (*computeSmall) (uint64_t*,uint64_t const*), void (*computeLarge) (uint64_t*,uint64_t const*)) {
    bool success = true;
    
    try {
            PublicData pdata;
            SecretData sdata;
        
            zksnark::Verifier<FieldT> verifier(statement, irs, t, l, generator, stringGenerator, mtParams, pdata, sdata, computeSmall, computeLarge);
            // verifier.verifyConsolidatedZKSnarkFile(statement, "testVerifier");
    } 
    catch(internalError& e) {
        DBG("Internal Error: " << e.what());
        success = false;
    }
    catch(failedTest& e) {
        DBG("Failed Test: " << e.what());
        success = false;
    }

    return success;
}

template<typename FieldT>
void testProver(Statements statement, zksnark::InterleavedRSCode& irs, size_t l, size_t t, const hash::HashIntegerSeedGeneratorF<FieldT>& generator, const hash::HashStringSeedGeneratorF& stringGenerator, ligero::MTParameters<FieldT>& mtParams, PublicData& pdata, void (*computeSmall) (uint64_t*,uint64_t const*), void (*computeLarge) (uint64_t*,uint64_t const*)) {
    try {
            SecretData sdata;

            const SocketId _socketId = boost::uuids::random_generator()();
            ligero::ZeroMQClientTransport transport(full_addr, _socketId);

            zksnark::Prover<FieldT> prover(statement, irs, t, l, generator, stringGenerator, mtParams, pdata, sdata, computeSmall, computeLarge, transport);

            prover.produceArgumentOfKnowledge(1, ip, "5555");
    } 
    catch(std::exception& e) {
        DBG("Aborting: " << e.what());
    }
}

template<typename FieldT>
void testGatherer(size_t numberParties, size_t numberProofs) {
    ZeroMQCoordinatorTransport transport(full_addr, 1);

    GatheringCoordinator::gatherAndSaveArgumentsOfKnowledge<FieldT>(transport, "Test", numberProofs, true);
}

template<typename FieldT>
bool testZKcomponent(size_t numberProofs, Statements st, size_t l, size_t k, size_t n, size_t t, void (*computeSmall) (uint64_t*,uint64_t const*), void (*computeLarge) (uint64_t*,uint64_t const*)) {
        // Setting Protocol Parameters
        // ===========================

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
            const hash::HashIntegerSeedGeneratorF<FieldT>& randomGenerator = hash::blake2bIntegerSeedGenerator<FieldT>;

            // Empty Public Data
            PublicData pdata;

        // Running the test
        // ==============================
            DBG("Running component " << componentsDescription[static_cast<uint64_t>(st)]);

            // Generate the proof
            std::thread g(testGatherer<FieldT>, numberProofs);
            testProver(st, irs, l, t, randomGenerator, mtParams, pdata, computeSmall, computeLarge);
            g.join();

            // Verify
            return(testVerifier(st, irs, l, t, randomGenerator, mtParams, computeSmall, computeLarge));
}