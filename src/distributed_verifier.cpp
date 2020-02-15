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
#include "protocol/ExpressNPStatement.hpp"

#include "SecretSharingNTT.hpp"
#include "FiniteFields.hpp"

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

using namespace ligero;

// Parameters Zero-Knowledge
// ============================

// Settings Zero-Knowledge
RegularInteger p_zk(10514644122014433281ULL);
typedef Fpp_Fixed<p_zk> FieldT;

zksnark::InterleavedRSCode irs;                           // Reed-Solomon Interleaved Code
ligero::MTParameters<FieldT> mtParams;                    // Merkle Tree Parameters
constexpr size_t l = 4096;                                // Bloc size
constexpr size_t k = 8192;                                // Degree
constexpr size_t n = 16384;                               // # of Shares

constexpr size_t t = 256;                                 // Columns to Open Up (By default half the columns)
constexpr auto ptNumberEvaluation = 3;

// Random Generator
const hash::HashIntegerSeedGeneratorF<FieldT>& fieldGenerator = hash::blake2bIntegerSeedGenerator<FieldT>;
const hash::HashStringSeedGeneratorF& stringGenerator = hash::blake2bStringHash;

// Threading
constexpr size_t maxNumberThreads = 40; // [m5.metal has 96 threads over 48 cores]


int main(int argc, char** argv)
{
    if (argc<3) {
        std::cout << "USAGE: <Prod/Test> <partyId>" << std::endl;
    }

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

    FieldT::modulusIdx_ = 0;

    // Verifying
    // ==============================================================================
    size_t aggregateSuccess = 0;

    Statements statement = Statements::RSACeremony;
    Statements localStatement = statement;

    // Primitive Roots of Unity for All Moduli
    nfl::poly_p<uint64_t, degree, q> roots({0});
    for (size_t idx = 0; idx < q; idx++) {roots(idx,1) = 1ull;}
    roots.ntt_pow_phi();

    ligero::zksnark::FullFSTranscript<FieldT> transcript;

    SecretData sdata;
    PublicData pdata;
    
    std::string fileRoot(argv[1]);
    std::string stagPublicData("_PublicData_");
    std::string stagProof("_Proof_");
        
    std::string ssocketId(argv[2]);

    // Import Public Data
    std::string publicDataFile = fileRoot + stagPublicData + ssocketId;

    // Wait for Public Data to be Available 
    while (!doesTouchFileExist(publicDataFile)) {
        usleep(10);
    }

    fetchFromFile<PublicData>(pdata, publicDataFile);

    for (size_t modulusIdx = 0; modulusIdx< 21; modulusIdx++) { // TODO hardcoded

        std::string argumentFile = fileRoot + stagProof + ssocketId + std::string("_") + std::to_string(modulusIdx+1);

        while (!doesTouchFileExist(argumentFile)) {
            usleep(10);
        }

        fetchFromFile<ligero::zksnark::FullFSTranscript<FieldT>>(transcript, argumentFile);

        bool success = true;

        std::cout << "Modulus idx:" << modulusIdx << std::endl;

        // Zero-Knowledge
        p_zk = nfl::params<uint64_t>::P[modulusIdx];
        FieldT::modulusIdx_ = modulusIdx;
        FieldT::preprocessing(log2(n) + 1);

        pdata.modulusIdx = modulusIdx;

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
            std::cout << "Verified Argument Of Knowledge for Party " << ssocketId << ", modulus " << modulusIdx << std::endl;
            aggregateSuccess++;
        } else std::cout << "!!Argument Of Knowledge Failed for Party " << ssocketId << ", modulus " << modulusIdx << std::endl;
    }

    std::cout << "Verified All Arguments Of Knowledge for Party " << ssocketId << " " << aggregateSuccess << "/" << 21 << std::endl;
    return 0;
}
