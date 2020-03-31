// Boost
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// STL
#include <chrono>
#include <condition_variable>
#include <vector>
#include <thread>
#include <future>
#include <mutex>

// Ligero
#include "Common.hpp"

#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"

#include "ZkArgument.hpp"
#include "protocol/ExpressNPStatement.hpp"

#include "SecretSharingNTT.hpp"
#include "FiniteFields.hpp"
#include "Hash.hpp"

using namespace ligero;
using namespace std::chrono_literals;

// ===== Thread synchronization ===== //
std::mutex mutex;
std::condition_variable cv;
bool thread_aborted = false;

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

constexpr auto receive_timeout = 20min;

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
constexpr size_t maxNumberThreads = 50; // [m5.metal has 96 threads over 48 cores]

// Selecting a port for distributed verifying
constexpr char* port = "5556";

// Depositing the data here
PublicData pdata;
SigmaProtocolPublicData spdata;
std::vector<ligero::zksnark::FullFSTranscript<FieldT>> transcripts;
size_t availableTranscriptsIndex = 0;
bool firstIteration = true;

std::string localIpAddr = "UNKNOWN";

// Stitching Proofs
constexpr int nbIterations     = 2;
constexpr int nbVarsStitched   = 4;
constexpr int nbBlocksPerVar   = degree/l;


// Data
std::string hashedPublicData;

void resetGlobalVariables() {
    pdata = PublicData{};
    transcripts.clear();
    availableTranscriptsIndex = 0;
    firstIteration = true;
    thread_aborted = false;
    hashedPublicData.clear();
}

// Concerned that otherwise in the event of verifier failure and retry
// We might fall behind on receiving public data & proof because of processing delays
void collectData(ZeroMQClientTransport& transport) {
    DBG("Verifier " << transport.getSocketId() << " begining data acquisition.");

    // Tell the gathering coodrinator we are ready
    if (firstIteration) transport.send(MessageType::VERIFIER_READY, localIpAddr);

    // Receive Sigma Protocol Public Data 
    auto maybe_sigma = transport.awaitReply<SigmaProtocolPublicData>(MessageType::GATHER_SIGMA_DATA, true);

    if (hasError(maybe_sigma)) {
        LOG(ERROR) << "Kill/Restart received";
        thread_aborted = true;
        cv.notify_one();
        return;
    } else {
        spdata = getResult(maybe_sigma);
        // ligero::dumpSigmaProtocolPublicData(spdata);
    }

    // Receive public data
    auto maybe_pdata = transport.awaitReply<PublicData>(MessageType::GATHER_PUBLIC_DATA, true);

    if (hasError(maybe_pdata)) {
        LOG(ERROR) << "Kill/Restart received";
        thread_aborted = true;
        cv.notify_one();
        return;
    } else {
        pdata = getResult(maybe_pdata);
        DBG("Received and deserialized public data.");

        // Hashing Public Data
        pdata.modulusIdx = 0;
        pdata.roots.clear();
        pdata.ptNumberEvaluation = ptNumberEvaluation;
        // ligero::dumpPublicData(pdata);
        hashedPublicData = hash::hashPublicData(stringGenerator, pdata, spdata);
    }

    // Verify Sigma Protocol
   
    LOG(INFO) << "Verifying Sigma protocol proof";

    auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3210, primesGCD,175);
    std::vector<mpz_class> coefs;
    std::vector<mpz_class> nfl_primes;
    std::vector<mpz_class> nfl_primes_full;
    for(int i = 0; i < 9;i ++) nfl_primes.emplace_back(mpz_class(nfl::params<uint64_t>::P[i]));
    for(int i = 0; i < NbPrimesQ;i ++) nfl_primes_full.emplace_back(mpz_class(nfl::params<uint64_t>::P[i]));
    
    // Recomputing the final candidate from the CRT factors
    std::vector<mpz_class> finalModuli_GCD;
    finalModuli_GCD.resize(alphasGCD.size());
    for(int alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++)
    {
	    std::vector<mpz_class> prime_decomp_f;
	    for(int i = 0; i < NbPrimesQ; i ++)
	    {
		    prime_decomp_f.emplace_back(mpz_class(pdata.finalModuli_GCD[alphaGCDidx * NbPrimesQ + i]));
			    	    
	    }
	    finalModuli_GCD[alphaGCDidx] = math::crt_reconstruct(prime_decomp_f,coefs,nfl_primes_full);
    }
    std::vector<mpz_class> candidates;
    candidates.emplace_back(math::crt_reconstruct(finalModuli_GCD,coefs,alphasGCD));

    // Recomputing the gamma values 
    mpz_class gammaSeed = pdata.gammaSeed;
    std::vector<mpz_class> gammaValues;
    gammaValues.resize(candidates.size() * 128);
    std::vector<mpz_class> engine = math::generateRandomVector(gammaSeed,
                                    kJacobiNumberOfRandomValues, 2048);
        mpz_t one;
        mpz_init(one);
        mpz_set_ui(one, 1);
    size_t engineCount = 0;

    {
        mpz_t gcdResult;
        mpz_init(gcdResult);

        for (int i = 0; i < candidates.size(); ++i) {
            const mpz_class &N = candidates[i];

            DBG(" N % 2 = 1; N = " << N);
            // We expect that N is 3 (mod 4), and therefore must be odd.
            assert(N % 2 == 1);

            int foundGamma = 0;
            while (foundGamma < 128) {
                mpz_class g = engine[engineCount];
                ++engineCount;
                if (engineCount >= engine.size()) {
                    engineCount = 0;
                    DBG("Regenerating random values for Jacobi (GCD) gamma values");
                    engine = math::generateRandomVector(engine[0], kJacobiNumberOfRandomValues,
                                                        2048);
                }
                g = g % N;
                assert(g != 0);

                DBG("g = " << g);
                if (mpz_jacobi(g.get_mpz_t(), N.get_mpz_t()) == 1) {
                    gammaValues[i * 128 + foundGamma] = g;
                    foundGamma++;
                }
            }

            // Sanity check, assert if gcd(gammaValues[i], N) == 1
            for (int j = 0; j < 128; ++j) {
                mpz_gcd(gcdResult, gammaValues[i * 128 + j].get_mpz_t(),
                        N.get_mpz_t());
                assert(mpz_cmp(gcdResult, one) == 0);
            }
        }

        mpz_clear(one);
        mpz_clear(gcdResult);
    }
    for(int j = 0; j < 128; j++)
    {
	    std::vector<mpz_class> sigma_a_GCD;
	    sigma_a_GCD.resize(alphasGCD.size());
	    std::vector<mpz_class> sigma_e_GCD;
	    sigma_e_GCD.resize(alphasGCD.size());
	    std::vector<mpz_class> sigma_z_GCD;
	    sigma_z_GCD.resize(alphasGCD.size());
	    std::vector<mpz_class> sigma_g_GCD;
	    sigma_g_GCD.resize(alphasGCD.size());
	    for(int alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++)
	    {
		    std::vector<mpz_class> prime_decomp_a;
		    std::vector<mpz_class> prime_decomp_e;
		    std::vector<mpz_class> prime_decomp_z;
		    std::vector<mpz_class> prime_decomp_g;
		    for(int i = 0; i < 9; i ++)
		    {
			    prime_decomp_a.emplace_back(mpz_class(spdata.sigmaaGCD[i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx]));
			    prime_decomp_e.emplace_back(mpz_class(spdata.sigmaeGCD[i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx]));
			    prime_decomp_z.emplace_back(mpz_class(spdata.sigmazGCD[i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx]));
			    prime_decomp_g.emplace_back(mpz_class(spdata.sigmagGCD[i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx]));
		    }
		    sigma_a_GCD[alphaGCDidx] = math::crt_reconstruct(prime_decomp_a,coefs,nfl_primes);
		    sigma_e_GCD[alphaGCDidx] = math::crt_reconstruct(prime_decomp_e,coefs,nfl_primes);
		    sigma_z_GCD[alphaGCDidx] = math::crt_reconstruct(prime_decomp_z,coefs,nfl_primes);
		    sigma_g_GCD[alphaGCDidx] = math::crt_reconstruct(prime_decomp_g,coefs,nfl_primes);
	    }
	    mpz_class sigmaa = math::crt_reconstruct(sigma_a_GCD,coefs,alphasGCD);
	    mpz_class sigmae = math::crt_reconstruct(sigma_e_GCD,coefs,alphasGCD);
	    mpz_class sigmaz = math::crt_reconstruct(sigma_z_GCD,coefs,alphasGCD);
	    mpz_class sigmag = math::crt_reconstruct(sigma_g_GCD,coefs,alphasGCD);
            mpz_class sigma_lhs = math::powm(gammaValues[j], sigmaz, candidates[0]);
            mpz_class sigma_rhs = math::mod(sigmaa * math::powm(sigmag, sigmae, candidates[0]), candidates[0]);

            assert(sigma_lhs == sigma_rhs);

    }

    LOG(INFO) << "Verified Sigma protocol proof";


    for (size_t modulusIdx = 0; modulusIdx < q; modulusIdx++) {
        // Receive proofs
        auto maybe_proof = transport.awaitReply<ligero::zksnark::FullFSTranscript<FieldT>>(MessageType::GATHER_PROOFS, true);

        if (hasError(maybe_proof)) {
            LOG(ERROR) << "Kill/Restart received";
            thread_aborted = true;
            cv.notify_one();
            return;
        } else {
            // get lock and notify verifier
            std::unique_lock<std::mutex> lock(mutex);

            transcripts.push_back(getResult(maybe_proof));
            availableTranscriptsIndex++;

            lock.unlock();
            cv.notify_one();
        }
    }

    LOG(INFO) << "Got all proofs";

    return;
}

expected<Unit> verify(ZeroMQClientTransport& transport, timers& myTimers) {

    // Reset all global variables
    resetGlobalVariables();

    // Verifying
    // ==============================================================================
    FieldT::modulusIdx_ = 0;
    size_t aggregateSuccess = 0;

    Statements statement = Statements::RSACeremony;
    Statements localStatement = statement;

    // Primitive Roots of Unity for All Moduli
    nfl::poly_p<uint64_t, degree, q> roots({0});
    for (size_t idx = 0; idx < q; idx++) {roots(idx,1) = 1ull;}
    roots.ntt_pow_phi();

    ligero::zksnark::FullFSTranscript<FieldT> transcript;

    SecretData sdata;

    // Dimensioning ai scalars
    sdata.ai.resize(nbIterations);
    for (size_t iter = 0; iter < nbIterations; iter++) {
        sdata.ai[iter].resize(nbVarsStitched);
        for (size_t var = 0; var < nbVarsStitched; var++) {
            sdata.ai[iter][var].resize(nbBlocksPerVar);
            for (size_t block = 0; block < nbBlocksPerVar; block++) {
                sdata.ai[iter][var][block].resize(l);
            }
        }
    }

    // Storing the results of the verification per modulus
    std::vector<int> report(q, 0); 

    // Import Public Data
    std::thread* dataCollection = new std::thread(collectData, std::ref(transport));

    // CRT representation to check integrity fofthe proofs
    std::vector<mpz_class> stitchingCRT;
    mpz_class stitching0;

    for (size_t modulusIdx = 0; modulusIdx< q; modulusIdx++) {

        std::string timerVerificationIdle = std::string("6.a. verification, idle, modulus ") + std::to_string(modulusIdx);
        myTimers.begin(1,timerVerificationIdle.c_str(),0);

        {
            // get lock and wait for complete
            std::unique_lock<std::mutex> lock(mutex);
            auto received = cv.wait_for(lock, receive_timeout, [&]{ return thread_aborted || (availableTranscriptsIndex > modulusIdx); });

            // check if thread failed or timeout
            // we won't restart this verifier because the corresponding party will be kicked out
            if (thread_aborted || !received) {
                LOG(INFO) << "Working thread failed or timed out, quit verifier";
                dataCollection->join();
                delete dataCollection;
                return Error::TIMED_OUT;
            }

            lock.unlock();
        }

        myTimers.end(1,timerVerificationIdle.c_str(),0, true);

        std::string timerVerificationCompute = std::string("6.a. verification, compute, modulus ") + std::to_string(modulusIdx);
        myTimers.begin(1,timerVerificationCompute.c_str(),0);

        auto& transcript = transcripts[modulusIdx];
        bool success = true;

        LOG(INFO) << "Verifying modulus idx:" << modulusIdx;

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

        zksnark::Verifier<FieldT> verifier(localStatement, irs, t, l, fieldGenerator, stringGenerator, mtParams, pdata, spdata, sdata, permut<k>::compute, permut<n>::compute, hashedPublicData);

        try {
            verifier.verifyConsolidatedZKSnarkFile(transcript);
            
            LOG(INFO) << "Verifier mod " << modulusIdx << ":" << verifier._stitching;
            if (modulusIdx==0) {
                stitching0 = verifier._stitching;
            } else {
                stitchingCRT.push_back(verifier._stitching);
            }
        }
        catch(failedTest& e) {
            DBG("Failed Test: " << e.what());
            success = false;
        }

        if (success) {
            LOG(INFO) << "Verification for modulus idx:" << modulusIdx << " succeeded";
            report[modulusIdx] = 1;
            aggregateSuccess++;
        } else {
            LOG(INFO) << "Verification for modulus idx:" << modulusIdx << " failed";
            report[modulusIdx] = 0;
        }

        myTimers.end(1,timerVerificationCompute.c_str(),0, true);

    }

    // Cross-Moduli Verification
    {
        bool success = false;
        std::vector<mpz_class> mods, coefs;
        for (size_t i = 1; i < q; i++) {mods.emplace_back(params<uint64_t>::P[i]);}

        success = (math::mod(stitching0 - math::crt_reconstruct(stitchingCRT,coefs,mods), mpz_class(params<uint64_t>::P[0])) == mpz_class(0));

        if (success) {
            LOG(INFO) << "Verification for cross-moduli integrity succeeded";
        } else {
            LOG(INFO) << "Verification for cross-moduli integrity failed";
            for (size_t idx = 0; idx < q; idx++) {report[idx] = 0;}
            aggregateSuccess = 0;
        }
    }

    dataCollection->join();
    delete dataCollection;

    // Wait for the signal from the gathering coordinator
    auto maybe_response_reports = transport.awaitReply();

    if (hasError(maybe_response_reports)) {
        if (getError(maybe_response_reports) == Error::RESTART) {
            LOG(INFO) << "Restart verifier";
            myTimers.reset();
            return verify(transport, myTimers);
        }
        else {
            LOG(INFO) << "Timed out";
            return Unit{};
        }
    } else {
        if (getResult(maybe_response_reports) != MessageType::GATHER_REPORTS) {
            LOG(ERROR) << "Out of sequence message received when expecting gathering reports, aborting.";
            return Error::OUT_OF_SYNC;
        }
    }

    if (aggregateSuccess == q) { 
        LOG(INFO) << "All proofs successfully verified.";
    } else {
        LOG(INFO) << "Failed to verify some of the proofs.";
    }

    // Send the verifier report to the gathering coordinator
    transport.send<std::vector<int>>(MessageType::VERIFIER_REPORT,report);

    // Finally wait for restart message if any verification failed
    auto maybeRestart = transport.awaitReply();

    if (hasError(maybeRestart) && getError(maybeRestart) == Error::RESTART) {
        LOG(INFO) << "Restart message received, restarting verifier";
        myTimers.reset();
        return verify(transport, myTimers);
    }
    else if (getResult(maybeRestart) == MessageType::END) {
        LOG(INFO) << "All done.";
        return Unit{};
    }
    else {
        LOG(INFO) << "Unknown error happened!";
        return Error::UNKNOWN_ERROR;
    }
}

int main(int argc, char** argv)
{
    std::string coordinatorIP("127.0.0.1");

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("localIPaddress", po::value<std::string>(), "sets the local IP address")
        ("ip", po::value<std::string>(), "coordinator IP address");


    po::variables_map vm;
    po::positional_options_description p;
    p.add("coordinator-ip", -1);

    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);

    po::notify(vm);

    if (vm.count("help")) {
        LOG(INFO) << desc << "\n";
        LOG(INFO) << "USAGE: --localIPaddress <IP Distributed Verifier> <IP Coordinator>";
        return 0;
    }

    if (vm.count("localIPaddress")) {
        localIpAddr = vm["localIPaddress"].as<std::string>();
        LOG(INFO) << "Local IP Address " << localIpAddr;
    }

    if (vm.count("ip")) {
        coordinatorIP = vm["ip"].as<std::string>();
        LOG(INFO) << "Gathering Coordinator IP " << coordinatorIP;
    }

    SocketId socketId = boost::uuids::random_generator()();
    ZeroMQClientTransport transport(std::string("tcp://") + coordinatorIP + std::string(":") + std::string(port), socketId, receive_timeout);

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
    configureEasyLogging("distributed_verifier");

    // Add timers
    timers myTimers;
    myTimers.initialize_timer();

    auto result = verify(transport, myTimers);

    if (hasError(result)) {
        return EXIT_FAILURE;
    }

    return 0;
}
