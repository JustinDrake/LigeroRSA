#include <boost/program_options.hpp>
#include <thread>
#include <chrono>

#include "EncryptedCoordinator.hpp"
#include "LatticeEncryption.hpp"
#include "Ligero.hpp"

#include "Transport.hpp"
#include "ZkArgument.hpp"
#include "GatheringCoordinator.hpp"
#include "SecretSharingNTT.hpp"

#define NDEBUG 1

namespace po = boost::program_options;

//==============================================================================
/** A small helper function for wrapping the construction of a client and
 *  contributing to the computation.
 *
 *  This is primarily to ensure the socket construction happens on the same thread
 *  as the socket usage.
 */

using namespace ligero;
using namespace ligero::lattice;
using namespace std::chrono_literals;

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

// Parameters Zero-Knowledge
// ============================

Statements localStatement = Statements::RSACeremony;

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
const hash::HashIntegerSeedGeneratorF<FieldT>& randomGenerator = hash::blake2bIntegerSeedGenerator<FieldT>;
const hash::HashStringSeedGeneratorF& stringGenerator = hash::blake2bStringHash;

// Threading
constexpr size_t maxNumberThreads = 95; // [m5.metal has 96 threads over 48 cores]

//bool& success, 
void verifyProof(uint8_t& success, PublicData& pdata, ligero::zksnark::FullFSTranscript<FieldT>& transcript) {

        // Verifying
        // ==============================================================================
        success = 1;

        // For portions 7 to 12, only run the first 9 moduli
        if ((localStatement == Statements::RSACeremony)&&(pdata.modulusIdx>8)) localStatement = Statements::RSACeremony_Rounds3to6;

        SecretData sdata;
        zksnark::Verifier<FieldT> verifier(localStatement, irs, t, l, randomGenerator, stringGenerator, mtParams, pdata, sdata, permut<k>::compute, permut<n>::compute);

        try {
                verifier.verifyConsolidatedZKSnarkFile(transcript);
            }
            catch(failedTest& e) {
                DBG("Failed Test: " << e.what());
                success = 0;
            }
            catch(internalError& e) {
                DBG("Internal Error: " << e.what());
                success = 0;
            }
            catch(...) {
                DBG("Other Error");
                success = 0;
            }

        return;
}

/** An equivalent helper function for the coordinator side of things. */
void coordinate (int parties, std::string ipAddr, int numCandidates, bool passive) {
       
	try {
	    ProtocolConfig<uint64_t> config(parties, p, q, degree, sigma, lambda, tau_limit_bit, pbs);

        // Generate all proofs
        // =====================================================================================================================
        auto coordinator = EncryptedCoordinator<uint64_t, degree, p, q>();
        ZeroMQCoordinatorTransport transport(std::string("tcp://") + ipAddr + std::string(":5555"), config.numParties());
        coordinator.host(transport, config, {duration, wait_time, cleanup_time, data_size}, throughput_cutoff);

        if (!passive) {
            GatheringCoordinator::gatherAndSaveArgumentsOfKnowledge<FieldT>(transport, "Prod");
            LOG(INFO) << "Gathered all proofs";
        }
    }
    catch (zmq::error_t& e) {
        LOG(ERROR) << "Caught an error trying to coordinate... " << e.what();
    }
}

//==============================================================================
int main(int argc, char** argv)
{
    std::string ipAddr("127.0.0.1");
    std::string logfile("default");

    int numBits         = 100;
    int numParties      = 2;
    int numCandidates   = 2048;
    double chiStd       = 3.2;

    bool passive = false;

    // Set Zero-Knowledge Parameters 
    irs.k                       = k;       // Degree of the Polynomial
    irs.n                       = n;       // Number of Shares

     mtParams.digestLength       = 32;
     mtParams.zkByteSize         = 8;
     mtParams.hashInnerNodeF     = hash::blake2bTwoToOneHash;
     mtParams.hashLeafContentF   = hash::blake2bFieldElementHash<FieldT>;
     mtParams.hashZKRandomnessF  = hash::blake2bZKElementHash;
     mtParams.leavesNumber       = irs.n;

    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("version", "show version")
            ("parties", po::value<int>(), "set the number of participating clients")
            ("logfile", po::value<std::string>(), "sets the logfile name")
            ("ip", po::value<std::string>(), "set the ip address for the coordinator")
            ("bitsize", po::value<int>(), "sets the bit size")
            ("chiStd", po::value<double>(), "sets the std distribution randomness")
            ("numcandidates", po::value<int>(), "sets the number of shares attempted per participant")
            ("passive", "do not include a zero-knowledge argument");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            LOG(TRACE) << desc << "\n";
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "coordinator version " << std::hex << PROTOCOL_VERSION_NUMBER << std::endl;
            return 0;
        }

        if (vm.count("passive")) {
            passive = true;
        }

        if (vm.count("parties")) {
            numParties = vm["parties"].as<int>();
            LOG(TRACE) << "nb parties = " << numParties;
        }

        if (vm.count("logfile")) {
            logfile = vm["logfile"].as<std::string>();
            LOG(TRACE) << "logfile = " << logfile;
        }

       if (vm.count("ip")) {
            ipAddr = vm["ip"].as<std::string>();
            LOG(TRACE) << "Coordinator resides at " << ipAddr;
        }

       if (vm.count("numcandidates")) {
            numCandidates = (int)vm["numcandidates"].as<int>();
            LOG(TRACE) << "numCandidates = " << numCandidates;
        }

       if (vm.count("chiStd")) {
            chiStd = (double)vm["chiStd"].as<double>();
            LOG(TRACE) << "chiStd = " << chiStd;
            std::cout << chiStd;
        }

        const int packingFactor = numCandidates;

        // Configuring Easylogging
        if (doesFileExist("logging")) {
            el::Configurations conf("logging");
            conf.set(el::Level::Global,el::ConfigurationType::Filename, logfile);
            el::Loggers::reconfigureAllLoggers(conf);
        } else {
            el::Configurations conf;
            conf.setToDefault();
            conf.parseFromText("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s} %msg\n TO_FILE = false\n TO_STANDARD_OUTPUT = true\n*DEBUG:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = true\n*INFO:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = true");
            el::Loggers::reconfigureAllLoggers(conf);
        }

        LOG(TRACE) << std::string(120, '=');
        LOG(TRACE) << "Coordinator started.";

        coordinate(numParties, ipAddr, numCandidates, passive);
    }
    catch(std::exception& e) {
        LOG(FATAL) << "Error: " << e.what();
    }
    catch(...) {
        LOG(FATAL) << "Caught an exception of unknown type.";
    }

    return 0;
}
