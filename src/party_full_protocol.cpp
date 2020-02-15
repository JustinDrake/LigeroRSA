#include <boost/program_options.hpp>
#include <thread>

#include "EncryptedClient.hpp"
#include "ZkArgument.hpp"
#include "SecretSharingNTT.hpp"

#include "Ligero.hpp"

#define NDEBUG 1

namespace po = boost::program_options;

using namespace ligero;
using namespace ligero::lattice;

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

// Parameters Zero-Knowledge
// ============================

// Settings Zero-Knowledge
RegularInteger p_zk(10514644122014433281ULL);
typedef Fpp_Fixed<p_zk> FieldT;

zksnark::InterleavedRSCode irs;                     // Reed-Solomon Interleaved Code
ligero::MTParameters<FieldT> mtParams;              // Merkle Tree Parameters
constexpr size_t l = 4096;                                // Bloc size
constexpr size_t k = 8192;                                // Degree
constexpr size_t n = 16384;                               // # of Shares

Statements statement = Statements::RSACeremony;

constexpr size_t t = 256;                                 // Columns to Open Up (By default half the columns)
constexpr auto ptNumberEvaluation = 3;

//==============================================================================
/** A small helper function for wrapping the construction of a client and
 *  contributing to the computation.
 *
 *  This is primarily to ensure the socket construction happens on the same thread
 *  as the socket usage.
 */
void participate (std::string ipAddr, std::string localIpAddr, int numBits, int numCandidates, bool passive) {
    try {

        // RSA Ceremony
        SocketId socketId = boost::uuids::random_generator()();
 
        auto client = new EncryptedClient<FieldT, uint64_t, degree, p, q>(socketId);
        ZeroMQClientTransport transport(std::string("tcp://") + ipAddr + std::string(":5555"), socketId);

        client->start(transport, wait_timeout, timeout, data_size, nb_max_send, localIpAddr);

        if (!passive) {
            // Zero-Knowledge Argument

            // Setting Protocol Parameters
            // ===========================

            // Reed-Solomon Interleaved Code
            zksnark::InterleavedRSCode irs;

            irs.k       = k;            // Degree of the Polynomial
            irs.n       = n;            // Number of Shares

            // Merkle Tree Parameters
            ligero::MTParameters<FieldT> mtParams;

            mtParams.digestLength       = 32;
            mtParams.zkByteSize         = 8;
            mtParams.hashInnerNodeF     = hash::blake2bTwoToOneHash;
            mtParams.hashLeafContentF   = hash::blake2bFieldElementHash<FieldT>;
            mtParams.hashZKRandomnessF  = hash::blake2bZKElementHash;
            mtParams.leavesNumber       = irs.n;

            // Fiat-Shamir Seed Generators
            const hash::HashIntegerSeedGeneratorF<FieldT>& fieldGenerator = hash::blake2bIntegerSeedGenerator<FieldT>;
            const hash::HashStringSeedGeneratorF& stringGenerator = hash::blake2bStringHash;

            PublicData pdata;
            SecretData sdata;

            // Set all parameters for active and passive protocols
            client->gatherData(pdata, sdata);

            delete client;

            // core
            pdata.degree = degree;
            pdata.p = p;
            pdata.q = q; 
            pdata.sigma = sigma;
            pdata.lambda = lambda;
            pdata.tau_limit_bit = tau_limit_bit; 
            pdata.tau = tau;

            // Primitive Roots of Unity for All Moduli
            nfl::poly_p<uint64_t, degree, q> roots({0});
            for (size_t idx = 0; idx < q; idx++) {roots(idx,1) = 1ull;}
            roots.ntt_pow_phi();

            LOG(INFO) << "Generating Zero-Knowledge Argument."; 

            for (pdata.modulusIdx = 0; pdata.modulusIdx < q; pdata.modulusIdx++) {
                Statements localStatement = statement;
                LOG(INFO) << "Modulus idx:" << pdata.modulusIdx; 

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
        }
    }
    catch (zmq::error_t& e) {
        LOG(ERROR) << "Caught an error trying to participate... " << e.what();
    }

}

//==============================================================================
int main(int argc, char** argv)
{
    std::string ipAddr("127.0.0.1");
    std::string localIpAddr("");
    int bitsize        = 100;
    int numCandidates  = 2048;

    bool passive = false;

    // Configuring Easylogging
    if (doesFileExist("logging")) {
        el::Configurations conf("logging");
        el::Loggers::reconfigureAllLoggers(conf);
    } else {
        el::Configurations conf;
        conf.setToDefault();
        conf.parseFromText("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s} %msg\n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = false\n*INFO:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = true\n*ERROR:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = true");
        el::Loggers::reconfigureAllLoggers(conf);
    }

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("version", "show version")
        ("ip", po::value<std::string>(), "set the ip address for the coordinator")
        ("bitsize", po::value<int>(), "sets the bit size")
        ("localIPaddress", po::value<std::string>(), "sets the local IP address")
        ("numcandidates", po::value<int>(), "sets the number of shares attempted per participant")
        ("passive", "do not include a zero-knowledge argument");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        LOG(INFO) << desc << "\n";
        LOG(INFO) << "bitsize 8 16 32 64 128 256";
        return 0;
    }

    if (vm.count("passive")) {
        passive = true;
    }

    if (vm.count("version")) {
        std::cout << "party version " << std::hex << PROTOCOL_VERSION_NUMBER << std::endl;
        return 0;
    }

    if (vm.count("ip")) {
        ipAddr = vm["ip"].as<std::string>();
        LOG(INFO) << "Coordinator resides at " << ipAddr;
    }

    if (vm.count("bitsize")) {
        bitsize = vm["bitsize"].as<int>();
        LOG(INFO) << "bitsize = " << bitsize;
    }

    if (vm.count("numcandidates")) {
        numCandidates = vm["numcandidates"].as<int>();
        LOG(INFO) << "Num candidates" << numCandidates;
    }

    if (vm.count("localIPaddress")) {
        localIpAddr = vm["localIPaddress"].as<std::string>();
        LOG(INFO) << "Local IP Address" << localIpAddr;
    }

    participate(ipAddr, localIpAddr, bitsize, numCandidates, passive);
    return 0;
}
