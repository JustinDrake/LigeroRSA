/** Boost */
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/** STL */
#include <chrono>
#include <condition_variable>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

/** Ligero */
#include "Common.hpp"

#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"

#include "ZkArgument.hpp"
#include "protocol/ExpressNPStatement.hpp"

#include "FiniteFields.hpp"
#include "Hash.hpp"
#include "SecretSharingNTT.hpp"

using namespace ligero;
using namespace std::chrono_literals;

/** Parameters RSA Ceremony */
/** ======================= */
constexpr auto degree = 1 << 16;
constexpr auto p = 9;
constexpr auto q = 21;
constexpr auto sigma = 8;
constexpr auto lambda = 128;
constexpr auto tau_limit_bit = 175;
constexpr auto tau = 1000;
constexpr auto numCandidates = 2048;

/** throughput parties */
constexpr auto wait_timeout = 300 * 1000; /* wait 5 min when start */
constexpr auto timeout = 25 * 1000;       /* a little longer than test */
constexpr auto data_size = 32 * 1024;     /* 32 KB */
constexpr auto nb_max_send = 1024;
constexpr auto parties = 1;

/** throughtput coordinator */
constexpr auto pbs = 1000;
constexpr auto duration = 20s;
constexpr auto wait_time = 120s;
constexpr auto cleanup_time = 5s;
constexpr auto throughput_cutoff = 100 * 1024; /* 200 KB/s */

constexpr auto receive_timeout = 20min;

/** Parameters Zero-Knowledge */
/** ========================= */

/** Settings Zero-Knowledge */
RegularInteger p_zk(10514644122014433281ULL);
typedef Fpp_Fixed<p_zk> FieldT;

zksnark::InterleavedRSCode irs;        /** Reed-Solomon Interleaved Code */
ligero::MTParameters<FieldT> mtParams; /** Merkle Tree Parameters */
constexpr size_t l = 4096;             /** Bloc size */
constexpr size_t k = 8192;             /** Degree */
constexpr size_t n = 16384;            /** # of Shares */

/** Columns to Open Up (By default half the columns) */
constexpr size_t t = 256;
constexpr auto ptNumberEvaluation = 3;

/** Random Generator */
const hash::HashIntegerSeedGeneratorF<FieldT>& fieldGenerator =
    hash::blake2bIntegerSeedGenerator<FieldT>;
const hash::HashStringSeedGeneratorF& stringGenerator = hash::blake2bStringHash;

/** Threading */
constexpr size_t maxNumberThreads = 50;
/** [ a potential target machine, m5.metal has 96 threads over 48 cores] */

/** Selecting a port for distributed verifying */
constexpr char* port = "5556";

/** Depositing the data here */
PublicData* pd = new PublicData;
PublicData& pdata = *pd;
SigmaProtocolPublicData spdata;
std::vector<ligero::zksnark::FullFSTranscript<FieldT>> transcripts;
size_t availableTranscriptsIndex = 0;
bool firstIteration = true;

std::string localIpAddr = "UNKNOWN";

/** Stitching Proofs */
constexpr int nbIterations = 3;
constexpr int nbVarsStitched = 4;
constexpr int nbBlocksPerVar = degree / l;

/** Data */
std::string hashedPublicData;

expected<Unit> verify(std::string basename) {
  /** Verifying */
  /** ====================================================================== */
  FieldT::modulusIdx_ = 0;
  size_t aggregateSuccess = 0;

  Statements statement = Statements::RSACeremony;
  Statements localStatement = statement;

  LOG(INFO) << "Fetching Sigma Data";
  fetchFromFile<SigmaProtocolPublicData>(spdata,
                                         basename + std::string(".sigmadata"));

  LOG(INFO) << "Fetching Public Data";
  fetchFromFile<PublicData>(pdata, basename + std::string(".publicdata"));

  LOG(INFO) << "Processing";
  /** Hashing Public Data */
  pdata.modulusIdx = 0;
  pdata.roots.clear();
  pdata.ptNumberEvaluation = ptNumberEvaluation;
  hashedPublicData = hash::hashPublicData(stringGenerator, pdata, spdata);

  /** Primitive Roots of Unity for All Moduli */
  nfl::poly_p<uint64_t, degree, q> roots({0});
  for (size_t idx = 0; idx < q; idx++) {
    roots(idx, 1) = 1ull;
  }
  roots.ntt_pow_phi();

  SecretData sdata;

  LOG(INFO) << "Preparing Scalars";
  /** Dimensioning ai scalars */
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

  /** Storing the results of the verification per modulus */
  std::vector<int> report(q, 0);

  /** CRT representation to check integrity for the proofs */
  std::vector<mpz_class> stitchingCRT;
  mpz_class stitching0;

  for (size_t modulusIdx = 0; modulusIdx < q; modulusIdx++) {
    ligero::zksnark::FullFSTranscript<FieldT> transcript;

    LOG(INFO) << "Fetching Proof " << modulusIdx;
    fetchFromFile<ligero::zksnark::FullFSTranscript<FieldT>>(
        transcript,
        basename + std::string(".proof." + std::to_string(modulusIdx)));

    bool success = true;

    LOG(INFO) << "Verifying modulus idx:" << modulusIdx;

    // Zero-Knowledge
    p_zk = nfl::params<uint64_t>::P[modulusIdx];
    FieldT::modulusIdx_ = modulusIdx;
    FieldT::preprocessing(log2(n) + 1);

    pdata.modulusIdx = modulusIdx;

    /** For portions 7 to 12, only run the first 9 moduli */
    if ((localStatement == Statements::RSACeremony) && (pdata.modulusIdx > 8))
      localStatement = Statements::RSACeremony_Rounds3to6;

    /** Roots of Unity */
    pdata.roots.assign(
        roots.poly_obj().data() + pdata.modulusIdx * (size_t)roots.degree,
        roots.poly_obj().data() +
            (pdata.modulusIdx + 1) * (size_t)roots.degree);
    pdata.ptNumberEvaluation = ptNumberEvaluation;

    zksnark::Verifier<FieldT> verifier(
        localStatement, irs, t, l, fieldGenerator, stringGenerator, mtParams,
        pdata, spdata, sdata, permut<k>::compute, permut<n>::compute,
        hashedPublicData);

    try {
      verifier.verifyProof(transcript);

      LOG(INFO) << "Verifier mod " << modulusIdx << ":" << verifier._stitching;
      if (modulusIdx == 0) {
        stitching0 = verifier._stitching;
      } else {
        stitchingCRT.push_back(verifier._stitching);
      }
    } catch (failedTest& e) {
      DBG("Failed Test: " << e.what());
      success = false;
    }

    if (success) {
      LOG(INFO) << "Verification for modulus idx:" << modulusIdx
                << " succeeded";
      report[modulusIdx] = 1;
      aggregateSuccess++;
    } else {
      LOG(INFO) << "Verification for modulus idx:" << modulusIdx << " failed";
      report[modulusIdx] = 0;
    }
  }

  /** Cross-Moduli Verification */
  {
    bool success = false;
    std::vector<mpz_class> mods, coefs;
    for (size_t i = 1; i < q; i++) {
      mods.emplace_back(params<uint64_t>::P[i]);
    }

    success = (math::mod(stitching0 -
                             math::crt_reconstruct(stitchingCRT, coefs, mods),
                         mpz_class(params<uint64_t>::P[0])) == mpz_class(0));

    if (success) {
      LOG(INFO) << "Verification for cross-moduli integrity succeeded";
    } else {
      LOG(INFO) << "Verification for cross-moduli integrity failed";
      for (size_t idx = 0; idx < q; idx++) {
        report[idx] = 0;
      }
      aggregateSuccess = 0;
    }
  }

  if (aggregateSuccess == q) {
    LOG(INFO) << "All proofs successfully verified.";
  } else {
    LOG(INFO) << "Failed to verify some of the proofs.";
  }
}

int main(int argc, char** argv) {
  std::string proof;

  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "proof", po::value<std::string>(), "base filename for the proof");

  po::variables_map vm;
  po::positional_options_description p;

  po::store(
      po::command_line_parser(argc, argv).options(desc).positional(p).run(),
      vm);

  po::notify(vm);

  if (vm.count("help")) {
    LOG(INFO) << desc << "\n";
    LOG(INFO) << "USAGE: --proof <Proof basename>";
    return 0;
  }

  if (vm.count("proof")) {
    proof = vm["proof"].as<std::string>();
  }

  /** Set Zero-Knowledge Parameters */
  irs.k = k; /** Degree of the Polynomial */
  irs.n = n; /** Number of Shares */

  mtParams.digestLength = 32;
  mtParams.zkByteSize = 8;
  mtParams.hashInnerNodeF = hash::blake2bTwoToOneHash;
  mtParams.hashLeafContentF = hash::blake2bFieldElementHash<FieldT>;
  mtParams.hashZKRandomnessF = hash::blake2bZKElementHash;
  mtParams.leavesNumber = irs.n;

  //** Configure easylogging */
  configureEasyLogging("manual_verifier");

  if (hasError(verify(proof))) {
    LOG(TRACE) << "Failed to process the proofs, exiting now.";
    return EXIT_FAILURE;
  }

  return 0;
}
