#include <boost/program_options.hpp>
#include <chrono>
#include <thread>

#include "EncryptedClient.hpp"
#include "SecretSharingNTT.hpp"
#include "ZkArgument.hpp"

#include "Ligero.hpp"

#include "Common.hpp"
#include "Hash.hpp"

namespace po = boost::program_options;

using namespace ligero;
using namespace ligero::lattice;
using namespace std::chrono_literals;

/** Parameters RSA Ceremony */
constexpr auto degree = 1 << 16;
constexpr auto p = 9;
constexpr auto q = 21;
constexpr auto sigma = 8;
constexpr auto lambda = 128;
constexpr auto tau_limit_bit = 175;
constexpr auto tau = 1000;
constexpr auto numCandidates = 2048;

/** throughput parties */
constexpr auto wait_timeout = 1200 * 1000; /** wait 20 min when start */
constexpr auto timeout = 25 * 1000;        /** a little longer than test */
constexpr auto data_size = 32 * 1024;      /** 32 KB */
constexpr auto nb_max_send = 1024;

/** general timeout */
constexpr auto receive_timeout = 60min;

/** Parameters Zero-Knowledge */

/** Settings Zero-Knowledge */
RegularInteger p_zk(10514644122014433281ULL);
typedef Fpp_Fixed<p_zk> FieldT;

zksnark::InterleavedRSCode irs;        /** Reed-Solomon Interleaved Code */
ligero::MTParameters<FieldT> mtParams; /** Merkle Tree Parameters */
constexpr size_t l = 8192;             // 4096;             /** Bloc size */
constexpr size_t k = 16384;            // 8192;             /** Degree */
constexpr size_t n = 65536;            // 16384;            /** # of Shares */

Statements statement = Statements::RSACeremony;

/** Columns to Open Up (By default half the columns) */
constexpr size_t t = 189;
constexpr auto ptNumberEvaluation = 3;

/** Stitching Proofs */
constexpr int nbIterations = 3;
constexpr int nbVarsStitched = 4;
constexpr int nbBlocksPerVar = degree / l;

void compromise(uint64_t& u) {
  int r = rand() % 10;

  if (r == 3) {
    u += ligero::math::generateRandomValue(mpz_class(std::random_device()()))
             .get_ui();
  }
};

void compromise(std::vector<uint64_t>& u) {
  int r = rand() % 10;

  if (r == 3) {
    if (u.size() != 0) {
      size_t pos =
          ligero::math::generateRandomValue(mpz_class(std::random_device()()))
              .get_ui() %
          u.size();

      u[pos] +=
          ligero::math::generateRandomValue(mpz_class(std::random_device()()))
              .get_ui();
    }
  }
};

void compromise(mpz_class& u) {
  int r = rand() % 10;

  if (r == 3) {
    u += ligero::math::generateRandomValue(mpz_class(std::random_device()()));
  }
};

/** This function handles the ceremony itself, whether it is performed with
 * passive or active security.
 * @param transport transport object to communicate with coordinator
 * @param socketId unique identifier for this party
 * @param localIpAddr local ip address for this party
 * @param passive boolean flag indicating whether to generate a ZK proof or
 * not
 * @return Unit if successful
 */
expected<Unit> participateHelper(ZeroMQClientTransport& transport,
                                 SocketId socketId, std::string localIpAddr,
                                 bool passive) {
  auto client = new EncryptedClient<FieldT, uint64_t, degree, p, q>(socketId);
  expected<Unit> result = Unit{};

  if ((!passive) && (statement == Statements::RSACeremony)) {
    Statements localStatement = statement;

    /** Zero-Knowledge Argument */

    /** Setting Protocol Parameters */
    /** =========================== */

    /** Reed-Solomon Interleaved Code */
    zksnark::InterleavedRSCode irs;

    irs.k = k;  // Degree of the Polynomial
    irs.n = n;  // Number of Shares

    /** Merkle Tree Parameters */
    ligero::MTParameters<FieldT> mtParams;

    mtParams.digestLength = 32;
    mtParams.zkByteSize = 8;
    mtParams.hashInnerNodeF = hash::blake2bTwoToOneHash;
    mtParams.hashLeafContentF = hash::blake2bFieldElementHash<FieldT>;
    mtParams.hashZKRandomnessF = hash::blake2bZKElementHash;
    mtParams.leavesNumber = irs.n;

    /** Fiat-Shamir Seed Generators */
    const hash::HashIntegerSeedGeneratorF<FieldT>& fieldGenerator =
        hash::blake2bIntegerSeedGenerator<FieldT>;
    const hash::HashStringSeedGeneratorF& stringGenerator =
        hash::blake2bStringHash;

    PublicData pdata;
    SigmaProtocolPublicData spdata;
    SecretData sdata;

    /** No Fiat-Shamir in Early Witness */
    std::string hashedPublicData("");

    /** Core */
    pdata.degree = degree;
    pdata.p = p;
    pdata.q = q;
    pdata.sigma = sigma;
    pdata.lambda = lambda;
    pdata.tau_limit_bit = tau_limit_bit;
    pdata.tau = tau;

    zksnark::Prover<FieldT> prover(
        localStatement, irs, t, l, fieldGenerator, stringGenerator, mtParams,
        pdata, spdata, sdata, permut<k>::compute, permut<n>::compute, transport,
        hashedPublicData);

    result = client->start(transport, sigma, lambda, localIpAddr, &prover);
  } else {
    result = client->start(transport, sigma, lambda, localIpAddr);
  }

  if (hasError(result)) {
    LOG(INFO) << "Party failed with error message: "
              << showError(getError(result));
    exit(EXIT_FAILURE);
  }

  if (passive) {
    delete client;
    return Unit{};
  } else {
    LOG(INFO) << "Generating Zero-Knowledge Argument.";
    transport.myTimers.begin(1, "Active Security", 1);

    /** Zero-Knowledge Argument */

    /** Setting Protocol Parameters */
    /** =========================== */

    /** Reed-Solomon Interleaved Code */
    zksnark::InterleavedRSCode irs;

    irs.k = k; /** Degree of the Polynomial */
    irs.n = n; /** Number of Shares */

    /** Merkle Tree Parameters */
    ligero::MTParameters<FieldT> mtParams;

    mtParams.digestLength = 32;
    mtParams.zkByteSize = 8;
    mtParams.hashInnerNodeF = hash::blake2bTwoToOneHash;
    mtParams.hashLeafContentF = hash::blake2bFieldElementHash<FieldT>;
    mtParams.hashZKRandomnessF = hash::blake2bZKElementHash;
    mtParams.leavesNumber = irs.n;

    /** Fiat-Shamir Seed Generators */
    const hash::HashIntegerSeedGeneratorF<FieldT>& fieldGenerator =
        hash::blake2bIntegerSeedGenerator<FieldT>;
    const hash::HashStringSeedGeneratorF& stringGenerator =
        hash::blake2bStringHash;

    PublicData pdata;
    SigmaProtocolPublicData spdata;
    SecretData sdata;

    client->gatherData(pdata, spdata, sdata);
    transport.send(MessageType::GATHER_SIGMA_DATA, spdata);

    /** Set all parameters for active and passive protocols */

    delete client;

    /** core */
    pdata.degree = degree;
    pdata.p = p;
    pdata.q = q;
    pdata.sigma = sigma;
    pdata.lambda = lambda;
    pdata.tau_limit_bit = tau_limit_bit;
    pdata.tau = tau;

    /** Hashing Public Data */
    pdata.modulusIdx = 0;
    pdata.roots.clear();
    pdata.ptNumberEvaluation = ptNumberEvaluation;

    /** Compromise some of the public data **/
    compromise(pdata.ps);
    compromise(pdata.A);
    compromise(pdata.bi);
    compromise(pdata.b);
    compromise(pdata.ci_1);
    compromise(pdata.ci_2);
    compromise(pdata.ci_1_prime);
    compromise(pdata.ci_2_prime);
    compromise(pdata.c_1);
    compromise(pdata.c_2);
    compromise(pdata.di);
    compromise(pdata.c_1_prime);
    compromise(pdata.c_2_prime);
    compromise(pdata.ax_shares);
    compromise(pdata.by_shares);
    compromise(pdata.fixedCoefs);
    compromise(pdata.coefsCAN);
    compromise(pdata.cans);
    compromise(pdata.prodcans);
    compromise(pdata.ax);
    compromise(pdata.by);
    compromise(pdata.axby);
    compromise(pdata.ax_shares_GCD);
    compromise(pdata.by_shares_GCD);
    compromise(pdata.coefsGCD);
    compromise(pdata.gcds);
    compromise(pdata.prodgcds);
    compromise(pdata.axGCD);
    compromise(pdata.byGCD);
    compromise(pdata.axbyGCD);
    compromise(pdata.finalModuli_GCD);
    compromise(pdata.roots);
    compromise(pdata.sigmaeGCD);
    compromise(pdata.sigmazGCD);
    compromise(pdata.gammaSeed);
    compromise(pdata.log2Ds);
    compromise(pdata.Cs);
    compromise(pdata.Ds);

    /** Compromise some of the sigma public data **/
    compromise(spdata.sigmaeGCD);
    compromise(spdata.sigmazGCD);
    compromise(spdata.sigmaaGCD);
    compromise(spdata.sigmagGCD);

    /** Compromise some of the secret data **/
    compromise(sdata.siP);
    compromise(sdata.eiP);
    compromise(sdata.q_r3);
    compromise(sdata.xi);
    compromise(sdata.ux);
    compromise(sdata.vx);
    compromise(sdata.wx);
    compromise(sdata.vxP);
    compromise(sdata.wxP);
    compromise(sdata.q_r4_1);
    compromise(sdata.q_r4_2);
    compromise(sdata.zp);
    compromise(sdata.yi);
    compromise(sdata.zi);
    compromise(sdata.uz);
    compromise(sdata.vz);
    compromise(sdata.vzP);
    compromise(sdata.wz);
    compromise(sdata.wzP);
    compromise(sdata.q_r5_1);
    compromise(sdata.q_r5_2);
    compromise(sdata.rNoise);
    compromise(sdata.q_r6);
    compromise(sdata.x_sharesCAN);
    compromise(sdata.x_sharesPS);
    compromise(sdata.y_sharesCAN);
    compromise(sdata.y_sharesPS);
    compromise(sdata.z_sharesCAN);
    compromise(sdata.q_p_prod_r7);
    compromise(sdata.q_p_r7);
    compromise(sdata.q_q_prod_r7);
    compromise(sdata.q_q_r7);
    compromise(sdata.q_r8);
    compromise(sdata.x_sharesGCD);
    compromise(sdata.y_sharesGCD);
    compromise(sdata.z_sharesGCD);
    compromise(sdata.q_p_prod_r11);
    compromise(sdata.q_q_prod_r11);
    compromise(sdata.q_pq_r11);
    compromise(sdata.r_CRTs);
    compromise(sdata.q_r_r11);
    compromise(sdata.ss_GCD);
    compromise(sdata.q_r12);
    compromise(sdata.alpha);
    compromise(sdata.W1);
    compromise(sdata.W2);
    compromise(sdata.XplusCs);
    compromise(sdata.XplusDs);
    compromise(sdata.Xs);
    compromise(sdata.sigmarGCD);
    compromise(sdata.sigmaxGCD);
    compromise(sdata.sigmaqGCD);
    compromise(sdata.exp_q);
    compromise(sdata.expqGCD);

    std::string hashedPublicData =
        hash::hashPublicData(stringGenerator, pdata, spdata);

    /** Primitive Roots of Unity for All Moduli */
    nfl::poly_p<uint64_t, degree, q> roots({0});
    for (size_t idx = 0; idx < q; idx++) {
      roots(idx, 1) = 1ull;
    }
    roots.ntt_pow_phi();

    LOG(INFO) << "Generating Zero-Knowledge Argument.";

    /** Compute Stitching */
    /** Variables we store across proofs */
    mpz_class alpha =
        math::generateRandomValue(mpz_class(std::random_device()()), 210);
    uint64_t W1 = ligero::math::scaleToPrime(
        math::generateRandomValue(mpz_class(std::random_device()()), 62)
            .get_ui(),
        0);
    uint64_t W2 = ligero::math::scaleToPrime(
        math::generateRandomValue(mpz_class(std::random_device()()), 62)
            .get_ui(),
        0);
    std::vector<uint64_t> rs(4 * degree);

    /** Preliminary checks */
    {
      /** Compute Stitching */
      std::vector<mpz_class> mods, alpha_CRT;
      for (size_t i = 1; i < q; i++) {
        mods.emplace_back(params<uint64_t>::P[i]);
      }

      /** CRT deconstruction */
      alpha_CRT = math::crt_deconstruct(
          alpha * mpz_class(params<uint64_t>::P[0]), mods);
      LOG(INFO) << "CRT 1:" << alpha_CRT[0];
      LOG(INFO) << "CRT 2:" << alpha_CRT[1];

      /** CRT reconstruction */
      std::vector<mpz_class> coefs;
      mpz_class reconstructed_alpha =
          math::crt_reconstruct(alpha_CRT, coefs, mods);

      assert(math::mod(reconstructed_alpha -
                           alpha * mpz_class(params<uint64_t>::P[0]),
                       mpz_class(params<uint64_t>::P[0])) == mpz_class(0));
      LOG(INFO) << math::mod(
          reconstructed_alpha - alpha * mpz_class(params<uint64_t>::P[0]),
          mpz_class(params<uint64_t>::P[0]));
      LOG(INFO) << "CRT Check Successful";
    }

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

    for (pdata.modulusIdx = 0; pdata.modulusIdx < q; pdata.modulusIdx++) {
      Statements localStatement = statement;
      LOG(INFO) << "Generating proof for modulus idx:" << pdata.modulusIdx;

      /** Roots of Unity */
      pdata.roots.assign(
          roots.poly_obj().data() + pdata.modulusIdx * (size_t)roots.degree,
          roots.poly_obj().data() +
              (pdata.modulusIdx + 1) * (size_t)roots.degree);

      /** Zero-Knowledge */
      p_zk = nfl::params<uint64_t>::P[pdata.modulusIdx];
      FieldT::modulusIdx_ = pdata.modulusIdx;
      FieldT::preprocessing(log2(n) + 1);

      /** Cross-Moduli Integrity */
      sdata.alpha = math::mod(alpha, mpz_class(p_zk)).get_ui();
      sdata.W1 = math::mod(W1, mpz_class(p_zk)).get_ui();
      sdata.W2 = math::mod(W2, mpz_class(p_zk)).get_ui();

      /** For portions 7 to 12, only run the first 9 moduli */
      if ((statement == Statements::RSACeremony) && (pdata.modulusIdx > 8))
        localStatement = Statements::RSACeremony_Rounds3to6;

      zksnark::Prover<FieldT> prover(
          localStatement, irs, t, l, fieldGenerator, stringGenerator, mtParams,
          pdata, spdata, sdata, permut<k>::compute, permut<n>::compute,
          transport, hashedPublicData);

      prover.produceArgumentOfKnowledge(1, false);
    }

    LOG(INFO) << "Waiting for possible restart";
    auto command = transport.awaitReply();

    if (!hasError(command) && (getResult(command) == MessageType::END)) {
      transport.myTimers.end(2, "1.b. Overall speed", 1, true);
      transport.myTimers.end(1, "Active Security", 1, true);
      return Unit{};
    } else if (hasError(command) && (getError(command) == Error::RESTART)) {
      LOG(INFO) << "Restart protocol";
      transport.myTimers.reset();
      return participateHelper(transport, socketId, localIpAddr, passive);
    } else {
      LOG(INFO) << "Error happen, quit";
      exit(EXIT_FAILURE);
    }
  }

  return Unit{};
}

/** A small helper function for setting up the ceremony, handling optional
 * preliminary steps such as throughput test. This is primarily to ensure the
 * socket construction happens on the same thread as the socket usage.
 *
 * @param ipAddr ip address of the coordinator for the ceremony
 * @param locaIpAddr ip address of the party joining the ceremony, for
 * logging/references only
 * @param passive instructs whether or not to generate a ZKP
 * @param speedtest instructs whether or not to run a throughput test ahead of
 * the ceremony
 */
expected<Unit> participate(std::string ipAddr, std::string localIpAddr,
                           bool passive, bool speedtest) {
  try {
    // RSA Ceremony
    SocketId socketId = boost::uuids::random_generator()();
    ZeroMQClientTransport transport(
        std::string("tcp://") + ipAddr + std::string(":5555"), socketId,
        receive_timeout);

    if (speedtest) {
      auto proceed = transport.joinThroughputTest(wait_timeout, timeout,
                                                  data_size, nb_max_send);

      if (hasError(proceed)) {
        LOG(INFO) << "Kicked out due to poor network connection";
        LOG(INFO) << "Error message is: " << showError(getError(proceed));
        return Error::KILLED_BY_COORDINATOR;
      }
    }

    participateHelper(transport, socketId, localIpAddr, passive);
  } catch (zmq::error_t& e) {
    LOG(ERROR) << "Caught an error trying to participate... " << e.what();
  }
  return Unit{};
}

/**
 * =================================================================================
 */
int main(int argc, char** argv) {
  std::string ipAddr("127.0.0.1");
  std::string localIpAddr("");

  bool passive = false;
  bool speedtest = false;

  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")("version", "show version")(
      "ip", po::value<std::string>(), "set the ip address for the coordinator")(
      "localIPaddress", po::value<std::string>(), "sets the local IP address")(
      "passive", "do not include a zero-knowledge argument")(
      "speedtest", "perform throughput test before kick off the ceremony");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    LOG(INFO) << desc << "\n";
    return 0;
  }

  if (vm.count("passive")) {
    passive = true;
  }

  if (vm.count("speedtest")) {
    speedtest = true;
    LOG(TRACE) << "Running with throughput test, make sure running the "
                  "coordinator with `--speedtest` also";
  }

  if (vm.count("version")) {
    std::cout << "party version " << std::hex << PROTOCOL_VERSION_NUMBER
              << std::endl;
    return 0;
  }

  if (vm.count("ip")) {
    ipAddr = vm["ip"].as<std::string>();
    LOG(INFO) << "Coordinator resides at " << ipAddr;
  }

  if (vm.count("localIPaddress")) {
    localIpAddr = vm["localIPaddress"].as<std::string>();
    LOG(INFO) << "Local IP Address" << localIpAddr;
  }

  /** Configure easylogging */
  configureEasyLogging("malicious_party");

  auto result = participate(ipAddr, localIpAddr, passive, speedtest);
  if (hasError(result)) {
    return 1;
  }
  return 0;
}
