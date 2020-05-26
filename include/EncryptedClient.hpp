#pragma once

#include <vector>

#include "Transport.hpp"
#include "Math.hpp"
#include "Common.hpp"
#include "LatticeEncryption.hpp"
#include <stdlib.h>
#include <optional>
#include "FiniteFields.hpp"
#include "Hash.hpp"
#include "ZkArgument.hpp"
 
#include "eccrypto.h"
#include "osrng.h"
#include "oids.h"

extern RegularInteger p_zk;
using namespace CryptoPP;

namespace ligero {

template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
class EncryptedClient {
    public:
        EncryptedClient(SocketId& sid) : _socketId(sid) {};
        using Q = nfl::poly_p<T, Degree, NbPrimesQ>;
        using P = nfl::poly_p<T, Degree, NbPrimesP>;
        using cipherText = typename
                           lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>::cipherText;

        expected<lattice::key_pair<T, Degree, NbPrimesQ>> generateKeyPair(
                    ZeroMQClientTransport &trans, nfl::gaussian<uint16_t, T, 2> &chi);
        std::vector<std::vector<mpz_class>> pruneAndReorderShares(
                                             const std::vector<mpz_class> &as,
                                             const std::vector<mpz_class> &bs,
                                             const std::vector<int> &flags,
                                             int numAlphas,
                                             std::vector<size_t> bucketSize
                                         );

        expected<ProtocolConfig<T>> registerAndAwaitConfiguration(ZeroMQClientTransport &transport, std::string ipAddress, nfl::gaussian<uint16_t, T, 2>& chi, zksnark::Prover<FieldT>* prover);

        expected<std::tuple<
        std::vector<mpz_class>, // xcan
            std::vector<mpz_class>, // ycan
            std::vector<mpz_class>, // zcan
            std::vector<mpz_class>, // xgcd
            std::vector<mpz_class>, // ygcd
            std::vector<mpz_class>, // zgcd
            std::vector<mpz_class> // both_shares
            >>
            generatePreSievingShares(
                lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ> &e,
                ZeroMQClientTransport &transport,
                const lattice::pub<T, Degree, NbPrimesQ> &publicKey,
                const lattice::sec<T, Degree, NbPrimesQ> &secretKey,
                bool special,
                ProtocolConfig<T> &config
            );

        expected<std::tuple<std::vector<mpz_class>, std::vector<mpz_class>, std::vector<mpz_class>>>
        performModulusGeneration(
            ZeroMQClientTransport &transport,
            const std::vector<mpz_class> xcan,
            const std::vector<mpz_class> ycan,
            const std::vector<mpz_class> zcan,
            const std::vector<mpz_class> &both_shares,
            bool special,
            ProtocolConfig<T> &config
        );

        expected<int> performJacobiTest(
            ZeroMQClientTransport &transport,
            const std::vector<mpz_class> &candidates,
            const std::vector<mpz_class> &p_shares,
            const std::vector<mpz_class> &q_shares,
            bool special
        );

        expected<int> performGCDandJacobiTest(
            ZeroMQClientTransport &transport,
            const std::vector<mpz_class> &candidates,
            const std::vector<mpz_class> &p_shares,
            const std::vector<mpz_class> &q_shares,
            const std::vector<mpz_class> &x,
            const std::vector<mpz_class> &y,
            const std::vector<mpz_class> &z,
            bool special,
            ProtocolConfig<T> &config
        );

        /* ========== getters and setters ========== */
        SocketId &socketId() {
            return _socketId;
        }
        bool &speciald() {
            return _special;
        }
        const bool &speciald() const {
            return _special;
        }
        Q &siP() {
            return _siP;
        }
        const Q &siP() const {
            return _siP;
        }
        Q &eiP() {
            return _eiP;
        }
        const Q &eiP() const {
            return _eiP;
        }
        Q &si() {
            return _si;
        }
        const Q &si() const {
            return _si;
        }
        Q &ei() {
            return _ei;
        }
        const Q &ei() const {
            return _ei;
        }
        Q &bi() {
            return _bi;
        }
        const Q &bi() const {
            return _bi;
        }
        Q &A() {
            return _A;
        }
        const Q &A() const {
            return _A;
        }
        Q &ai() {
            return _ai;
        }
        const Q &ai() const {
            return _ai;
        }
        mpz_class &jacobiSeedShares() {
            return _jacobiSeedShares;
        }
        const mpz_class &jacobiSeedShares() const {
            return _jacobiSeedShares;
        }
        mpz_class &jacobiAndGCDSeedShares() {
            return _jacobiAndGCDSeedShares;
        }
        const mpz_class &jacobiAndGCDSeedShares() const {
            return _jacobiAndGCDSeedShares;
        }

        Q &b() {
            return _b;
        }
        const Q &b() const {
            return _b;
        }
				std::vector<mpz_class> &ss() {
            return _ss;
        }
        const std::vector<mpz_class> &ss() const {
            return _ss;
        }


				std::vector<mpz_class> &ssGCD() {
            return _ssGCD;
        }
        const std::vector<mpz_class> &ssGCD() const {
            return _ssGCD;
        }

        mpz_class &gammaSeed_g() {
            return _gammaSeed;
        }
        const mpz_class &gammaSeed_g() const {
            return _gammaSeed;
        }

        std::vector<mpz_class> &moduli() {
            return _moduli;
        }
        const std::vector<mpz_class> &moduli() const {
            return _moduli;
        }

        std::vector<mpz_class> &sigma_r_GCD() {
            return _sigma_r_GCD;
        }
        const std::vector<mpz_class> &sigma_r_GCD() const {
            return _sigma_r_GCD;
        }

        std::vector<mpz_class> &sigma_x_GCD() {
            return _sigma_x_GCD;
        }
        const std::vector<mpz_class> &sigma_x_GCD() const {
            return _sigma_x_GCD;
        }

        std::vector<mpz_class> &exp_q_GCD() {
            return _exp_q_GCD;
        }
        const std::vector<mpz_class> &exp_q_GCD() const {
            return _exp_q_GCD;
        }

        std::vector<mpz_class> &sigma_e_GCD() {
            return _sigma_e_GCD;
        }
        const std::vector<mpz_class> &sigma_e_GCD() const {
            return _sigma_e_GCD;
        }

        std::vector<mpz_class> &sigma_a_GCD() {
            return _sigma_a_GCD;
        }
        const std::vector<mpz_class> &sigma_a_GCD() const {
            return _sigma_a_GCD;
        }

        std::vector<mpz_class> &sigma_g_GCD() {
            return _sigma_g_GCD;
        }
        const std::vector<mpz_class> &sigma_g_GCD() const {
            return _sigma_g_GCD;
        }

        std::vector<mpz_class> &sigma_z_GCD() {
            return _sigma_z_GCD;
        }
        const std::vector<mpz_class> &sigma_z_GCD() const {
            return _sigma_z_GCD;
        }

        std::vector<mpz_class> &sigma_q_GCD() {
            return _sigma_q_GCD;
        }
        const std::vector<mpz_class> &sigma_q_GCD() const {
            return _sigma_q_GCD;
        }

        std::vector<mpz_class> &sigmas_r() {
            return _sigmas_r;
        }
        const std::vector<mpz_class> &sigmas_r() const {
            return _sigmas_r;
        }

        std::vector<std::vector<mpz_class>> &coefs() {
            return _coefs;
        }
        const std::vector<std::vector<mpz_class>> &coefs() const {
            return _coefs;
        }

        std::vector<size_t> &index_candidates() {
            return _index_candidates;
        }
        const std::vector<size_t> &index_candidates() const {
            return _index_candidates;
        }

        uint32_t &final_index() {
            return _finalindex;
        }
        const uint32_t &final_index() const {
            return _finalindex;
        }

        const mpz_class &finalModuli() const {
            return _finalModuli;
        }
        
        mpz_class &finalModuli() {
            return _finalModuli;
        }

        const mpz_class &gammaExponent() const {
            return _gammaExponent;
        }
        mpz_class &gammaExponent() {
            return _gammaExponent;
        }

        std::vector<mpz_class> &candidatesCAN() {
            return _candidatesCAN;
        }
        const std::vector<mpz_class> &candidatesCAN() const {
            return _candidatesCAN;
        }

        std::vector<mpz_class> &x_shares() {
            return _x_shares;
        }
        const std::vector<mpz_class> &x_shares() const {
            return _x_shares;
        }
        std::vector<mpz_class> &y_shares() {
            return _y_shares;
        }
        const std::vector<mpz_class> &y_shares() const {
            return _y_shares;
        }
        std::vector<mpz_class> &z_shares() {
            return _z_shares;
        }
        const std::vector<mpz_class> &z_shares() const {
            return _z_shares;
        }
        std::vector<mpz_class> &a_shares() {
            return _a_shares;
        }
        const std::vector<mpz_class> &a_shares() const {
            return _a_shares;
        }
        std::vector<mpz_class> &b_shares() {
            return _b_shares;
        }
        const std::vector<mpz_class> &b_shares() const {
            return _b_shares;
        }

        std::vector<mpz_class> &p_sharesCRTs() {
            return _p_sharesCRTs;
        }
        const std::vector<mpz_class> &p_sharesCRTs() const {
            return _p_sharesCRTs;
        }
        std::vector<mpz_class> &q_sharesCRTs() {
            return _q_sharesCRTs;
        }
        const std::vector<mpz_class> &q_sharesCRTs() const {
            return _q_sharesCRTs;
        }

        std::vector<mpz_class> &ax() {
            return _ax;
        }
        const std::vector<mpz_class> &ax() const {
            return _ax;
        }

        std::vector<mpz_class> &by() {
            return _by;
        }
        const std::vector<mpz_class> &by() const {
            return _by;
        }

        std::vector<mpz_class> &axby() {
            return _axby;
        }
        const std::vector<mpz_class> &axby() const {
            return _axby;
        }

        std::vector<mpz_class> &rCRTs() {
            return _rCRTs;
        }
        const std::vector<mpz_class> &rCRTs() const {
            return _rCRTs;
        }

        mpz_class &gcdRX() {
            return _gcdRX;
        }
        const mpz_class &gcdRX() const {
            return _gcdRX;
        }

        std::vector<mpz_class> &axGCD() {
            return _axGCD;
        }
        const std::vector<mpz_class> &axGCD() const {
            return _axGCD;
        }

        std::vector<mpz_class> &byGCD() {
            return _byGCD;
        }
        const std::vector<mpz_class> &byGCD() const {
            return _byGCD;
        }

        std::vector<mpz_class> &axbyGCD() {
            return _axbyGCD;
        }
        const std::vector<mpz_class> &axbyGCD() const {
            return _axbyGCD;
        }

        std::vector<mpz_class> &gBeforeExp() {
            return _gBeforeExp;
        }
        const std::vector<mpz_class> &gBeforeExp() const {
            return _gBeforeExp;
        }

        std::vector<mpz_class> &gAfterExp() {
            return _gAfterExp;
        }
        const std::vector<mpz_class> &gAfterExp() const {
            return _gAfterExp;
        }

        Q &ux() {
            return _ux;
        }
        const Q &ux() const {
            return _ux;
        }
        Q &vx() {
            return _vx;
        }
        const Q &vx() const {
            return _vx;
        }
        Q &wx() {
            return _wx;
        }
        const Q &wx() const {
            return _wx;
        }

        lattice::cipher<T, Degree, NbPrimesQ> &enc_z_shares() {
            return _enc_z_shares;
        }
        const lattice::cipher<T, Degree, NbPrimesQ> &enc_z_shares() const {
            return _enc_z_shares;
        }

        Q &g_r() {
            return _g_r;
        }
        const Q &g_r() const {
            return _g_r;
        }

        P &p_y_shares() {
            return _p_y_shares;
        }
        const P &p_y_shares() const {
            return _p_y_shares;
        }

        Q &uz() {
            return _uz;
        }
        const Q &uz() const {
            return _uz;
        }
        Q &vz() {
            return _vz;
        }
        const Q &vz() const {
            return _vz;
        }
        Q &wz() {
            return _wz;
        }
        const Q &wz() const {
            return _wz;
        }

        Q &partial_xyz_shares() {
            return _partial_xyz_shares;
        }
        const Q &partial_xyz_shares() const {
            return _partial_xyz_shares;
        }

        CandidateIndices &candidateIndices() {
            return _candidateIndices;
        }
        const CandidateIndices &candidateIndices() const {
            return _candidateIndices;
        }

        lattice::cipher<T, Degree, NbPrimesQ> &xyz_sum() {
            return _xyz_sum;
        }
        const lattice::cipher<T, Degree, NbPrimesQ> &xyz_sum() const {
            return _xyz_sum;
        }

        nfl::poly_p<T, Degree, NbPrimesQ> &zp() {
            return _zp;
        }
        const nfl::poly_p<T, Degree, NbPrimesQ> &zp() const {
            return _zp;
        }

        lattice::cipher<T, Degree, NbPrimesQ> &xsum_first() {
            return _xsum_first;
        }
        const lattice::cipher<T, Degree, NbPrimesQ> &xsum_first() const {
            return _xsum_first;
        }
        lattice::cipher<T, Degree, NbPrimesQ> &xsum_final() {
            return _xsum_final;
        }
        const lattice::cipher<T, Degree, NbPrimesQ> &xsum_final() const {
            return _xsum_final;
        }

        Q &up() {
            return _up;
        }
        const Q &up() const {
            return _up;
        }
        Q &vp() {
            return _vp;
        }
        const Q &vp() const {
            return _vp;
        }
        Q &wp() {
            return _wp;
        }
        const Q &wp() const {
            return _wp;
        }

        lattice::cipher<T, Degree, NbPrimesQ> &g_enc_x_shares() {
            return _enc_x_shares;
        }
        const lattice::cipher<T, Degree, NbPrimesQ> &g_enc_x_shares() const {
            return _enc_x_shares;
        }

        cipherText &as() {
            return _as;
        }
        const cipherText &as() const {
            return _as;
        }
        cipherText &es() {
            return _es;
        }
        const cipherText &es() const {
            return _es;
        }
        std::vector<cipherText> &as_vec() {
            return _as_vec;
        }
        const std::vector<cipherText> &as_vec() const {
            return _as_vec;
        }
        std::vector<cipherText> &es_vec() {
            return _es_vec;
        }
        const std::vector<cipherText> &es_vec() const {
            return _es_vec;
        }
        Q &partial_e_shares() {
            return _partial_e_shares;
        }
        const Q &partial_e_shares() const {
            return _partial_e_shares;
        }
        std::vector<int> &sievingFlags() {
            return _sievingFlags;
        }
        const std::vector<int> &sievingFlags() const {
            return _sievingFlags;
        }
        ProtocolConfig<T> &config() {
            return _config;
        }
        const ProtocolConfig<T> &config() const {
            return _config;
        }

        /** Clear proofs data between protocol restarts */
        void clearProofs() {
            {
                Q tmp = Q{};
                si() = tmp;
                ei() = tmp;
                A() = tmp;
                b() = tmp;
                partial_e_shares() = tmp;
                as().first = tmp;
                as().second = tmp;
                es().first = tmp;
                es().second = tmp;
            }

            a_shares().clear();
            b_shares().clear();
            as_vec().clear();
            es_vec().clear();
            sievingFlags().clear();
        }

        /** Protocol start helper
         * @param transport communication helper
         * @param needRegistration is true when registration needed
         * @param localUri coordinator's IP address
         * @param chi protocol parameter
         * @param prover zero-knowledge prover
         */
        expected<Unit> startHelper(ZeroMQClientTransport& transport, 
            bool needRegistration,
            const std::string& localUri,
            std::unique_ptr<nfl::gaussian<uint16_t, T, 2>> &chi,
            zksnark::Prover<FieldT>* prover
        ) {
            if (needRegistration) {
                auto conf = registerAndAwaitConfiguration(transport, localUri, *chi, prover);
                if (hasError(conf)) {
                    LOG(INFO) << "Failed in registration, the error message is: " << showError(getError(conf));
                    exit(EXIT_FAILURE);
                }
                config() = getResult(conf);
                // set coordinators public key
                transport.coordinatorPublicKey = config().publicKey();
            }

            auto enc = lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>(config());

            auto result = start_rsa_ceremony(transport, config(), enc);

            if (!hasError(result)) {
                return Unit{};
            }

            auto err = getError(result);
            if (err == Error::MODULUS_NOT_FOUND) {
                // restart without registration
                clearProofs();
                return startHelper(transport, false, localUri, chi, prover);
            }
            else if (err == Error::RESTART) {
                clearProofs();
                return startHelper(transport, true, localUri, chi, prover);
            }
            return err;
        }

        expected<Unit> start(
            ZeroMQClientTransport& transport,
            int sigma,
            int lambda,
            const std::string &localUri = "n/a",
            zksnark::Prover<FieldT>* prover = static_cast<zksnark::Prover<FieldT>*>(NULL)
        ) {
            _fg = std::make_unique<nfl::FastGaussianNoise<uint16_t, T, 2>>(sigma, lambda, Degree);
            _chi = std::make_unique<nfl::gaussian<uint16_t, T, 2>>(_fg.get());
            _e = std::make_unique<lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>>(config(), *_chi);

            return startHelper(transport, true, localUri, _chi, prover);
        }

        expected<int> start_rsa_ceremony(
            ZeroMQClientTransport &transport,
            ProtocolConfig<T> &config,
            lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ> &e
        ) {

            transport.myTimers.begin(2,"1.b. Overall speed", 1);
            auto maybe_keys = generateKeyPair(transport, e.chi());
            if (hasError(maybe_keys)) {
                LOG(ERROR) << "Kill/Restart received when keygen";
                return getError(maybe_keys);
            }
            auto [publicKey, secretKey] = getResult(maybe_keys);

            // Await designation
            DBG("Awaiting designation.");

            auto maybe_assignment = transport.awaitReply();
            if (hasError(maybe_assignment)) {
                LOG(ERROR) << "Kill/Restart when wait assignment";
                return getError(maybe_assignment);
            }
            MessageType assignment = getResult(maybe_assignment);
            bool special = assignment == MessageType::ASSIGNMENT_P1;
            speciald() = special;
            DBG("Participant (" << socketId() << ") got my designation.");

            LOG(INFO) << "Connected to coordinator.";

            bool foundModuli = false;

            LOG(INFO) << "Generating shares for pre-sieving.";
            auto maybe_shares = generatePreSievingShares(e, transport, publicKey, secretKey,
                                special, config);
            if (hasError(maybe_shares)) {
                LOG(ERROR) << "Kill/Restart received when generate presieve shares";
                return getError(maybe_shares);
            }
            auto [xcan, ycan, zcan, xgcd, ygcd, zgcd, both_shares] = getResult(maybe_shares);

            LOG(INFO) << "Generated shares for " << both_shares.size() / 2 <<
                        " candidates.";

            LOG(INFO) << "Using pre-sieved shares for candidate generation.";
            auto maybe_modulus = performModulusGeneration(transport, xcan, ycan, zcan,
                                    both_shares, special, config);
            if (hasError(maybe_modulus)) {
                LOG(ERROR) << "Kill/Restart received when perform modulus generation";
                return getError(maybe_modulus);
            }
            auto [candidates, p_shares, q_shares] = getResult(maybe_modulus);

            {
                for (int i = 0; i < p_shares.size(); ++i) {
                    if (special) {
                        if (p_shares[i] % 4 != 3 || q_shares[i] % 4 != 3) {
                            DBG("special");
                            DBG("i = " << i);
                            DBG("p_shares[i] " << p_shares[i]);
                            DBG("p_shares[i] " << q_shares[i]);
                            assert(false);
                        }
                    } else {
                        if (p_shares[i] % 4 != 0 || q_shares[i] % 4 != 0) {
                            DBG("ordinary");
                            DBG("i = " << i);
                            DBG("p_shares[i] " << p_shares[i]);
                            DBG("p_shares[i] " << q_shares[i]);
                            assert(false);
                        }
                    }
                }
            }

            auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>
                                        (MessageType::POST_SIEVE);
            if (hasError(maybe_discardFlags)) {
                LOG(ERROR) << "Kill/Restart when wait discard flags";
                return getError(maybe_discardFlags);
            }
            boost::dynamic_bitset<> discardFlags = getResult(maybe_discardFlags);
            candidates = discardCandidates(candidates, discardFlags);
            p_shares = discardCandidates(p_shares, discardFlags);
            q_shares = discardCandidates(q_shares, discardFlags);

            // NOTE: here we need to adjust
            //indicesArray = discardCandidates(indicesArray, discardFlags);

            LOG(INFO) << "Running distributed primality test on the surviving " <<
                        candidates.size() << " candidates.";

            //Finally, we move onto the biprimality test, which we implement here
            //in a loop because the coordinator may decide to run each test a variable
            //number of times.
            bool stillRunningTests = true;
            bool ranJacobi = false;

            while (stillRunningTests) {
                auto success = transport.awaitReply();
                if (hasError(success)) {
                    LOG(ERROR) << "Kill/Restart received";
                    return getError(success);
                }
                switch (getResult(success)) {
                    case MessageType::GAMMA_SHARES: {
                        if (!ranJacobi) {
                            LOG(INFO) << "Running Jacobi test on candidates.";
                            ranJacobi = true;
                        }
                        auto maybe_JacobiResult = performJacobiTest(transport, candidates, p_shares, q_shares, special);
                        if (hasError(maybe_JacobiResult)) {
                            return getError(maybe_JacobiResult);
                        }

                        auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>(MessageType::DISCARD_FLAGS);
                        if (hasError(maybe_discardFlags)) {
                            LOG(ERROR) << "Kill/Restart received during wait discard flags";
                            return getError(maybe_discardFlags);
                        }
                        boost::dynamic_bitset<> discardFlags = getResult(maybe_discardFlags);

                        candidates = discardCandidates(candidates, discardFlags);
                        p_shares = discardCandidates(p_shares, discardFlags);
                        q_shares = discardCandidates(q_shares, discardFlags);

                        break;
                    }
                    case MessageType::GCD_RAND_SHARES: {
                        LOG(INFO) << "Running GCD test on candidates.";

                        if (candidates.size() > 1) {
                            candidates.resize(1);
                            p_shares.resize(1);
                            q_shares.resize(1);
                        }

                        auto maybe_JacobiGCDResult = performGCDandJacobiTest(
                                transport,
                                candidates,
                                p_shares,
                                q_shares,
                                xgcd,
                                ygcd,
                                zgcd,
                                special,
                                config
                                );

                        if (hasError(maybe_JacobiGCDResult)) {
                            return getError(maybe_JacobiGCDResult);
                        }

                        auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>(MessageType::DISCARD_FLAGS);
                        if (hasError(maybe_discardFlags)) {
                            LOG(ERROR) << "Kill/Restart received during wait discard flags";
                            return getError(maybe_discardFlags);
                        }
                        boost::dynamic_bitset<> discardFlags = getResult(maybe_discardFlags);

                        candidates = discardCandidates(candidates, discardFlags);
                        p_shares = discardCandidates(p_shares, discardFlags);
                        q_shares = discardCandidates(q_shares, discardFlags);

                        break;
                    }
                    case MessageType::FOUND_MODULI:
                        foundModuli = true;
                        stillRunningTests = false;
                        break;
                    case MessageType::NO_MODULI: {
                        LOG(INFO) << "Did not find a modulus. Restarting";
                        stillRunningTests = false;
                        transport.send(MessageType::NO_MODULI_VERIFICATION_SHARES, std::make_pair(x_shares(), std::make_pair(y_shares(), z_shares())));
                        auto maybeNoModuliVerification = transport.awaitReply<bool>(MessageType::NO_MODULI_VERIFICATION_FINISHED);
                        if (hasError(maybe_discardFlags)) {
                            LOG(ERROR) << "Kill/Restart received on no moduli verification";
                            return getError(maybe_discardFlags);
                        }

                        bool noModuliVerification = getResult(maybeNoModuliVerification);
                        if (!noModuliVerification) {
                            transport.send(MessageType::BAD_EVENT_DATA_SHARES, BadEventData<T, Degree, NbPrimesQ>(jacobiSeedShares(), jacobiAndGCDSeedShares(), si(), ei(), bi(), gcdRX(), ss()));

                            transport.awaitReply<bool>(MessageType::BAD_EVENT_DATA_RESPONSE);
                        }
                        return Error::MODULUS_NOT_FOUND;

                    }
                    default:
                        return Error::OUT_OF_SYNC;
                }
            }

            if (foundModuli) {
                LOG(INFO) << "Found " << candidates.size() << " valid moduli:";

                for (int i = 0; i < candidates.size(); ++i) {
                    LOG(INFO) << candidates[i];
                }
                size_t idx = 0;
                bool rec = false;
                for (; idx < candidatesCAN().size(); idx++) {
                    if (candidatesCAN()[idx] == candidates[0]) {
                        rec = true;
                        break;
                    }
                }

                if (!rec) {
                    // handle error
                }

                candidateIndices() = getCandidateIndices(idx, index_candidates());
                final_index() = idx;
                finalModuli() = candidates[0];
                gammaExponent() = ((speciald() ? mpz_class(candidates[0]+1) : mpz_class(0)) - p_shares[0] - q_shares[0])/4;

                if (config.protocolMode() == ProtocolMode::RECORD || config.protocolMode() == ProtocolMode::REPLAY) {
                    // The following code transmits the p shares and q shares in the clear to the coordinator.
                    // This is ONLY for checking the result. We need to remove this part finally.
                    // DBG("Participant (" << socketId() <<
                    //     ") Transmitting the p and q values in the clear for debugging purposes only.");
                    transport.send(MessageType::P_CLEAR_DEBUG, std::pair{p_shares, q_shares});
                    transport.awaitReply<int>(MessageType::P_CLEAR_DEBUG);
                }

            }

            return 0;
        }

        /** Gathers and updates public and secrete data passed as params
         * @param d public data to update
         * @param s sigma protocol public data to update
         * @param e secret data to update
         */
        void gatherData(PublicData &d, SigmaProtocolPublicData &s, SecretData &e) {

            // Shift a and bi to poly representation
            A().invntt_pow_invphi();
            bi().invntt_pow_invphi();

            // Compute all q polys
            // ========================================================================================
            using Qwide = nfl::poly_p<T, 2 * Degree, NbPrimesQ>;

            // transform Q/P to FFT polynomial
            mpz_class q_over_p = 1;
            for (auto i = NbPrimesP; i < NbPrimesQ; i++) {
                q_over_p *= nfl::params<T>::P[i];
            }

            Q qdivp{q_over_p};
            qdivp.ntt_pow_phi();

            // Promoting to a poly of degree 2*d for poly multiplications
            auto promote = [&](Q & source, Qwide & target) -> void {
                for (size_t primeIdx = 0; primeIdx < NbPrimesQ; primeIdx++) {
                    for (size_t idx = 0; idx < Degree; idx++) {
                        target.poly_obj()(primeIdx, idx) = source.poly_obj()(primeIdx, idx);
                    }

                    for (size_t idx = Degree; idx < 2 * Degree; idx++) {
                        target.poly_obj()(primeIdx, idx) = 0;
                    }
                }
            };
            auto checkZero = [&](Q & target) -> void {
                for (size_t primeIdx = 0; primeIdx < NbPrimesQ; primeIdx++) {
                    for (size_t coefIdx = 0; coefIdx < Degree; coefIdx++) {
                        if (target(primeIdx, coefIdx) != 0) {
                            std::cout << "primeIdx = " << primeIdx << " and coefIdx = " << coefIdx <<
                                      std::endl;
                            std::cout << "target(" << primeIdx << "," << coefIdx << ") = " << target(
                                          primeIdx, coefIdx) << std::endl;
                            std::cout << "Prime = " << Q::get_modulus(primeIdx) << std::endl;
                        }
                        assert(target(primeIdx, coefIdx) == 0);
                    }
                }
            };

            auto check = [&](Qwide & target) -> void {
                for (size_t primeIdx = 0; primeIdx < NbPrimesQ; primeIdx++) {
                    for (size_t coefIdx = 0; coefIdx < Degree; coefIdx++) {
                        assert(target(primeIdx, coefIdx) == target(primeIdx, coefIdx + Degree));
                    }
                }
            };

            auto deriveQ = [&](Q & q, Qwide & left) -> void {
                for (size_t primeIdx = 0; primeIdx < NbPrimesQ; primeIdx++) {
                    for (size_t coef = 0; coef < Degree; coef++) {
                        q(primeIdx, coef) = left(primeIdx, coef);
                    }
                }
            };

            auto convert_P = [&](P & source, Q & target) -> void {
                source.invntt_pow_invphi();
                auto s_coeff = source.poly2mpz();
                target.set_mpz(s_coeff);
                std::for_each(s_coeff.begin(), s_coeff.end(), mpz_clear);
            };

            auto convert_mpz = [&](std::vector<mpz_class> &source, Q & target) -> void {
                // preparing the message polynomial
                auto *ptr = new std::array<mpz_class, Degree>{};
                std::copy(source.begin(), source.end(), ptr->begin());
                target.set_mpz(*ptr);
                target = target * qdivp;
                delete (ptr);
            };

            auto convert_noscale_mpz_P = [&](std::vector<mpz_class> &source,
            P & target) -> void {
                // preparing the message polynomial
                auto *ptr = new std::array<mpz_class, Degree>{};
                std::copy(source.begin(), source.end(), ptr->begin());
                target.set_mpz(*ptr);
                delete (ptr);
            };

            auto convert_noscale_mpz = [&](std::vector<mpz_class> &source,
            Q & target) -> void {
                // preparing the message polynomial
                auto *ptr = new std::array<mpz_class, Degree>{};
                std::copy(source.begin(), source.end(), ptr->begin());
                target.set_mpz(*ptr);
                delete (ptr);
            };

            auto modifiedProduct = [&](std::vector<mpz_class> &a,
                                       std::vector<size_t> &range1, std::vector<mpz_class> &b,
            std::vector<size_t> &range2) -> mpz_class {
                mpz_class tmp = 0;
                size_t idx1 = 0;

                for (size_t idx2 : range2) {
                    tmp = tmp + a[range1[idx1++]] * b[idx2];
                }

                tmp += (speciald() ? a.back() * 3 : mpz_class(0));
                return tmp;
            };

            auto positiveRemainder = [&](mpz_class & a, mpz_class & b) -> mpz_class {
                return ((a % b + b) % b);
            };

            Q q_r3, q_r4_1, q_r4_2, q_r5_1, q_r5_2, q_r6;

            std::vector<mpz_class> q_p_r7, q_q_r7, q_r8;
            std::vector<mpz_class> q_p_prod_r7, q_q_prod_r7;
            std::vector<mpz_class> q_p_prod_r11, q_q_prod_r11;
            std::vector<mpz_class> q_pq_r11, q_r_r11,r_CRTs;
						std::vector<mpz_class> q_r12;

            // Round 3
            {
                Qwide A_wide, siP_wide, eiP_wide, bi_wide;

                // Must promote the ring elements to degree 2*n prevent overflow
                promote(A(), A_wide);
                promote(siP(), siP_wide);
                promote(eiP(), eiP_wide);
                promote(bi(), bi_wide);

                // Compute left handside of eq
                // A*si + ei - bi = (x^n+1) * q
                A_wide.ntt_pow_phi();
                siP_wide.ntt_pow_phi();
                Qwide prod = A_wide * siP_wide;
                prod.invntt_pow_invphi();
                Qwide left = prod + eiP_wide - bi_wide;

                // Derive q
                check(left);                // Make sure that poly left is indeed a ring element
                deriveQ(q_r3, left);        // Now solve for q coefs
            }

            // Round 4
            // Shift a and bi to poly representation
            // check equation on first coef
            b().invntt_pow_invphi();
            ux().invntt_pow_invphi();
            vx().invntt_pow_invphi();
            wx().invntt_pow_invphi();
            g_enc_x_shares().first.invntt_pow_invphi();
            g_enc_x_shares().second.invntt_pow_invphi();

            Q x_prime;
            convert_mpz(x_shares(), x_prime);

            Q copy = x_prime;
            x_prime.invntt_pow_invphi();

            {
                Qwide ux_wide, A_wide, vx_wide, ci_1_wide;

                // Must promote the ring elements to degree 2*n prevent overflow
                promote(ux(), ux_wide);
                promote(A(), A_wide);
                promote(vx(), vx_wide);
                promote(g_enc_x_shares().first, ci_1_wide);

                // Compute left handside of eq
                // u * a + v - ci_1 = (x^n+1) * q
                A_wide.ntt_pow_phi();
                ux_wide.ntt_pow_phi();
                Qwide prod = A_wide * ux_wide;
                prod.invntt_pow_invphi();

                Qwide left = prod + vx_wide - ci_1_wide;

                // Derive q
                check(left);                    // Make sure that poly left is indeed a ring element
                deriveQ(q_r4_1, left);          // Now solve for q coefs
            }

            {
                Qwide ux_wide, b_wide, wx_wide, x_prime_wide, ci_2_wide;

                // Must promote the ring elements to degree 2*n prevent overflow
                promote(ux(), ux_wide);
                promote(b(), b_wide);
                promote(wx(), wx_wide);
                promote(x_prime, x_prime_wide);
                promote(g_enc_x_shares().second, ci_2_wide);

                // Compute left handside of eq
                // u * b + w + (Q/P) m' - ci_2 = (x^n+1) * q
                ux_wide.ntt_pow_phi();
                b_wide.ntt_pow_phi();
                Qwide prod = ux_wide * b_wide;
                prod.invntt_pow_invphi();

                Qwide left = prod + wx_wide + x_prime_wide - ci_2_wide;

                // Derive q
                check(left);                    // Make sure that poly left is indeed a ring element
                deriveQ(q_r4_2, left);          // Now solve for q coefs
            }

            // Round 5
            Q y_prime;
            P p_y_prime;
            convert_noscale_mpz_P(y_shares(), p_y_prime);

            p_y_prime.invntt_pow_invphi();
            auto s_coeff = p_y_prime.poly2mpz();
            y_prime.set_mpz(s_coeff);
            std::for_each(s_coeff.begin(), s_coeff.end(), mpz_clear);

            Q z_prime;
            convert_mpz(z_shares(), z_prime);
            z_prime.invntt_pow_invphi();

            enc_z_shares().first.invntt_pow_invphi();
            xsum_final().first.invntt_pow_invphi();         // ci'_1
            xsum_final().second.invntt_pow_invphi();
            xsum_first().first.invntt_pow_invphi();         // c_1
            xsum_first().second.invntt_pow_invphi();
            uz().invntt_pow_invphi();
            vz().invntt_pow_invphi();
            wz().invntt_pow_invphi();
            zp().invntt_pow_invphi();

            {
                // First Constraint: ci_1_prime = yi*c_1 - u*a - v
                Qwide uz_wide, vz_wide, A_wide, ci_1_prime_wide, c_1_wide, yi_wide;
                Qwide enc_z_wide;

                // Must promote the ring elements to degree 2*n prevent overflow
                promote(xsum_final().first, ci_1_prime_wide);
                promote(xsum_first().first, c_1_wide);
                promote(uz(), uz_wide);
                promote(vz(), vz_wide);
                promote(A(), A_wide);

                promote(y_prime, yi_wide);
                promote(enc_z_shares().first, enc_z_wide);

                // Compute left handside of eq
                yi_wide.ntt_pow_phi();
                c_1_wide.ntt_pow_phi();
                Qwide prod1 = yi_wide * c_1_wide; // yi_wide* ;
                prod1.invntt_pow_invphi();

                uz_wide.ntt_pow_phi();
                A_wide.ntt_pow_phi();
                Qwide prod2 = uz_wide * A_wide;
                prod2.invntt_pow_invphi();

                Qwide left = prod1 - prod2 - vz_wide - ci_1_prime_wide;

                // Derive q
                check(left);                    // Make sure that poly left is indeed a ring element
                deriveQ(q_r5_1, left);          // Now solve for q coefs
            }

            {
                // Second Constraint: ci_2_prime = yi*c_2 + zp - uz*b - wz - (Q/P) m'
                Qwide ci_2_prime_wide, c_2_wide, yi_wide, zp_wide, uz_wide, b_wide, wz_wide,
                      z_prime_wide;

                // Must promote the ring elements to degree 2*n prevent overflow
                promote(xsum_final().second, ci_2_prime_wide);
                promote(xsum_first().second, c_2_wide);
                promote(y_prime, yi_wide);
                promote(zp(), zp_wide);
                promote(uz(), uz_wide);
                promote(b(), b_wide);
                promote(wz(), wz_wide);
                promote(z_prime, z_prime_wide);

                // Compute left handside of eq
                yi_wide.ntt_pow_phi();
                c_2_wide.ntt_pow_phi();
                Qwide prod1 = yi_wide * c_2_wide;
                prod1.invntt_pow_invphi();

                uz_wide.ntt_pow_phi();
                b_wide.ntt_pow_phi();
                Qwide prod2 = uz_wide * b_wide;
                prod2.invntt_pow_invphi();

                Qwide left = prod1 + zp_wide - prod2 - wz_wide - z_prime_wide - ci_2_prime_wide;

                // Derive q
                check(left);                    // Make sure that poly left is indeed a ring element
                deriveQ(q_r5_2, left);          // Now solve for q coefs
            }

            g_r().invntt_pow_invphi();
            partial_xyz_shares().invntt_pow_invphi();
            xyz_sum().first.invntt_pow_invphi();
            xyz_sum().second.invntt_pow_invphi();

            // Round 6
            {
                // di = c'_2 - c'_1*s + r;
                Qwide di_wide, c_prime_2_wide, c_prime_1_wide, siP_wide, r_wide;

                // Must promote the ring elements to degree 2*n prevent overflow
                promote(g_r(), r_wide);
                promote(partial_xyz_shares(), di_wide);
                promote(xyz_sum().first, c_prime_1_wide);
                promote(xyz_sum().second, c_prime_2_wide);
                promote(siP(), siP_wide);

                // Compute left handside of eq
                c_prime_1_wide.ntt_pow_phi();
                siP_wide.ntt_pow_phi();
                Qwide prod = c_prime_1_wide * siP_wide;
                prod.invntt_pow_invphi();
                Qwide left = r_wide - di_wide - prod;

                if (speciald()) {
                    left = left + c_prime_2_wide;
                }

                // Derive q
                check(left);                // Make sure that poly left is indeed a ring element
                deriveQ(q_r6, left);        // Now solve for q coefs
            }

            std::vector<mpz_class> ax_shares;
            std::vector<mpz_class> by_shares;
            std::vector<size_t> coefIndices = {0, 1, 2, 3, 4, 5};
            std::vector<size_t> axbyindices;

            // These are common for Rounds 7,8,11,12
            auto [alphasCAN, _bsz] = math::fixed_bucket_n_primes(1048, primesCAN, 175);
            std::vector<mpz_class> alphasPSprod_CAN, moduli_CAN;
            auto [alphasPS, bucketSizePS] = math::balanced_bucket_n_primes(1000, primesPS, 175, 1);
            mpz_class alphasPSprod = mpz_class(4);
            
            for (size_t j = 0; j < alphasPS.size(); ++j) {
                alphasPSprod *= alphasPS[j];
            }

            mpz_class recomputed_pshare;
            mpz_class recomputed_qshare;

            // Round 7 and 8
            {
                size_t alphaCANidx = 0;

                for (auto j : candidateIndices().can) {
                    std::vector<mpz_class> moduliCAN(moduli().size());
                    mpz_class alphasPSprodCAN = alphasPSprod % alphasCAN[alphaCANidx];
                    for (size_t i = 0; i < moduli().size(); i ++) {
                        moduliCAN[i] = moduli()[i] % alphasCAN[alphaCANidx];
                        moduli_CAN.push_back(moduliCAN[i]);
                    }

                    alphasPSprod_CAN.push_back(alphasPSprodCAN);

                    // pushing pshare items
                    recomputed_pshare = modifiedProduct(moduli(), coefIndices, x_shares(), candidateIndices().ps);
                    mpz_class pshare_quotient = (recomputed_pshare / alphasPSprod) % alphasCAN[alphaCANidx];

                    mpz_class recomputed_pshare_CAN = modifiedProduct(moduliCAN, coefIndices, x_shares(), candidateIndices().ps);
                    mpz_class axshares = math::mod(p_sharesCRTs()[final_index() * alphasCAN.size() + alphaCANidx] - x_shares()[j], alphasCAN[alphaCANidx]);
                    mpz_class pdividend = axshares + x_shares()[j] - (recomputed_pshare_CAN - pshare_quotient * alphasPSprodCAN);
                    mpz_class p_shareCRT_quotient  = (pdividend / alphasCAN[alphaCANidx]);
                    assert(pdividend - p_shareCRT_quotient * alphasCAN[alphaCANidx] == 0);

                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        q_p_prod_r7.push_back(positiveRemainder(pshare_quotient, prime));
                        q_p_r7.push_back(positiveRemainder(p_shareCRT_quotient, prime));
                        ax_shares.push_back(positiveRemainder(axshares, prime));
                    }

                    //pushing qshare items
                    recomputed_qshare = modifiedProduct(moduli(), coefIndices, y_shares(), candidateIndices().ps);
                    mpz_class qshare_quotient = (recomputed_qshare / alphasPSprod) % alphasCAN[alphaCANidx];

                    mpz_class recomputed_qshare_CAN = modifiedProduct(moduliCAN, coefIndices, y_shares(), candidateIndices().ps);
                    mpz_class byshares = math::mod(q_sharesCRTs()[final_index() * alphasCAN.size() + alphaCANidx] - y_shares()[j], alphasCAN[alphaCANidx]);
                    mpz_class qdividend = byshares + y_shares()[j] - (recomputed_qshare_CAN - qshare_quotient * alphasPSprodCAN);
                    mpz_class q_shareCRT_quotient  = (qdividend / alphasCAN[alphaCANidx]);
                    assert(qdividend - q_shareCRT_quotient * alphasCAN[alphaCANidx] == 0);

                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        q_q_prod_r7.push_back(positiveRemainder(qshare_quotient, prime));
                        q_q_r7.push_back(positiveRemainder(q_shareCRT_quotient, prime));
                        by_shares.push_back(positiveRemainder(byshares, prime)); //mpz_class(nfl::params<uint64_t>::P[i])));
                    }
                    size_t pos = alphasCAN.size() * final_index() + alphaCANidx;
                    axbyindices.push_back(pos);
                    mpz_class dividend = (ax()[pos] * (recomputed_qshare_CAN - qshare_quotient *
                                                       alphasPSprodCAN) + by()[pos] * (recomputed_pshare_CAN - pshare_quotient *
                                                               alphasPSprodCAN) + z_shares()[j] - axby()[pos]);
                    mpz_class quotient = dividend / alphasCAN[alphaCANidx];
                    assert(dividend - quotient * alphasCAN[alphaCANidx] == 0);
                    
                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        q_r8.push_back(positiveRemainder(quotient, prime));
                    }

                    alphaCANidx++;
                }
            }

            auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3210, primesGCD,175);
            std::vector<mpz_class> alphasPSprod_GCD, moduli_GCD;
            std::vector<mpz_class> ax_shares_GCD;
            std::vector<mpz_class> by_shares_GCD;
            std::vector<mpz_class> ss_GCD;
            std::vector<mpz_class> finalModuli_GCD;

            std::vector<mpz_class> expsharesGCD;

            // Round 11
            // and Round 12
            {
                size_t alphaGCDidx = 0;
                const mpz_class maybeOne = speciald() ? mpz_class(1) : mpz_class(0);
                for (auto j : candidateIndices().gcd) {
                    std::vector<mpz_class> moduliGCD(moduli().size());
                    mpz_class alphasPSprodGCD = alphasPSprod % alphasGCD[alphaGCDidx];
                    for (size_t i = 0; i < moduli().size(); i ++) {
                        moduliGCD[i] = moduli()[i] % alphasGCD[alphaGCDidx];
                        moduli_GCD.push_back(moduliGCD[i]);
                    }

                    alphasPSprod_GCD.push_back(alphasPSprodGCD);

                    mpz_class pshare_quotient = (recomputed_pshare / alphasPSprod) % alphasGCD[alphaGCDidx];
                    mpz_class qshare_quotient = (recomputed_qshare / alphasPSprod) % alphasGCD[alphaGCDidx];
                    mpz_class recomputed_pshare_GCD = modifiedProduct(moduliGCD, coefIndices, x_shares(), candidateIndices().ps);
                    mpz_class recomputed_qshare_GCD = modifiedProduct(moduliGCD, coefIndices, y_shares(), candidateIndices().ps);
                    expsharesGCD.push_back(recomputed_pshare_GCD + recomputed_qshare_GCD - qshare_quotient * alphasPSprodGCD - pshare_quotient * alphasPSprodGCD);
                    mpz_class bysharesGCD = math::mod(recomputed_pshare_GCD + recomputed_qshare_GCD - maybeOne - qshare_quotient * alphasPSprodGCD - pshare_quotient * alphasPSprodGCD - y_shares()[j] , alphasGCD[alphaGCDidx]) ;
                    mpz_class pqdividendGCD = bysharesGCD + y_shares()[j] - (recomputed_pshare_GCD + recomputed_qshare_GCD - maybeOne - pshare_quotient * alphasPSprodGCD - qshare_quotient * alphasPSprodGCD);
                    mpz_class pqshare_quotient  = (pqdividendGCD/alphasGCD[alphaGCDidx]);
                    assert(pqdividendGCD - pqshare_quotient * alphasGCD[alphaGCDidx] == 0);

                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        q_p_prod_r11.push_back(positiveRemainder(pshare_quotient, prime));
                        q_q_prod_r11.push_back(positiveRemainder(qshare_quotient, prime));
                        q_pq_r11.push_back(positiveRemainder(pqshare_quotient, prime));
                        by_shares_GCD.push_back(positiveRemainder(bysharesGCD, prime));
                    }

                    mpz_class axsharesGCD = math::mod(rCRTs()[alphaGCDidx] - x_shares()[j], alphasGCD[alphaGCDidx]);
                    mpz_class rdividendGCD = axsharesGCD + x_shares()[j] - rCRTs()[alphaGCDidx];
                    mpz_class rquotientGCD = rdividendGCD/alphasGCD[alphaGCDidx];
                    assert(rdividendGCD - rquotientGCD * alphasGCD[alphaGCDidx] == 0);

                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        r_CRTs.push_back(positiveRemainder(rCRTs()[alphaGCDidx], prime));
                        q_r_r11.push_back(positiveRemainder(rquotientGCD, prime));
                        ax_shares_GCD.push_back(positiveRemainder(axsharesGCD, prime));
                    }

                    mpz_class dividend = axGCD()[alphaGCDidx] * ( 
                        recomputed_pshare_GCD + recomputed_qshare_GCD - maybeOne - qshare_quotient * alphasPSprodGCD - pshare_quotient * alphasPSprodGCD 
                        ) + byGCD()[alphaGCDidx] * (rCRTs()[alphaGCDidx]) + z_shares()[j] + math::mod(ss()[0],alphasGCD[alphaGCDidx])*math::mod(finalModuli(),alphasGCD[alphaGCDidx]) - axbyGCD()[alphaGCDidx];

					assert(dividend % alphasGCD[alphaGCDidx] == 0);
                    mpz_class quotient = dividend / alphasGCD[alphaGCDidx];
                    assert(dividend - quotient * alphasGCD[alphaGCDidx] == 0);
                    
                    for (size_t i = 0; i < NbPrimesQ; i++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        mpz_class out = math::mod(ss()[0],alphasGCD[alphaGCDidx]);
                        ss_GCD.push_back(positiveRemainder(out,prime));
                        out = math::mod(finalModuli(),alphasGCD[alphaGCDidx]);
                        finalModuli_GCD.push_back(positiveRemainder(out,prime));
                        q_r12.push_back(positiveRemainder(quotient, prime));
                    }

                    alphaGCDidx++;
                }
            }
            // Data Extraction
            // ========================================================================================

            auto dataSub = [&](std::vector<uint64_t> &target, std::vector<uint64_t> &source, std::vector<size_t> &idxs) -> void {
                target.reserve(idxs.size());
                for (auto idx : idxs) {
                    target.emplace_back(source[idx]);
                }
            };

            auto dataExtractMatrix = [&](std::vector<std::vector<uint64_t>> &assignee,
            std::vector<std::vector<mpz_class>> &target) -> void {
                assignee.resize(target.size());
                for (size_t row = 0; row < target.size(); row++) {
                    assignee[row].resize(target[row].size());
                    for (size_t col = 0; col < target.size(); col++) {
                        assignee[row][col] = mpz_get_ui(target[row][col].get_mpz_t());
                    }
                }
            };

            auto pickSubVector = [&](std::vector<mpz_class> &assignee, std::vector<mpz_class> &target, std::vector<size_t> &idxs) -> void {
                std::vector<mpz_class> tmp(idxs.size());

                for (size_t idx = 0; idx < idxs.size(); idx++) {
                        tmp[idx] = target[idxs[idx]];
                }

                if (assignee.size()>0) assignee.insert(assignee.end(), idxs.begin(), idxs.end());
                else assignee.assign(idxs.begin(), idxs.end());
            };

	    // Extract Gamm Seed for Sigma Protocol
	    d.gammaSeed = gammaSeed_g();

            // Randomness attached to surviving candidates
            std::vector<size_t> indices(candidateIndices().can);
            indices.insert(indices.end(), candidateIndices().ps.begin(), candidateIndices().ps.end());
            indices.insert(indices.end(), candidateIndices().gcd.begin(), candidateIndices().gcd.end());

            // Alphas
            dataExtractVectorCRT<NbPrimesQ>(d.ps,      alphasPS);

            // Keygen
            dataExtractQ(d.A,                A());
            dataExtractQ(d.bi,               bi());
            dataExtractQ(e.eiP,              eiP());
            dataExtractQ(e.siP,              siP());
            dataExtractQ(e.q_r3,             q_r3);

            // Round4
            dataExtractQ(d.ci_1,             g_enc_x_shares().first);
            dataExtractQ(d.ci_2,             g_enc_x_shares().second);
            dataExtractQ(d.b,                b());
            dataExtractQ(e.ux,               ux());
            dataExtractQ(e.vx,               vx());
            dataExtractQ(e.wx,               wx());
            dataExtractQ(e.xi,               x_prime);
            dataExtractQ(e.q_r4_1,           q_r4_1);
            dataExtractQ(e.q_r4_2,           q_r4_2);

            dataSub(e.vxP,                   e.vx, indices);
            dataSub(e.wxP,                   e.wx, indices);

            // Round5
            dataExtractQ(d.ci_1_prime,       xsum_final().first);
            dataExtractQ(d.ci_2_prime,       xsum_final().second);
            dataExtractQ(d.c_1,              xsum_first().first);
            dataExtractQ(d.c_2,              xsum_first().second);
            dataExtractQ(e.yi,               y_prime);
            dataExtractQ(e.zi,               z_prime);
            dataExtractQ(e.zp,               zp());
            dataExtractQ(e.uz,               uz());
            dataExtractQ(e.vz,               vz());
            dataExtractQ(e.wz,               wz());
            dataExtractQ(e.q_r5_1,           q_r5_1);
            dataExtractQ(e.q_r5_2,           q_r5_2);

            dataSub(e.vzP,                   e.vz, indices);
            dataSub(e.wzP,                   e.wz, indices);

            // Round6
            d.special = speciald();
            dataExtractQ(d.di,               partial_xyz_shares());
            dataExtractQ(d.c_1_prime,        xyz_sum().first);
            dataExtractQ(d.c_2_prime,        xyz_sum().second);
            dataExtractQ(e.rNoise,           g_r());
            dataExtractQ(e.q_r6,             q_r6);

            // Round7
            dataExtractSubVector<NbPrimesQ>(e.y_sharesCAN, y_shares(), candidateIndices().can);
            dataExtractSubVector<NbPrimesQ>(e.y_sharesPS,  y_shares(), candidateIndices().ps);
            dataExtractSubVector<NbPrimesQ>(e.x_sharesCAN, x_shares(), candidateIndices().can);
            dataExtractSubVector<NbPrimesQ>(e.x_sharesPS,  x_shares(), candidateIndices().ps);

            dataExtractVectorCRT<NbPrimesQ>(d.coefsCAN,	    moduli_CAN);
            dataExtractVectorCRT<NbPrimesQ>(d.cans,         alphasCAN);
            dataExtractVectorCRT<NbPrimesQ>(d.prodcans,     alphasPSprod_CAN);

            d.indicesPS.assign(candidateIndices().ps.begin(),   candidateIndices().ps.end());
            d.indicesCAN.assign(candidateIndices().can.begin(), candidateIndices().can.end());
            d.indicesGCD.assign(candidateIndices().gcd.begin(), candidateIndices().gcd.end());

            dataExtractVector(d.by_shares,      by_shares);
            dataExtractVector(d.ax_shares,      ax_shares);
            dataExtractVector(e.q_p_prod_r7,    q_p_prod_r7);
            dataExtractVector(e.q_p_r7,         q_p_r7);
            dataExtractVector(e.q_q_prod_r7,    q_q_prod_r7);
            dataExtractVector(e.q_q_r7,         q_q_r7);

            // Round8
            dataExtractSubVector<NbPrimesQ>(e.z_sharesCAN, z_shares(), candidateIndices().can);
            dataExtractSubVector<NbPrimesQ>(d.ax,		    ax(), axbyindices);
            dataExtractSubVector<NbPrimesQ>(d.by,		    by(), axbyindices);
            dataExtractSubVector<NbPrimesQ>(d.axby,	    axby(), axbyindices);
            dataExtractVector(e.q_r8,           q_r8);

            // Bounding x_shares, y_shares and z_shares
            std::vector<mpz_class> vars;

            pickSubVector(vars, x_shares(), indices);
            pickSubVector(vars, y_shares(), indices);
            pickSubVector(vars, z_shares(), candidateIndices().can);
            
            std::vector<size_t> idxs = {0};
            pickSubVector(vars, ss(), idxs);

            // Bounding

            // x, y and z
            std::vector<mpz_class> cans(alphasCAN);
            cans.insert(cans.end(), alphasPS.begin(), alphasPS.end());
            cans.insert(cans.end(), alphasGCD.begin(), alphasGCD.end());

            // Z and R
            std::vector<mpz_class> mpz_z, mpz_r;
            std::array<mpz_t, Degree> exp_z = zp().poly2mpz();
            std::array<mpz_t, Degree> exp_r = g_r().poly2mpz();
            
            std::transform(exp_z.begin(), exp_z.end(), back_inserter(mpz_z), [](mpz_t in)->mpz_class {return mpz_class(in);});
            std::transform(exp_r.begin(), exp_r.end(), back_inserter(mpz_r), [](mpz_t in)->mpz_class {return mpz_class(in);});

            // pickSubVector(vars, mpz_z, indices);
            pickSubVector(vars, mpz_r, indices);
            for (auto& it : vars) { it = mpz_class(0);}

            std::vector<mpz_class> zbounds;

            for (auto tau : cans) {
                /** verifying the ranging from [±Z], where Z = τ * N * (2 ^ λ) */
                zbounds.emplace_back(tau * tau * config().numParties() * mpz_class(2)^mpz_class(config().lambda()));
            }

            /** verifying the ranging from [±R], where R = (2 ^ λ) * (2 * σ * P * N^1.5 * n) */
            mpz_class pown = config().numParties();
            ligero::lattice::mpz_pow_d(pown.get_mpz_t(), pown.get_mpz_t(), 1.5);

            mpz_class rbound = mpz_class(1) << config().lambda();
            rbound *= mpz_class(2) * mpz_class(config().sigma()) * mpz_class(pown) *
                    mpz_class(Degree) *
                    mpz_class(nfl::poly_p<T, Degree, NbPrimesP>::moduli_product());

            size_t size = alphasCAN.size() + alphasPS.size() + alphasGCD.size();
            std::vector<mpz_class> rbounds(size, rbound);

            // auto params = {cans, cans, alphasCAN, {mpz_class(2)^mpz_class(1234)}, zbounds, rbounds};
            auto params = {cans, cans, alphasCAN, {mpz_class(2)^mpz_class(1234)}, rbounds};

            std::vector<mpz_class> Cs;

            e.bitDecompositions.clear();
            e.XplusCs.clear();
            e.XplusDs.clear();
            size_t idx = 0;

            for (auto param : params) {
                for (auto bound : param) {

                    auto [dd, C, D, log2D] = rephraseAs2Nct(bound);
                    d.log2Ds.emplace_back(log2D);

                    for (size_t modulusIdx = 0; modulusIdx < NbPrimesQ; modulusIdx++) {
                        mpz_class prime(nfl::params<uint64_t>::P[modulusIdx]);
                        mpz_class calc = math::mod(vars[idx] + C, prime);
                        mpz_class calc2 = math::mod(vars[idx] + dd, prime);
                        mpz_class calc3 = math::mod(vars[idx], prime);
                        mpz_class calc4 = math::mod(C, prime);
                        mpz_class calc5 = math::mod(dd, prime);

                        e.XplusCs.emplace_back(calc.get_ui());
                        e.XplusDs.emplace_back(calc2.get_ui());
                        e.Xs.emplace_back(calc3.get_ui());
                        d.Cs.emplace_back(calc4.get_ui());
                        d.Ds.emplace_back(calc5.get_ui());
                    }

                    mpz_class powerOfTwo(1);

                    std::vector<bool> bC, bD;

                    for (size_t bitIndex = 0; bitIndex < log2D; bitIndex++) {
                        bool bitPlusC = false;
                        bool bitPlusD = false;

                        bitPlusC = (((vars[idx] + C) & powerOfTwo) > 0);
                        bitPlusD = (((vars[idx] + dd) & powerOfTwo) > 0);

                        // interleaving bit decompositions as well
                        bC.emplace_back(bitPlusC);
                        bD.emplace_back(bitPlusD);

                        powerOfTwo = powerOfTwo * mpz_class(2);
                    }

                    e.bitDecompositions.emplace_back(bC);
                    e.bitDecompositions.emplace_back(bD);
                    idx++;
                }
            }

            // Rounds 11 & 12

            // Data Extraction
            // ===================================================================================================================
            dataExtractVector(d.by_shares_GCD,  by_shares_GCD);
            dataExtractVector(d.ax_shares_GCD,  ax_shares_GCD);
            dataExtractVector(e.q_p_prod_r11,   q_p_prod_r11);
            dataExtractVector(e.q_pq_r11,       q_pq_r11);
            dataExtractVector(e.q_q_prod_r11,   q_q_prod_r11);
            dataExtractVector(e.q_r_r11,        q_r_r11);
            dataExtractVector(e.r_CRTs,         r_CRTs);

            dataExtractVector(e.q_r12,              q_r12);
            dataExtractVector(d.finalModuli_GCD,    finalModuli_GCD);
            dataExtractVector(e.ss_GCD,	            ss_GCD);

            dataExtractSubVector<NbPrimesQ>(e.y_sharesGCD,  y_shares(), candidateIndices().gcd);
            dataExtractSubVector<NbPrimesQ>(e.x_sharesGCD,  x_shares(), candidateIndices().gcd);
            dataExtractSubVector<NbPrimesQ>(e.z_sharesGCD,  z_shares(), candidateIndices().gcd);

            dataExtractVectorCRT<NbPrimesQ>(d.coefsGCD,	    moduli_GCD);
            dataExtractVectorCRT<NbPrimesQ>(d.gcds,         alphasGCD);
            dataExtractVectorCRT<NbPrimesQ>(d.prodgcds,     alphasPSprod_GCD);
            dataExtractVectorCRT<NbPrimesQ>(d.axGCD,        axGCD());
            dataExtractVectorCRT<NbPrimesQ>(d.byGCD,        byGCD());
            dataExtractVectorCRT<NbPrimesQ>(d.axbyGCD,	    axbyGCD());

            // Computing sigma_x and exp_q
            // ===================================================================================================================
            sigma_x_GCD()  = std::vector<mpz_class>(NbPrimesP * alphasGCD.size());
            exp_q_GCD()    = std::vector<mpz_class>(129 * NbPrimesP * alphasGCD.size());;

            std::vector<mpz_class> expGCD;

            for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                mpz_class exponentGCD = math::mod(gammaExponent(), alphasGCD[alphaGCDidx]);
                mpz_class exponentGCDquotient = (4 * exponentGCD - ((speciald() ? math::mod(finalModuli(), alphasGCD[alphaGCDidx]) + 1 : mpz_class(0)) - expsharesGCD[alphaGCDidx]));

                assert(math::mod(exponentGCDquotient, alphasGCD[alphaGCDidx]) == 0);
                exponentGCDquotient /= alphasGCD[alphaGCDidx];

                for (size_t i = 0; i < NbPrimesP; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        sigma_x_GCD()[i * alphasGCD.size() + alphaGCDidx] = positiveRemainder(exponentGCD, prime);
                        expGCD.push_back(sigma_x_GCD()[i * alphasGCD.size() + alphaGCDidx]);
                        exp_q_GCD()[i * alphasGCD.size() + alphaGCDidx] = positiveRemainder(exponentGCDquotient, prime);
                }
            }

            // Computing sigma_e, sigma_z and sigma_q
            // ===================================================================================================================
            std::vector<mpz_class> sigmas_e(128), sigmas_z(128);
            sigma_a_GCD() = std::vector<mpz_class> (128 * alphasGCD.size() * NbPrimesP);
            sigma_g_GCD() = std::vector<mpz_class> (128 * alphasGCD.size() * NbPrimesP);
            sigma_e_GCD() = std::vector<mpz_class> (128 * alphasGCD.size() * NbPrimesP);
            sigma_z_GCD() = std::vector<mpz_class> (128 * alphasGCD.size() * NbPrimesP);
            sigma_q_GCD() = std::vector<mpz_class> (128 * alphasGCD.size() * NbPrimesP);
            
            std::vector<mpz_class> sigmas_a(128);

            // Generate the Seed for Fiat-Shamir Randomness Generation
            std::string msg1 = serialize<std::vector<mpz_class>>(sigma_r_GCD());
            std::string msg2 = serialize<PublicData>(d);
            std::string msg = msg1 + msg2;

            std::vector<unsigned char> seed = hash::blake2bStringHash(msg.c_str(), msg.size());
            uint64_t randomness[128*2];
            randombytes_buf_deterministic((void *)&randomness[0], 128*2*sizeof(uint64_t), &seed[0]);
            mpz_class higherBitsScalar = mpz_class(2)^mpz_class(64);

            // First completing the sigma ZK protocol and then verifying the main ZK proof
            for (size_t j = 0 ; j < 128; j++) {

                mpz_class sigma_a = math::powm(gBeforeExp()[j], sigmas_r()[j], finalModuli());

                sigmas_a[j] = sigma_a;
                sigmas_e[j] = mpz_class(randomness[j*2])+(higherBitsScalar*mpz_class(randomness[j*2+1]));

                sigmas_z[j] = sigmas_r()[j] + sigmas_e[j] * gammaExponent();
		if(j==0)
		{
			LOG(INFO) << "sigma a = " << sigmas_a[0];
			LOG(INFO) << "sigma e = " << sigmas_e[0];
			LOG(INFO) << "sigma z = " << sigmas_z[0];
		}
                for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                    mpz_class sigma_gGCD = math::mod(gAfterExp()[j], alphasGCD[alphaGCDidx]);                    
                    mpz_class sigma_aGCD = math::mod(sigmas_a[j], alphasGCD[alphaGCDidx]);
                    mpz_class sigma_xGCD = math::mod(gammaExponent(), alphasGCD[alphaGCDidx]);
                    mpz_class sigma_rGCD = math::mod(sigmas_r()[j], alphasGCD[alphaGCDidx]);
                    mpz_class sigma_eGCD = math::mod(sigmas_e[j], alphasGCD[alphaGCDidx]);
                    mpz_class sigma_zGCD = math::mod(sigmas_z[j], alphasGCD[alphaGCDidx]);

                    mpz_class dividend = sigma_zGCD - (sigma_rGCD + sigma_eGCD * sigma_xGCD);

                    assert(dividend % alphasGCD[alphaGCDidx] == 0);
                    mpz_class quotient = dividend / alphasGCD[alphaGCDidx];

                    for (size_t i = 0; i < NbPrimesP; i++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        assert(positiveRemainder(sigma_xGCD,prime) == sigma_x_GCD()[i * alphasGCD.size() + alphaGCDidx]); 
                        sigma_g_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(sigma_gGCD,prime);
                        sigma_a_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(sigma_aGCD,prime);
                        sigma_e_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(sigma_eGCD,prime);
                        sigma_q_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(quotient,prime);
                        sigma_z_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(sigma_zGCD,prime);

                        assert( 
                            math::mod(
                                sigma_x_GCD()[i * alphasGCD.size() + alphaGCDidx]
                                    * sigma_e_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size()+ alphaGCDidx]
                                + sigma_r_GCD()[i * alphasGCD.size() * 128 + j * alphasGCD.size() + alphaGCDidx] 
                                + math::mod(alphasGCD[alphaGCDidx],prime)
                                    * sigma_q_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx]
                                - sigma_z_GCD() [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx]
                            , prime) == 0);
                    }

                }
            }

            // Extracting Variables for Sigma Protocol
            // ===================================================================================================================
            dataExtractVector(e.sigmarGCD,      sigma_r_GCD());
            dataExtractVector(e.sigmaxGCD,      sigma_x_GCD());
            dataExtractVector(e.expqGCD,        exp_q_GCD());
            dataExtractVector(e.sigmaqGCD,      sigma_q_GCD());

            // sigma protocol
            dataExtractVector(s.sigmazGCD,      sigma_z_GCD());
            dataExtractVector(s.sigmaeGCD,      sigma_e_GCD());
            dataExtractVector(s.sigmaaGCD,      sigma_a_GCD());
            dataExtractVector(s.sigmagGCD,      sigma_g_GCD());

            e.isAvailable = true;

            // Checking Extractions [Mirrored by ZK tests]
            // ===================================================================================================================

            // Checking extraction for Rounds 7 & 8
            for (size_t p = 0; p < NbPrimesQ; p++) {
                for (size_t alphaCANidx = 0; alphaCANidx < alphasCAN.size(); alphaCANidx++) {

                    size_t i = alphaCANidx * NbPrimesQ;
                    mpz_class rem4 = mpz_class(0);
                    if (speciald()) {
                        rem4 = mpz_class(3);
                    }

                    // Round 7
                    // Checking p equation
                    {
                        mpz_class rem = mpz_class(d.ax_shares[i + p]) 
                                        + mpz_class(e.x_sharesCAN[i + p]) 
                                        - (
                                            mpz_class(e.x_sharesPS[p]) * mpz_class(d.coefsCAN[7 * i + p])
                                            + mpz_class(e.x_sharesPS[NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + NbPrimesQ + p])
                                            + mpz_class(e.x_sharesPS[2 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 2 * NbPrimesQ + p])
                                            + mpz_class(e.x_sharesPS[3 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 3 * NbPrimesQ + p])
                                            + mpz_class(e.x_sharesPS[4 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 4 * NbPrimesQ + p])
                                            + mpz_class(e.x_sharesPS[5 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 5 * NbPrimesQ + p])
                                            + rem4 * mpz_class(d.coefsCAN[7 * i + 6 * NbPrimesQ + p])
                                            - mpz_class(e.q_p_prod_r7[i + p]) * mpz_class(d.prodcans[i + p])
                                            )
                                        - mpz_class(e.q_p_r7[i + p]) * mpz_class(d.cans[i + p]);
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                    //Checking q equation
                    {
                        mpz_class rem = mpz_class(d.by_shares[i + p]) + mpz_class(e.y_sharesCAN[i + p]) -
                                        (mpz_class(e.y_sharesPS[p]) * mpz_class(d.coefsCAN[7 * i + p])
                                         + mpz_class(e.y_sharesPS[NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + NbPrimesQ + p])
                                         + mpz_class(e.y_sharesPS[2 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 2 * NbPrimesQ + p])
                                         + mpz_class(e.y_sharesPS[3 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 3 * NbPrimesQ + p])
                                         + mpz_class(e.y_sharesPS[4 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 4 * NbPrimesQ + p])
                                         + mpz_class(e.y_sharesPS[5 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 5 * NbPrimesQ + p])
                                         + rem4 * mpz_class(d.coefsCAN[7 * i + 6 * NbPrimesQ + p])
                                         - mpz_class(e.q_q_prod_r7[i + p]) * mpz_class(d.prodcans[i + p]))
                                         - mpz_class(e.q_q_r7[i + p]) * mpz_class(d.cans[i + p]);
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                    //Round 8
                    {
                        mpz_class rem =
                            mpz_class(d.ax[i + p]) *
                            (mpz_class(e.y_sharesPS[p])                   * mpz_class(d.coefsCAN[7 * i + 0 * NbPrimesQ + p])
                             + mpz_class(e.y_sharesPS[NbPrimesQ + p])     * mpz_class(d.coefsCAN[7 * i + 1 * NbPrimesQ + p])
                             + mpz_class(e.y_sharesPS[2 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 2 * NbPrimesQ + p])
                             + mpz_class(e.y_sharesPS[3 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 3 * NbPrimesQ + p])
                             + mpz_class(e.y_sharesPS[4 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 4 * NbPrimesQ + p])
                             + mpz_class(e.y_sharesPS[5 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 5 * NbPrimesQ + p])
                             + rem4                                       * mpz_class(d.coefsCAN[7 * i + 6 * NbPrimesQ + p])
                             - mpz_class(e.q_q_prod_r7[i + p]) * mpz_class(d.prodcans[i + p]))
                             + mpz_class(d.by[i + p]) * 
                             (mpz_class(e.x_sharesPS[p])                  * mpz_class(d.coefsCAN[7 * i + 0 * NbPrimesQ + p])
                             + mpz_class(e.x_sharesPS[NbPrimesQ + p])     * mpz_class(d.coefsCAN[7 * i + 1 * NbPrimesQ + p])
                             + mpz_class(e.x_sharesPS[2 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 2 * NbPrimesQ + p])
                             + mpz_class(e.x_sharesPS[3 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 3 * NbPrimesQ + p])
                             + mpz_class(e.x_sharesPS[4 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 4 * NbPrimesQ + p])
                             + mpz_class(e.x_sharesPS[5 * NbPrimesQ + p]) * mpz_class(d.coefsCAN[7 * i + 5 * NbPrimesQ + p])
                             + rem4                                       * mpz_class(d.coefsCAN[7 * i + 6 * NbPrimesQ + p])
                             - mpz_class(e.q_p_prod_r7[i + p]) * mpz_class(d.prodcans[i + p]))
                            + mpz_class(e.z_sharesCAN[i + p]) - mpz_class(d.axby[i + p])
                            - mpz_class(e.q_r8[i + p]) * mpz_class(d.cans[i + p]);
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }
				}
			}

            // Round 11 and 12
            for (size_t p = 0; p < NbPrimesP; p++) {
                for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                    size_t i = alphaGCDidx * NbPrimesQ;
                    mpz_class rem4 = mpz_class(0);

                    if (speciald()) {
                                    rem4 = mpz_class(3);
                    }

                    {
                        // Constraint for (p + q - 1)/4 = x 
                        // x : sigma_x_GCD 
                        // p : recomputed_p_share_GCD
                        // q : recomputed_q_share_GCD
                        mpz_class prime                 = mpz_class(nfl::params<uint64_t>::P[p]);
                        mpz_class finalMod              = math::mod(finalModuli(),alphasGCD[alphaGCDidx]);
                        mpz_class out                   = positiveRemainder(finalMod,prime);
                        const mpz_class maybeNplusOne   = d.special ? mpz_class(d.finalModuli_GCD[i+p]) + mpz_class(1) : mpz_class(0);

                        mpz_class rem = 
                                    (
                                            // recomputed_p_share_GCD
                                        (         mpz_class(e.x_sharesPS[0 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 0 * NbPrimesQ + p])
                                                + mpz_class(e.x_sharesPS[1 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 1 * NbPrimesQ + p])
                                                + mpz_class(e.x_sharesPS[2 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 2 * NbPrimesQ + p])
                                                + mpz_class(e.x_sharesPS[3 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 3 * NbPrimesQ + p])
                                                + mpz_class(e.x_sharesPS[4 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 4 * NbPrimesQ + p])
                                                + mpz_class(e.x_sharesPS[5 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 5 * NbPrimesQ + p])
                                                + rem4                                       * mpz_class(d.coefsGCD[7 * i + 6 * NbPrimesQ + p])
                                                - mpz_class(e.q_p_prod_r11[i + p]) * mpz_class(d.prodgcds[i + p]))
                                            // recomputed_q_share_GCD
                                        + (       mpz_class(e.y_sharesPS[0 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 0 * NbPrimesQ + p])
                                                + mpz_class(e.y_sharesPS[1 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 1 * NbPrimesQ + p])
                                                + mpz_class(e.y_sharesPS[2 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 2 * NbPrimesQ + p])
                                                + mpz_class(e.y_sharesPS[3 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 3 * NbPrimesQ + p])
                                                + mpz_class(e.y_sharesPS[4 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 4 * NbPrimesQ + p])
                                                + mpz_class(e.y_sharesPS[5 * NbPrimesQ + p]) * mpz_class(d.coefsGCD[7 * i + 5 * NbPrimesQ + p])
                                                + rem4                                       * mpz_class(d.coefsGCD[7 * i + 6 * NbPrimesQ + p])
                                                - mpz_class(e.q_q_prod_r11[i + p]) * mpz_class(d.prodgcds[i + p]))
                                                - maybeNplusOne
                                            ) 
                                        - exp_q_GCD()[p*alphasGCD.size()+alphaGCDidx]  * d.gcds[i+p]
                                        + mpz_class(4*sigma_x_GCD()[p * alphasGCD.size() + alphaGCDidx]);

                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                    LOG(INFO) << "Verified sigma" ;
  
                    {
                        const mpz_class maybeOne = d.special ? mpz_class(1) : mpz_class(0);
                        mpz_class rem = mpz_class(d.by_shares_GCD[i + p]) + mpz_class(e.y_sharesGCD[i + p]) 
                                        - (mpz_class(e.x_sharesPS[0 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 0 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[1 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 1 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[2 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 2 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[3 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 3 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[4 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 4 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[5 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 5 * NbPrimesQ + p])
                                        + rem4 					                    * mpz_class(d.coefsGCD[7 * i + 6 * NbPrimesQ + p])
                                        - mpz_class(e.q_p_prod_r11[i + p]) * mpz_class(d.prodgcds[i + p])
                                        + mpz_class(e.y_sharesPS[0 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 0 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[1 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 1 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[2 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 2 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[3 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 3 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[4 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 4 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[5 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 5 * NbPrimesQ + p])
                                        + rem4 					                    * mpz_class(d.coefsGCD[7 * i + 6 * NbPrimesQ + p])
                                        - mpz_class(e.q_q_prod_r11[i + p]) * mpz_class(d.prodgcds[i + p])
                                        - maybeOne)
                                        - mpz_class(e.q_pq_r11[i + p]) * mpz_class(d.gcds[i + p]);
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                    LOG(INFO) << "Verified round11" ;

                    {
                        mpz_class rem = mpz_class(d.ax_shares_GCD[i + p])
                                        + mpz_class(e.x_sharesGCD[i + p])
                                        - mpz_class(e.r_CRTs[i + p])
                                        - mpz_class(e.q_r_r11[i + p]) * mpz_class(d.gcds[i + p]);
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                    {
                        const mpz_class maybeOne = d.special ? mpz_class(1) : mpz_class(0);
                        mpz_class rem = mpz_class(d.axbyGCD[i + p]) 
                                - mpz_class(d.axGCD[i + p]) * 
                                        (
                                        mpz_class(e.x_sharesPS[0 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 0 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[1 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 1 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[2 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 2 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[3 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 3 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[4 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 4 * NbPrimesQ + p])
                                        + mpz_class(e.x_sharesPS[5 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 5 * NbPrimesQ + p])
                                        + rem4 					                    * mpz_class(d.coefsGCD[7 * i + 6 * NbPrimesQ + p])
                                        - mpz_class(e.q_p_prod_r11[i + p]) * mpz_class(d.prodgcds[i + p])
                                        + mpz_class(e.y_sharesPS[0 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 0 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[1 * NbPrimesQ + p])   * mpz_class(d.coefsGCD[7 * i + 1 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[2 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 2 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[3 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 3 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[4 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 4 * NbPrimesQ + p])
                                        + mpz_class(e.y_sharesPS[5 * NbPrimesQ + p])	* mpz_class(d.coefsGCD[7 * i + 5 * NbPrimesQ + p])
                                        + rem4 					                    * mpz_class(d.coefsGCD[7 * i + 6 * NbPrimesQ + p])
                                        - mpz_class(e.q_q_prod_r11[i + p]) * mpz_class(d.prodgcds[i + p])
                                        - maybeOne)
                                                - mpz_class(d.byGCD[i + p]) * mpz_class(e.r_CRTs[i + p])
                                                - mpz_class(e.z_sharesGCD[i + p]) 
                                                - mpz_class(e.ss_GCD[i + p]) * mpz_class(d.finalModuli_GCD[i + p]) 
                                                + mpz_class(e.q_r12[i + p]) * mpz_class(d.gcds[i + p]);

                                        LOG(INFO) << "i = " << i << "and p = " << p;
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

        		    LOG(INFO) << "Verified round12";
                }
            }
            
	        LOG(INFO) << "finished gathering" ;
        }

        /* ========== variables used for zero knowledge ========== */
    protected:
        // variables from `generateKeyPair`
        Q _siP, _eiP, _si, _ei, _bi, _A, _ai, _b;
        mpz_class _jacobiSeedShares, _jacobiAndGCDSeedShares;

        std::vector<mpz_class> _ss;

        bool _special;
        // variables from `generatePreSievingShares`
        //
        Q _zp;
        P _p_y_shares;
        Q _ux, _vx, _wx;
        Q _uz, _vz, _wz;
        Q _up, _vp, _wp;
        Q _partial_xyz_shares;
        Q _g_r;
        CandidateIndices _candidateIndices;
        std::vector<mpz_class> _moduli;
        std::vector<std::vector<mpz_class>> _coefs;
        std::vector<size_t> _index_candidates;
        std::vector<mpz_class> _sigma_r_GCD;
        std::vector<mpz_class> _sigma_x_GCD;
        std::vector<mpz_class> _sigma_a_GCD;
        std::vector<mpz_class> _sigma_g_GCD;
        std::vector<mpz_class> _sigma_e_GCD;
        std::vector<mpz_class> _sigma_z_GCD;
        std::vector<mpz_class> _sigma_q_GCD;
        std::vector<mpz_class> _exp_q_GCD;
        std::vector<mpz_class> _sigmas_r;

        uint32_t _finalindex;

        mpz_class _finalModuli,_gammaExponent;


        std::vector<mpz_class> _candidatesCAN;
        std::vector<mpz_class> _candidatesGCD;

        std::vector<mpz_class> _x_shares;
        std::vector<mpz_class> _y_shares;
        std::vector<mpz_class> _z_shares;
        lattice::cipher<T, Degree, NbPrimesQ>  _xsum_first, _xsum_final, _enc_z_shares, _xyz_sum, _enc_x_shares;

        std::vector<mpz_class> _p_sharesCRTs;
        std::vector<mpz_class> _q_sharesCRTs;

        //Round 8
        std::vector<mpz_class> _ax;
        std::vector<mpz_class> _by;
        std::vector<mpz_class> _axby;

        //Round 11 and 12
        std::vector<mpz_class> _rCRTs;
        std::vector<mpz_class> _axGCD;
        std::vector<mpz_class> _byGCD;
        std::vector<mpz_class> _axbyGCD;
        std::vector<mpz_class> _ssGCD;

        mpz_class               _gcdRX;
        mpz_class               _gammaSeed;

        // sigma
        
        std::vector<mpz_class> _gBeforeExp;
        std::vector<mpz_class> _gAfterExp;

        // legacy
        std::vector<mpz_class>  _a_shares, _b_shares;
        cipherText              _as, _es;
        std::vector<cipherText> _as_vec, _es_vec;
        Q                       _partial_e_shares;
        std::vector<int>        _sievingFlags;

        // others
        SocketId&               _socketId;
        ProtocolConfig<T>       _config;

        // passive protocol
        std::unique_ptr<nfl::FastGaussianNoise<uint16_t, T, 2>> _fg;
        std::unique_ptr<nfl::gaussian<uint16_t, T, 2>> _chi;
        std::unique_ptr<lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>> _e;
};

/** Registers party and awaits for configuration from the protocol coordinator
 * @param transport communication helper
 * @param ipAddress coordinator's IP address
 * @param chi protocol configuration parameter
 * @param prover zero-knowledge prover
 * @return protocol configuration
 */
template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
expected<ProtocolConfig<T>> EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::registerAndAwaitConfiguration(
    ZeroMQClientTransport &transport, std::string ipAddress, nfl::gaussian<uint16_t, T, 2>& chi, zksnark::Prover<FieldT>* prover) {

    DBG("Participant (" << transport.getSocketId() << "," << ipAddress << ") registering with coordinator.");

    // gaussian dist
    siP() = Q(chi);
    eiP() = Q(chi);

    auto positiveRemainder = [&](mpz_class & a, mpz_class & b) -> mpz_class {
        return ((a % b + b) % b);
    };

    auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3210, primesGCD,175);
    std::vector<mpz_class> alphasPSprod_GCD, moduli_GCD;
    std::vector<mpz_class> ax_shares_GCD;
    std::vector<mpz_class> by_shares_GCD;
    std::vector<mpz_class> expsharesGCD;

    // Sigma Protocol
    // Step 1: Commit to r and x before executing the sigma protocol
    
    sigmas_r()     = std::vector<mpz_class>(128);
    sigma_r_GCD()  = std::vector<mpz_class>(129 * NbPrimesP * alphasGCD.size());

    // assuming that second jacobi test only has one moduli and security parameter is 128.
    for (size_t j = 0; j < 128; j++) {

        mpz_class sigma_r = math::generateRandomValue(mpz_class(std::random_device()()), 2048 + 128);
        sigmas_r()[j] = sigma_r;

        for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
            mpz_class sigma_rGCD = math::mod(sigma_r, alphasGCD[alphaGCDidx]);
            
            for (size_t i = 0; i < NbPrimesP; i ++) {
                mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                sigma_r_GCD()[i * alphasGCD.size() * 128 + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(sigma_rGCD, prime);
            }
        }
    }

    // Primitive Roots of Unity for All Moduli
    nfl::poly_p<uint64_t, 65536, 21> roots({0});
    for (size_t idx = 0; idx < 21; idx++) {roots(idx,1) = 1ull;}
    roots.ntt_pow_phi();
    std::vector<std::string> earlyWitnessHash;
    
    for (size_t modulusIdx = 0; modulusIdx < NbPrimesP; modulusIdx++) {
        // Roots of Unity
        prover->_publicData.modulusIdx = modulusIdx;
        prover->_publicData.roots.assign(roots.poly_obj().data()+prover->_publicData.modulusIdx*(size_t)roots.degree,roots.poly_obj().data()+(prover->_publicData.modulusIdx+1)*(size_t)roots.degree);

        // Zero-Knowledge
        p_zk = nfl::params<uint64_t>::P[prover->_publicData.modulusIdx];
        FieldT::modulusIdx_ = prover->_publicData.modulusIdx;

        dataExtractQ(prover->_secretData.eiP, eiP());
        dataExtractQ(prover->_secretData.siP, siP());
        dataExtractVector(prover->_secretData.sigmarGCD, sigma_r_GCD());

        prover->produceArgumentOfKnowledge(1, true);

        hash::digest& d = prover->_transcript[0].proverCommitment_early;
        earlyWitnessHash.emplace_back(d.begin(), d.end());
    }

    // Generate all values we need to commit to coordinator
    // first uniformly choose a
    ai() = Q(nfl::uniform{});
    ai().ntt_pow_phi();
    // generate seed1 for jacobi
    jacobiSeedShares() = math::generateRandomValue(mpz_class(std::random_device()()));
    // generate seed2 for jacobiAndGCD
    jacobiAndGCDSeedShares() = math::generateRandomValue(mpz_class(std::random_device()()));

    transport.send(MessageType::ID_PARTY, 
            PartyCommitment(
                ipAddress,
                earlyWitnessHash,
                math::computeHash(ai()),
                math::computeHash(jacobiSeedShares()),
                math::computeHash(jacobiAndGCDSeedShares()),
                transport.publicKey
                ));

    auto success = transport.awaitReply<ProtocolConfig<T>> (MessageType::PROTOCOL_CONFIG);
    if (hasError(success)) { return getError(success); }

    return getResult(success);
}

/** Generate key pair in a distributed way.
 * @param trans communication helper
 * @param chi protocol parameter
 * @return a key pair on success and error otherwise
 */
template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
expected<lattice::key_pair<T, Degree, NbPrimesQ>> EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::generateKeyPair(ZeroMQClientTransport& trans, nfl::gaussian<uint16_t, T, 2>& chi) {

    DBG("Generating Keys");
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

    // NOTE: ai is now initialized in registerAndAwaitConfiguration
    
    DBG("Sending A shares");
    trans.send<Q>(MessageType::PUBLIC_KEY_A_SHARES, ai());
    auto maybe_a = trans.awaitReply<Q>(MessageType::PUBLIC_KEY_A_VALUE);
    if (hasError(maybe_a)) {
        LOG(ERROR) << "Kill/Restart received during keygen A";
        return getError(maybe_a);
    }
    A() = getResult(maybe_a);
    DBG("Received A");

    si() = siP();
    ei() = eiP();

    si().ntt_pow_phi();
    ei().ntt_pow_phi();

    bi() = si() * A() + ei();

    DBG("Sending B shares");
    trans.send<Q>(MessageType::PUBLIC_KEY_B_SHARES, bi());
    auto maybe_b = trans.awaitReply<Q>(MessageType::PUBLIC_KEY_B_VALUE);
    if (hasError(maybe_b)) {
        LOG(ERROR) << "Kill/Restart received during keygen B";
        return getError(maybe_b);
    }
    DBG("Received B");
    b() = getResult(maybe_b);

    // return (public key, private key)
    return std::make_pair(std::make_pair(A(), b()), si());
}

/** Helper function to prune and reorder party's shares.
 * @param as shares
 * @param bs shares
 * @param flags used to prune as and bs
 * @param numAlphas number of alpha's in the buckets
 * @param bucketSize array of bucket sizes
 * @return a matrix with valid shares only
 */
template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::vector<std::vector<mpz_class>>
                                 EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::pruneAndReorderShares(
                                     const std::vector<mpz_class> &as,
                                     const std::vector<mpz_class> &bs,
                                     const std::vector<int> &flags,
                                     int numAlphas,
std::vector<size_t> bucketSize) {

    assert(as.size() == flags.size());
    assert(bs.size() == flags.size());

    // Now we can construct a result matrix
    std::vector<std::vector<mpz_class>> result;

    // And finally we fill the result matrix compactly with valid shares from `as`
    int k = 0;
    for (int c = 0; c < numAlphas; ++c) {
        int j = 0;

        std::vector<mpz_class> col;
        for (int r = 0; r < bucketSize[c]; ++r) {
            if (flags[k] == 1) {
                col.push_back(as[k]);
                col.push_back(bs[k]);
                index_candidates()[k] = j++;
            } else {
                index_candidates()[k] = -1;
            }
            k++;
        }
        result.push_back(col);
    }

    return result;
}


/** Runs pre-sieving rounds of the protocol
 * @param e encryption helper
 * @param transport communication helper
 * @param publicKey used for encryption
 * @param secretKey used for partial decryption
 * @param special flag for a special party
 * @param config protocol configuration
 */
template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
expected<std::tuple<
std::vector<mpz_class>, // xcan
    std::vector<mpz_class>, // ycan
    std::vector<mpz_class>, // zcan
    std::vector<mpz_class>, // xgcd
    std::vector<mpz_class>, // ygcd
    std::vector<mpz_class>, // zgcd
    std::vector<mpz_class> // both_shares
>>
    EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::generatePreSievingShares(
        lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ> &e,
        ZeroMQClientTransport &transport,
        const lattice::pub<T, Degree, NbPrimesQ> &publicKey,
        const lattice::sec<T, Degree, NbPrimesQ> &secretKey,
        bool special,
        ProtocolConfig<T> &config
) {

    assert(primesPS + primesCAN + primesGCD <= Degree);

    using P = nfl::poly_p<T, Degree, NbPrimesP>;
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

    // For the pre-sieving part of our candidate generation, we need prime
    // buckets s.t. the product of buckets has a bit-length greater than 1024.
    auto [alphasPS, bucketSizePS] = math::balanced_bucket_n_primes(config.pbs(), primesPS, config.tauLimitBit(), 1);
    auto [alphasCAN, bucketSizeCAN] = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());
    auto [alphasGCD, bucketSizeGCD] = math::fixed_bucket_n_primes(3*config.pbs()+210, primesGCD, config.tauLimitBit());

    // set equal buckets for bucketSizeCan and bucketSizeGCD
    const int bucketSizeCAN_value = lrint(floor(double(primesCAN) / double(
            alphasCAN.size())));
    for (int i = 0; i < bucketSizeCAN.size(); ++i) {
        bucketSizeCAN[i] = bucketSizeCAN_value;
    }
    const int bucketSizeGCD_value = lrint(floor(double(primesGCD) / double(
            alphasGCD.size())));
    for (int i = 0; i < bucketSizeGCD.size(); ++i) {
        bucketSizeGCD[i] = bucketSizeGCD_value;
    }



    auto [alphasdummy, bucketSizedummy] = math::balanced_bucket_n_primes(
            config.pbs(), Degree - primesTOTAL, config.tauLimitBit(), 1);


    a_shares().resize(primesPS);
    b_shares().resize(primesPS);

    LOG(INFO) << "Check protocol mode";
    if (config.protocolMode() == ProtocolMode::NORMAL
        || config.protocolMode() == ProtocolMode::RECORD) {
        x_shares() = math::generateRandomVector(std::random_device()(), Degree,
                config.tauLimitBit());
        y_shares() = math::generateRandomVector(std::random_device()(), Degree,
                config.tauLimitBit());
        z_shares() = math::generateRandomVector(std::random_device()(), Degree,
                                            config.tauLimitBit());

        if (config.protocolMode() == ProtocolMode::RECORD) {

            LOG(INFO) << "Protocol mode record";
            transport.send(MessageType::RECORD_PROTOCOL_SHARES, std::make_pair(x_shares(), std::make_pair(y_shares(), z_shares())));
            LOG(INFO) << "Sent mode record";

            auto maybe_rep = transport.awaitReply<int>(MessageType::RECORD_PROTOCOL_RESPONSE);
            LOG(INFO) << "Continued after record";

        }
    } else if (config.protocolMode() == ProtocolMode::REPLAY) {
        LOG(INFO) << "Protocol replay";

        transport.send(MessageType::REPLAY_PROTOCOL_SHARES, int(1));

        LOG(INFO) << "Sent out that we are ready for replay";
        auto maybe_replay = transport.awaitReply<TripleVector>(MessageType::REPLAY_PROTOCOL_RESPONSE);
        if (hasError(maybe_replay)) {
            LOG(ERROR) << "Kill/Restart received during replay";
            return getError(maybe_replay);
        }
        auto [_xs, _yszs] = getResult(maybe_replay);
        x_shares() = _xs;
        std::tie(y_shares(), z_shares()) = _yszs;
    } else {
        throw std::runtime_error("Unknown Protocol execution mode.");
    }

    index_candidates() = std::vector<size_t>(Degree);
    std::vector<mpz_class> xcan(primesCAN);
    std::vector<mpz_class> ycan(primesCAN);
    std::vector<mpz_class> zcan(primesCAN);
    std::vector<mpz_class> xgcd(primesGCD);
    std::vector<mpz_class> ygcd(primesGCD);
    std::vector<mpz_class> zgcd(primesGCD);

    std::array<mpz_class, Degree> tauVector;

    // so here we will pack values as following:
    {
        gmp_randclass grc(gmp_randinit_default);
        grc.seed(std::random_device()());
        int k = 0;
        // pre-sieving
        for (size_t j = 0; j < alphasPS.size(); ++j) {
            for (size_t i = 0; i < bucketSizePS[j]; ++i) {
                x_shares()[k] = x_shares()[k] % alphasPS[j];
                y_shares()[k] = y_shares()[k] % alphasPS[j];
                z_shares()[k] = 0;
                a_shares()[k] = x_shares()[k];
                b_shares()[k] = y_shares()[k];
                tauVector[k] = alphasPS[j];
                k++;
            }
        }

        // candidate generation
        int can_i = 0;
        for (size_t i = 0; i < bucketSizeCAN_value; ++i) {
            for (size_t j = 0; j < alphasCAN.size(); ++j) {
                index_candidates()[k] = i;
                x_shares()[k] = x_shares()[k] % alphasCAN[j];
                y_shares()[k] = y_shares()[k] % alphasCAN[j];
                z_shares()[k] = z_shares()[k] % alphasCAN[j];
                tauVector[k] = alphasCAN[j];
                xcan[can_i] = x_shares()[k];
                ycan[can_i] = y_shares()[k];
                zcan[can_i] = z_shares()[k];
                can_i++;
                k++;
            }
        }

        // GCD
        int gcd_i = 0;
        for (size_t i = 0; i < bucketSizeGCD_value; ++i) {
            for (size_t j = 0; j < alphasGCD.size(); ++j) {
                if (i < 1) {
                    index_candidates()[k] = -2;
                } else {
                    index_candidates()[k] = -3;
                }
                x_shares()[k] = x_shares()[k] % alphasGCD[j];
                y_shares()[k] = y_shares()[k] % alphasGCD[j];
                z_shares()[k] = z_shares()[k] % alphasGCD[j];
                tauVector[k] = alphasGCD[j];
                xgcd[gcd_i] = x_shares()[k];
                ygcd[gcd_i] = y_shares()[k];
                zgcd[gcd_i] = z_shares()[k];
                gcd_i++;
                k++;
            }
        }

        // setting up tauVector for dummy
        for (size_t j = 0; j < alphasdummy.size(); ++j) {
            for (size_t i = 0; i < bucketSizedummy[j]; ++i) {
                y_shares()[k] = 0;
                tauVector[k] = alphasdummy[j];
                k++;
            }
        }
    }

    // Pi sends Enc(xi) to coordinator. Coordinator aggregates all encryptions and returns Enc(\sum x_i) = Enc(x)
    DBG("Participant (" << transport.getSocketId() << ") sending `Enc(xi)` shares");

    auto [enc_x_shares, uxi, vxi, wxi] = e.encrypt(x_shares(), publicKey);
    ux() = uxi;
    vx() = vxi;
    wx() = wxi;
    g_enc_x_shares().first = enc_x_shares.first;
    g_enc_x_shares().second = enc_x_shares.second;

    transport.send(MessageType::ENCRYPTED_X_SHARES, enc_x_shares);

    auto maybe_xsum = transport.awaitReply<lattice::cipher<T, Degree, NbPrimesQ>>
                      (MessageType::ENCRYPTED_X_VALUE);
    if (hasError(maybe_xsum)) {
        LOG(ERROR) << "Kill/Restart received during presieve";
        return getError(maybe_xsum);
    }
    auto xsum = getResult(maybe_xsum);
    xsum_first() = xsum;
    DBG("Participant (" << transport.getSocketId() << ") received `Enc(a)`");

    auto [enc_z_s, uzi, vzi, wzi] = e.encrypt(z_shares(), publicKey);
    enc_z_shares() = enc_z_s;
    uz() = uzi;
    vz() = vzi;
    wz() = wzi;

    {
        auto *ptr = new std::array<mpz_class, Degree> {};
        std::transform(y_shares().begin(), y_shares().end(),
        ptr->begin(), [](const auto & x) {
            return x;
        });
        p_y_shares().set_mpz(*ptr);
        delete (ptr);

        auto [ep, zpi] = e.product(xsum, p_y_shares(), tauVector, publicKey);
        zp() = zpi;

        xsum = ep;
        xsum = lattice::pair_sub(xsum, enc_z_shares());
        xsum_final() = xsum;

        DBG("Participant (" << transport.getSocketId() <<
            ") sending `Enc(x) * y_i - Enc(zi)` shares");
        transport.send(MessageType::ENCRYPTED_XY_PLUS_Z_SHARES, xsum);
    }

    auto maybe_xyz_sum =
        transport.awaitReply<lattice::cipher<T, Degree, NbPrimesQ>>
        (MessageType::ENCRYPTED_XY_PLUS_Z_VALUE);
    if (hasError(maybe_xyz_sum)) {
        LOG(ERROR) << "Kill/Restart received (xyz_sum)";
        return getError(maybe_xyz_sum);
    }
    xyz_sum() = getResult(maybe_xyz_sum);
    DBG("Participant (" << transport.getSocketId() <<
        ") received `Enc(x)*y - Enc(z)`");

    if (special) {
        auto [ pxyzs, r ] = e.partial_decrypt(xyz_sum(), secretKey);
        partial_xyz_shares() = pxyzs;
        g_r() = r;
    } else {
        auto [ pxyzs, r ] = e.partial_decrypt({xyz_sum().first, 0}, secretKey);
        partial_xyz_shares() = pxyzs;
        g_r() = r;
    }

    transport.send(MessageType::PARTIAL_XY_MINUS_Z_SHARES, partial_xyz_shares());
    DBG("Participant (" << transport.getSocketId() <<
        ") sending `Dec_si(xyz)` shares");


    // Now the coordinator can evaluate which of the `a_shares` we sent are valid. Many of them
    // will be eliminated, and the coordinator will communicate which values to eliminate via
    // a binary matrix, which we then use to shuffle our shares around before the CRT reconstruction step.
    auto maybe_sievingFlags = transport.awaitReply<std::vector<int>>
                              (MessageType::PS_SIEVING_FLAGS);
    if (hasError(maybe_sievingFlags)) {
        LOG(ERROR) << "Kill/Received (sievingFlags)";
        return getError(maybe_sievingFlags);
    }
    sievingFlags() = getResult(maybe_sievingFlags);
    std::vector<std::vector<mpz_class>> valid_shares = pruneAndReorderShares(
                                         a_shares(), b_shares(), sievingFlags(), alphasPS.size(), bucketSizePS);

    // Finally we can use CRT reconstruction across each row of the surviving a_shares to
    // construct our p_i,q_i shares for the candidate generation phase, taking the opportunity
    // to satisfy the 3 (mod 4) constraint at the same time.
    std::vector<mpz_class> prime_shares;

    // To do the 3 (mod 4) piece, we add an extra bucket to the alphas vector with a value of
    // four. The P1 party will use a corresponding coefficient of 3 while the others use 0.
    std::vector<mpz_class> _alphas(alphasPS.size() + 1);
    std::vector<mpz_class> _coeffs(_alphas.size());
    for (size_t i = 0; i < alphasPS.size(); ++i) {
        _alphas[i] = alphasPS[i];
    }
    _alphas[_alphas.size() - 1] = mpz_class(4);

    // Compute min row size

    size_t minRowSize = valid_shares[0].size();
    DBG("alphas[" << 0 << "] " << alphasPS[0]);
    DBG("valid_shares.size()= " << valid_shares[0].size());
    for (size_t i = 1; i < valid_shares.size(); ++i) {
        DBG("alphas[" << i << "] " << alphasPS[i]);
        DBG("valid_shares.size()= " << valid_shares[i].size());
        minRowSize = std::min(valid_shares[i].size(), minRowSize);
    }

    moduli().resize(_alphas.size());

    // NOTE: Here we prune x_index PS
    for (size_t i = 0; i < minRowSize; ++i) {
        for (size_t j = 0; j < alphasPS.size(); ++j) {
            _coeffs[j] = valid_shares[j][i];
        }
        _coeffs[_alphas.size() - 1] = (special ? mpz_class(3) : mpz_class(0));

        coefs().push_back(_coeffs);
        prime_shares.push_back(math::crt_reconstruct(_coeffs, moduli(), _alphas));
    }

    DBG("prime_shares.size() = " << prime_shares.size());
    return std::make_tuple(xcan, ycan, zcan, xgcd, ygcd, zgcd, prime_shares);
}

/**
 * Generates modulus candidates.
 * @param transport
 * @param xcan the part of x array from the triples dedicated to the candidate generation
 * @param ycan the part of x array from the triples dedicated to the candidate generation
 * @param zcan the part of x array from the triples dedicated to the candidate generation
 * @param both_shares candidate shares from which we derive both p_i and q_i
 * @param special flag for a special party
 * @param config protocol configuration
 * @return a tuple of three arrays: candidates, p_shares, q_shares
 */
template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
expected<std::tuple<std::vector<mpz_class>, std::vector<mpz_class>, std::vector<mpz_class>>>
EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::performModulusGeneration(
    ZeroMQClientTransport &transport,
    const std::vector<mpz_class> xcan,
    const std::vector<mpz_class> ycan,
    const std::vector<mpz_class> zcan,
    const std::vector<mpz_class> &both_shares,
    bool special,
    ProtocolConfig<T> &config
) {
    auto [alphasCAN, _bsz] = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());

    // set equal buckets for bucketSizeCan and bucketSizeGCD
    const auto bucketSize = lrint(floor(double(primesCAN) / double(
                                            alphasCAN.size())));
    const auto pq_size = std::min(bucketSize,
                                  lrint(floor(double(both_shares.size()) / double(2))));

    DBG("bucketSize = " << bucketSize);
    DBG("pq_size = " << pq_size);
    // Step one: generate some uniformly random numbers contributing
    // to the generation of the prime candidates.
    std::vector<mpz_class> p_shares(bucketSize);
    std::vector<mpz_class> q_shares(bucketSize);

    {
        for (int i = 0; i < pq_size; ++i) {
            p_shares[i] = both_shares[2 * i];
            q_shares[i] = both_shares[2 * i + 1];
        }
    }

    DBG("Doing CRT");
    std::vector<mpz_class> ax_shares(p_shares.size() * alphasCAN.size());
    std::vector<mpz_class> by_shares(p_shares.size() * alphasCAN.size());

    p_sharesCRTs().clear();
    q_sharesCRTs().clear();

    size_t imi = 0;
    for (int i = 0; i < p_shares.size(); ++i) {
        auto pCRT = math::crt_deconstruct(p_shares[i], alphasCAN);
        auto qCRT = math::crt_deconstruct(q_shares[i], alphasCAN);
        for (int j = 0; j < alphasCAN.size(); ++j) {
            p_sharesCRTs().push_back(pCRT[j]);
            q_sharesCRTs().push_back(qCRT[j]);
            ++imi;
        }
    }

    DBG("Done CRT");

    // compute ax/by
    for (int i = 0; i < ax_shares.size(); ++i) {
        ax_shares[i] = p_sharesCRTs()[i] - xcan[i];
        by_shares[i] = q_sharesCRTs()[i] - ycan[i];
    }

    {
        int k = 0;
        for (int i = 0; i < bucketSize; ++i) {
            for (int j = 0; j < alphasCAN.size(); ++j) {
                ax_shares[k] = math::mod(ax_shares[k], alphasCAN[j]);
                by_shares[k] = math::mod(by_shares[k], alphasCAN[j]);
                ++k;
            }
        }
    }

    DBG("Participant (" << transport.getSocketId() <<
        ") sending `ax and by` shares");
    transport.send(MessageType::AX_BY_SHARES, std::pair{ax_shares, by_shares});

    auto ax_by =
        transport.awaitReply<std::pair<std::vector<mpz_class>, std::vector<mpz_class>>>
        (MessageType::AX_BY_VALUE);
    if (hasError(ax_by)) {
        LOG(ERROR) << "Kill/Restart received (modulus generation, ax_by";
        return getError(ax_by);
    }
    auto [tax, tby] = getResult(ax_by);
    ax() = tax;
    by() = tby;
    DBG("Participant (" << transport.getSocketId() <<
        ") received `ax and by` values");

    // candidate generation
    axby().resize(ax().size());
    {
        int k = 0;
        for (size_t i = 0; i < bucketSize; ++i) {
            for (size_t j = 0; j < alphasCAN.size(); ++j) {
                axby()[k] = math::mod((ax()[k] * q_sharesCRTs()[k] + by()[k] * p_sharesCRTs()[k]
                                       + zcan[k]), alphasCAN[j]);
                ++k;
            }
        }
    }

    transport.send(MessageType::AXB_MINUS_BYA_SHARES, axby());
    DBG("Participant (" << transport.getSocketId() <<
        ") received `axb - bya` shares");

    auto maybe_candidates = transport.awaitReply<std::vector<mpz_class>>
                            (MessageType::MODULUS_CANDIDATE);
    if (hasError(maybe_candidates)) {
        LOG(ERROR) << "Kill/Restart received (modulus generation, candidates";
        return getError(maybe_candidates);
    }
    auto candidates = getResult(maybe_candidates);
    DBG("Participant (" << transport.getSocketId() << ") received `N` candidates");
    // In return the coordinator will send us a plaintext vector of modulus candidates.

    candidates.resize(pq_size);
    p_shares.resize(pq_size);
    q_shares.resize(pq_size);

    std::vector<int> syncF;
    syncF.push_back(1);
    transport.send(MessageType::SYNCHRONIZE_NOW, int(1));

    candidatesCAN() = candidates;
    return std::make_tuple(candidates, p_shares, q_shares);
}


/** Runs Jacobi biprimality test to prune candidates
 * @param transport communication helper
 * @param candidates array is a multiplication of sum(p_i) * sum(q_i) where i runs for all parties
 * @param p_shares party's private array of p_shares
 * @param q_shares party's private array of q_shares
 * @param special flag for a special party
 * @return 0 and success or error otherwise
 */
template<typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
expected<int>
EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::performJacobiTest(
    ZeroMQClientTransport &transport,
    const std::vector<mpz_class> &candidates,
    const std::vector<mpz_class> &p_shares,
    const std::vector<mpz_class> &q_shares,
    bool special
) {

    DBG("Going to send gamma seed share: " << jacobiSeedShares() << " for participant " << transport.getSocketId());
    transport.send(MessageType::GAMMA_RANDOM_SEED_SHARES, jacobiSeedShares());
    DBG("Participant (" << transport.getSocketId() << ") sending gamma shares");

    auto maybe_gammaSeed = transport.awaitReply<mpz_class>
                           (MessageType::GAMMA_RANDOM_SEED_VALUE);
    if (hasError(maybe_gammaSeed)) {
        LOG(ERROR) << "Kill/Restart received (performJacobiTest, gammaSeed";
        return getError(maybe_gammaSeed);
    }
    mpz_class gammaSeed = getResult(maybe_gammaSeed);
    // We find a set of valid gamma values by their jacobi symbol...
    std::vector<mpz_class> gammaValues(candidates.size());
    DBG("gammaSeed = " << gammaSeed);
    std::vector<mpz_class> engine = math::generateRandomVector(gammaSeed, kJacobiNumberOfRandomValues, 2048);
    size_t engineCount = 0;

    {
        mpz_t one;
        mpz_init(one);
        mpz_set_ui(one, 1);

        mpz_t gcdResult;
        mpz_init(gcdResult);

        for (int i = 0; i < candidates.size(); ++i) {
            const mpz_class &N = candidates[i];

            DBG(" N % 2 = 1; N = " << N);
            // We expect that N is 3 (mod 4), and therefore must be odd.
            assert(N % 2 == 1);

            while (1) {
                mpz_class g = engine[engineCount];
                ++engineCount;
                if (engineCount >= engine.size()) {
                    engineCount = 0;
                    DBG("Regenerating random values for Jacobi gamma values");
                    engine = math::generateRandomVector(engine[0], kJacobiNumberOfRandomValues,2048);
                }

                g = g % N;
                assert(g != 0);

                DBG("g = " << g);
                if (mpz_jacobi(g.get_mpz_t(), N.get_mpz_t()) == 1) {
                    gammaValues[i] = g;
                    break;
                }
            }

            // Sanity check, assert if gcd(gammaValues[i], N) == 1
            mpz_gcd(gcdResult, gammaValues[i].get_mpz_t(), N.get_mpz_t());
            assert(mpz_cmp(gcdResult, one) == 0);

            // Export 
        }

        mpz_clear(one);
        mpz_clear(gcdResult);
    }

    for (int i = 0; i < gammaValues.size(); ++i) {
        const mpz_class &N = candidates[i];
        mpz_class &g = gammaValues[i];

        mpz_class exp = (special ? mpz_class(N + 1) : mpz_class(
                             0)) - p_shares[i] - q_shares[i];

        if (special) {
            if (p_shares[i] % 4 != 3 || q_shares[i] % 4 != 3) {
                DBG("special");
                DBG("i = " << i);
                DBG("p_shares[i] " << p_shares[i]);
                DBG("p_shares[i] " << q_shares[i]);
                assert(false);
            }
        } else {
            if (p_shares[i] % 4 != 0 || q_shares[i] % 4 != 0) {
                DBG("ordinary");
                DBG("i = " << i);
                DBG("p_shares[i] " << p_shares[i]);
                DBG("p_shares[i] " << q_shares[i]);
                assert(false);
            }
        }

        // Sanity check before the division operator
        assert(exp % 4 == 0);
        exp /= 4;
        g = math::powm(g, exp, N);
    }

    // Finally we send our gamma values back and this concludes the test
    // from the client perspective. We will repeat the same test if the coordinator
    // asks.
    transport.send(MessageType::EXPONENTIATED_GAMMA_VALUE, gammaValues);
    DBG("Participant (" << transport.getSocketId() <<
        ") sending exponentiated gamma shares");
    return 0;
}


/** Runs GCD and Jacobi biprimality tests to prune candidates
 * @param transport communication helper
 * @param candidates array is a multiplication of sum(p_i) * sum(q_i) where i runs for all parties
 * @param p_shares party's private array of p_shares
 * @param q_shares party's private array of q_shares
 * @param x the part of x array from the triples dedicated to the gcd rounds
 * @param y the part of y array from the triples dedicated to the gcd rounds
 * @param z the part of z array from the triples dedicated to the gcd rounds
 * @param special flag for a special party
 * @param config protocol configuration
 * @return 0 on success and error otherwise
 */
template<typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
expected<int>
EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::performGCDandJacobiTest(
    ZeroMQClientTransport &transport,
    const std::vector<mpz_class> &candidates,
    const std::vector<mpz_class> &p_shares,
    const std::vector<mpz_class> &q_shares,
    const std::vector<mpz_class> &x,
    const std::vector<mpz_class> &y,
    const std::vector<mpz_class> &z,
    bool special,
    ProtocolConfig<T> &config
) {

    DBG("performGCDTest candidates.size() = " << candidates.size());
    assert(candidates.size() == p_shares.size());
    assert(candidates.size() == q_shares.size());

    // The first step of the GCD test is to generate shares of a random `r` value,
    // use CRT to deconstruct this value, encrypt deconstructed representation
    // and send them to the coordinator
    DBG("Initialize alphasGCD");
    auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3*config.pbs() + 210, primesGCD, config.tauLimitBit());
    const int bucketSize = lrint(floor(double(primesGCD) / double(alphasGCD.size())));

    DBG("bucketSize = " << bucketSize);
    DBG("alphasGCD.size() = " << alphasGCD.size());

    DBG("Initialize ax_shares and bx_shares");
    std::vector<mpz_class> ax_shares(primesGCD);
    std::vector<mpz_class> by_shares(primesGCD);

    DBG("p_shares.size() = " << p_shares.size());
    assert(x.size() >= (p_shares.size() * alphasGCD.size()));
    assert(y.size() >= (p_shares.size() * alphasGCD.size()));
    assert(z.size() >= (p_shares.size() * alphasGCD.size()));

    DBG("Initialize random");
    gmp_randclass grc(gmp_randinit_default);
    grc.seed(std::random_device()());

    DBG("Initialize a and b");
    // a crts
    rCRTs().resize(p_shares.size() * alphasGCD.size());

    // b crts
    std::vector<mpz_class> p_plus_qCRTs(p_shares.size() * alphasGCD.size());

    const mpz_class maybeOne = special ? 1 : 0;
    {
        int rc = 0;
        for (int i = 0; i < candidates.size(); ++i) {

            // a = random number
            DBG("a = random ");
            gcdRX() = grc.get_z_bits(2 * config.pbs() + 48); 
            DBG("doing CRT deconstruct ");
            auto aCRT = math::crt_deconstruct(gcdRX(), alphasGCD);
            assert(aCRT.size() == alphasGCD.size());

            DBG("b = p_i ");
            // b = p_i + q_i - maybeOne
            auto bCRT = math::crt_deconstruct(p_shares[i] + q_shares[i] - maybeOne,
                                              alphasGCD);
            assert(bCRT.size() == alphasGCD.size());

            DBG("copy numbers");
            for (int j = 0; j < aCRT.size(); ++j) {
                DBG("ax");
                rCRTs()[rc] = aCRT[j];
                ax_shares[rc] = rCRTs()[rc] - x[rc];

                DBG("by 1");
                p_plus_qCRTs[rc] = bCRT[j];
                DBG("by 2");
                by_shares[rc] = p_plus_qCRTs[rc] - y[rc];
                DBG("by 3");
                ++rc;
            }
        }
    }

    DBG("Initialize ax and bx");
    {
        DBG("ax_shares.size()" << ax_shares.size());
        int k = 0;
        for (int i = 0; i < bucketSize; ++i) {
            for (int j = 0; j < alphasGCD.size(); ++j) {
                DBG("k " << k);
                ax_shares[k] = math::mod(ax_shares[k], alphasGCD[j]);
                by_shares[k] = math::mod(by_shares[k], alphasGCD[j]);
                ++k;
            }
        }
    }

    // Jacobi generates gamma seed share

    DBG("Participant (" << transport.getSocketId() << ") sending `ax and by` shares");
    transport.send(MessageType::GCD_AX_BY_SHARES, std::pair{std::pair{ax_shares, by_shares}, jacobiAndGCDSeedShares()});

    auto maybe_axby =
        transport.awaitReply<std::pair<std::pair<std::vector<mpz_class>,
        std::vector<mpz_class>>,
        mpz_class>>(MessageType::AX_BY_VALUE);
    if (hasError(maybe_axby)) {
        LOG(ERROR) << "Kill/Restarted received (GCDandJacobiTest, axby)";
        return getError(maybe_axby);
    }
    auto [axby_pair, gammaSeed] = getResult(maybe_axby);
    
    gammaSeed_g() = gammaSeed;
    auto [ax_GCD, by_GCD] = axby_pair;
    axGCD() = ax_GCD;
    byGCD() = by_GCD;
    DBG("Participant (" << transport.getSocketId() <<
        ") received `ax and by` values");


    // Jacobi: We find a set of valid gamma values by their jacobi symbol...
    std::vector<mpz_class> gammaValues(config.lambda() * candidates.size());

    DBG("gammaSeed = " << gammaSeed);
    std::vector<mpz_class> engine = math::generateRandomVector(gammaSeed,
                                    kJacobiNumberOfRandomValues, 2048);
    size_t engineCount = 0;

    {
        mpz_t one;
        mpz_init(one);
        mpz_set_ui(one, 1);

        mpz_t gcdResult;
        mpz_init(gcdResult);

        for (int i = 0; i < candidates.size(); ++i) {
            const mpz_class &N = candidates[i];

            DBG(" N % 2 = 1; N = " << N);
            // We expect that N is 3 (mod 4), and therefore must be odd.
            assert(N % 2 == 1);

            int foundGamma = 0;
            while (foundGamma < config.lambda()) {
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
                    gammaValues[i * config.lambda() + foundGamma] = g;
                    foundGamma++;
                }
            }

            // Sanity check, assert if gcd(gammaValues[i], N) == 1
            for (int j = 0; j < config.lambda(); ++j) {
                mpz_gcd(gcdResult, gammaValues[i * config.lambda() + j].get_mpz_t(),
                        N.get_mpz_t());
                assert(mpz_cmp(gcdResult, one) == 0);
            }
        }

        mpz_clear(one);
        mpz_clear(gcdResult);
    }

    // copy and store g before exp
    gBeforeExp() = gammaValues;
    gAfterExp().resize(gammaValues.size());

    for (int i = 0; i < candidates.size(); ++i) {
        for (int j = 0; j < config.lambda(); j++) {
            const mpz_class &N = candidates[i];
            mpz_class &g = gammaValues[i * config.lambda() + j];

            mpz_class exp = (special ? mpz_class(N + 1) : mpz_class(0)) - p_shares[i] - q_shares[i];

            if (special) {
                if (p_shares[i] % 4 != 3 || q_shares[i] % 4 != 3) {
                    DBG("special");
                    DBG("i = " << i);
                    DBG("p_shares[i] " << p_shares[i]);
                    DBG("p_shares[i] " << q_shares[i]);
                    assert(false);
                }
            } else {
                if (p_shares[i] % 4 != 0 || q_shares[i] % 4 != 0) {
                    DBG("ordinary");
                    DBG("i = " << i);
                    DBG("p_shares[i] " << p_shares[i]);
                    DBG("p_shares[i] " << q_shares[i]);
                    assert(false);
                }
            }

            // Sanity check before the division operator
            assert(exp % 4 == 0);
            exp /= 4;
            mpz_class orig_g = mpz_class(g);
            g = math::powm(g, exp, N);

            gAfterExp()[i * config.lambda() + j] = mpz_class(g);

        }
    }


    // GCD: Compute  axby
    axbyGCD().resize(axGCD().size());
    {
        int k = 0;
        ss().resize(bucketSize);
        ssGCD().resize(ss().size());
        LOG(INFO) << "candidates.size() " << candidates.size();
        LOG(INFO) << "bucketSize " << bucketSize;
        for (size_t i = 0; i < bucketSize; ++i) {
            ss()[i] = grc.get_z_bits(1024 + 210);
            ssGCD()[i] = ss()[i] * candidates[i];
            for (size_t j = 0; j < alphasGCD.size(); ++j) {
                axbyGCD()[k] = math::mod((axGCD()[k] * p_plus_qCRTs[k] + byGCD()[k] * rCRTs()[k]
                                          + z[k] + ssGCD()[i]), alphasGCD[j]);
                ++k;
            }
        }
    }

    for (int i = 0; i < gammaValues.size(); ++i) {
        if (gammaValues[i] == 0) {
            DBG("Found 0 gamma value at i = " << i);
            assert(false);
        }
    }
    transport.send(MessageType::AXB_MINUS_BYA_SHARES, std::pair{axbyGCD(), gammaValues});
    DBG("Participant (" << transport.getSocketId() <<
        ") received `axb - bya` shares");
    return 0;
}


} // namespace ligero
