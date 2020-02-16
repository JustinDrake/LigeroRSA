#pragma once

#include <vector>

#include "Transport.hpp"
#include "Math.hpp"
#include "Common.hpp"
#include "LatticeEncryption.hpp"
#include <stdlib.h>
#include "FiniteFields.hpp"

uint64_t ppp = 4611686018326724609;

namespace ligero {

template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
class EncryptedClient {
    public:
        EncryptedClient(SocketId& sid) : _socketId(sid) {};
        using Q = nfl::poly_p<T, Degree, NbPrimesQ>;
        using P = nfl::poly_p<T, Degree, NbPrimesP>;
        using cipherText = typename
                           lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>::cipherText;

        CandidateIndices getCandidateIndices(size_t idx);

        std::optional<lattice::key_pair<T, Degree, NbPrimesQ>> generateKeyPair(
                    ZeroMQClientTransport &trans, nfl::gaussian<uint16_t, T, 2> &chi);
        std::vector<std::vector<mpz_class>> pruneAndReorderShares(
                                             const std::vector<mpz_class> &as,
                                             const std::vector<mpz_class> &bs,
                                             const std::vector<int> &flags,
                                             int numAlphas,
                                             std::vector<size_t> bucketSize
                                         );

        ProtocolConfig<T> registerAndAwaitConfiguration(ZeroMQClientTransport &transport, std::string ipAddress);

        std::optional<std::tuple<
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

        std::optional<std::tuple<std::vector<mpz_class>, std::vector<mpz_class>, std::vector<mpz_class>>>
        performModulusGeneration(
            ZeroMQClientTransport &transport,
            const std::vector<mpz_class> xcan,
            const std::vector<mpz_class> ycan,
            const std::vector<mpz_class> zcan,
            const std::vector<mpz_class> &both_shares,
            bool special,
            ProtocolConfig<T> &config
        );

        std::optional<int> performJacobiTest(
            ZeroMQClientTransport &transport,
            const std::vector<mpz_class> &candidates,
            const std::vector<mpz_class> &p_shares,
            const std::vector<mpz_class> &q_shares,
            bool special
        );

        std::optional<int> performGCDandJacobiTest(
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

        std::vector<mpz_class> &moduli() {
            return _moduli;
        }
        const std::vector<mpz_class> &moduli() const {
            return _moduli;
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

        /* ========== ========== */
        void start(ZeroMQClientTransport& transport, size_t wait_timeout, size_t test_timeout,
                   size_t data_size, size_t nb_max_send, const std::string &localUri = "n/a") {

            auto proceed = transport.joinThroughputTest(wait_timeout, test_timeout,
                           data_size, nb_max_send);

            if (!proceed) {
                LOG(INFO) << "Kicked out due to poor network connection";
                return;
            }

            auto config = registerAndAwaitConfiguration(transport, localUri);
            auto e = lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>(config);

            std::optional<int> success;
            try {
                success = start_rsa_ceremony(transport, config, e);
            } catch (const std::runtime_error &e) {
                LOG(ERROR) << e.what();
                LOG(INFO) << "Killed by coordinator, aborting";
                exit(EXIT_FAILURE);
            }
            while (!success) {
                //LOG(ERROR) <<
                //           "First run terminated by coordinator, restarting and try again...";
                LOG(INFO) << "Run failed, restarting";
                // cleanup and restart
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

		
                LOG(INFO) << "========== starting next run ==========";
                try {
                    success = start_rsa_ceremony(transport, config, e);
                } catch (const std::runtime_error &e) {
                    LOG(ERROR) << e.what();
                    LOG(INFO) << "Killed by coordinator, aborting";
                    exit(EXIT_FAILURE);
                }
                //if (!success) {
                //    LOG(INFO) << "Run failed, restarting";
                    //LOG(FATAL) << "Second run failed, aborting";
                //}
            }
        }

        std::optional<int> start_rsa_ceremony(
            ZeroMQClientTransport &transport,
            ProtocolConfig<T> &config,
            lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ> &e
        ) {
            auto maybe_keys = generateKeyPair(transport, e.chi());
            if (!maybe_keys) {
                LOG(ERROR) << "Kill/Restart received when keygen";
                return std::nullopt;
            }
            auto [publicKey, secretKey] = *maybe_keys;

            // Await designation
            DBG("Awaiting designation.");

            auto maybe_assignment = transport.awaitReply();
            if (!maybe_assignment) {
                LOG(ERROR) << "Kill/Restart when wait assignment";
                return std::nullopt;
            }
            MessageType assignment = *maybe_assignment;
            bool special = assignment == MessageType::ASSIGNMENT_P1;
            speciald() = special;
            DBG("Participant (" << socketId() << ") got my designation.");

            LOG(INFO) << "Connected to coordinator.";

            bool foundModuli = false;
            int numberOfTrials = 0;
            while (!foundModuli && (numberOfTrials < kTrialsThreshold)) {
                numberOfTrials++;

                LOG(INFO) << "Generating shares for pre-sieving.";
                auto maybe_shares = generatePreSievingShares(e, transport, publicKey, secretKey,
                                    special, config);
                if (!maybe_shares) {
                    LOG(ERROR) << "Kill/Restart received when generate presieve shares";
                    return std::nullopt;
                }
                auto [xcan, ycan, zcan, xgcd, ygcd, zgcd, both_shares] = *maybe_shares;

                LOG(INFO) << "Generated shares for " << both_shares.size() / 2 <<
                          " candidates.";

                LOG(INFO) << "Using pre-sieved shares for candidate generation.";
                auto maybe_modulus = performModulusGeneration(transport, xcan, ycan, zcan,
                                     both_shares, special, config);
                if (!maybe_modulus) {
                    LOG(ERROR) << "Kill/Restart received when perform modulus generation";
                    return std::nullopt;
                }
                auto [candidates, p_shares, q_shares] = *maybe_modulus;

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
                if (!maybe_discardFlags) {
                    LOG(ERROR) << "Kill/Restart when wait discard flags";
                    return std::nullopt;
                }
                boost::dynamic_bitset<> discardFlags = *maybe_discardFlags;
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
                    if (!success) {
                        LOG(ERROR) << "Kill/Restart received";
                        return std::nullopt;
                    }
                    switch (*success) {
                        case MessageType::GAMMA_SHARES: {
                            if (!ranJacobi) {
                                LOG(INFO) << "Running Jacobi test on candidates.";
                                ranJacobi = true;
                            }
                            auto maybe_JacobiResult = performJacobiTest(transport, candidates, p_shares, q_shares, special);
                            if (!maybe_JacobiResult) {
                                return std::nullopt;
                            }

                            auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>(MessageType::DISCARD_FLAGS);
                            if (!maybe_discardFlags) {
                                LOG(ERROR) << "Kill/Restart received during wait discard flags";
                                return std::nullopt;
                            }
                            boost::dynamic_bitset<> discardFlags = *maybe_discardFlags;

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

                            if (!maybe_JacobiGCDResult) {
                                return std::nullopt;
                            }

                            auto maybe_discardFlags = transport.awaitReply<boost::dynamic_bitset<>>(MessageType::DISCARD_FLAGS);
                            if (!maybe_discardFlags) {
                                LOG(ERROR) << "Kill/Restart received during wait discard flags";
                                return std::nullopt;
                            }
                            boost::dynamic_bitset<> discardFlags = *maybe_discardFlags;

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
                            break;
                        }
                        default:
                            throw std::runtime_error("Received a message type out of order.");
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

                    candidateIndices() = getCandidateIndices(idx);
                    final_index() = idx;
                    finalModuli() = candidates[0];
                    gammaExponent() = ((speciald() ? mpz_class(candidates[0]+1) : mpz_class(0)) - p_shares[0] - q_shares[0])/4;

                    //std::cout<< speciald() << ": p_shares[0] = " << p_shares[0] << std::endl;

                    // The following code transmits the p shares and q shares in the clear to the coordinator.
                    // This is ONLY for checking the result. We need to remove this part finally.
                    DBG("Participant (" << socketId() <<
                        ") Transmitting the p and q values in the clear for debugging purposes only.");

                    transport.send(MessageType::MUTHU_ACK, int(1));
                    DBG("Participant (" << socketId() << ") done.");

                }
            }

            return 0;
        }

        /* ========== public data, will eventually be removed ========== */
        /*
         * Gathers and returns all public data.
         */
        void gatherData(PublicData &d, SecretData &e) {

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

            //y_prime.invntt_pow_invphi();

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
            auto [alphasPS, bucketSizePS] = math::balanced_bucket_n_primes(1000, primesPS,
                                            175, 1);
            mpz_class alphasPSprod = mpz_class(4);
            for (size_t j = 0; j < alphasPS.size(); ++j) {
                alphasPSprod *= alphasPS[j];
            }
            mpz_class recomputed_pshare;
            mpz_class recomputed_qshare;

            // Round 7
            // and Round 8
            // axby[k] = math::mod((ax()[k] * q_sharesCRTs()[k] + by()[k] * p_sharesCRTs()[k] + zcan[k]), alphasCAN[j]);
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
                    recomputed_pshare = modifiedProduct(moduli(), coefIndices, x_shares(),
                                                        candidateIndices().ps);
                    mpz_class pshare_quotient = (recomputed_pshare / alphasPSprod) % alphasCAN[alphaCANidx];

                    mpz_class recomputed_pshare_CAN = modifiedProduct(moduliCAN, coefIndices,
                                                      x_shares(), candidateIndices().ps);
                    mpz_class axshares = math::mod(p_sharesCRTs()[final_index() * alphasCAN.size() +
                                                        alphaCANidx] - x_shares()[j], alphasCAN[alphaCANidx]);
                    mpz_class pdividend = axshares + x_shares()[j] - (recomputed_pshare_CAN -
                                          pshare_quotient * alphasPSprodCAN);
                    mpz_class p_shareCRT_quotient  = (pdividend / alphasCAN[alphaCANidx]);
                    assert(pdividend - p_shareCRT_quotient * alphasCAN[alphaCANidx] == 0);

                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        q_p_prod_r7.push_back(positiveRemainder(pshare_quotient, prime));
                        q_p_r7.push_back(positiveRemainder(p_shareCRT_quotient, prime));
                        ax_shares.push_back(positiveRemainder(axshares, prime));
                    }

                    //pushing qshare items
                    recomputed_qshare = modifiedProduct(moduli(), coefIndices, y_shares(),
                                                        candidateIndices().ps);
                    mpz_class qshare_quotient = (recomputed_qshare / alphasPSprod) % alphasCAN[alphaCANidx];

                    mpz_class recomputed_qshare_CAN = modifiedProduct(moduliCAN, coefIndices,
                                                      y_shares(), candidateIndices().ps);
                    mpz_class byshares = math::mod(q_sharesCRTs()[final_index() * alphasCAN.size() +
                                                        alphaCANidx] - y_shares()[j], alphasCAN[alphaCANidx]);
                    mpz_class qdividend = byshares + y_shares()[j] - (recomputed_qshare_CAN -
                                          qshare_quotient * alphasPSprodCAN);
                    mpz_class q_shareCRT_quotient  = (qdividend / alphasCAN[alphaCANidx]);
                    assert(qdividend - q_shareCRT_quotient * alphasCAN[alphaCANidx] == 0);

                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        q_q_prod_r7.push_back(positiveRemainder(qshare_quotient, prime));
                        q_q_r7.push_back(positiveRemainder(q_shareCRT_quotient, prime));
                        by_shares.push_back(positiveRemainder(byshares,
                                                              prime)); //mpz_class(nfl::params<uint64_t>::P[i])));
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

            auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3096, primesGCD,175);
            std::vector<mpz_class> alphasPSprod_GCD, moduli_GCD;
            std::vector<mpz_class> ax_shares_GCD;
            std::vector<mpz_class> by_shares_GCD;

            std::vector<mpz_class> expsharesGCD;
            // Round 11
            // ax_shares[rc] = rCRTs()[rc] - x[rc];
            // auto bCRT = math::crt_deconstruct(p_shares[i] + q_shares[i] - maybeOne, alphasGCD);
            // by_shares[rc] = p_plus_qCRTs[rc] - y[rc];
            // p_plus_qCRTs[rc] = bCRT[j];
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

		
										//axbyGCD()[k] = math::mod((axGCD()[k] * p_plus_qCRTs[k] + byGCD()[k] * rCRTs()[k]
                    //                      + z[k]), alphasGCD[j]);

                    mpz_class dividend = axGCD()[alphaGCDidx] * ( recomputed_pshare_GCD + recomputed_qshare_GCD - maybeOne - qshare_quotient * alphasPSprodGCD - pshare_quotient * alphasPSprodGCD ) 
																		+ byGCD()[alphaGCDidx] * (rCRTs()[alphaGCDidx]) + z_shares()[j] - axbyGCD()[alphaGCDidx];
										assert(dividend % alphasGCD[alphaGCDidx] == 0);
                    mpz_class quotient = dividend / alphasGCD[alphaGCDidx];
                    assert(dividend - quotient * alphasGCD[alphaGCDidx] == 0);
                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        q_r12.push_back(positiveRemainder(quotient, prime));
                    }

                    alphaGCDidx++;
                }
            }



        // Sigma Protocol
        //
        // Step 1: Commit to r and x before executing the sigma protocol
         std::vector<mpz_class> sigmas_r(128);
         std::vector<mpz_class> sigma_r_GCD(129 * NbPrimesP * alphasGCD.size()),
         sigma_x_GCD(NbPrimesP * alphasGCD.size()), exp_q_GCD (129 * NbPrimesP * alphasGCD.size());;
         for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
            mpz_class exponentGCD = math::mod(gammaExponent(), alphasGCD[alphaGCDidx]);
            
            mpz_class exponentGCDquotient = (4 * exponentGCD - ((speciald() ? math::mod(finalModuli(),
                                                                          alphasGCD[alphaGCDidx]) + 1 : mpz_class(0)) - expsharesGCD[alphaGCDidx]));
            assert(math::mod(exponentGCDquotient, alphasGCD[alphaGCDidx]) == 0);
            exponentGCDquotient /= alphasGCD[alphaGCDidx];
             for (size_t i = 0; i < NbPrimesP; i ++) {
                mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                sigma_x_GCD[i * alphasGCD.size() + alphaGCDidx] = positiveRemainder(exponentGCD, prime);
                exp_q_GCD[i * alphasGCD.size() + alphaGCDidx] 
												                      = positiveRemainder(exponentGCDquotient, prime);
            }
        }
        // assuming that second jacobi test only has one moduli and security parameter is 128.
        for (size_t j = 0; j < 128; j++) {
						mpz_class sigma_r = math::generateRandomValue(mpz_class(std::random_device()()), 2048 + 128);
						sigmas_r[j] = sigma_r;
            for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                mpz_class sigma_rGCD = math::mod(sigma_r, alphasGCD[alphaGCDidx]);
                for (size_t i = 0; i < NbPrimesP; i ++) {
                    mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                    sigma_r_GCD[i * alphasGCD.size() * 128 + j * alphasGCD.size() + alphaGCDidx] 
														                                    = positiveRemainder(sigma_rGCD, prime);
                }
                
            }
            
        }

        // Create Merkle tree for committing to a block that contain the rAndX = concatenation of exp_x_GCD and exp_r_GCD
        // Step 2: Executing the Sigma Protocol
        std::vector<mpz_class> sigmas_e(128), sigmas_z(128);
        std::vector<mpz_class> sigma_e_GCD(128 * alphasGCD.size() * NbPrimesP);
        std::vector<mpz_class> sigma_z_GCD(128 * alphasGCD.size() * NbPrimesP);
        std::vector<mpz_class> sigma_q_GCD(128 * alphasGCD.size() * NbPrimesP);
        std::vector<mpz_class> sigmas_a(128);
        // First completing the sigma ZK protocol and then verifying the main ZK proof
        for (size_t j = 0 ; j < 128; j++) {

            //mpz_class sigma_r = math::generateRandomValue(mpz_class(std::random_device()()), 2 * config.pbs() + 48 + config.lambda());
            mpz_class sigma_a = math::powm(gBeforeExp()[j], sigmas_r[j], finalModuli());
            sigmas_a[j] = sigma_a;
            sigmas_e[j] = math::generateRandomValue(mpz_class(std::random_device()()),
                    128); //obtain from Fiat-Shamir
            sigmas_z[j] = sigmas_r[j] + sigmas_e[j] * gammaExponent();
            for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                mpz_class sigma_xGCD = math::mod(gammaExponent(), alphasGCD[alphaGCDidx]);
                mpz_class sigma_rGCD = math::mod(sigmas_r[j], alphasGCD[alphaGCDidx]);
                mpz_class sigma_eGCD = math::mod(sigmas_e[j], alphasGCD[alphaGCDidx]);
                mpz_class sigma_zGCD = math::mod(sigmas_z[j], alphasGCD[alphaGCDidx]);

                mpz_class dividend = sigma_zGCD - (sigma_rGCD + sigma_eGCD * sigma_xGCD);

                assert(dividend % alphasGCD[alphaGCDidx] == 0);
                mpz_class quotient = dividend / alphasGCD[alphaGCDidx];
                for (size_t i = 0; i < NbPrimesP; i++) {
                    mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                    sigma_e_GCD [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(sigma_eGCD,prime);
                    sigma_q_GCD [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(quotient,prime);
                    sigma_z_GCD [i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx] = positiveRemainder(sigma_zGCD,prime);
                }
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

            auto dataExtractQ = [&](std::vector<uint64_t> &assignee, Q & target) -> void {
                assignee.assign(
                    target.poly_obj().data(),
                    target.poly_obj().data() + (size_t)target.nmoduli * (size_t)target.degree);
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

            auto dataExtractVectorCRT = [&](std::vector<uint64_t> &assignee,
            std::vector<mpz_class> &target) -> void {
                assignee.resize(target.size()* NbPrimesQ);
                size_t ctr = 0;
                for (size_t idx = 0; idx < target.size(); idx++) {
                    for (size_t pidx = 0; pidx < NbPrimesQ; pidx++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[pidx]);
                        mpz_class u = positiveRemainder(target[idx], prime);
                        assignee[ctr++] = mpz_get_ui(u.get_mpz_t());
                    }
                }
            };

            auto dataExtractVector = [&](std::vector<uint64_t> &assignee,
            std::vector<mpz_class> &target) -> void {
                assignee.resize(target.size());
                for (size_t idx = 0; idx < target.size(); idx++) {
                    assignee[idx] = mpz_get_ui(target[idx].get_mpz_t());
                }
            };

            auto dataExtractSubVector = [&](std::vector<uint64_t> &assignee,
            std::vector<mpz_class> &target, std::vector<size_t> &idxs) -> void {
                assignee.resize(target.size()*NbPrimesQ);
                int ctr = 0;
                for (size_t idx = 0; idx < idxs.size(); idx++) {
                    for (size_t pidx = 0; pidx < NbPrimesQ; pidx++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[pidx]);
                        mpz_class u = positiveRemainder(target[idxs[idx]], prime);
                        assignee[ctr++] = mpz_get_ui(u.get_mpz_t());
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

            // Randomness attached to surviving candidates
            std::vector<size_t> indices(candidateIndices().can);
            indices.insert(indices.end(), candidateIndices().ps.begin(), candidateIndices().ps.end());
            indices.insert(indices.end(), candidateIndices().gcd.begin(), candidateIndices().gcd.end());

            // Alphas
            dataExtractVectorCRT(d.ps,       alphasPS);

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
            dataExtractVector(d.by_shares,      by_shares);
            dataExtractSubVector(e.y_sharesCAN, y_shares(), candidateIndices().can);
            dataExtractSubVector(e.y_sharesPS,  y_shares(), candidateIndices().ps);
            dataExtractVector(d.ax_shares,      ax_shares);
            dataExtractSubVector(e.x_sharesCAN, x_shares(), candidateIndices().can);
            dataExtractSubVector(e.x_sharesPS,  x_shares(), candidateIndices().ps);
            d.indicesPS.assign(candidateIndices().ps.begin(), candidateIndices().ps.end());
            d.indicesCAN.assign(candidateIndices().can.begin(),candidateIndices().can.end());
            d.indicesGCD.assign(candidateIndices().gcd.begin(),candidateIndices().gcd.end());
            dataExtractVectorCRT(d.coefsCAN,	moduli_CAN);
            dataExtractVectorCRT(d.cans,        alphasCAN);
            dataExtractVectorCRT(d.prodcans,    alphasPSprod_CAN);
            dataExtractVector(e.q_p_prod_r7,    q_p_prod_r7);
            dataExtractVector(e.q_p_r7,         q_p_r7);
            dataExtractVector(e.q_q_prod_r7,    q_q_prod_r7);
            dataExtractVector(e.q_q_r7,         q_q_r7);

            // Round8
            dataExtractSubVector(e.z_sharesCAN, z_shares(), candidateIndices().can);
            dataExtractSubVector(d.ax,          ax(), axbyindices);
            dataExtractSubVector(d.by,          by(), axbyindices);
            dataExtractSubVector(d.axby,        axby(), axbyindices);
            dataExtractVector(e.q_r8,           q_r8);

            // Bounding x_shares, y_shares and z_shares
            std::vector<mpz_class> vars;

            pickSubVector(vars, x_shares(), indices);
            pickSubVector(vars, y_shares(), indices);
            pickSubVector(vars, z_shares(), candidateIndices().can);

            std::vector<mpz_class> cans(alphasCAN);
            cans.insert(cans.end(), alphasPS.begin(), alphasPS.end());
            cans.insert(cans.end(), alphasGCD.begin(), alphasGCD.end());

            auto params = {cans, cans, alphasCAN};
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
            dataExtractVector(d.by_shares_GCD,  by_shares_GCD);
            dataExtractSubVector(e.y_sharesGCD, y_shares(), candidateIndices().gcd);
            dataExtractVector(d.ax_shares_GCD,  ax_shares_GCD);
            dataExtractSubVector(e.x_sharesGCD, x_shares(), candidateIndices().gcd);
            dataExtractVectorCRT(d.coefsGCD,    moduli_GCD);
            dataExtractVectorCRT(d.gcds,        alphasGCD);
            dataExtractVectorCRT(d.prodgcds,    alphasPSprod_GCD);
            dataExtractVector(e.q_p_prod_r11,   q_p_prod_r11);
            dataExtractVector(e.q_pq_r11,       q_pq_r11);
            dataExtractVector(e.q_q_prod_r11,   q_q_prod_r11);
            dataExtractVector(e.q_r_r11,        q_r_r11);
            dataExtractVector(e.r_CRTs,         r_CRTs);
            dataExtractVector(d.ax_shares_GCD,  ax_shares_GCD);

            dataExtractSubVector(e.z_sharesGCD, z_shares(), candidateIndices().gcd);
            dataExtractVectorCRT(d.axGCD,       axGCD());
            dataExtractVectorCRT(d.byGCD,       byGCD());
            dataExtractVectorCRT(d.axbyGCD,     axbyGCD());
            dataExtractVector(e.q_r12,          q_r12);
    	    // Sigma Protocol - this is not quite extraction. We are preparing variables to
        // perform the sigma protocol. In particular, we need to store the r and x in the
        // witness of the main ZK proof
        
        // e.sigmar = sigmas_r;
        // d.sigmaa = sigmas_a;
        // d.sigmae = sigmas_e;
        // d.sigmaz = sigmas_z;
        dataExtractVector(e.sigmarGCD,      sigma_r_GCD);
        dataExtractVector(e.sigmaxGCD,      sigma_x_GCD);
        dataExtractVector(d.sigmaeGCD,      sigma_e_GCD);
        dataExtractVector(d.sigmazGCD,      sigma_z_GCD);
        dataExtractVector(e.exp_q,          exp_q_GCD);
        dataExtractVector(e.sigmaqGCD,      sigma_q_GCD);
        //dataExtractVector(e.sigma_blind,    sigma_blind);
 
   
            // debug: test evaluation of x_shares on PS, CAN and GCD
            // x_prime.ntt_pow_phi();

            // std::cout << e.x_sharesPS[0] << std::endl;
            // std::cout << copy.poly_obj().data()[candidateIndices().ps[0]] << std::endl;
            // std::cout << x_prime.poly_obj().data()[candidateIndices().ps[0]] << std::endl;
            // std::cout << qdivp.poly_obj().data()[candidateIndices().ps[0]] << std::endl;
            // std::cout << nfl::params<uint64_t>::P[0] << std::endl;

            // for (size_t idx =0; idx < 6; idx++) {
            //     uint64_t mod = static_cast<__uint64_t>(static_cast<__uint128_t>(e.x_sharesPS[idx*NbPrimesQ])*static_cast<__uint128_t>(qdivp.poly_obj().data()[candidateIndices().ps[idx]]) %
            //     static_cast<__uint128_t>(nfl::params<uint64_t>::P[0]));

            //     std::cout << mod << std::endl;

            //     assert(mod == x_prime.poly_obj().data()[candidateIndices().ps[idx]]);
            // }
            // signal that we did gather all secret data
            e.isAvailable = true;

            // Checking extraction for Rounds 7 & 8
            //

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
                        //std::cout << speciald() << ": rem = " << rem % mpz_class(nfl::params<uint64_t>::P[p]) <<  std::endl;
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
                        //std::cout << speciald() << ": rem = " << rem % mpz_class(nfl::params<uint64_t>::P[p]) <<  std::endl;
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                    //Round 8
                    //Checking dividend - quotient * alphasCAN[alphaCANidx] == 0 in
                    //	 mpz_class dividend = (ax()[pos]*(recomputed_qshare_CAN - qshare_quotient*alphasPSprodCAN)+ by()[pos]*(recomputed_pshare_CAN - pshare_quotient*alphasPSprodCAN)+ z_shares()[j] - axby()[pos]);
                    //   mpz_class quotient = dividend/alphasCAN[alphaCANidx];
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
                        // std::cout << speciald() << ": rem = " << rem % mpz_class(nfl::params<uint64_t>::P[p]) <<  std::endl;
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                    // Round 11 and 12
                    // 
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
                        //std::cout << speciald() << ": rem = " << rem % mpz_class(nfl::params<uint64_t>::P[p]) <<  std::endl;
                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }
                    {
                        mpz_class rem = mpz_class(d.ax_shares_GCD[i + p]) 
                                        + mpz_class(e.x_sharesGCD[i + p])
                                        - mpz_class(e.r_CRTs[i + p])
                                        - mpz_class(e.q_r_r11[i + p]) * mpz_class(d.gcds[i + p]);
                        //std::cout << speciald() << ": rem = " << rem % mpz_class(nfl::params<uint64_t>::P[p]) <<  std::endl;
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
																				+ mpz_class(e.q_r12[i + p]) * mpz_class(d.gcds[i + p]);

                        assert(rem % mpz_class(nfl::params<uint64_t>::P[p]) == 0);
                    }

                }
            }
            for (size_t p = 0; p < NbPrimesQ; p++) {
                for (size_t alphaCANidx = 0; alphaCANidx < alphasCAN.size(); alphaCANidx++) {

                    size_t i = alphaCANidx * NbPrimesQ;
								}
						}
/*  
					// Mock Sigma protocol test
                
                std::vector<mpz_class> sigma_e(128), sigma_z(128);
                std::vector<mpz_class> sigma_e_GCD(alphasGCD.size()*NbPrimesP);
                std::vector<mpz_class> weight_per_gamma(128);
                std::vector<mpz_class> sigma_z_GCD(alphasGCD.size());
                {
                    // In the U matrix, the prover needs to include e.sigamrGCD and e.sigmaxGCD
                    // Then, concatenate the root hash of the Merkle tree with e.sigmar and then do Fiat-Shamir for the
                    // sigma_e and the degree tests.
                    // Note that sigmarGCD is basically 18 values obtained by computing the exponent used in computing sigmar modulo the 18 elements in alphasGCD
                    // sigmar is a vector of 128 elements. sigmargcd is vector of 128*18*9, where 9 is the number of zk prime moduli
                    
                    
                    // First completing the sigma ZK protocol and then verifying the main ZK proof
                    for (size_t j = 0 ; j < 128; j++) {
                        sigma_e[j] = math::generateRandomValue(mpz_class(std::random_device()()),
                                                               128); //obtain from Fiat-Shamir
                        weight_per_gamma[j] = math::generateRandomValue(mpz_class(std::random_device()()),
                                                                        64); //obtain from Fiat-Shamir
                        sigma_z[j] = e.sigmar[j] + sigma_e[j] * gammaExponent();
                        // Verifier's action for the sigma protocol
                        mpz_class sigma_lhs = math::powm(gBeforeExp()[j], sigma_z[j], finalModuli());
                        mpz_class sigma_rhs = math::mod(d.sigmaa[j] * math::powm(gAfterExp()[j], sigma_e[j], finalModuli()), finalModuli());
                        
                        assert(sigma_lhs == sigma_rhs);
                        
                        for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                            sigma_z_GCD[alphaGCDidx] += math::mod(weight_per_gamma[j] * sigma_z[j], alphasGCD[alphaGCDidx]);
                            
                            for (size_t i = 0; i < NbPrimesP; i ++) {
                                mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                                size_t pos = i * 128 * alphasGCD.size() + j * alphasGCD.size() + alphaGCDidx;
                                mpz_class sigmaeGCD = math::mod(sigma_e[j], alphasGCD[alphaGCDidx]);
                                mpz_class weight = math::mod(weight_per_gamma[j], alphasGCD[alphaGCDidx]);
                                sigma_e_GCD[i * alphasGCD.size() + alphaGCDidx] += positiveRemainder(sigmaeGCD,
                                                                                                     prime) * positiveRemainder(weight, prime);
                                //sigma_z_GCD[pos] = positiveRemainder(sigmazGCD,  prime);
                                //mpz_class sigmarGCD = math::mod(e.sigmar[j],alphasGCD[alphaGCDidx]);
                                //mpz_class sigmaxGCD = math::mod(gammaExponent(),alphasGCD[alphaGCDidx]);
                                //mpz_class temp = sigma_r_GCD[pos]+sigma_e_GCD[pos]*sigma_x_GCD[i*alphasGCD.size() + alphaGCDidx];
                                }
                        }
                    }
                    
                    std::vector<mpz_class> lintest(129 * alphasGCD.size()*NbPrimesP);
                    
                    for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                        for (size_t i = 0; i < NbPrimesP; i ++) {
                            sigma_e_GCD[i * alphasGCD.size() + alphaGCDidx] = math::mod(sigma_e_GCD[i * alphasGCD.size() +
                                                                                                    alphaGCDidx], mpz_class(nfl::params<uint64_t>::P[i]));
                            lintest[i * alphasGCD.size() * 129 + alphaGCDidx] = sigma_blind[i * alphasGCD.size() * 129 +
                                                                                            alphaGCDidx]
                                                                                + sigma_x_GCD[i * alphasGCD.size() + alphaGCDidx]
                                                                                * sigma_e_GCD[i * alphasGCD.size() + alphaGCDidx];
                                                                                
                        }
                    }
                    
                    for (size_t j = 0; j < 128; j++) {
                        for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                            for (size_t i = 0; i < NbPrimesP; i ++) {
                                mpz_class weight = math::mod(weight_per_gamma[j], alphasGCD[alphaGCDidx]);
                                lintest[i * alphasGCD.size() * 129 + (j + 1)*alphasGCD.size() + alphaGCDidx]
                                    = sigma_blind[i * alphasGCD.size() * 129 + (j + 1) * alphasGCD.size() + alphaGCDidx]
                                      + weight * sigma_r_GCD[i * alphasGCD.size() * 128 + j * alphasGCD.size() + alphaGCDidx];
                            }
                        }
                    }
                    
                    // Checking the mock linear test
                    //
                    std::vector<mpz_class> crt_primes(NbPrimesP);
                    
                    for (size_t alphaGCDidx = 0; alphaGCDidx < alphasGCD.size(); alphaGCDidx++) {
                        std::vector<mpz_class> crt_vec(NbPrimesP);
                        
                        for (size_t i = 0; i < NbPrimesP; i ++) {
                            crt_primes[i] = mpz_class(nfl::params<uint64_t>::P[i]);
                            
                            for (size_t j = 0; j < 129; j++) {
                                crt_vec[i] += lintest[i * alphasGCD.size() * 129 + j * alphasGCD.size() + alphaGCDidx];
                            }
                        }
                        
                        for (size_t i = 0; i < NbPrimesP; i ++) {
                            crt_vec[i] = math::mod(crt_vec[i], crt_primes[i]);
                        }
                        
                        std::vector<mpz_class> coeffs_temp;
                        mpz_class crt_result = math::crt_reconstruct(crt_vec, coeffs_temp, crt_primes);
                        //std::cout<< speciald() << ": CRT result = " << math::mod(crt_result, alphasGCD[alphaGCDidx])<< std::endl;
                        //std::cout<< speciald() << ": Z mod alpha[" << alphaGCDidx << "] = " << math::mod(sigma_z_GCD[alphaGCDidx], alphasGCD[alphaGCDidx])<< std::endl;
                        assert(math::mod(crt_result, alphasGCD[alphaGCDidx]) == math::mod(sigma_z_GCD[alphaGCDidx],
                                                                                          alphasGCD[alphaGCDidx]));
                                                                                          
                    }
                }
				*/
        }

        /* ========== variables used for zero knowledge ========== */
    protected:
        // variables from `generateKeyPair`
        Q _siP, _eiP, _si, _ei, _bi, _A, _ai, _b;
        mpz_class _jacobiSeedShares, _jacobiAndGCDSeedShares;

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
        uint32_t _finalindex;

        mpz_class _finalModuli,_gammaExponent;


        std::vector<mpz_class> _candidatesCAN;
        std::vector<mpz_class> _candidatesGCD;

        std::vector<mpz_class> _x_shares;
        std::vector<mpz_class> _y_shares;
        std::vector<mpz_class> _z_shares;
        lattice::cipher<T, Degree, NbPrimesQ>  _xsum_first, _xsum_final, _enc_z_shares,
                _xyz_sum, _enc_x_shares;

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

        // sigma
        
        std::vector<mpz_class> _gBeforeExp;
        std::vector<mpz_class> _gAfterExp;

        // legacy
        std::vector<mpz_class> _a_shares, _b_shares;
        cipherText _as, _es;
        std::vector<cipherText> _as_vec, _es_vec;
        Q _partial_e_shares;
        std::vector<int> _sievingFlags;

        // others
        SocketId& _socketId;
};

template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
ProtocolConfig<T> EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::registerAndAwaitConfiguration(
    ZeroMQClientTransport &transport, std::string ipAddress) {
    DBG("Participant (" << transport.getSocketId() << "," << ipAddress <<
        ") registering with coordinator.");


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
                math::computeHash(ai()),
                math::computeHash(jacobiSeedShares()),
                math::computeHash(jacobiAndGCDSeedShares())
                ));

    auto success = transport.awaitReply<ProtocolConfig<T>> (MessageType::PROTOCOL_CONFIG);
    if (!success) {
        throw std::runtime_error("Kill/Restart message received during registration!");
    }
    return *success;
}

/*
 * Generate key pair in a distributed way.
 * return a public key (a, b) and a private key si
 */
template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::optional<lattice::key_pair<T, Degree, NbPrimesQ>> EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::generateKeyPair(ZeroMQClientTransport& trans, nfl::gaussian<uint16_t, T, 2>& chi) {

    DBG("Generating Keys");
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

    // NOTE: ai is now initialized in registerAndAwaitConfiguration
    
    DBG("Sending A shares");
    trans.send<Q>(MessageType::PUBLIC_KEY_A_SHARES, ai());
    auto maybe_a = trans.awaitReply<Q>(MessageType::PUBLIC_KEY_A_VALUE);
    if (!maybe_a) {
        LOG(ERROR) << "Kill/Restart received during keygen A";
        return std::nullopt;
    }
    A() = *maybe_a;
    DBG("Received A");

    // gaussian dist
    siP() = Q(chi);
    eiP() = Q(chi);

    // keep everything in FFT domain
    si() = siP();
    ei() = eiP();

    si().ntt_pow_phi();
    ei().ntt_pow_phi();

    bi() = si() * A() + ei();

    DBG("Sending B shares");
    trans.send<Q>(MessageType::PUBLIC_KEY_B_SHARES, bi());
    auto maybe_b = trans.awaitReply<Q>(MessageType::PUBLIC_KEY_B_VALUE);
    if (!maybe_b) {
        LOG(ERROR) << "Kill/Restart received during keygen B";
        return std::nullopt;
    }
    DBG("Received B");
    b() = *maybe_b;

    // return (public key, private key)
    return {{{A(), b()}, si()}};
}

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


template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
CandidateIndices
EncryptedClient<FieldT, T, Degree, NbPrimesP, NbPrimesQ>::getCandidateIndices(
    size_t idx) {

    std::vector<size_t> ps;
    size_t i;
    for (i = 0; i < primesPS; ++i) {
        if (index_candidates()[i] == idx) {
            ps.push_back(i);
        }
    }

    std::vector<size_t> can;
    for (; i < primesPS + primesCAN; ++i) {
        if (index_candidates()[i] == idx) {
            can.push_back(i);
        }
    }

    std::vector<size_t> gcd;
    for (; i < primesPS + primesCAN + primesGCD; ++i) {
        if (index_candidates()[i] ==
                -2) { // -2 is a special value used only for GCD. Only one candidate uses it
            gcd.push_back(i);
        }
    }

    return CandidateIndices({ps, can, gcd});
}

template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::optional<std::tuple<
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
    using P = nfl::poly_p<T, Degree, NbPrimesP>;
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;


    // For the pre-sieving part of our candidate generation, we need prime
    // buckets s.t. the product of buckets has a bit-length greater than 1024.
    auto [alphasPS, bucketSizePS] = math::balanced_bucket_n_primes(config.pbs(), primesPS, config.tauLimitBit(), 1);
    auto [alphasCAN, bucketSizeCAN] = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());
    auto [alphasGCD, bucketSizeGCD] = math::fixed_bucket_n_primes(3*config.pbs()+96, primesGCD, config.tauLimitBit());

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


    // set bucketSizeGCD same size = lrint(floor(double(primesCAN)/alphasCAN.size()));

    auto [alphasdummy, bucketSizedummy] = math::balanced_bucket_n_primes(
            config.pbs(), Degree - primesTOTAL, config.tauLimitBit(), 1);

    // Now each party picks `degree` shares of random numbers `a` from
    // the set {0, ..., alpha_t - 1} for each alpha_t.

    a_shares().resize(primesPS);
    b_shares().resize(primesPS);

    x_shares() = math::generateRandomVector(std::random_device()(), Degree,
                                            config.tauLimitBit());
    y_shares() = math::generateRandomVector(std::random_device()(), Degree,
                                            config.tauLimitBit());
    z_shares() = math::generateRandomVector(std::random_device()(), Degree,
                                            config.tauLimitBit());

    index_candidates() = std::vector<size_t>(Degree);

    // i from 0 to degree - 1
    // get(i) -> candidate_number from 0 to 3600
    // i from 0 to primesPS - 1 -> pre-sieving
    // i from primesPS to primesPS +primesCAN - 1 -> candidate generation
    // i from primesPS + primesCAN  to degree - 1 -> GCD
    // x = [ 1 1 1 2 2 2.. .... -1 -1 -1 ... ... ]
    // y = [ 1 1 1 2 2 2.. .... ... ... ]
    // z = [ 1 1 1 2 2 2.. .... ... ... ]
    // only one array needed because they are the same
    // -1 not leaded to any candidate
    //
    // if I take candidate i = 1, then it will point exactly 6 positions in [0, primesPS - 1],
    // and exactly to 6 positions in [primesPS, primesPS + primesCAN]
    //
    // given a candidate give me a map of these 6 indeces from first and 6 from second so it is 12 indices
    //
    //                                               pre-sieving         candidateGen        GCD
    // giveCandidateIndeces(MPInt: candidate) -> indices({{1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5,6, 7, ..., 18 }})
    //
    //
    //
    // Structure indices
    // {
    //    std::vector preSieving;
    //    std::vector candidateGeneration;
    //    std::vector* GCD; // NULL
    // }
    //
    //
    // Create an array [ 1 2 3 4 ...]
    // When you prune discard prune this array the same way
    //

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

    // Store x, y, z, for ZK

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
    if (!maybe_xsum) {
        LOG(ERROR) << "Kill/Restart received during presieve";
        return std::nullopt;
    }
    auto xsum = *maybe_xsum;
    xsum_first() = xsum;
    DBG("Participant (" << transport.getSocketId() << ") received `Enc(a)`");

    auto [enc_z_s, uzi, vzi, wzi] = e.encrypt(z_shares(), publicKey);
    enc_z_shares() = enc_z_s;
    uz() = uzi;
    vz() = vzi;
    wz() = wzi;

    // std::cout << "sanity check:" << enc_z_shares().first.poly_obj().data()[0] - (uz().poly_obj().data()[0]*publicKey.first.poly_obj().data()[0] + wz().poly_obj().data()[0]) << std::endl;
    // std::cout << "sanity check:" << A().poly_obj().data()[0] << std::endl;
    // std::cout << "sanity check:" << publicKey.first.poly_obj().data()[0] << std::endl;

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
    if (!maybe_xyz_sum) {
        LOG(ERROR) << "Kill/Restart received (xyz_sum)";
        return std::nullopt;
    }
    xyz_sum() = *maybe_xyz_sum;
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


    /*
    // Next for debugging purposes sending the shares in the clear
    P p_b_shares,p_a_shares;
    {
            auto *ptr = new std::array<mpz_class, Degree>{};
            std::transform(b_shares().begin(), b_shares().end(), ptr->begin(), [](const auto& x) { return x; });
            p_b_shares.set_mpz(*ptr);
            delete(ptr);
    }

    {
            auto *ptr = new std::array<mpz_class, Degree>{};
            std::transform(a_shares().begin(), a_shares().end(), ptr->begin(), [](const auto& x) { return x; });
            p_a_shares.set_mpz(*ptr);
            delete(ptr);
    }

    transport.send(MessageType::P_CLEAR_DEBUG, p_a_shares);
    DBG("Participant (" << transport.getSocketId() << ") sending `a` shares in the clear");

    uint64_t flag = transport.awaitReply<uint64_t>(MessageType::MUTHU_ACK);
    DBG("Participant (" << transport.getSocketId() << ") received 0 as an ACK");
    transport.send(MessageType::Q_CLEAR_DEBUG, p_b_shares);
    DBG("Participant (" << transport.getSocketId() << ") sending `b` shares in the clear");

    */

    // Now the coordinator can evaluate which of the `a_shares` we sent are valid. Many of them
    // will be eliminated, and the coordinator will communicate which values to eliminate via
    // a binary matrix, which we then use to shuffle our shares around before the CRT reconstruction step.
    auto maybe_sievingFlags = transport.awaitReply<std::vector<int>>
                              (MessageType::PS_SIEVING_FLAGS);
    if (!maybe_sievingFlags) {
        LOG(ERROR) << "Kill/Received (sievingFlags)";
        return std::nullopt;
    }
    sievingFlags() = *maybe_sievingFlags;
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

template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::optional<std::tuple<std::vector<mpz_class>, std::vector<mpz_class>, std::vector<mpz_class>>>
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
    if (!ax_by) {
        LOG(ERROR) << "Kill/Restart received (modulus generation, ax_by";
        return std::nullopt;
    }
    auto [tax, tby] = *ax_by;
    ax() = tax;
    by() = tby;
    DBG("Participant (" << transport.getSocketId() <<
        ") received `ax and by` values");

    // candidate generation
    //std::vector<mpz_class> axby(ax().size());
    axby().resize(ax().size());
    {
        int k = 0;
        for (size_t i = 0; i < bucketSize; ++i) {
            for (size_t j = 0; j < alphasCAN.size(); ++j) {
                axby()[k] = math::mod((ax()[k] * q_sharesCRTs()[k] + by()[k] * p_sharesCRTs()[k]
                                       + zcan[k]), alphasCAN[j]);
                /*if (special) {*/
                //axby[k] = math::mod(axby[k] - ax[k]*by[k], alphasCAN[j]);
                /*}*/
                ++k;
            }
        }
    }

    transport.send(MessageType::AXB_MINUS_BYA_SHARES, axby());
    DBG("Participant (" << transport.getSocketId() <<
        ") received `axb - bya` shares");

    auto maybe_candidates = transport.awaitReply<std::vector<mpz_class>>
                            (MessageType::MODULUS_CANDIDATE);
    if (!maybe_candidates) {
        LOG(ERROR) << "Kill/Restart received (modulus generation, candidates";
        return std::nullopt;
    }
    auto candidates = *maybe_candidates;
    DBG("Participant (" << transport.getSocketId() << ") received `N` candidates");
    // In return the coordinator will send us a plaintext vector of modulus candidates.

    candidates.resize(pq_size);
    p_shares.resize(pq_size);
    q_shares.resize(pq_size);

    std::vector<int> syncF;
    syncF.push_back(1);
    transport.send(MessageType::MUTHU_ACK, int(1));

    /*
    transport.send(MessageType::MUTHU_ACK, xcan);
    auto syncF = transport.awaitReply<int>(MessageType::MUTHU_ACK);
    transport.send(MessageType::MUTHU_ACK, ycan);
    syncF = transport.awaitReply<int>(MessageType::MUTHU_ACK);
    transport.send(MessageType::MUTHU_ACK, zcan);
    syncF = transport.awaitReply<int>(MessageType::MUTHU_ACK);
    */
    //transport.send(MessageType::MUTHU_ACK, std::pair{p_shares, q_shares});
    //auto syncF = transport.awaitReply<int>(MessageType::MUTHU_ACK);


    candidatesCAN() = candidates;
    return std::make_tuple(candidates, p_shares, q_shares);
}


/** Runs the client protocol for performing the first of the biprimality tests. */
template<typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::optional<int>
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
    if (!maybe_gammaSeed) {
        LOG(ERROR) << "Kill/Restart received (performJacobiTest, gammaSeed";
        return std::nullopt;
    }
    auto gammaSeed = MPInt(maybe_gammaSeed->get_mpz_t());
    // We find a set of valid gamma values by their jacobi symbol...
    std::vector<mpz_class> gammaValues(candidates.size());

    DBG("gammaSeed = " << gammaSeed);
    std::vector<mpz_class> engine = math::generateRandomVector(*maybe_gammaSeed,
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


/** Runs the client protocol for performing the second of the biprimality tests. */
template<typename FieldT, typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::optional<int>
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
    auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3*config.pbs() + 96, primesGCD, config.tauLimitBit());
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
            mpz_class rx = grc.get_z_bits(2 * config.pbs() + 48);
            DBG("doing CRT deconstruct ");
            auto aCRT = math::crt_deconstruct(rx, alphasGCD);
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
    if (!maybe_axby) {
        LOG(ERROR) << "Kill/Restarted received (GCDandJacobiTest, axby)";
        return std::nullopt;
    }
    auto [axby_pair, gammaSeed] = *maybe_axby;
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
            mpz_class orig_g = mpz_class(g);
            g = math::powm(g, exp, N);

            gAfterExp()[i * config.lambda() + j] = mpz_class(g);

            mpz_class sigma_r = math::generateRandomValue(mpz_class(std::random_device()()), 2 * config.pbs() + 48 + config.lambda());
            mpz_class sigma_a = math::powm(orig_g, sigma_r, N);

            mpz_class sigma_e = math::generateRandomValue(mpz_class(std::random_device()()), config.pbs() + 24 + config.lambda());
            mpz_class sigma_z = sigma_r + sigma_e * exp;

            mpz_class sigma_lhs = math::powm(orig_g, sigma_z, N);
            mpz_class sigma_rhs = math::mod(sigma_a * math::powm(g, sigma_e, N), N);

            if (sigma_lhs != sigma_rhs) {
                DBG("g = " << g);
                DBG("orig_g = " << orig_g);
                DBG("sigma_r = " << sigma_r);
                DBG("sigma_a = " << sigma_a);
                DBG("sigma_e = " << sigma_e);
                DBG("sigma_z = " << sigma_z);

                DBG("sigma_lhs = " << sigma_lhs);
                DBG("sigma_rhs = " << sigma_rhs);
            }
            assert(sigma_lhs == sigma_rhs);

        }
    }


    // GCD: Compute  axby
    axbyGCD().resize(axGCD().size());
    {
        int k = 0;
        for (size_t i = 0; i < bucketSize; ++i) {
            for (size_t j = 0; j < alphasGCD.size(); ++j) {
                axbyGCD()[k] = math::mod((axGCD()[k] * p_plus_qCRTs[k] + byGCD()[k] * rCRTs()[k]
                                          + z[k]), alphasGCD[j]);
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
    //DBG("Client gammaValues = " << gammaValues);
    transport.send(MessageType::AXB_MINUS_BYA_SHARES, std::pair{axbyGCD(), gammaValues});
    DBG("Participant (" << transport.getSocketId() <<
        ") received `axb - bya` shares");
    return 0;
}


} // namespace ligero
