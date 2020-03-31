#pragma once


#include "Transport.hpp"
#include "Common.hpp"
#include "Math.hpp"
#include <boost/multiprecision/miller_rabin.hpp>
#include "LatticeEncryption.hpp"
#include "Factoring.hpp"
#include <functional>
#include <unordered_set>
#include <boost/functional/hash.hpp>

namespace ligero
{


//==============================================================================
template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
class EncryptedCoordinator
{
public:

    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

    ~EncryptedCoordinator() {
        saveProtocolRecord();
    }

    EncryptedCoordinator(ProtocolConfig<T>& _config) : config(_config) {}

    /*
     * Host key generation in the coordinator side.
     */
        expected<int> hostGenerateKeyPair(ZeroMQCoordinatorTransport& trans) {
            DBG("Starting keygen");
            // sum up shares of a_i
            DBG("Waiting for a shares");
            auto maybeA = trans.awaitAggregateVectorInput<Q>(MessageType::PUBLIC_KEY_A_SHARES, trans.ids, std::plus<Q>(), 0);
            if (hasError(maybeA)) { return getError(maybeA); }

            std::tie(A, std::ignore) = getResult(maybeA);
            trans.broadcast(MessageType::PUBLIC_KEY_A_VALUE, trans.ids, A);
            A.invntt_pow_invphi();

            DBG("Received and broadcasted a, waiting for b");
            auto maybeB = trans.awaitAggregateVectorInput<Q>(MessageType::PUBLIC_KEY_B_SHARES, trans.ids, std::plus<Q>(), 0);
            if (hasError(maybeB)) { return getError(maybeB); }

            std::tie(B, bi) = getResult(maybeB);
            trans.broadcast(MessageType::PUBLIC_KEY_B_VALUE, trans.ids, B);
            B.invntt_pow_invphi();

            DBG("Received and broadcasted b");
            return 0;
        }

        std::vector<std::vector<mpz_class>>
            pruneAndReorderShares(
                const std::vector<mpz_class> &as,
                const std::vector<mpz_class> &bs,
                const std::vector<int> &flags,
                int numAlphas,
                std::vector<size_t> bucketSize) {

            assert(as.size() == flags.size());
            assert(bs.size() == flags.size());

            std::vector<std::vector<mpz_class>> result;

            int k = 0;
            for (int c = 0; c < numAlphas; ++c) {
                int j = 0;

                std::vector<mpz_class> col;
                for (int r = 0; r < bucketSize[c]; ++r) {
                    if (flags[k] == 1) {
                        col.push_back(as[k]);
                        col.push_back(bs[k]);
                    } else {
                    }
                    k++;
                }
                result.push_back(col);
            }

            return result;
        }

        bool checkGCD(bool special,
                      const SocketId &sid, 
                      const std::vector<mpz_class>& p_shares, 
                      const std::vector<mpz_class>& q_shares, 
                      const std::vector<mpz_class>& candidates,
                      const std::vector<mpz_class> &x,
                      const std::vector<mpz_class> &y,
                      const std::vector<mpz_class> &z,
                      mpz_class gcdRX,
                      const std::vector<mpz_class> ss)
        {
            auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3*config.pbs() + 210, primesGCD, config.tauLimitBit());
            const int bucketSize = lrint(floor(double(primesGCD) / double(alphasGCD.size())));
            DBG("bucketSize = " << bucketSize);
            DBG("alphasGCD.size() = " << alphasGCD.size());

            std::vector<mpz_class> ax_shares(primesGCD);
            std::vector<mpz_class> by_shares(primesGCD);

            DBG("p_shares.size() = " << p_shares.size());
            assert(x.size() >= (p_shares.size() * alphasGCD.size()));
            assert(y.size() >= (p_shares.size() * alphasGCD.size()));
            assert(z.size() >= (p_shares.size() * alphasGCD.size()));

            std::vector<mpz_class> rCRTs;
            rCRTs.resize(p_shares.size() * alphasGCD.size());
            std::vector<mpz_class> p_plus_qCRTs(p_shares.size() * alphasGCD.size());

            const mpz_class maybeOne = special ? 1 : 0;
            {
                int rc = 0;
                for (int i = 0; i < candidates.size(); ++i) {

                    // a = random number
                    auto aCRT = math::crt_deconstruct(gcdRX, alphasGCD);
                    assert(aCRT.size() == alphasGCD.size());

                    // b = p_i + q_i - maybeOne
                    auto bCRT = math::crt_deconstruct(p_shares[i] + q_shares[i] - maybeOne, alphasGCD);
                    assert(bCRT.size() == alphasGCD.size());

                    for (int j = 0; j < aCRT.size(); ++j) {
                        rCRTs[rc] = aCRT[j];
                        ax_shares[rc] = rCRTs[rc] - x[rc];
                        p_plus_qCRTs[rc] = bCRT[j];
                        by_shares[rc] = p_plus_qCRTs[rc] - y[rc];
                        ++rc;
                    }
                }
            }

            {
                int k = 0;
                for (int i = 0; i < bucketSize; ++i) {
                    for (int j = 0; j < alphasGCD.size(); ++j) {
                        ax_shares[k] = math::mod(ax_shares[k], alphasGCD[j]);
                        by_shares[k] = math::mod(by_shares[k], alphasGCD[j]);
                        ++k;
                    }
                }
            }

            // here check if that is what party sent
            {
                auto axbyShareItem = axby_sharesGCD.find(sid);
                if (axbyShareItem == axby_sharesGCD.end()) {
                    LOG(INFO) << "Failed to find axby GCD share for " << sid;
                    return false;
                }
                auto [axFromParty, byFromParty] = axbyShareItem->second.first;
                if ((ax_shares.size() != axFromParty.size()) ||
                    (by_shares.size() != byFromParty.size())) {
                    LOG(INFO) << "ax/by computed size " << ax_shares.size() << " does not match ax/by size from party " << sid;
                    return false;
                }

                for (size_t i = 0; i < ax_shares.size(); i++) {
                    if (ax_shares[i] != axFromParty[i]) {
                        LOG(INFO) << "Computed ax[" << i << "] = " << ax_shares[i] << " does not match ax from party " << axFromParty[i] << " for party " << sid;
                        return false;
                    }
                    if (by_shares[i] != byFromParty[i]) {
                        LOG(INFO) << "Computed by[" << i << "] = " << by_shares[i] << " does not match by from party " << byFromParty[i] << " for party " << sid;
                        return false;
                    }
                }
            }


            // now lets check axbyGCD
            std::vector<mpz_class> axbyGCD;
            axbyGCD.resize(axGCD.size());
            {
                std::vector<mpz_class> ssGCD;
                int k = 0;
                ssGCD.resize(ss.size());
                LOG(INFO) << "candidates.size() " << candidates.size();
                LOG(INFO) << "bucketSize " << bucketSize;
                for (size_t i = 0; i < bucketSize; ++i) {
                    ssGCD[i] = ss[i] * candidates[i];
                    for (size_t j = 0; j < alphasGCD.size(); ++j) {
                        axbyGCD[k] = math::mod((axGCD[k] * p_plus_qCRTs[k] + byGCD[k] * rCRTs[k] + z[k] + ssGCD[i]), alphasGCD[j]);
                        ++k;
                    }
                }
            }

            auto axbyGCD_share = zcrt_shares.find(sid);
            if (axbyGCD_share == zcrt_shares.end()) {
                LOG(INFO) << "Could not find encrypted axbyGCD share for " << sid;
                return false;
            }
            std::vector<mpz_class> axbyGCDFromParty = axbyGCD_share->second.first;

            if (axbyGCDFromParty.size() != axbyGCD.size()) {
                return false;
            }

            for (size_t i = 0; i < axbyGCD.size(); i++) {
                if (axbyGCD[i] != axbyGCDFromParty[i]) {
                    return false;
                }
            }

            return true;
        }


        bool checkJacobi(bool special,
                         const SocketId &sid, 
                         const std::vector<mpz_class>& p_shares, 
                         const std::vector<mpz_class>& q_shares, 
                         const std::vector<mpz_class>& candidates, 
                         mpz_class gammaSeedValue,
                         const std::vector<mpz_class>& gvShare)
        {
            std::vector<mpz_class> engine = math::generateRandomVector(gammaSeedValue, kJacobiNumberOfRandomValues, 2048);
            std::vector<mpz_class> gammaValues(p_shares.size());
            size_t engineCount = 0;

            {
                mpz_t one;
                mpz_init(one);
                mpz_set_ui(one, 1);

                mpz_t gcdResult;
                mpz_init(gcdResult);

                for (int i = 0; i < p_shares.size(); ++i) {
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

            // Now Verify that gammaValues are equal
            if (gvShare.size() != gammaValues.size()) {
                LOG(INFO) << "Party " << sid << " gamma values size " << gvShare.size() 
                    << " did not match computed gamma values size " << gammaValues.size();
                return false;
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

                if (gammaValues[i] != gvShare[i]) {
                    LOG(INFO) << "Compute gamma value " << gammaValues[i] 
                        << " did not match shared gamma value " << gvShare[i]
                        << " for party " << sid;
                    return false;
                }
            }
            return true;
        }


        bool verifyNoModuli(ZeroMQCoordinatorTransport& transport,
                            lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& e) {
            LOG(INFO) << "Verifying party randomness";

            bool passed = true;
            std::vector<mpz_class> tauVector(Degree);
            TripleVector _tv;
            auto maybeXYZShares = transport.awaitAggregateVectorInput<TripleVector>
                (
                 MessageType::NO_MODULI_VERIFICATION_SHARES,
                 ZeroMQCoordinatorTransport::ids,
                 [](const TripleVector& a, const TripleVector& b) {
                 return a;
                 },
                 _tv
                );

            if (hasError(maybeXYZShares)) {
                DBG("Failed to receive randomness shares");
                transport.broadcast(MessageType::NO_MODULI_VERIFICATION_FINISHED, ZeroMQCoordinatorTransport::ids, false);
                return false;
            }


            std::unordered_map<SocketId, TripleVector, boost::hash<SocketId>> xyzShares;

            std::tie(std::ignore, xyzShares) = getResult(maybeXYZShares);

            // Verify PS

            std::vector<mpz_class> ashares(primesPS);
            std::vector<mpz_class> bshares(primesPS);

            std::vector<mpz_class> prime_shares;
            // extract x, y, z
            std::vector<mpz_class> x_shares(Degree, mpz_class(0)),
                                   y_shares(Degree, mpz_class(0)),
                                   z_shares(Degree, mpz_class(0));


            std::unordered_map<SocketId, std::vector<mpz_class>, boost::hash<SocketId>> primeSharesMap;

            if (passed) {

                auto [alphasdummy, bucketSizedummy] = math::balanced_bucket_n_primes(
                        config.pbs(), Degree - primesTOTAL, config.tauLimitBit(), 1);

                for (auto it = xyzShares.begin(); it != xyzShares.end(); ++it ) {
                    auto [xi, yzi] = it->second;
                    auto [yi, zi] = yzi;

                    std::vector<mpz_class> ai_shares(primesPS);
                    std::vector<mpz_class> bi_shares(primesPS);

                    int k = 0;
                    for (size_t j = 0; j < alphasPS.size(); ++j) {
                        for (size_t i = 0; i < bucketSizePS[j]; ++i) {
                            x_shares[k] += xi[k];
                            y_shares[k] += yi[k];
                            z_shares[k] += zi[k];
                            ai_shares[k] = xi[k] % alphasPS[j];
                            bi_shares[k] = yi[k] % alphasPS[j];
                            ashares[k] = (ashares[k] + ai_shares[k]) % alphasPS[j];
                            bshares[k] = (bshares[k] + bi_shares[k]) % alphasPS[j];
                            tauVector[k] = alphasPS[j];
                            k++;
                        }
                    }

                    // candidate generation
                    for (size_t i = 0; i < bucketSizeCAN_value; ++i) {
                        for (size_t j = 0; j < alphasCAN.size(); ++j) {
                            x_shares[k] += xi[k];
                            y_shares[k] += yi[k];
                            z_shares[k] += zi[k];
                            tauVector[k] = alphasCAN[j];
                            k++;
                        }
                    }

                    // GCD
                    for (size_t i = 0; i < bucketSizeGCD_value; ++i) {
                        for (size_t j = 0; j < alphasGCD.size(); ++j) {
                            x_shares[k] += xi[k];
                            y_shares[k] += yi[k];
                            z_shares[k] += zi[k];
                            tauVector[k] = alphasGCD[j];
                            k++;
                        }
                    }

                    // setting up tauVector for dummy
                    for (size_t j = 0; j < alphasdummy.size(); ++j) {
                        for (size_t i = 0; i < bucketSizedummy[j]; ++i) {
                            tauVector[k] = alphasdummy[j];
                            k++;
                        }
                    }



                    bool special = it->first == ZeroMQCoordinatorTransport::ids[0];
                    std::vector<std::vector<mpz_class>> valid_shares = pruneAndReorderShares(ai_shares, bi_shares, flags_can, alphasPS.size(), bucketSizePS);
                    std::vector<mpz_class> prime_shares_i;
                    std::vector<mpz_class> _alphas(alphasPS.size() + 1);
                    std::vector<mpz_class> _coeffs(_alphas.size());
                    for (size_t i = 0; i < alphasPS.size(); ++i) {
                        _alphas[i] = alphasPS[i];
                    }
                    _alphas[_alphas.size() - 1] = mpz_class(4);


                    std::vector<mpz_class> moduli(_alphas.size());
                    size_t minRowSize = valid_shares[0].size();
                    for (size_t i = 1; i < valid_shares.size(); ++i) {
                        minRowSize = std::min(valid_shares[i].size(), minRowSize);
                    }

                    for (size_t i = 0; i < minRowSize; ++i) {
                        for (size_t j = 0; j < alphasPS.size(); ++j) {
                            _coeffs[j] = valid_shares[j][i];
                        }
                        _coeffs[_alphas.size() - 1] = (special ? mpz_class(3) : mpz_class(0));

                        prime_shares_i.push_back(math::crt_reconstruct(_coeffs, moduli, _alphas));
                    }

                    primeSharesMap[it->first] = prime_shares_i;

                    if (prime_shares.size() == 0) {
                        prime_shares = prime_shares_i;
                    } else {
                        for (size_t cpi = 0; cpi < prime_shares_i.size(); cpi++) {
                            prime_shares[cpi] += prime_shares_i[cpi];
                        }
                    }


                }

                // Check: eit[i] = ashares[i] * bshares[i] % tauvector[i]
                {
                    int k = 0;
                    for (size_t j = 0; j < alphasPS.size() && passed; ++j) {
                        for (size_t i = 0; i < bucketSizePS[j]; ++i) {
                            if ((ashares[k] * bshares[k]) % tauVector[k] != eitPS[k]) {
                                LOG(INFO) <<  "Failed verification for pre-sieving";
                                DBG("ashares[k] = " << ashares[k]);
                                DBG("bshares[k] = " << ashares[k]);
                                DBG("ashares[k] * bshares[k] \% tauVector[k] = " << ashares[k] * bshares[k] % tauVector[k]);
                                DBG("eitPS[k] = " << eitPS[k]);
                                passed = false;
                                break;
                            }
                            k++;
                        }
                    }
                }
            }

            const auto bucketSize = lrint(floor(double(primesCAN) / double(alphasCAN.size())));
            const auto pq_size = std::min(bucketSize, lrint(floor(double(prime_shares.size()) / double(2))));

            // Verify candidate generation
            if (passed) {
                
                // Verify candidates equal to p_i * q_i
                for (int i = 0; i < pq_size; ++i) {
                    if (candidatesCAN[i] != prime_shares[2 * i] * prime_shares[2 * i + 1]) {
                        LOG(INFO) << "Failed verification for candidates modulus generation";
                        DBG("candidatesCAN[" << i << "] = " << candidatesCAN[i]);
                        DBG("prime_shares[2 * i] * prime_shares[2 * i + 1] = " << prime_shares[2 * i] * prime_shares[2 * i + 1]);
                        for (size_t j = 0; j < alphasPS.size(); j++) {
                            DBG("candidatesCAN[" << i << "] \% alphasPS[0] = " << candidatesCAN[i] % alphasPS[j]);
                            DBG("prime_shares[2 * i] * prime_shares[2 * i + 1] \% alphasPS[0] = " << (prime_shares[2 * i] * prime_shares[2 * i + 1]) % alphasPS[j]);
                        }
                        passed = false;
                        break;
                    }
                }
                if (passed) {
                    // Run candidates through miller rabin
                    for (int i = 0; i < pq_size; ++i) {
                        if (boost::multiprecision::miller_rabin_test(MPInt(prime_shares[2 * i].get_mpz_t()), 4) &&
                           boost::multiprecision::miller_rabin_test(MPInt(prime_shares[2 * i + 1].get_mpz_t()), 4)) {

                            LOG(INFO) << "Verification for Miller-Rabin test failed (both are primes according to Miller-Rabin test).";
                            DBG("prime_shares[" << 2 * i << "] = " << prime_shares[2 * i]);
                            DBG("prime_shares[" << 2 * i + 1 << "] = " << prime_shares[2 * i + 1]);
                            passed = false;
                            break;
                        }

                    }
                }
            }


            if (!passed) {

                std::unordered_set<SocketId, boost::hash<SocketId>> cheatingParties;

                auto maybeBadEventData = transport.awaitAggregateVectorInput<BadEventData<T, Degree, NbPrimesQ>>
                    (
                     MessageType::BAD_EVENT_DATA_SHARES,
                     ZeroMQCoordinatorTransport::ids,
                     [](const BadEventData<T, Degree, NbPrimesQ> & a, const BadEventData<T, Degree, NbPrimesQ>& b) {
                     return a;
                     },
                     BadEventData<T, Degree, NbPrimesQ>(mpz_class(0), mpz_class(0), Q{0}, Q{0}, Q{0}, mpz_class(0), std::vector<mpz_class>())
                    );

                if (hasError(maybeBadEventData)) {
                    LOG(INFO) << "Failed to receive data shares from parties";
                    transport.broadcast(MessageType::BAD_EVENT_DATA_RESPONSE, ZeroMQCoordinatorTransport::ids, false);
                    return false;
                }

                std::unordered_map<SocketId, BadEventData<T, Degree, NbPrimesQ>, boost::hash<SocketId>> badEventDataShares;
                std::tie(std::ignore, badEventDataShares) = getResult(maybeBadEventData);

                Q SiSum = Q{0};

                for (auto it = badEventDataShares.begin(); it != badEventDataShares.end(); ++it) {

                    SiSum = SiSum + it->second.si;

                    // check ei bounding constaints, all elements should belong to (-80; 80)
                    bool passedEIBoundingConstraintVerification = true;
                    auto ei_coeffs = it->second.ei.poly2mpz();
                    for (auto& eic : ei_coeffs) {
                        mpz_class eicv = mpz_class(eic);
                        if ((eicv > mpz_class(80)) ||
                            (eicv < mpz_class(-80))) {

                            LOG(INFO) << "Failed check ei bounding constraints, should belong to (-80;80) for party " << it->first;
                            cheatingParties.insert(it->first);
                            passedEIBoundingConstraintVerification = false;
                        }
                    }
                    std::for_each(ei_coeffs.begin(), ei_coeffs.end(), mpz_clear);

                    if (passedEIBoundingConstraintVerification) {
                        // check bi == si * A + ei
                        if (it->second.bi != it->second.si * A + it->second.ei) {
                            LOG(INFO) << "Failed check bi != si * A + ei for parry " << it->first;
                            cheatingParties.insert(it->first);
                        }
                    }
                }


                // Check that xshares match after decryption
                {

                    for (auto& sid : ZeroMQCoordinatorTransport::ids) {
                        auto g_enc_x_s = g_enc_x_shares.find(sid);
                        if (g_enc_x_s == g_enc_x_shares.end()) {
                            LOG(INFO) << "Could not find encrypted x share for " << sid;
                            cheatingParties.insert(sid);
                        }
                        std::pair<Q, Q> cypherText = g_enc_x_s->second;
                        Q encX = cypherText.second - cypherText.first * SiSum;
                        auto decX = e.eval_poly(encX, tauVector);
                        auto xyzS = xyzShares.find(sid);
                        if (xyzS == xyzShares.end()) {
                            LOG(INFO) << "Could not find x,y,z shares for " << sid;
                            cheatingParties.insert(sid);
                        }
                        auto [xi, yzi] = xyzS->second;

                        // check if xi share for a given sid matches computed one
                        if (decX.size() != xi.size()) {
                            LOG(INFO) << "Decrypted x share size " << decX.size() << " does not match expected x share size " << xi.size() << " for party " << sid;
                            cheatingParties.insert(sid);
                        }

                        for (size_t i = 0; i < decX.size(); i++) {
                            if (mpz_class(decX[i]) != xi[i]) {
                                LOG(INFO) << "decrypted x[" << i << "] = " << decX[i] << " does not match original share x[" << i << "] = " << xi[i] << " for party " << sid;
                                cheatingParties.insert(sid);
                            }
                        }
                        std::for_each(decX.begin(), decX.end(), mpz_clear);
                    }
                }

                // check jacobi round and then jacobi and gcd round
                for (auto sid: ZeroMQCoordinatorTransport::ids) {

                    bool special = sid == ZeroMQCoordinatorTransport::ids[0];

                    // check Jacobi round
                    auto bothSharesItem = primeSharesMap.find(sid);
                    if (bothSharesItem == primeSharesMap.end()) {
                        LOG(INFO) << "Failed to find primeShare for party " << sid;
                        cheatingParties.insert(sid);
                        continue;
                    }

                    const auto bucketSize = lrint(floor(double(primesCAN) / double(alphasCAN.size())));
                    const auto pq_size = std::min(bucketSize, lrint(floor(double(bothSharesItem->second.size()) / double(2))));

                    std::vector<mpz_class> p_shares(pq_size);
                    std::vector<mpz_class> q_shares(pq_size);
                    for (int i = 0; i < pq_size; ++i) {
                        p_shares[i] = bothSharesItem->second[2 * i];
                        q_shares[i] = bothSharesItem->second[2 * i + 1];
                    }

                    // check Jacobi
                    {
                        p_shares = discardCandidates(p_shares, discardFlagsPostSieve);
                        q_shares = discardCandidates(q_shares, discardFlagsPostSieve);

                        auto ggsShareItem = ggsShares.find(sid);
                        if (ggsShareItem == ggsShares.end()) {
                            LOG(INFO) << "Failed to find ggs share for party " << sid;
                            return false;
                        }
                        std::vector<mpz_class> gvShare = ggsShareItem->second;
                        if (!checkJacobi(special, sid, p_shares, q_shares, candidatesPostSieve, gammaSeedAccum, gvShare)) {
                            cheatingParties.insert(sid);
                            continue;
                        }
                    }

                    // check Jacobi and GCD Round
                    {
                        p_shares = discardCandidates(p_shares, discardFlagsJacobi);
                        q_shares = discardCandidates(q_shares, discardFlagsJacobi);

                        auto gammaValuesGCD_share = zcrt_shares.find(sid);
                        if (gammaValuesGCD_share == zcrt_shares.end()) {
                            LOG(INFO) << "Could not find encrypted gammaValuesGCD share for " << sid;
                            cheatingParties.insert(sid);
                            continue;
                        }

                        std::vector<mpz_class> gammaValuesGCDFromParty = gammaValuesGCD_share->second.second;

                        if (!checkJacobi(special, sid, p_shares, q_shares, candidatesJacobi, gammaSeedJacobiAndGCD, gammaValuesGCDFromParty)) {
                            cheatingParties.insert(sid);
                            continue;
                        }

                        auto badEventDataItem = badEventDataShares.find(sid);
                        if (badEventDataItem == badEventDataShares.end()) {
                            cheatingParties.insert(sid);
                            continue;
                        }

                        auto xyzS = xyzShares.find(sid);
                        if (xyzS == xyzShares.end()) {
                            LOG(INFO) << "Could not find x,y,z shares for " << sid;
                            cheatingParties.insert(sid);
                            continue;
                        }
                        auto [xi, yzi] = xyzS->second;
                        auto [yi, zi] = yzi;
                        std::vector<mpz_class> xgcd(primesGCD);
                        std::vector<mpz_class> ygcd(primesGCD);
                        std::vector<mpz_class> zgcd(primesGCD);
                        {
                            int k = primesPS + primesCAN;
                            // GCD
                            int gcd_i = 0;
                            for (size_t i = 0; i < bucketSizeGCD_value; ++i) {
                                for (size_t j = 0; j < alphasGCD.size(); ++j) {
                                    xgcd[gcd_i] = xi[k];
                                    ygcd[gcd_i] = yi[k];
                                    zgcd[gcd_i] = zi[k];
                                    gcd_i++;
                                    k++;
                                }
                            }
                        }
                        if (!checkGCD(special,
                                    sid,
                                    p_shares,
                                    q_shares,
                                    candidatesJacobi,
                                    xgcd,
                                    ygcd,
                                    zgcd,
                                    badEventDataItem->second.gcdRX,
                                    badEventDataItem->second.gcdSS)) {
                            cheatingParties.insert(sid);
                        }
                    }
                }

                transport.broadcast(MessageType::BAD_EVENT_DATA_RESPONSE, ZeroMQCoordinatorTransport::ids, false);
                LOG(INFO) << "Party randomness verification failed";

                // kickout cheating parties
                std::vector<SocketId> restarts;
                for (auto sid : ZeroMQCoordinatorTransport::ids) {
                    if (cheatingParties.find(sid) == cheatingParties.end()) {
                        restarts.push_back(sid);
                    }
                }
                auto success = transport.update_ids(restarts);
                if (hasError(success)) {
                    throw std::runtime_error("Sometimes went wrong. Failed to kickout cheating parties");
                }
                else {
                    LOG(INFO) << "Update ids with " << getResult(success) << " parties";
                }

                return false;
            }

            transport.broadcast(MessageType::NO_MODULI_VERIFICATION_FINISHED, ZeroMQCoordinatorTransport::ids, passed);
            LOG(INFO) << "Party randomness verification passed";
            return true;
        }

        std::tuple<std::vector<mpz_class>,
                   std::vector<int>, 
                   std::vector<mpz_class>> computeEITAndFlags(std::array<mpz_t, Degree> c_mpz_t) {
                       
            std::vector<mpz_class> eit(primesPS);

            int pq_size = bucketSizePS[0];
            std::vector<int> flags(primesPS);
            assert(eit.size() == flags.size());
            std::vector<mpz_class> c(c_mpz_t.size());

            for (size_t i = 0; i < c.size(); ++i) {
                c[i] = mpz_class(c_mpz_t[i]);
            }


            mpz_t one;
            mpz_init(one);
            mpz_set_ui(one, 1);

            mpz_t gcdResult;
            mpz_init(gcdResult);


            int k = 0;
            for (int j = 0; j < alphasPS.size(); ++j) {

                int ick = 0;
                int sum = 0;
                for (int i = 0; i < bucketSizePS[j]; ++i) {

                    mpz_gcd(gcdResult, alphasPS[j].get_mpz_t(), c_mpz_t[k]);
                    eit[k] = c[k];

                    flags[k] = mpz_cmp(gcdResult, one) == 0;
                    sum += flags[k];
                    if (flags[k] == 1) {
                        index_candidates[k] = ick++;
                    } else {
                        index_candidates[k] = -1;
                    }

                    k++;
                }

                if (sum < pq_size) {
                    pq_size = sum;
                }
            }

            // free mem
            mpz_clear(one);
            mpz_clear(gcdResult);
            std::for_each(c_mpz_t.begin(), c_mpz_t.end(), mpz_clear);

            return {eit, flags, c};
        }

    expected<std::tuple< std::vector<std::vector<mpz_class>>,  // eit_pruned
             std::vector<mpz_class>,               // alphasPS_tick
             std::vector<mpz_class>,               // c_can
             std::vector<mpz_class>                // c_gcd
            >>
    hostPreSieving (
            ZeroMQCoordinatorTransport& transport,
            lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& e) {
        using P = nfl::poly_p<T, Degree, NbPrimesP>;
        using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

        std::tie(alphasPS, bucketSizePS) = math::balanced_bucket_n_primes(config.pbs(), primesPS, config.tauLimitBit(), 1);

        std::tie(alphasCAN, bucketSizeCAN) = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());

        std::vector<size_t> bucketSizeGCD;
        std::tie(alphasGCD, bucketSizeGCD) = math::fixed_bucket_n_primes(3*config.pbs()+210, primesGCD, config.tauLimitBit());

        // set equal buckets for bucketSizeCan and bucketSizeGCD
        bucketSizeCAN_value = lrint(floor(double(primesCAN) / double(alphasCAN.size())));
        for (int i = 0; i < bucketSizeCAN.size(); ++i) {
            bucketSizeCAN[i] = bucketSizeCAN_value;
        }
        bucketSizeGCD_value = lrint(floor(double(primesGCD) / double(alphasGCD.size())));
        for (int i = 0; i < bucketSizeGCD.size(); ++i) {
            bucketSizeGCD[i] = bucketSizeGCD_value;
        }

        if (config.protocolMode() == ProtocolMode::RECORD) {
            // get from parties

            TripleVector _tv;
            auto maybe_record = transport.awaitAggregateVectorInput<TripleVector>
                (
                 MessageType::RECORD_PROTOCOL_SHARES,
                 ZeroMQCoordinatorTransport::ids,
                 [](const TripleVector& a, const TripleVector& b) {
                 return a;
                 },
                 _tv
                );

            if (hasError(maybe_record)) {
                DBG("Failed to receive recorded shares");
                return getError(maybe_record);
            }
            std::tie(std::ignore, record_protocol_shares) = getResult(maybe_record);

            transport.broadcast(MessageType::RECORD_PROTOCOL_RESPONSE, ZeroMQCoordinatorTransport::ids, int(1));
        } else if (config.protocolMode() == ProtocolMode::REPLAY) {

            loadProtocolRecord();
            LOG(INFO) << ">>>>>>>>  Replaying experiment";
            auto maybe_replay = transport.awaitAggregateVectorInput<int>
                (
                 MessageType::REPLAY_PROTOCOL_SHARES,
                 ZeroMQCoordinatorTransport::ids,
                 [](const int a, const int& b) {
                 return a;
                 },
                 0
                );

            LOG(INFO) << "Received from all replaying parties...";
            if (hasError(maybe_replay)) {
                LOG(INFO) << "Failed to receive all replays";
                return getError(maybe_replay);
            }

            LOG(INFO) << "Sending out shares... " << record_protocol_shares.size();
            size_t socketIndex = 0;
            for ( auto it = record_protocol_shares.begin(); it != record_protocol_shares.end(); ++it ) {
                LOG(INFO) << "Sending shares " << socketIndex << " to party " << ZeroMQCoordinatorTransport::ids[socketIndex];
                transport.send(MessageType::REPLAY_PROTOCOL_RESPONSE, ZeroMQCoordinatorTransport::ids[socketIndex], it->second);
                socketIndex++;
            }
        }

        // Aggregate a matrix of `a` shares and broadcast the sum
        std::pair<Q, Q> acc = {0, 0};
        auto maybeAS = transport.awaitAggregateVectorInput<std::pair<Q, Q>>(
            MessageType::ENCRYPTED_X_SHARES,
            ZeroMQCoordinatorTransport::ids,
            lattice::pair_add<Q>,
            acc
        );

        if (hasError(maybeAS)) { return getError(maybeAS); }

        std::tie(xsum_first, g_enc_x_shares) = getResult(maybeAS);

        DBG("Received Enc(x) matrix.");
        DBG("Communication cost before = " << transport.communicationCost);
        transport.broadcast(MessageType::ENCRYPTED_X_VALUE, ZeroMQCoordinatorTransport::ids, xsum_first);
        DBG("Communication cost after = " << transport.communicationCost);
        xsum_first.first.invntt_pow_invphi();
        xsum_first.second.invntt_pow_invphi();
        // The parties will take the previous sum and multiply in their shares of
        // another random number `b_j,i,t`. These form shares of `e_i,t` which we
        // sum then broadcast for partial decryption
        auto maybeES = transport.awaitAggregateVectorInput<std::pair<Q, Q>>(
            MessageType::ENCRYPTED_XY_PLUS_Z_SHARES,
            ZeroMQCoordinatorTransport::ids,
            lattice::pair_add<Q>,
            std::pair<Q, Q>{0, 0}
        );

        if (hasError(maybeES)) { return getError(maybeES); }

        std::tie(xyz_sum, xsum_final_shares) = getResult(maybeES);

        DBG("Received Enc(e_i,t) matrix.");
        transport.broadcast(MessageType::ENCRYPTED_XY_PLUS_Z_VALUE, ZeroMQCoordinatorTransport::ids, xyz_sum);

        xyz_sum.first.invntt_pow_invphi();
        xyz_sum.second.invntt_pow_invphi();
        // Now we finish the decryption, giving us a set of e_i,t values that
        // we can use for validating the shares by gcd.
        auto maybeEit_poly = transport.awaitAggregateVectorInput<Q>(
            MessageType::PARTIAL_XY_MINUS_Z_SHARES,
            ZeroMQCoordinatorTransport::ids,
            std::plus<Q>(),
            Q{0}
        );

        if (hasError(maybeEit_poly)) { return getError(maybeEit_poly); }

        Q eit_poly;
        std::tie(eit_poly, partial_xyz_shares) = getResult(maybeEit_poly);

        // convert eit_poly to eit
        std::vector<mpz_class> tauVector(Degree);
        {
            int tvk = 0;


            // pre-sieving
            for (size_t i = 0; i < alphasPS.size(); ++i) {
                for (size_t j = 0; j < bucketSizePS[i]; ++j) {
                    tauVector[tvk] = alphasPS[i];
                    tvk++;
                }
            }

            DBG("tvk << " << tvk);
            DBG("primesPS << " << primesPS);
            assert(tvk == primesPS);
            tvk = primesPS;
            // candidate generation
            index_candidates.resize(Degree);
            for (size_t j = 0; j < bucketSizeCAN_value; ++j) {
                for (size_t i = 0; i < alphasCAN.size(); ++i) {
                    index_candidates[tvk] = j;
                    tauVector[tvk] = alphasCAN[i];
                    tvk++;
                }
            }

            // GCD
            for (size_t j = 0; j < bucketSizeGCD_value; ++j) {
                for (size_t i = 0; i < alphasGCD.size(); ++i) {
                    if (j < 1) {
                        index_candidates[tvk] = -2;
                    } else {
                        index_candidates[tvk] = -3;
                    }
                    tauVector[tvk] = alphasGCD[i];
                    tvk++;
                }
            }

            // fill the rest with dummy values from pre-sieving
            for (; tvk < Degree; ++tvk) {
                    tauVector[tvk] = alphasPS[0];
            }
        }

        DBG("before eval poly");

        auto c_mpz_t = e.eval_poly(eit_poly, tauVector);

        DBG("after eval poly");
        
        // Need to compute presieve flags
        std::vector<mpz_class> c;
        std::tie(eitPS, flags_can, c) = computeEITAndFlags(c_mpz_t);
        
        transport.broadcast(MessageType::PS_SIEVING_FLAGS, ZeroMQCoordinatorTransport::ids, flags_can);
        transport.writeToTranscript(MessageType::PS_SIEVING_FLAGS, flags_can);

        std::vector<mpz_class> c_can(primesCAN);

        eit_prunedPS = pruneAndReorderEIT(eitPS, flags_can, alphasPS.size());
        auto alphasPS_tick = alphasPS;

        alphasPS_tick.insert(alphasPS_tick.begin(), mpz_class(4));
        for (int i = primesPS; i < primesPS + primesCAN; ++i) {
            c_can[i - primesPS] = c[i];
        }

        std::vector<mpz_class> c_gcd(primesGCD);
        for (int i = primesPS + primesCAN; i < primesPS + primesCAN + primesGCD; ++i) {
            c_gcd[i - primesPS - primesCAN] = c[i];
        }

        return std::make_tuple(eit_prunedPS, alphasPS_tick, c_can, c_gcd);
    }

    expected<std::vector<mpz_class>>
    hostModulusCandidates (
        ZeroMQCoordinatorTransport& transport,
        std::vector<std::vector<mpz_class>> eit_pruned,
        std::vector<mpz_class> alphasPS_tick,
        std::vector<mpz_class> c_can) {

        DBG("Hosting modulus candidates");
        //auto [alphasCAN, _bsz] = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());
        auto bucketSize = lrint(floor(double(primesCAN) / double(alphasCAN.size())));
        if (bucketSize > eit_pruned.size()) {
            bucketSize = eit_pruned.size();
        }

        // combine alphas into alphas_combined
        std::vector<mpz_class> alphas_combined(alphasCAN.size() + alphasPS_tick.size());
        int aci = 0;
        for (size_t i = 0; i < alphasCAN.size(); ++i) {
            alphas_combined[aci] = alphasCAN[i];
            ++aci;
        }
        for (size_t i = 0; i < alphasPS_tick.size(); ++i) {
            alphas_combined[aci] = alphasPS_tick[i];
            ++aci;
        }

        DBG("Aggregate ax and by shares");
        // Step. 5 Parties send pair (ai-xi) , (bi-yi) to the coordinator and get it aggregated, so they learn a-x and b-y
        using PairOfVec = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>;
        auto acc_pair = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>(std::vector<mpz_class>(c_can.size()), std::vector<mpz_class>(c_can.size()));
        auto maybe_ax_by_sum = transport.awaitAggregateVectorInput<PairOfVec>
            (
            MessageType::AX_BY_SHARES,
            ZeroMQCoordinatorTransport::ids,
            [](const PairOfVec& a, const PairOfVec& b) {
            assert(a.first.size() == b.second.size());
            assert(a.second.size() == b.first.size());

            std::vector<mpz_class> rf(a.first.size());
            std::vector<mpz_class> rs(a.second.size());

            for(size_t i = 0; i < a.first.size(); ++i) {
                rf[i] = a.first[i] + b.first[i];
                rs[i] = a.second[i] + b.second[i];
            }
            return std::pair{rf, rs};
            },
            acc_pair
        );

        if (hasError(maybe_ax_by_sum)) {
            DBG("Failed to receive ax and by");
            return getError(maybe_ax_by_sum);
        }
        PairOfVec ax_by_sum;
        std::tie(ax_by_sum, ax_by_sum_shares) = getResult(maybe_ax_by_sum);

        {
            int k = 0;
            for (int j = 0; j < bucketSize; ++j) {
                for (int i = 0; i < alphasCAN.size(); ++i) {
                    ax_by_sum.first[k] = ax_by_sum.first[k] % alphasCAN[i];
                    ax_by_sum.second[k] = ax_by_sum.second[k] % alphasCAN[i];
                    ++k;
                }
            }
        }
        ax_modulus = ax_by_sum.first;
        by_modulus = ax_by_sum.second;
        DBG("Broadcast ax and by values");
        transport.broadcast(MessageType::AX_BY_VALUE, ZeroMQCoordinatorTransport::ids, ax_by_sum);


        // Step 6. Parties send (a-x)bi + (b-y)ai  - zi to the coordinator (special party sends a-x)bi+(b-y)ai - (a-x)(b-y) - zi - c) and get it aggregated.
        auto maybe_ab = transport.awaitAggregateVectorInput<std::vector<mpz_class>>
            (
             MessageType::AXB_MINUS_BYA_SHARES,
             ZeroMQCoordinatorTransport::ids,
             [](const std::vector<mpz_class>& a, const std::vector<mpz_class>& b) {
             assert(a.size() == b.size());

             std::vector<mpz_class> result(a.size());

             for(size_t i = 0; i < a.size(); ++i) {
             result[i] = a[i] + b[i];
             }
             return result;
             },
             std::vector<mpz_class>(c_can.size())
         );

        if (hasError(maybe_ab)) {
            DBG("Failed to receive ab");
            return getError(maybe_ab);
        }
        std::vector<mpz_class> ab;
        std::tie(ab, axby_modulus) = getResult(maybe_ab);

        DBG("Aggregate ab values");

        {
            int k = 0;
            for (int j = 0; j < bucketSize; ++j) {
                for (int i = 0; i < alphasCAN.size(); ++i) {
                    ab[k] = math::mod(ab[k] + c_can[k] - ax_by_sum.first[k] * ax_by_sum.second[k], alphasCAN[i]);
                    ++k;
                }
            }
        }

        DBG("Reconstruct candidates");
        std::vector<mpz_class> candidates_raw(bucketSize);
        DBG("bucketSize = " << bucketSize);
        DBG("eit_pruned.size = " << eit_pruned.size());

        int k = 0;
        auto pq_size = bucketSize;
        if (pq_size > eit_pruned.size()) {
            pq_size = eit_pruned.size();
        }
        for (int i = 0; i < pq_size; ++i) {

            std::vector<mpz_class> x(alphasCAN.size() + alphasPS_tick.size());

            for (int zz = 0; zz < alphasCAN.size(); ++zz) {
                x[zz] = mpz_class(ab[i * alphasCAN.size() + zz]);
            }
            for (int zz = 0; zz < alphasPS_tick.size(); ++zz) {
                x[zz + alphasCAN.size()] = mpz_class(eit_pruned[i][zz]);
            }

            std::vector<mpz_class> coefs;
            candidates_raw[k] = math::crt_reconstruct(x, coefs, alphas_combined);
            for(int zz=1; zz < 127; ++zz)
            {
                mpz_t result;
                mpz_init(result);
                mpz_class prime = boost::math::prime(zz);
                mpz_gcd(result, candidates_raw[k].get_mpz_t(),prime.get_mpz_t());
                if(mpz_cmp_ui(result,1) != 0)
                {
                    LOG(INFO) << "Candidate[" << k << "] = " << candidates_raw[k] << " is divisible by " << boost::math::prime(zz);
                    assert(false);
                }
                mpz_clear(result);
            }
            k++;
        }

        DBG("k = " << k);

        // Compute pq_size
        std::vector<mpz_class> candidates(k);
        for (int i = 0; i < k; ++i) {
            candidates[i] = candidates_raw[i];
        }

        // Then we broadcast the candidates to the clients
        DBG("Sending modulus candidates...");
        transport.broadcast(MessageType::MODULUS_CANDIDATE, ZeroMQCoordinatorTransport::ids, candidates);
        transport.writeToTranscript(MessageType::MODULUS_CANDIDATE, candidates);
        candidatesCAN = candidates;
        int number_responded = 0;
        transport.awaitAggregateVectorInput<int>(MessageType::SYNCHRONIZE_NOW
                , ZeroMQCoordinatorTransport::ids
                , std::plus<int>()
                , number_responded);
        return candidates;
    }

    expected<Unit> hostRegistration (ZeroMQCoordinatorTransport& transport) {
        auto result = transport.awaitRegistration();

        if (hasError(result)) { return getError(result); }

        transport.broadcast(
            MessageType::PROTOCOL_CONFIG,
            ZeroMQCoordinatorTransport::ids,
            config
        );
        return Unit{};
    }

    std::vector<std::vector<mpz_class>> pruneAndReorderEIT (
            const std::vector<mpz_class>& eit,
            const std::vector<int>& flags,
            int numAlphas) {

        assert (eit.size() == flags.size());

        // Now we can construct a result matrix
        std::vector<std::vector<mpz_class>> result;

        // And finally we fill the result matrix compactly with valid shares from `as`
        int k = 0;
        size_t min_size = 0;
        for (size_t c = 0; c < numAlphas; ++c) {

            std::vector<mpz_class> col;
            for (int r = 0; r < bucketSizePS[c]; ++r) {
                if (flags[k] == 1) {
                    col.push_back(eit[k]);
                }
                k++;
            }
            result.push_back(col);
            if (min_size == 0) { min_size = col.size(); }
            if (min_size > col.size()) { min_size = col.size(); }
        }

        std::vector<mpz_class> col(bucketSizePS[0], mpz_class(1));
        result.insert(result.begin(), col);

        minRowSize = min_size;
        // resize and transpose
        std::vector<std::vector<mpz_class>> transposed(min_size);

        for (size_t i = 0; i < result.size(); ++i) {
            for (size_t j = 0; j < min_size; ++j) {
                transposed[j].push_back(result[i][j]);
            }
        }

        return transposed;
    }


    expected<boost::dynamic_bitset<>> hostGCDandJacobiTest (
        ZeroMQCoordinatorTransport& transport,
        const std::vector<mpz_class>& c_gcd,
        const std::vector<mpz_class>& candidates) {


        DBG("hostGCDTest candidates.size() = " << candidates.size());
        const int bucketSize = lrint(floor(double(primesGCD) / double(alphasGCD.size())));

        DBG("Broadcast begin GCD" << candidates.size());
        DBG("c_gcd.size()" << c_gcd.size());
        // First, we broadcast an ask for shares of a random value
        transport.broadcast(MessageType::GCD_RAND_SHARES, ZeroMQCoordinatorTransport::ids);


        DBG("Aggregate ax and by shares");
        // Step. 5 Parties send pair (ai-xi) , (bi-yi) to the coordinator and get it aggregated, so they learn a-x and b-y

        auto acc_pair = std::pair{std::pair{std::vector<mpz_class>(c_gcd.size()),
                                            std::vector<mpz_class>(c_gcd.size())},
                                  mpz_class(0)};

        auto maybe_result = transport.awaitAggregateVectorInput<PairOfVecGCD>
            (
             MessageType::GCD_AX_BY_SHARES,
             ZeroMQCoordinatorTransport::ids,
             [](const PairOfVecGCD& a, const PairOfVecGCD& b) {
             auto [pair_a, at] = a;
             auto [pair_b, bt] = b;
             auto [af, as] = pair_a;
             auto [bf, bs] = pair_b;
             assert(af.size() == bf.size());
             assert(as.size() == bs.size());

             std::vector<mpz_class> rf(af.size());
             std::vector<mpz_class> rs(as.size());

             for(size_t i = 0; i < af.size(); ++i) {
                 rf[i] = af[i] + bf[i];
                 rs[i] = as[i] + bs[i];
             }

             mpz_class rt = at ^ bt;

             return std::pair{std::pair{rf, rs}, rt};
             },
             acc_pair
            );

        if (hasError(maybe_result)) {
            return getError(maybe_result);
        }
        PairOfVecGCD result_value;
        std::tie(result_value, axby_sharesGCD) = getResult(maybe_result);

        auto [pair_axby, g] = result_value;
        gammaSeed = g;

        auto [ax_sum, by_sum] = pair_axby;

        {
            int k = 0;
            for (int j = 0; j < bucketSize; ++j) {
                for (int i = 0; i < alphasGCD.size(); ++i) {
                    ax_sum[k] = ax_sum[k] % alphasGCD[i];
                    by_sum[k] = by_sum[k] % alphasGCD[i];
                    ++k;
                }
            }
        }

        gammaSeedJacobiAndGCD = gammaSeed;
        transport.broadcast(MessageType::AX_BY_VALUE, ZeroMQCoordinatorTransport::ids, std::pair{std::pair{ax_sum, by_sum}, gammaSeed});
        axGCD = ax_sum;
        byGCD = by_sum;

        // Step 6. Parties send (a-x)bi + (b-y)ai  - zi to the coordinator (special party sends a-x)bi+(b-y)ai - (a-x)(b-y) - zi - c) and get it aggregated.
        auto gcd_jacobi_pair = std::pair{std::vector<mpz_class>(c_gcd.size()),
                                         std::vector<mpz_class>(config.lambda() * candidates.size(), mpz_class(1))};
        auto maybe_zcrtsggs = transport.awaitAggregateVectorInput<std::pair<std::vector<mpz_class>, std::vector<mpz_class>>>
            (
             MessageType::AXB_MINUS_BYA_SHARES,
             ZeroMQCoordinatorTransport::ids,
             [](const std::pair<std::vector<mpz_class>, std::vector<mpz_class>>& pa,
                const std::pair<std::vector<mpz_class>, std::vector<mpz_class>>& pb) {
             auto [a, expGa] = pa;
             auto [b, expGb] = pb;

             assert(expGa.size() == expGb.size());
             std::vector<mpz_class> expResult(expGb.size());
             for (size_t i = 0; i < expResult.size(); ++i) {
                expResult[i] = expGa[i] * expGb[i];
             }
             assert(a.size() == b.size());

             std::vector<mpz_class> result(a.size());

             for(size_t i = 0; i < a.size(); ++i) {
             result[i] = a[i] + b[i];
             }
             return std::pair{result, expResult};
             },
             gcd_jacobi_pair
         );

        if (hasError(maybe_zcrtsggs)) { return getError(maybe_zcrtsggs); }

        std::pair<std::vector<mpz_class>, std::vector<mpz_class>> zcrt_value;
        std::tie(zcrt_value, zcrt_shares) = getResult(maybe_zcrtsggs);
        std::vector<mpz_class> zCRTs;

        std::tie(zCRTs, ggsGCD) = zcrt_value;

        DBG("Aggregate ab values");

        {
            int k = 0;
            for (int j = 0; j < bucketSize; ++j) {
                for (int i = 0; i < alphasGCD.size(); ++i) {
                    zCRTs[k] = math::mod(zCRTs[k] + c_gcd[k] - ax_sum[k] * by_sum[k], alphasGCD[i]);
                    ++k;
                }
            }
        }
        // reconstruct z from CRT representation
        std::vector<mpz_class> zs(candidates.size());
        int zc = 0;
        DBG("candidates.size() = " << candidates.size());
        for (int i = 0; i < candidates.size(); ++i) {

            // TODO: extract into a method
            std::vector<mpz_class> x(alphasGCD.size());
            for (int zz = 0; zz < alphasGCD.size(); ++zz) {
                x[zz] = zCRTs[zc];
                ++zc;
            }

            std::vector<mpz_class> coefs;
            zs[i] = math::crt_reconstruct(x, coefs, alphasGCD);
        }


        // GCD: compute discard vector
        DBG("Received `z`, continuing GCD test.");
        DBG("zs.size()" << zs.size());
        boost::dynamic_bitset<> discardGCD (zs.size());
        {
            {
                mpz_t one;
                mpz_init(one);
                mpz_set_ui(one, 1);

                mpz_t gcdResult;
                mpz_init(gcdResult);
                for (int i = 0; i < zs.size(); ++i) {
                    DBG("i = " << i);
                    const mpz_class& N = candidates[i];
                    DBG("N = " << N);

                    const mpz_class z = zs[i] % N;

                    DBG("z = " << z);
                    DBG("before mpz_gcd");
                    mpz_gcd(gcdResult, N.get_mpz_t(), z.get_mpz_t());

                    if (mpz_cmp(gcdResult, one) != 0) {
                        discardGCD[i] = 1;
                    }

                }
                mpz_clear(gcdResult);
                mpz_clear(one);
            }
        }

        assert(zs.size() == candidates.size());

        // Jacobi: compute discard vector
        boost::dynamic_bitset<> discardJacobi (candidates.size());
        {

            DBG("Taking gamma values modulo N");
            for (int i = 0; i < candidates.size(); ++i) {
                for (int j = 0; j < config.lambda(); ++j) {
                    const mpz_class& N = candidates[i];
                    mpz_class& gg = ggsGCD[i * config.lambda() + j];
                    mpz_class x = gg % N;

                    DBG("x = " << x);
                    DBG("N = " << N);
                    DBG("ggs[" << i << "] =" << ggsGCD[i * config.lambda() + j]);

                    // Elimination
                    if (x != 1 && x != (N - 1)) {
                        discardJacobi[i] = 1;
                    }
                }
            }
        }

        // Final step: merge discardGCD and discardJacobi
        boost::dynamic_bitset<> discard (candidates.size());
        DBG("hostGCDandJacobiTest candidates.size() = " << candidates.size());
        for (int i = 0; i < candidates.size(); ++i) {
            DBG("discardGCD[" << i << "] = " << discardGCD[i]);
            DBG("discardJacobi[" << i << "] = " << discardJacobi[i]);
            discard[i] = discardGCD[i] | discardJacobi[i];
        }

        return discard;
    }


    expected<boost::dynamic_bitset<>> hostJacobiProtocol (
        ZeroMQCoordinatorTransport& transport,
        const std::vector<mpz_class>& candidates
    ) {
        // First, we broadcast an ask for shares of a gamma value
        transport.broadcast(MessageType::GAMMA_SHARES, ZeroMQCoordinatorTransport::ids);

        // Aggregate and broadcast gamma seed
        auto maybeGammaSeed = transport.awaitAggregateVectorInput<mpz_class>(
            MessageType::GAMMA_RANDOM_SEED_SHARES,
            ZeroMQCoordinatorTransport::ids,
            [](const mpz_class& a,
               const mpz_class& b) { return a ^ b; },
            mpz_class(0)
        );

        if (hasError(maybeGammaSeed)) { return getError(maybeGammaSeed); }

        
        std::tie(gammaSeedAccum, std::ignore) = getResult(maybeGammaSeed);

        transport.broadcast(MessageType::GAMMA_RANDOM_SEED_VALUE, ZeroMQCoordinatorTransport::ids, gammaSeedAccum);

        // Now the input we'll receive is a multiplied value into which
        // we'll finally multiply g^(N+1) for each gamma/modulus pair.
        // If the result is not +/-1 we discard the candidate and move on.
        auto maybe_ggs = transport.awaitAggregateVectorInput<std::vector<mpz_class>>(
            MessageType::EXPONENTIATED_GAMMA_VALUE,
            ZeroMQCoordinatorTransport::ids,
            [](const std::vector<mpz_class>& a, const std::vector<mpz_class>& b) {
                std::vector<mpz_class> c(b.size());
                for (size_t i = 0; i < c.size(); ++i) {
                    c[i] = a[i] * b[i];
                }
                return c;
            },
            std::vector<mpz_class>(candidates.size(), mpz_class(1))
        );

        if (hasError(maybe_ggs)) { return getError(maybe_ggs); }

        std::vector<mpz_class> ggs;
        std::tie(ggs, ggsShares) = getResult(maybe_ggs);

        boost::dynamic_bitset<> discard (candidates.size());

        int discarded = 0;
        for (int i = 0; i < candidates.size(); ++i) {
            const mpz_class& N = candidates[i];
            mpz_class& gg = ggs[i];
            mpz_class x = gg % N;

            // Elimination
            if (x != 1 && x != (N - 1)) {
                discard[i] = 1;
                discarded++;
            }
        }
        LOG(INFO) << "After first Jacobi iteration " << candidates.size() - discarded << " candidates survived";

        return discard;
    }

    expected<Unit> hostThroughputTest(
        ZeroMQCoordinatorTransport& transport,
        throughput_test_config tconfig,
        size_t throughput_cutoff
    ) {
        DBG("Host throughput test");
        auto survivor = transport.hostThroughputTest(tconfig, [&](size_t t) { return t > throughput_cutoff; });

        // No party survive the test, quit without restarting
        if (survivor.size() == 0) {
            LOG(ERROR) << "No party survives, Quiting.";
            return Error::TOO_FEW_PARTIES;
        }
        // Some parties survive but still too few to restart the protocol
        else if (survivor.size() < kMinAmountOfPartiesRequired) {
            LOG(INFO) << "Too few survivors, killing the protocol and quiting.";
            // broadcast only to survivors because we don't start registration
            transport.broadcast(MessageType::PROTOCOL_KILLED, survivor);
            return Error::TOO_FEW_PARTIES;
        }

        // success
        return Unit{};
    }

    expected<Unit> hostHelper(ZeroMQCoordinatorTransport& transport,
        lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& enc,
        bool needRegistration,
        int *restartRemaining
    ) {
        // base case
        if (*restartRemaining <= 0) {
            return Error::TOO_MANY_RESTART;
        }

        auto result = registerAndHostRSACeremony(transport, enc, restartRemaining, needRegistration);
        if (!hasError(result)) {
            // No errors, we are good to go
            return Unit{};
        }

        // There are some errors...
        auto err = getError(result);
        if (err == Error::TOO_FEW_PARTIES) {
            transport.broadcast(MessageType::PROTOCOL_KILLED, ZeroMQCoordinatorTransport::ids);
            LOG(INFO) << "Too few parties left, killing the protocol and quiting.";
            return Error::TOO_FEW_PARTIES;
        }
        else if (err == Error::MODULUS_NOT_FOUND) {
            // Since there are no kick out, restart without re-registration
            (*restartRemaining)--;
            return hostHelper(transport, enc, false, restartRemaining);
        }
        else if (err == Error::TIMED_OUT 
                || err == Error::OUT_OF_SYNC 
                || err == Error::DESERIALIZE_FAIL 
                || err == Error::RESTART)
        {
            // Kick out some parties and restart
            auto parties = ZeroMQCoordinatorTransport::ids.size();
            LOG(INFO) << "Restarting the protocol with " << parties << " parties";
            clearPublicData();
            config.numParties() = parties;
            transport.parties() = parties;

            (*restartRemaining)--;

            transport.broadcast(MessageType::PROTOCOL_RESTART, ZeroMQCoordinatorTransport::ids);
            return hostHelper(transport, enc, true, restartRemaining);
        }
        else {
            LOG(INFO) << "Coordinator got impossible error: " << showError(err);
            return Error::UNKNOWN_ERROR;
        }

        LOG(INFO) << "Too many restart, coordinator aborting.";
        transport.broadcast(MessageType::PROTOCOL_KILLED, ZeroMQCoordinatorTransport::ids);
        return Error::TOO_MANY_RESTART;
    }

    /** Host the RSA MPC ceremony. */
    expected<Unit> host(
        ZeroMQCoordinatorTransport& transport,
        int *restartRemaining
    ) {
        // After running the throughput test,
        // now we have enough parties to start the protocol

        // Assumption:
        //    Whenever error happened the transport will update `ids` then proxy the error

        auto e = lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>(config);
        bool needRegistration = true;

        return hostHelper(transport, e, needRegistration, restartRemaining);
    }

    expected<Unit> registerAndHostRSACeremony(
        ZeroMQCoordinatorTransport& transport, 
        lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& enc,
        int *restartRemainging,
        bool registration = true
    ) {
        if (registration) {
            auto result = hostRegistration(transport);
            if (hasError(result)) { return getError(result); }
        }

        // Start timer
        transport.myTimers.initialize_timer();
        transport.myTimers.begin(1, "RSA Ceremony", transport.communicationCost);

        // Kickoff the protocol
        auto result = host_rsa_ceremony(transport, enc);
        if (!hasError(result)) {
            // RSA protocol success
            transport.myTimers.end(1,"RSA Ceremony",transport.communicationCost);
            return Unit{};
        }

        // Cleanup and proxy the error
        transport.myTimers.reset();
        return getError(result);
    }

    expected<int> host_rsa_ceremony(
        ZeroMQCoordinatorTransport& transport,
        lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& e
        ) {

        DBG("DBG001 GENERATE KEYS");
        transport.myTimers.begin(2,"Key Generation",transport.communicationCost);
        transport.myTimers.begin(2,"1.a. Overall speed", transport.communicationCost);
        auto keygen_success = hostGenerateKeyPair(transport);

        if (hasError(keygen_success)) {
            LOG(ERROR) << "Keygen timed out!";
            transport.myTimers.end(2,"Key Generation",transport.communicationCost);
            return getError(keygen_success);
        }

        DBG("DBG001 HOSTING");

        // Assign P1
        transport.broadcast(MessageType::ASSIGNMENT_P1, std::vector<SocketId>{ZeroMQCoordinatorTransport::ids[0]});
        DBG("Sending ASSIGNMENT_P1");

        // Assign Pn
        transport.broadcast(MessageType::ASSIGNMENT_PN, std::vector<SocketId>(ZeroMQCoordinatorTransport::ids.begin() + 1, ZeroMQCoordinatorTransport::ids.end()));
        DBG("Sending ASSIGNMENT_PN");
        transport.myTimers.end(2,"Key Generation",transport.communicationCost);


        bool foundModuli = false;

        // precompute bounds and M with min=4096 max=104729(10k th prime) and 2 tests
        auto [Bs, Ms] = ligero::math::compute_m_b_vec(4096, 104729, 2);
        auto nb_threads = 90;
        auto postSieve = std::bind(ligero::math::test_factorizable_threaded, std::placeholders::_1, Bs, Ms, nb_threads);

        transport.myTimers.begin(2,"Pre-Sieving",transport.communicationCost);
        auto presieve_success = hostPreSieving(transport, e);

        if (hasError(presieve_success)) {
            LOG(ERROR) << "Presieving timed out!";
            transport.myTimers.end(2,"Pre-Sieving",transport.communicationCost);
            return getError(presieve_success);
        }

        auto [eit_pruned, alphasPS_tick, c_can, c_gcd] = getResult(presieve_success);

        transport.myTimers.end(2,"Pre-Sieving",transport.communicationCost);

        DBG("Beginning Candidate generation");
        transport.myTimers.begin(2,"Generate Candidates",transport.communicationCost);
        auto candidates_success = hostModulusCandidates (transport, eit_pruned, alphasPS_tick, c_can);

        if (hasError(candidates_success)) {
            LOG(ERROR) << "Modulus candidate timed out!";
            transport.myTimers.end(2,"Generate Candidates",transport.communicationCost);
            return getError(candidates_success);
        }

        auto candidates = getResult(candidates_success);

        transport.myTimers.end(2,"Generate Candidates",transport.communicationCost);
        if (candidates.size() == 0) {
            LOG(INFO) << "No candidates found.";
            transport.broadcast(MessageType::NO_MODULI, ZeroMQCoordinatorTransport::ids);
            if (!verifyNoModuli(transport, e)) {
                LOG(INFO) << "Failed randomness verification";
            }
            return Error::MODULUS_NOT_FOUND;
        }
        // sanity check after modulus
        {
            DBG("Start sanity check");
            for (int i = 0; i < 100; ++i) {
                auto x = boost::math::prime(i);
                for (int j = 0; j < candidates.size(); ++j) {
                    if (candidates[j] % x == 0) {
                        DBG("candidates["<< j <<"] = " << candidates[j] << " x = " << x << "candidates[j] mod x = " << (candidates[j] % x));
                        assert(false);
                    }
                }
            }
            DBG("Passed sanity check");
        }

        boost::dynamic_bitset<> discardFlags;
        {
            DBG("Running postSieve");
            discardFlags = postSieve(std::vector<mpz_class>(candidates.data(), candidates.data() + candidates.size()));
            discardFlagsPostSieve = discardFlags;
            transport.broadcast(MessageType::POST_SIEVE, ZeroMQCoordinatorTransport::ids, discardFlags);
            transport.writeToTranscript(MessageType::POST_SIEVE, discardFlags);
            candidates = discardCandidates(candidates, discardFlags);
            DBG("Completed postSieve");
            candidatesPostSieve = candidates;
        }


        for (int i = 0; i < 2; ++i) {
            DBG("before sanity check");
            DBG("candidates[" << i <<  "] " << candidates[i]);
            assert(candidates[i] != 0);
        }
        // sanity check after postSieve
        {
            DBG("Start sanity check");
            for (int i = 0; i < 100; ++i) {
                auto x = boost::math::prime(i);
                for (int j = 0; j < candidates.size(); ++j) {
                    if (candidates[j] % x == 0) {
                        DBG("candidates["<<j<<"] = " << candidates[j] << " x = " << x << "candidates[j] mod x = " << (candidates[j] % x));
                        assert(false);
                    }
                }
            }
            DBG("Passed sanity check");
        }
        // Now we move on to our biprimality test. Here we run the
        // jacobi test multiple times, eliminating candidates after each trial

        {
            DBG("Beginning Jacobi");
            transport.myTimers.begin(2,"Jacobi Test",transport.communicationCost);
            auto jacobi_success = hostJacobiProtocol(transport, candidates);

            if (hasError(jacobi_success)) {
                LOG(ERROR) << "Jacobi timed out!";
                transport.myTimers.end(2,"Jacobi Test",transport.communicationCost);
                return getError(jacobi_success);
            }

            boost::dynamic_bitset<> discardFlags = getResult(jacobi_success);

            transport.broadcast(MessageType::DISCARD_FLAGS, ZeroMQCoordinatorTransport::ids, discardFlags);
            transport.writeToTranscript(MessageType::DISCARD_FLAGS, discardFlags);

            candidates = discardCandidates(candidates, discardFlags);
            candidatesJacobi = candidates;
            transport.myTimers.end(2,"Jacobi Test",transport.communicationCost);
        }

        if (candidates.size() == 0) {
            LOG(INFO) << "No candidates found.";
            transport.broadcast(MessageType::NO_MODULI, ZeroMQCoordinatorTransport::ids);
            if (!verifyNoModuli(transport, e)) {
                LOG(INFO) << "Failed randomness verification";
            }
            return Error::MODULUS_NOT_FOUND;
        }

        DBG("After Jacobi candidates.size() = " << candidates.size());
        //// Complete the biprimality test with the GCD test, discarding any
        //// other candidates that fail.
        transport.myTimers.begin(2,"GCD",transport.communicationCost);

        if (candidates.size() > 1) {
            candidates.resize(1);
        }
        auto gcd_jacobi_success = hostGCDandJacobiTest (transport, c_gcd, candidates);

        if (hasError(gcd_jacobi_success)) {
            LOG(ERROR) << "GCD and Jacobi timed out!";
            transport.myTimers.end(2,"GCD",transport.communicationCost);
            return getError(gcd_jacobi_success);
        }

        discardFlags = getResult(gcd_jacobi_success);
        discardFlagsJacobi = discardFlags;

        transport.broadcast(MessageType::DISCARD_FLAGS, ZeroMQCoordinatorTransport::ids, discardFlags);
        transport.writeToTranscript(MessageType::DISCARD_FLAGS, discardFlags);

        transport.myTimers.end(2,"GCD",transport.communicationCost);

        transport.myTimers.begin(2,"Discarding",transport.communicationCost);

        candidates = discardCandidates(candidates, discardFlags);

        // Here now, any modulus candidate remaining in our vector has passed
        // the test.
        if (candidates.size() == 0) {
            LOG(INFO) << "No candidates found.";
            transport.broadcast(MessageType::NO_MODULI, ZeroMQCoordinatorTransport::ids);
            if (!verifyNoModuli(transport, e)) {
                LOG(INFO) << "Failed randomness verification";
            }
            return Error::MODULUS_NOT_FOUND;
        }
        else {
            LOG(INFO) << "Found " << candidates.size() << " valid moduli:";

            for (int i = 0; i < candidates.size(); ++i){
                LOG(INFO) << candidates[i];
            }
            transport.broadcast(MessageType::FOUND_MODULI, ZeroMQCoordinatorTransport::ids);
            candidatesFinal = candidates;
            foundModuli = true;

            return 0;
        }
        transport.myTimers.end(2,"Discarding",transport.communicationCost);
        numCandidates = candidates.size();

        if (config.protocolMode() == ProtocolMode::RECORD || config.protocolMode() == ProtocolMode::REPLAY) {

            // First check that candidates match
            if (config.protocolMode() == ProtocolMode::REPLAY) {
                if (candidatesFinal[0] != candidatesToCheck[0]) {
                    LOG(FATAL) << "Recorded moduli" <<  candidatesToCheck[0] << " does not match computed moduli" << candidatesFinal[0];
                }
            }

            // Below we are going to check if our output is a product of two primes
            // There is no other way to check it than actually computing p and q as we
            // cannot possibly factor 2048 bit numbers!
            // TODO: This code has to be removed

            //transport.myTimers.begin(2,"Miller Rabin Tests",transport.communicationCost);
            using PairOfVec = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>;
            if (foundModuli) {
                auto maybe_pq = transport.awaitAggregateVectorInput<PairOfVec>(
                        MessageType::P_CLEAR_DEBUG,
                        ZeroMQCoordinatorTransport::ids,
                        [](const PairOfVec& a, const PairOfVec& b) {
                        assert(a.first.size() == b.first.size());
                        assert(a.second.size() == b.second.size());

                        std::vector<mpz_class> fsum(a.first.size());
                        std::vector<mpz_class> ssum(b.second.size());

                        for (size_t i = 0; i < a.first.size(); ++i) {
                        fsum[i] = a.first[i] + b.first[i];
                        ssum[i] = a.second[i] + b.second[i];
                        }
                        return std::pair{fsum, ssum};
                        },
                        std::pair{std::vector<mpz_class>(candidates.size(), mpz_class(0)),
                        std::vector<mpz_class>(candidates.size(), mpz_class(0))}
                        );

                if (hasError(maybe_pq)) {
                    LOG(ERROR) << "PQ timed out!";
                    return getError(maybe_pq);
                }
                auto [pq, pq_shares] = getResult(maybe_pq);

                transport.broadcast(MessageType::P_CLEAR_DEBUG, ZeroMQCoordinatorTransport::ids, int(1));

                // sanity check
                for(size_t i = 0; i < pq.first.size(); ++i)
                {
                    assert(boost::multiprecision::miller_rabin_test(MPInt(pq.first[i].get_mpz_t()), 4));
                    assert(boost::multiprecision::miller_rabin_test(MPInt(pq.second[i].get_mpz_t()), 4));
                }

                LOG(INFO) << "Prime factorization of the candidates are: ";
                if (config.protocolMode() == ProtocolMode::RECORD) {
                    recorded_pq_first.resize(pq.first.size());
                    recorded_pq_second.resize(pq.second.size());
                }
                for (size_t i = 0; i < pq.first.size(); ++i) {
                    LOG(INFO) << pq.first[i];
                    LOG(INFO) << pq.second[i];
                    if (config.protocolMode() == ProtocolMode::RECORD) {
                        recorded_pq_first[i] = pq.first[i];
                        recorded_pq_second[i] = pq.second[i];
                    }
                    if (config.protocolMode() == ProtocolMode::REPLAY) {
                        if (recorded_pq_first[i] != pq.first[i]) {
                            LOG(FATAL) << "Recorded prime factorization " << recorded_pq_first[i] 
                                << " does not match computed prime factorization " << pq.first[i];
                        }
                        if (recorded_pq_second[i] != pq.second[i]) {
                            LOG(FATAL) << "Recorded prime factorization " << recorded_pq_second[i] 
                                << " does not match computed prime factorization " << pq.second[i];
                        }
                    }
                }

                //transport.myTimers.end(2,"Miller Rabin Tests",transport.communicationCost);

                assert(pq.first.size() == candidates.size());

                for(size_t i = 0; i< pq.first.size(); ++i)
                {
                    assert(pq.first[i] * pq.second[i] == candidates[i]);
                }
                LOG(INFO) << "All candidate moduli checked.";
                DBG("Done.");

            }
        }

        //transport.myTimers.end(1,"RSA Ceremony",transport.communicationCost);
    }

    void updatePublicData(const SocketId& id, PublicData *pdata_p) {

        PublicData& pdata = *pdata_p;

        // Keygen
        {
            dataExtractQ(pdata.A, A);

            auto bi_share = bi.find(id);
            if (bi_share  == bi.end()) {
                LOG(FATAL) << "Could not find bi share for " << id;
            }
            Q bis = bi_share->second;
            bis.invntt_pow_invphi();
            dataExtractQ(pdata.bi, bis);
            dataExtractQ(pdata.b, B);
        }


        {
            auto g_enc_x_s = g_enc_x_shares.find(id);
            if (g_enc_x_s == g_enc_x_shares.end()) {
                LOG(FATAL) << "Could not find encrypted x share for " << id;
            }
            Q ci_1, ci_2;
            std::tie(ci_1, ci_2) = g_enc_x_s->second;
            ci_1.invntt_pow_invphi();
            ci_2.invntt_pow_invphi();
            dataExtractQ(pdata.ci_1, ci_1);
            dataExtractQ(pdata.ci_2, ci_2);

            dataExtractQ(pdata.c_1_prime, xyz_sum.first);
            dataExtractQ(pdata.c_2_prime, xyz_sum.second);

            auto di_share = partial_xyz_shares.find(id);
            if (di_share  == partial_xyz_shares.end()) {
                LOG(FATAL) << "Could not find partial xyz share for " << id;
            }
            Q dis = di_share->second;
            dis.invntt_pow_invphi();
            dataExtractQ(pdata.di, dis);
        }


        // Round 5
        auto xf_share = xsum_final_shares.find(id);
        if (xf_share == xsum_final_shares.end()) {
            LOG(FATAL) << "Could not find encrypted xsum_final share for " << id;
        }
        Q ci_1_prime, ci_2_prime;
        std::tie(ci_1_prime, ci_2_prime) = xf_share->second;
        ci_1_prime.invntt_pow_invphi();
        ci_2_prime.invntt_pow_invphi();
        dataExtractQ(pdata.ci_1_prime, ci_1_prime);
        dataExtractQ(pdata.ci_2_prime, ci_2_prime);

        dataExtractQ(pdata.c_1, xsum_first.first);
        dataExtractQ(pdata.c_2, xsum_first.second);

        // Round 11 & 12
        dataExtractVectorCRT<NbPrimesQ>(pdata.axGCD, axGCD);
        dataExtractVectorCRT<NbPrimesQ>(pdata.byGCD, byGCD);

        auto axbyGCD_share = zcrt_shares.find(id);
        if (axbyGCD_share == zcrt_shares.end()) {
            LOG(FATAL) << "Could not find encrypted axbyGCD share for " << id;
        }
        std::vector<mpz_class> axbyGCD = axbyGCD_share->second.first;
        dataExtractVectorCRT<NbPrimesQ>(pdata.axbyGCD, axbyGCD);

        // Misc
        dataExtractVectorCRT<NbPrimesQ>(pdata.ps, alphasPS);
        dataExtractVectorCRT<NbPrimesQ>(pdata.cans, alphasCAN);
        dataExtractVectorCRT<NbPrimesQ>(pdata.gcds, alphasGCD);

        pdata.degree = Degree;
        pdata.p = NbPrimesP;
        pdata.q = NbPrimesQ;
        pdata.sigma = config.sigma();
        pdata.lambda = config.lambda();
        pdata.tau_limit_bit = config.tauLimitBit();

        {
            std::vector<mpz_class> _alphas(alphasPS.size() + 1);
            std::vector<mpz_class> moduli(_alphas.size());

            for (size_t i = 0; i < alphasPS.size(); ++i) {
                _alphas[i] = alphasPS[i];
            }
            _alphas[_alphas.size() - 1] = mpz_class(4);

            for (size_t i = 0; i < minRowSize; ++i) {
                //prime_shares.push_back(math::crt_reconstruct(_coeffs, moduli(), _alphas));
                mpz_class p = 1;
                for (size_t i = 0; i < _alphas.size(); ++i) {
                    p *= _alphas[i];
                }

                for (size_t i = 0; i < _alphas.size(); ++i) {
                    mpz_class pa = p / _alphas[i];
                    mpz_class x = math::mod(pa, _alphas[i]);
                    moduli[i] = pa * math::mod_inverse(x, _alphas[i]);
                }
            }

            size_t idx = 0;
            for (; idx < candidatesCAN.size(); idx++) {
                if (candidatesCAN[idx] == candidatesFinal[0]) {
                    break;
                }
            }
            auto candidateIndices = getCandidateIndices(idx, index_candidates);

            std::vector<mpz_class> alphasPSprod_CAN, moduli_CAN;

            mpz_class alphasPSprod = mpz_class(4);
            for (size_t j = 0; j < alphasPS.size(); ++j) {
                alphasPSprod *= alphasPS[j];
            }

            auto ax_by_sum_share = ax_by_sum_shares.find(id);
            if (ax_by_sum_share == ax_by_sum_shares.end()) {
                LOG(FATAL) << "Could not find encrypted ax by share GCD share for " << id;
            }
            auto [ax_shareCAN, by_shareCAN] = ax_by_sum_share->second;

            {

                std::vector<size_t> axbyindices;
                std::vector<mpz_class> ax_shares, by_shares;
                size_t alphaCANidx = 0;
                for (auto j : candidateIndices.can) {
                    std::vector<mpz_class> moduliCAN(moduli.size());
                    mpz_class alphasPSprodCAN = alphasPSprod % alphasCAN[alphaCANidx];

                    for (size_t i = 0; i < moduli.size(); i++) {
                        moduliCAN[i] = moduli[i] % alphasCAN[alphaCANidx];
                        moduli_CAN.push_back(moduliCAN[i]);
                    }

                    mpz_class axshares = ax_shareCAN[idx * alphasCAN.size() + alphaCANidx];
                    mpz_class byshares = by_shareCAN[idx * alphasCAN.size() + alphaCANidx];

                    for (size_t i = 0; i < NbPrimesQ; i ++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                        ax_shares.push_back(positiveRemainder(axshares, prime));
                        by_shares.push_back(positiveRemainder(byshares, prime));
                    }

                    size_t pos = alphasCAN.size() * idx + alphaCANidx;
                    axbyindices.push_back(pos);
                    alphasPSprod_CAN.push_back(alphasPSprodCAN);
                    alphaCANidx++;
                }

                auto axby_modulus_share = axby_modulus.find(id);
                if (axby_modulus_share  == axby_modulus.end()) {
                    LOG(FATAL) << "Could not find axby share for " << id;
                }
                std::vector<mpz_class> axby_modulus_value = axby_modulus_share->second;

                // axby
                dataExtractSubVector<NbPrimesQ>(pdata.ax, ax_modulus, axbyindices);
                dataExtractSubVector<NbPrimesQ>(pdata.by, by_modulus, axbyindices);

                dataExtractSubVector<NbPrimesQ>(pdata.axby, axby_modulus_value, axbyindices);

                // Extract data for CAN
                dataExtractVectorCRT<NbPrimesQ>(pdata.prodcans, alphasPSprod_CAN);
                dataExtractVectorCRT<NbPrimesQ>(pdata.coefsCAN,	moduli_CAN);
                dataExtractVector(pdata.ax_shares,  ax_shares);
                dataExtractVector(pdata.by_shares,  by_shares);

                // Extract gamma for sigma protocol verification


            }

            auto ax_by_shareGCD = axby_sharesGCD.find(id);
            if (ax_by_shareGCD == axby_sharesGCD.end()) {
                LOG(FATAL) << "Could not find encrypted ax by share GCD share for " << id;
            }
            auto [pair_axby, gammaSeed] = ax_by_shareGCD->second;
            auto [axshare, byshare] = pair_axby;


            std::vector<mpz_class> alphasPSprod_GCD, moduli_GCD, ax_shares_GCD, by_shares_GCD;
            size_t alphaGCDidx = 0;
            for (auto j : candidateIndices.gcd) {
                std::vector<mpz_class> moduliGCD(moduli.size());
                mpz_class alphasPSprodGCD = alphasPSprod % alphasGCD[alphaGCDidx];
                mpz_class axsharesGCD = axshare[alphaGCDidx];
                mpz_class bysharesGCD = byshare[alphaGCDidx];

                for (size_t i = 0; i < moduli.size(); i++) {
                    moduliGCD[i] = moduli[i] % alphasGCD[alphaGCDidx];
                    moduli_GCD.push_back(moduliGCD[i]);
                }

                alphasPSprod_GCD.push_back(alphasPSprodGCD);

                for (size_t i = 0; i < NbPrimesQ; i++) {
                    mpz_class prime = mpz_class(nfl::params<uint64_t>::P[i]);
                    ax_shares_GCD.push_back(positiveRemainder(axsharesGCD, prime));
                    by_shares_GCD.push_back(positiveRemainder(bysharesGCD, prime));
                }
                alphaGCDidx++;

                std::vector<mpz_class> feed;
                for (size_t i = 0; i < alphasGCD.size(); i++) {
                    for (size_t j = 0; j < NbPrimesQ; j++) {
                        mpz_class prime = mpz_class(nfl::params<uint64_t>::P[j]);
                        mpz_class modGCD = math::mod(candidatesFinal[0], alphasGCD[i]);
                        feed.emplace_back(positiveRemainder(modGCD,prime));
                    }
                }

                dataExtractVector(pdata.finalModuli_GCD, feed);
            }

            // Extract data for GCD
            dataExtractVectorCRT<NbPrimesQ>(pdata.coefsGCD,	moduli_GCD);
            dataExtractVectorCRT<NbPrimesQ>(pdata.prodgcds, alphasPSprod_GCD);

            dataExtractVector(pdata.by_shares_GCD,  by_shares_GCD);
            dataExtractVector(pdata.ax_shares_GCD,  ax_shares_GCD);

            // indices
            pdata.indicesPS.assign(candidateIndices.ps.begin(),candidateIndices.ps.end());
            pdata.indicesCAN.assign(candidateIndices.can.begin(),candidateIndices.can.end());
            pdata.indicesGCD.assign(candidateIndices.gcd.begin(),candidateIndices.gcd.end());

            pdata.tau = 1000;
            // update special
            pdata.special = id == ZeroMQCoordinatorTransport::ids[0];

            // update bounding
            std::vector<mpz_class> cans(alphasCAN);
            cans.insert(cans.end(), alphasPS.begin(), alphasPS.end());
            cans.insert(cans.end(), alphasGCD.begin(), alphasGCD.end());

            std::vector<mpz_class> zbounds;

            for (auto tau : cans) {
                zbounds.emplace_back(tau * tau * config.numParties() * mpz_class(2)^mpz_class(config.lambda()));
            }

            mpz_class rbound = 2*config.lambda()*2*(config.sigma()*config.numParties()*mpz_class(2)^mpz_class(64*9)*Degree);                        

            size_t size = alphasCAN.size() + alphasPS.size() + alphasGCD.size();
            std::vector<mpz_class> rbounds(size, rbound);

            auto params = {cans, cans, alphasCAN, {mpz_class(2)^mpz_class(1234)}, zbounds, rbounds};

            std::vector<mpz_class> Cs;
            for (auto param : params) {
                for (auto bound : param) {
                    auto [dd, C, D, log2D] = rephraseAs2Nct(bound);
                    pdata.log2Ds.emplace_back(log2D);

                    for (size_t modulusIdx = 0; modulusIdx < NbPrimesQ; modulusIdx++) {
                        mpz_class prime(nfl::params<uint64_t>::P[modulusIdx]);
                        mpz_class calc4 = math::mod(C, prime);
                        mpz_class calc5 = math::mod(dd, prime);

                        pdata.Cs.emplace_back(calc4.get_ui());
                        pdata.Ds.emplace_back(calc5.get_ui());
                    }
                }
            }

            // Sigma protocol
            pdata.gammaSeed = gammaSeedJacobiAndGCD;
        }
    }

    void clearPublicData() {
        axGCD.clear();
        byGCD.clear();

        bucketSizeCAN.clear();

        partial_xyz_shares.clear();

        g_enc_x_shares.clear();
        xsum_final_shares.clear();
        zcrt_shares.clear();
        axby_sharesGCD.clear();
        ax_by_sum_shares.clear();

        axby_modulus.clear();

        alphasPS.clear();
        alphasCAN.clear();
        alphasGCD.clear();

        ax_modulus.clear();
        by_modulus.clear();

        candidatesCAN.clear();
        candidatesFinal.clear();

        index_candidates.clear();
    }

    void saveProtocolRecord() {
        LOG(INFO) << "Saving protocol recording into record.data file for " << record_protocol_shares.size() << " parties";
        boost::filesystem::ofstream record_ofs;
        record_ofs.open("record.data", std::ios::out | std::ios::binary);
        boost::archive::binary_oarchive recordOA(record_ofs, boost::archive::no_header);
        recordOA << record_protocol_shares;
        recordOA << candidatesFinal;
        recordOA << recorded_pq_first;
        recordOA << recorded_pq_second;
        record_ofs.close();
    }

    void loadProtocolRecord() {
        LOG(INFO) << "Loading protocol recording from record.data file" << record_protocol_shares.size() << " parties";
        boost::filesystem::ifstream record_ifs;
        record_ifs.open("record.data", std::ios::in | std::ios::binary);
        boost::archive::binary_iarchive recordIA(record_ifs, boost::archive::no_header);
        recordIA >> record_protocol_shares;
        recordIA >> candidatesToCheck;
        recordIA >> recorded_pq_first;
        recordIA >> recorded_pq_second;
        record_ifs.close();

        // Check if protocol recorded data matches current setup
        if (record_protocol_shares.size() != ZeroMQCoordinatorTransport::ids.size()) {
            LOG(FATAL) << "The number of parties in the recorded data " << record_protocol_shares.size() << " does not match current number of parties " << ZeroMQCoordinatorTransport::ids.size();
        }
        
    }

protected:
    size_t numCandidates;
    Q A, B;
    lattice::cipher<T, Degree, NbPrimesQ>  xsum_first, xyz_sum;
    ProtocolConfig<T>& config;
    std::vector<mpz_class> axGCD, byGCD;
    std::unordered_map<SocketId, Q, boost::hash<SocketId>> bi, partial_xyz_shares;
    std::unordered_map<SocketId, std::pair<Q, Q>, boost::hash<SocketId>> g_enc_x_shares, xsum_final_shares;
    std::unordered_map<SocketId, std::pair<std::vector<mpz_class>, std::vector<mpz_class>>, boost::hash<SocketId>> zcrt_shares;
    std::unordered_map<SocketId, PairOfVecGCD, boost::hash<SocketId>> axby_sharesGCD;
    std::unordered_map<SocketId, std::pair<std::vector<mpz_class>, std::vector<mpz_class>>, boost::hash<SocketId>> ax_by_sum_shares;

    std::unordered_map<SocketId, std::vector<mpz_class>, boost::hash<SocketId>> axby_modulus;
    std::unordered_map<SocketId, std::vector<mpz_class>, boost::hash<SocketId>> ggsShares;

    std::unordered_map<SocketId, TripleVector, boost::hash<SocketId>> record_protocol_shares;

    mpz_class gammaSeed;
    std::vector<mpz_class> alphasPS, alphasCAN, alphasGCD;
    std::vector<size_t> bucketSizePS;
    std::vector<size_t> bucketSizeCAN;
    std::vector<size_t> bucketSizeGCD;

    std::vector<mpz_class> ax_modulus, by_modulus;

    std::vector<std::vector<mpz_class>> eit_prunedPS;

    std::vector<mpz_class> candidatesCAN, candidatesFinal, candidatesToCheck, recorded_pq_first, recorded_pq_second;

    std::vector<mpz_class> ggsGCD;
    std::vector<mpz_class> candidatesPostSieve;
    std::vector<mpz_class> candidatesJacobi;
    mpz_class gammaSeedAccum;
    mpz_class gammaSeedJacobiAndGCD;

    std::vector<int> flags_can;
    int bucketSizeGCD_value;
    int bucketSizeCAN_value;
    std::vector<mpz_class> eitPS;
    std::vector<size_t> index_candidates;
    size_t minRowSize = -1;

    boost::dynamic_bitset<> discardFlagsPostSieve;
    boost::dynamic_bitset<> discardFlagsJacobi;
};


} // namespace ligero
