#pragma once


#include "Transport.hpp"
#include "Common.hpp"
#include "Math.hpp"
#include <boost/multiprecision/miller_rabin.hpp>
#include "LatticeEncryption.hpp"
#include "Factoring.hpp"
#include <functional>
#include <optional>

namespace ligero
{

template <typename NumberType>
void hostRegistration (ZeroMQCoordinatorTransport& transport, const ProtocolConfig<NumberType>& pc) {
    transport.awaitRegistration();

    transport.broadcast(
        MessageType::PROTOCOL_CONFIG,
        ZeroMQCoordinatorTransport::ids,
        pc
    );
}

/*
 * Host key generation in the coordinator side.
 */
template <typename T, size_t Degree, size_t NbPrimes>
std::optional<int> hostGenerateKeyPair(ZeroMQCoordinatorTransport& trans) {
    DBG("Starting keygen");
    using Q = nfl::poly_p<T, Degree, NbPrimes>;
    // sum up shares of a_i
    DBG("Waiting for a shares");
    auto a = trans.awaitAggregateVectorInput<Q>(MessageType::PUBLIC_KEY_A_SHARES, trans.ids, std::plus<Q>(), 0);
    if (!a) {
        return std::nullopt;
    }
    trans.broadcast(MessageType::PUBLIC_KEY_A_VALUE, trans.ids, *a);
    DBG("Received and broadcasted a, waiting for b");
    auto b = trans.awaitAggregateVectorInput<Q>(MessageType::PUBLIC_KEY_B_SHARES, trans.ids, std::plus<Q>(), 0);
    if (!b) {
        return std::nullopt;
    }
    trans.broadcast(MessageType::PUBLIC_KEY_B_VALUE, trans.ids, *b);
    DBG("Received and broadcasted b");
    return 0;
}


std::vector<std::vector<mpz_class>> pruneAndReorderEIT (
        const std::vector<mpz_class>& eit,
        const std::vector<int>& flags,
        int numAlphas,
        std::vector<size_t> bucketSize) {

    assert (eit.size() == flags.size());

    // Now we can construct a result matrix
    std::vector<std::vector<mpz_class>> result;

    // And finally we fill the result matrix compactly with valid shares from `as`
    int k = 0;
    size_t min_size = 0;
    for (size_t c = 0; c < numAlphas; ++c) {

        std::vector<mpz_class> col;
        for (int r = 0; r < bucketSize[c]; ++r) {
            if (flags[k] == 1) {
                col.push_back(eit[k]);
            }
            k++;
        }
        result.push_back(col);
        if (min_size == 0) { min_size = col.size(); }
        if (min_size > col.size()) { min_size = col.size(); }
    }

    std::vector<mpz_class> col(bucketSize[0], mpz_class(1));
    result.insert(result.begin(), col);

    // resize and transpose 
    std::vector<std::vector<mpz_class>> transposed(min_size);

    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = 0; j < min_size; ++j) {
            transposed[j].push_back(result[i][j]);
        }
    }

    return transposed;
}


template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::optional<std::tuple< std::vector<std::vector<mpz_class>>,  // eit_pruned
            std::vector<mpz_class>,               // alphasPS_tick
            std::vector<mpz_class>,               // c_can
            std::vector<mpz_class>                // c_gcd
          >>
hostPreSieving (
        ZeroMQCoordinatorTransport& transport,
        lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& e,
        ProtocolConfig<T>& config) {
    using P = nfl::poly_p<T, Degree, NbPrimesP>;
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;
   
    auto [alphasPS, bucketSizePS] = math::balanced_bucket_n_primes(config.pbs(), primesPS, config.tauLimitBit(), 1);
    auto [alphasCAN, bucketSizeCAN] = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());
    auto [alphasGCD, bucketSizeGCD] = math::fixed_bucket_n_primes(3*config.pbs()+96, primesGCD, config.tauLimitBit());
    //auto [alphasCAN, bucketSizeCAN] = math::balanced_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit(), 200);
    //auto [alphasGCD, bucketSizeGCD] = math::balanced_bucket_n_primes(3*config.pbs()+96, primesGCD, config.tauLimitBit(), 1);

    // set equal buckets for bucketSizeCan and bucketSizeGCD
    const int bucketSizeCAN_value = lrint(floor(double(primesCAN) / double(alphasCAN.size())));
    for (int i = 0; i < bucketSizeCAN.size(); ++i) {
        bucketSizeCAN[i] = bucketSizeCAN_value;
    }
    const int bucketSizeGCD_value = lrint(floor(double(primesGCD) / double(alphasGCD.size())));
    for (int i = 0; i < bucketSizeGCD.size(); ++i) {
        bucketSizeGCD[i] = bucketSizeGCD_value;
    }
    
    // Aggregate a matrix of `a` shares and broadcast the sum
    std::pair<Q, Q> acc = {0, 0};
    auto as = transport.awaitAggregateVectorInput<std::pair<Q, Q>>(
        MessageType::ENCRYPTED_X_SHARES,
        ZeroMQCoordinatorTransport::ids,
        lattice::pair_add<Q>,
        acc
    );

    if (!as) {
        return std::nullopt;
    }

    DBG("Received Enc(x) matrix.");
    DBG("Communication cost before = " << transport.communicationCost);
    transport.broadcast(MessageType::ENCRYPTED_X_VALUE, ZeroMQCoordinatorTransport::ids, *as);
    DBG("Communication cost after = " << transport.communicationCost);
    // The parties will take the previous sum and multiply in their shares of
    // another random number `b_j,i,t`. These form shares of `e_i,t` which we
    // sum then broadcast for partial decryption
    auto es = transport.awaitAggregateVectorInput<std::pair<Q, Q>>(
        MessageType::ENCRYPTED_XY_PLUS_Z_SHARES,
        ZeroMQCoordinatorTransport::ids,
        lattice::pair_add<Q>,
        std::pair<Q, Q>{0, 0}
    );

    if (!es) {
        return std::nullopt;
    }

    DBG("Received Enc(e_i,t) matrix.");
    transport.broadcast(MessageType::ENCRYPTED_XY_PLUS_Z_VALUE, ZeroMQCoordinatorTransport::ids, *es);

    // Now we finish the decryption, giving us a set of e_i,t values that
    // we can use for validating the shares by gcd.
    auto eit_poly = transport.awaitAggregateVectorInput<Q>(
        MessageType::PARTIAL_XY_MINUS_Z_SHARES,
        ZeroMQCoordinatorTransport::ids,
        std::plus<Q>(), 
        Q{0}
    );

    if (!eit_poly) {
        return std::nullopt;
    }

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
        for (size_t j = 0; j < bucketSizeCAN_value; ++j) { 
            for (size_t i = 0; i < alphasCAN.size(); ++i) {
                tauVector[tvk] = alphasCAN[i];
                tvk++;
            }
        }

        // GCD
        for (size_t j = 0; j < bucketSizeGCD_value; ++j) { 
            for (size_t i = 0; i < alphasGCD.size(); ++i) {
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

    auto c_mpz_t = e.eval_poly(*eit_poly, tauVector);
    std::vector<mpz_class> c(c_mpz_t.size());

    DBG("after eval poly");
    for (size_t i = 0; i < c.size(); ++i) {
        c[i] = mpz_class(c_mpz_t[i]);
    }

    // first eit is for ps 
    //transport.broadcast(MessageType::XY_MINUS_Z_VALUE, ids, c);
    // Need to compute presieve flags
    std::vector<mpz_class> eit (primesPS);

    std::vector<mpz_class> c_can(primesCAN);
    std::vector<mpz_class> c_gcd(primesGCD);

    int pq_size = bucketSizePS[0];
    std::vector<int> flags (primesPS);

    {
        mpz_t one;
        mpz_init(one);
        mpz_set_ui(one, 1);

        mpz_t gcdResult;
        mpz_init(gcdResult);


        //assert(eit.size() == alphas.size() * bucketSize);
        int k = 0;
        for (int j = 0; j < alphasPS.size(); ++j) {

            int sum = 0;
            for (int i = 0; i < bucketSizePS[j]; ++i) {

                mpz_gcd(gcdResult, alphasPS[j].get_mpz_t(), c_mpz_t[k]);
                eit[k] = c[k];

                flags[k] = mpz_cmp(gcdResult, one) == 0;
                sum += flags[k];

                mpz_clear(c_mpz_t[k]);
                k++;
            }

            if (sum < pq_size) {
                pq_size = sum;
            }
        }

        // free mem
        mpz_clear(one);
        mpz_clear(gcdResult);
    }

    transport.broadcast(MessageType::PS_SIEVING_FLAGS, ZeroMQCoordinatorTransport::ids, flags);


    auto eit_pruned = pruneAndReorderEIT(eit, flags, alphasPS.size(), bucketSizePS); 
    auto alphasPS_tick = alphasPS;
    alphasPS_tick.insert(alphasPS_tick.begin(), mpz_class(4));
    for (int i = primesPS; i < primesPS + primesCAN; ++i) {
        c_can[i - primesPS] = c[i];
    }

    for (int i = primesPS + primesCAN; i < primesPS + primesCAN + primesGCD; ++i) {
        c_gcd[i - primesPS - primesCAN] = c[i];
    }

    return std::make_tuple(eit_pruned, alphasPS_tick, c_can, c_gcd);
}


template <typename T>
std::optional<std::vector<mpz_class>>
hostModulusCandidates (
    ZeroMQCoordinatorTransport& transport,
    std::vector<std::vector<mpz_class>> eit_pruned, 
    std::vector<mpz_class> alphasPS_tick,
    std::vector<mpz_class> c_can,
    ProtocolConfig<T>& config
) {

    DBG("Hosting modulus candidates");
    //auto [alphasCAN, _bsz] = math::balanced_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit(), 200);
    auto [alphasCAN, _bsz] = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());
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

    if (!maybe_ax_by_sum) {
        DBG("Failed to receive ax and by");
        return std::nullopt;
    }
    auto ax_by_sum = *maybe_ax_by_sum;

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

    if (!maybe_ab) {
        DBG("Failed to receive ab");
        return std::nullopt;
    }
    auto ab = *maybe_ab;

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
        // copy from ps
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

    // Compute pq_size
    std::vector<mpz_class> candidates(k);
    for (int i = 0; i < k; ++i) {
        candidates[i] = candidates_raw[i];
    }

    // Then we broadcast the candidates to the clients
    transport.broadcast(MessageType::MODULUS_CANDIDATE, ZeroMQCoordinatorTransport::ids, candidates);
    DBG("Sending modulus candidates...");

    int number_responded = 0;
    transport.awaitAggregateVectorInput<int>(MessageType::MUTHU_ACK, ZeroMQCoordinatorTransport::ids,
         [](const int& a, const int& b) {
		    return a+b;
		    }
		    ,number_responded);
    return candidates;
}


template <typename T>
std::optional<boost::dynamic_bitset<>> hostGCDandJacobiTest (
    ZeroMQCoordinatorTransport& transport,
    const std::vector<mpz_class>& c_gcd,
    const std::vector<mpz_class>& candidates,
    ProtocolConfig<T>& config
) {


    DBG("hostGCDTest candidates.size() = " << candidates.size());
    //auto [alphasGCD, _drop_val] = math::balanced_bucket_n_primes(3*config.pbs() + 96, primesGCD, config.tauLimitBit(), 1);
    auto [alphasGCD, _drop_val] = math::fixed_bucket_n_primes(3*config.pbs() + 96, primesGCD, config.tauLimitBit());
    const int bucketSize = lrint(floor(double(primesGCD) / double(alphasGCD.size())));

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

    if (!maybe_result) {
        return std::nullopt;
    }
    auto [pair_axby, gammaSeed] = *maybe_result;

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
    DBG("Broadcast ax and by values");
    transport.broadcast(MessageType::AX_BY_VALUE, ZeroMQCoordinatorTransport::ids, std::pair{std::pair{ax_sum, by_sum}, gammaSeed});

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
    
    if (!maybe_zcrtsggs) {
        return std::nullopt;
    }
    auto [zCRTs, ggs] = *maybe_zcrtsggs;

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

        //DBG("my candidates = " << candidates);
        DBG("Taking gamma values modulo N");
        for (int i = 0; i < candidates.size(); ++i) {
            for (int j = 0; j < config.lambda(); ++j) {
                const mpz_class& N = candidates[i];
                mpz_class& gg = ggs[i * config.lambda() + j];
                mpz_class x = gg % N;

                DBG("x = " << x);
                DBG("N = " << N);
                DBG("ggs[" << i << "] =" << ggs[i * config.lambda() + j]);

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


std::optional<boost::dynamic_bitset<>> hostJacobiProtocol (
    ZeroMQCoordinatorTransport& transport,
    const std::vector<mpz_class>& candidates
) {
    // First, we broadcast an ask for shares of a gamma value
    transport.broadcast(MessageType::GAMMA_SHARES, ZeroMQCoordinatorTransport::ids);

    // Aggregate and broadcast gamma seed
    auto gammaSeed = transport.awaitAggregateVectorInput<mpz_class>(
        MessageType::GAMMA_RANDOM_SEED_SHARES,
        ZeroMQCoordinatorTransport::ids,
        [](const mpz_class& a,
           const mpz_class& b) { return a ^ b; },
        mpz_class(0)
    );

    if (!gammaSeed) {
        return std::nullopt;
    }

    transport.broadcast(MessageType::GAMMA_RANDOM_SEED_VALUE, ZeroMQCoordinatorTransport::ids, *gammaSeed);

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

    if (!maybe_ggs) {
        return std::nullopt;
    }
    auto ggs = *maybe_ggs;

    boost::dynamic_bitset<> discard (candidates.size());

    int discarded = 0;
    //DBG("my candidates = " << candidates);
    DBG("Taking gamma values modulo N");
    for (int i = 0; i < candidates.size(); ++i) {
        const mpz_class& N = candidates[i];
        mpz_class& gg = ggs[i];
        mpz_class x = gg % N;

	DBG("x = " << x);
	DBG("N = " << N);
	DBG("ggs[" << i << "] =" << ggs[i]);

        // Elimination
        if (x != 1 && x != (N - 1)) {
            discard[i] = 1;
	    discarded++;
        }
    }
    LOG(INFO) << "After first Jacobi iteration " << candidates.size() - discarded << " candidates survived";
    DBG("returning discard vector");

    return discard;
}

//==============================================================================
template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
class EncryptedCoordinator
{
public:

    /** Host the RSA MPC ceremony. */
    void host(
        ZeroMQCoordinatorTransport& transport, 
        ProtocolConfig<T>& config, 
        throughput_test_config tconfig,
        size_t throughput_cutoff
        ) {

        auto e = lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>(config);

        DBG("Initiate transport");

        DBG("Host throughput test");
        auto survivor = transport.hostThroughputTest(tconfig, [&](size_t t) { return t > throughput_cutoff; });
        if (survivor.size() < 2) {
            LOG(INFO) << "too few survivor, quit";
            exit(0);
        }
        DBG("Host registration");
        hostRegistration(transport, config);

        transport.myTimers.initialize_timer();
        transport.myTimers.begin(1,"RSA Ceremony",transport.communicationCost);

        auto success = host_rsa_ceremony(config, transport, e);
        while (!success) {
            // restart everything and try again from keygen
            //LOG(INFO) << "First run return with time out, try restarting...";
            //LOG(INFO) << "========= starting second run ==========";
            LOG(INFO) << "Run timed out, restarting";
            config.numParties() = ZeroMQCoordinatorTransport::ids.size();
            transport.parties() = ZeroMQCoordinatorTransport::ids.size();
            success = host_rsa_ceremony(config, transport, e);
            //if (!success) {
            //    LOG(INFO) << "Run timed out, restarting";
                //LOG(FATAL) << "Second run timed out, aborting";
            //}
        }
    }

    std::optional<int> host_rsa_ceremony(
        ProtocolConfig<T>& config,
        ZeroMQCoordinatorTransport& transport,
        lattice::LatticeEncryption<T, Degree, NbPrimesP, NbPrimesQ>& e
    ) {

        DBG("DBG001 GENERATE KEYS");
        transport.myTimers.begin(2,"Key Generation",transport.communicationCost);
        auto keygen_success = hostGenerateKeyPair<T, Degree, NbPrimesQ>(transport);

        if (!keygen_success) {
            LOG(ERROR) << "Keygen timed out!";
            transport.myTimers.end(2,"Key Generation",transport.communicationCost);
            return std::nullopt;
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
        int numberOfTrials = 0;

        // precompute bounds and M with min=4096 max=104729(10k th prime) and 2 tests
        auto [Bs, Ms] = ligero::math::compute_m_b_vec(4096, 104729, 2);
        auto nb_threads = 90;
        auto postSieve = std::bind(ligero::math::test_factorizable_threaded, std::placeholders::_1, Bs, Ms, nb_threads);

        while (!foundModuli && (numberOfTrials < kTrialsThreshold)) {
            numberOfTrials++;
            transport.myTimers.begin(2,"Pre-Sieving",transport.communicationCost);
            //int pq_size = hostPreSieving(transport, ids, e, config); 
            auto presieve_success = hostPreSieving(transport, e, config);

            if (!presieve_success) {
                LOG(ERROR) << "Presieving timed out!";
                transport.myTimers.end(2,"Pre-Sieving",transport.communicationCost);
                return std::nullopt;
            }

            auto [eit_pruned, alphasPS_tick, c_can, c_gcd] = *presieve_success;

            transport.myTimers.end(2,"Pre-Sieving",transport.communicationCost);

            DBG("Beginning Candidate generation");
            transport.myTimers.begin(2,"Generate Candidates",transport.communicationCost);
            auto candidates_success = hostModulusCandidates (transport, eit_pruned, alphasPS_tick, c_can, config);

            if (!candidates_success) {
                LOG(ERROR) << "Modulus candidate timed out!";
                transport.myTimers.end(2,"Generate Candidates",transport.communicationCost);
                return std::nullopt;
            }

            auto candidates = *candidates_success;

            //auto candidates = generateModulusCandidates (e, transport, ids, config);
            transport.myTimers.end(2,"Generate Candidates",transport.communicationCost);
            if (candidates.size() == 0) {
                LOG(INFO) << "No candidates found.";
                transport.broadcast(MessageType::NO_MODULI, ZeroMQCoordinatorTransport::ids);
                continue;
            }
            // sanity check after modulus
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

            boost::dynamic_bitset<> discardFlags; 
            {
                DBG("Running postSieve");
                discardFlags = postSieve(std::vector<mpz_class>(candidates.data(), candidates.data() + candidates.size()));
                transport.broadcast(MessageType::POST_SIEVE, ZeroMQCoordinatorTransport::ids, discardFlags);
                candidates = discardCandidates(candidates, discardFlags);
                DBG("Completed postSieve");
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

                if (!jacobi_success) {
                    LOG(ERROR) << "Jacobi timed out!";
                    transport.myTimers.end(2,"Jacobi Test",transport.communicationCost);
                    return std::nullopt;
                }

                boost::dynamic_bitset<> discardFlags = *jacobi_success;

                transport.broadcast(MessageType::DISCARD_FLAGS, ZeroMQCoordinatorTransport::ids, discardFlags);

                candidates = discardCandidates(candidates, discardFlags);
                transport.myTimers.end(2,"Jacobi Test",transport.communicationCost);
            }
            //LOG(INFO) << "Jacobi Test:" << transport.communicationCost;

            if (candidates.size() == 0) {
                LOG(INFO) << "No candidates found.";
                transport.broadcast(MessageType::NO_MODULI, ZeroMQCoordinatorTransport::ids);
                continue;
            }

            DBG("After Jacobi candidates.size() = " << candidates.size());
            //// Complete the biprimality test with the GCD test, discarding any
            //// other candidates that fail.
            transport.myTimers.begin(2,"GCD",transport.communicationCost);

            if (candidates.size() > 1) {
                candidates.resize(1);
            }
            auto gcd_jacobi_success = hostGCDandJacobiTest (transport, c_gcd, candidates, config);

            if (!gcd_jacobi_success) {
                LOG(ERROR) << "GCD and Jacobi timed out!";
                transport.myTimers.end(2,"GCD",transport.communicationCost);
                return std::nullopt;
            }

            discardFlags = *gcd_jacobi_success;

            transport.broadcast(MessageType::DISCARD_FLAGS, ZeroMQCoordinatorTransport::ids, discardFlags);

            transport.myTimers.end(2,"GCD",transport.communicationCost);

            transport.myTimers.begin(2,"Discarding",transport.communicationCost);

            candidates = discardCandidates(candidates, discardFlags);

            if (candidates.size() == 0) {
                LOG(INFO) << "No candidates found.";
                transport.broadcast(MessageType::NO_MODULI, ZeroMQCoordinatorTransport::ids);
                continue;
            }
            transport.myTimers.end(2,"Discarding",transport.communicationCost);

            //// Checking performance
            //PERFORMANCE_CHECKPOINT_WITH_ID(timerObj, "GCD Test");
            //LOG(INFO) << "GCD Test:" << transport.communicationCost;

            // Here now, any modulus candidate remaining in our vector has passed
            // the test.
            if (candidates.size() > 0) {
                LOG(INFO) << "Found " << candidates.size() << " valid moduli:";

                for (int i = 0; i < candidates.size(); ++i){
                    LOG(INFO) << candidates[i];
                }
                transport.broadcast(MessageType::FOUND_MODULI, ZeroMQCoordinatorTransport::ids);
                foundModuli = true;
            } else {
                LOG(INFO) << "Found no valid moduli restarting.";
                transport.broadcast(MessageType::NO_MODULI, ZeroMQCoordinatorTransport::ids);
            }

            // Below we are going to check if our output is a product of two primes
            // There is no other way to check it than actually computing p and q as we 
            // cannot possibly factor 2048 bit numbers! 

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
                
                if (!maybe_pq) {
                    LOG(ERROR) << "PQ timed out!";
                    return std::nullopt;
                }
                auto pq = *maybe_pq;

                // sanity check
                for(size_t i = 0; i < pq.first.size(); ++i)
                {
                    assert(boost::multiprecision::miller_rabin_test(MPInt(pq.first[i].get_mpz_t()), 4));
                    assert(boost::multiprecision::miller_rabin_test(MPInt(pq.second[i].get_mpz_t()), 4));
                }

                LOG(INFO) << "Prime factorization of the candidates are: ";
                for (size_t i = 0; i < pq.first.size(); ++i) {
                    LOG(INFO) << pq.first[i];
                    LOG(INFO) << pq.second[i];
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
            transport.myTimers.end(1,"RSA Ceremony",transport.communicationCost);
        }
        return 0;
    }

};

} // namespace ligero
