#include <sstream>
#include "Ligero.hpp"
#include "Math.hpp"
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/functional/hash.hpp>
#include <boost/uuid/uuid_io.hpp>


#include <boost/serialization/vector.hpp>
#include <fstream>
#include <iostream>

#include <boost/serialization/unordered_map.hpp>
#include "LatticeEncryption.hpp"
#include "Factoring.hpp"
#include "Common.hpp"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/archive/binary_iarchive.hpp>

using namespace ligero;

static ECDSA<ECP, SHA256>::PublicKey publicKey;
static ECDSA<ECP, SHA256>::PrivateKey privateKey;
AutoSeededRandomPool prng;

static std::vector<mpz_class> alphasPS, alphasCAN, alphasGCD;
static std::vector<size_t> bucketSizePS, bucketSizeCAN, bucketSizeGCD;

std::unordered_map<SocketId, ECDSA<ECP, SHA256>::PublicKey, boost::hash<SocketId>> partiesPublicKeys;

constexpr auto degree = 1 << 16;
constexpr auto p = 9;
constexpr auto q = 21;
constexpr auto sigma = 8;
constexpr auto lambda = 128;
constexpr auto tau_limit_bit = 175;
constexpr auto tau = 1000;
constexpr auto numCandidates = 2048;
constexpr auto pbs = 1000;

using Q = nfl::poly_p<uint64_t, degree, q>;

static boost::filesystem::ifstream ifs("script.data", std::ios::in | std::ios::binary);
static boost::archive::binary_iarchive ia(ifs, boost::archive::no_header);

static int numParties;
static int currentStep = 1;
static size_t minRowSize = -1;

/** Helper function to prune and reorder eit array */
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

/** Helper function to compute eit and flags arrays */
std::tuple<std::vector<mpz_class>, std::vector<int>, std::vector<mpz_class>> computeEITAndFlags(std::array<mpz_t, degree> c_mpz_t) {

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


    //assert(eit.size() == alphas.size() * bucketSize);
    int k = 0;
    for (int j = 0; j < alphasPS.size(); ++j) {

        int ick = 0;
        int sum = 0;
        for (int i = 0; i < bucketSizePS[j]; ++i) {

            mpz_gcd(gcdResult, alphasPS[j].get_mpz_t(), c_mpz_t[k]);
            eit[k] = c[k];

            flags[k] = mpz_cmp(gcdResult, one) == 0;
            sum += flags[k];
            /*
            if (flags[k] == 1) {
                index_candidates[k] = ick++;
            } else {
                index_candidates[k] = -1;
            }
            */

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


/** Helper to validate given round by comparing the expected result to 
 *  the computed result*/
template <typename T>
T validate(MessageType expectedMessageType, std::function<T(T, T)> op, T origAccumulator) {
    T computedAccumulator = origAccumulator;
    if (ifs.is_open()) {
        // Now read and validate each round
        MessageType t;
        std::cout << "Expected message type: " << std::endl; 
        std::cout << "READ for: " << msgs[1 + int(expectedMessageType) - int(MessageType::ID_PARTY)] << std::endl;
        ia >> t;

        int msgId = int(t);
        
        //std::cout << "int(t): " << msgId << std::endl;
        //std::cout << "t: " << (int)t << std::endl;
        //std::cout << "expectedMT: " << (int)expectedMessageType << std::endl;
        if (t != expectedMessageType) {
            std::cout << "Unexpected MessageType. Expected message type: " << msgs[1 + int(expectedMessageType) - int(MessageType::ID_PARTY)] << std::endl; 
            if (1 + int(t) - int(MessageType::ID_PARTY) < msgs.size()) {
                std::cout << "Got MessageType " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << std::endl;
            } else {
                std::cout << "Got MessageType" << 1 + int(t) - int(MessageType::ID_PARTY) << std::endl;
            }
            exit(1);
        } else {
            if (msgId > 0) msgId = 1 + msgId - int(MessageType::ID_PARTY);
            std::cout << currentStep++ << ". " << msgs[msgId] << std::endl;
        }

        bool initialized = false;
        for (size_t i = 0; i < numParties; ++i) {
            std::cout << "i = " << i << std::endl;
            if (ifs && ifs.peek() != EOF) {
                std::pair<MessageType, std::pair<SocketId, std::string>> zzz;
                ia >> zzz;
                auto [ti, idAndMsg] = zzz;
                auto [partyId, xWithSig] = idAndMsg;

                std::string signature = xWithSig.substr(0, 64);
                // verify signature

                std::string msgToVerify = xWithSig.substr(64, xWithSig.size() - 64);
                {
                    auto pcit = partiesPublicKeys.find(partyId);
                    if (pcit == partiesPublicKeys.end()) {
                        std::cout << "Could not find public key for party " << boost::uuids::to_string(partyId) << std::endl;
                        exit(1);
                    }

                    /*
                    bool valKeyResult = pcit->second.Validate(prng, 3);
                    if (!valKeyResult) {
                        std::cout << "Failed to validate public key for party " << boost::uuids::to_string(partyId) << std::endl;
                        exit(1);
                    } else {
                        std::cout << "Successfully validated public key for party " << boost::uuids::to_string(partyId) << std::endl;
                    }
                    */

                    ECDSA<ECP, SHA256>::Verifier verifier(pcit->second);

                    bool result = verifier.VerifyMessage((const byte*)&msgToVerify[0], msgToVerify.size(), (const byte*)&signature[0], signature.size());

                    // Verification failure?
                    if (!result) {
                        std::cout << "Failed signature verification for party " << boost::uuids::to_string(partyId) << std::endl;
                        exit(1);
                    } 
                }

                std::istringstream ssi(msgToVerify);
                boost::archive::binary_iarchive iassi(ssi);
                T x;
                iassi >> x;

                MessageType xt = MessageType(ti);
                if (xt != t) {
                    std::cout << "Unexpected MessageType" << std::endl;
                    std::cout << "Expected: " << ">>" << int(t) << "<<" << " -- "  << msgs[1 + int(t) - int(MessageType::ID_PARTY)];
                    std::cout << "Got: " << ">>" << int(xt) << "<<" << " -- "  << msgs[1 + int(xt) - int(MessageType::ID_PARTY)];
                    assert(false);
                }


                if (!initialized) {
                    initialized = true;
                    computedAccumulator = x;
                } else {
                    computedAccumulator = op(computedAccumulator, x);
                }
            } else {
                LOG(FATAL) << "Something went wrong, not enough data to read...";
            }
        }
        if (computedAccumulator == origAccumulator) {
            LOG(FATAL) << "Computed accumulator did not change, aborting...";
        }

        std::pair<MessageType, T> ea;
        ia >> ea;
        auto [at, expectedAccumulator] = ea;
        if (expectedMessageType != MessageType::SYNCHRONIZE_NOW) {
            if (at != expectedMessageType || expectedAccumulator != computedAccumulator) {
                LOG(FATAL) << "Accumulated result did NOT match expected value";
            }
        }
    } else {
        std::cout << "file is not good" << std::endl;
    }
    return computedAccumulator;
}



int main(int argc, char** argv)
{
    

    if (ifs.is_open()) {


        // Read the number of parties
        ia >> numParties;

        // Load coordinators public key from file
        {
            FileSource fs("publicKey.coordinator", true);
            publicKey.Load(fs);

            bool valKeyResult = publicKey.Validate(prng, 3);
            if (!valKeyResult) {
                std::cout << "Failed to validate public key for coordinator" << std::endl;
                exit(1);
            } else {
                std::cout << "Successfully validated public key for coordinator " << std::endl;
            }
        }
        // Load parties public keys
        {
            boost::filesystem::ifstream pcifs;
            pcifs.open("publicKey.parties", std::ios::in | std::ios::binary);
            boost::archive::binary_iarchive pcIArchive(pcifs);

            for (size_t i = 0; i < numParties; ++i) {
                std::pair<SocketId, std::string> x;
                pcIArchive >> x;

                ECDSA<ECP, SHA256>::PublicKey pk;
                StringSource pkss(x.second, true);
                pk.Load(pkss);

                // validate public key
                bool valKeyResult = pk.Validate(prng, 3);
                if (!valKeyResult) {
                    std::cout << "Failed to validate public key for party " << boost::uuids::to_string(x.first) << std::endl;
                    exit(1);
                } else {
                    std::cout << "Successfully validated public key for party " << to_string(x.first) << std::endl;
                }

                partiesPublicKeys[x.first] = std::move(pk);
            }

            pcifs.close();
        }


        std::vector<mpz_class> candidatesCAN, candidatesPostSieve, candidatesJacobi;

        ProtocolConfig<uint64_t> config = ProtocolConfig<uint64_t>(numParties, p, q, degree, sigma, lambda, tau_limit_bit, pbs, ProtocolMode::NORMAL, publicKey);
        auto e = lattice::LatticeEncryption<uint64_t, degree, p, q>(config);

        std::tie(alphasPS, bucketSizePS) = math::balanced_bucket_n_primes(config.pbs(), primesPS, config.tauLimitBit(), 1);
        std::tie(alphasCAN, bucketSizeCAN) = math::fixed_bucket_n_primes(config.pbs()+48, primesCAN, config.tauLimitBit());

        std::vector<size_t> bucketSizeGCD;
        std::tie(alphasGCD, bucketSizeGCD) = math::fixed_bucket_n_primes(3*config.pbs()+210, primesGCD, config.tauLimitBit());

        std::vector<std::vector<mpz_class>> eit_prunedPS;
        std::vector<mpz_class> c_can(primesCAN);
        auto alphasPS_tick = alphasPS;
        std::vector<mpz_class> c_gcd(primesGCD);


        // set equal buckets for bucketSizeCan and bucketSizeGCD
        auto bucketSizeCAN_value = lrint(floor(double(primesCAN) / double(alphasCAN.size())));
        for (int i = 0; i < bucketSizeCAN.size(); ++i) {
            bucketSizeCAN[i] = bucketSizeCAN_value;
        }
        auto bucketSizeGCD_value = lrint(floor(double(primesGCD) / double(alphasGCD.size())));
        for (int i = 0; i < bucketSizeGCD.size(); ++i) {
            bucketSizeGCD[i] = bucketSizeGCD_value;
        }


        // Validate Keygen
        {
            validate<Q>(MessageType::PUBLIC_KEY_A_SHARES, std::plus<Q>(), 0);
            validate<Q>(MessageType::PUBLIC_KEY_B_SHARES, std::plus<Q>(), 0);
            std::cout << "Keygen: [SUCCESS]" << std::endl;
        }

        // Validate pre-sieve
        {
            validate<std::pair<Q, Q>>(
                    MessageType::ENCRYPTED_X_SHARES,
                    lattice::pair_add<Q>, 
                    std::pair<Q, Q>{0, 0}
                    );

            validate<std::pair<Q, Q>>(
                    MessageType::ENCRYPTED_XY_PLUS_Z_SHARES, 
                    lattice::pair_add<Q>,
                    std::pair<Q, Q>{0, 0}
                    );

            {
                Q computedAccumulator = validate<Q>(
                    MessageType::PARTIAL_XY_MINUS_Z_SHARES,
                    std::plus<Q>(), 
                    Q{0}
                    );

                std::vector<mpz_class> eitPS;
                std::vector<int> flags_can;
                Q eit_poly = computedAccumulator;

                // convert eit_poly to eit
                std::vector<mpz_class> tauVector(degree);
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
                    //index_candidates.resize(Degree);
                    for (size_t j = 0; j < bucketSizeCAN_value; ++j) {
                        for (size_t i = 0; i < alphasCAN.size(); ++i) {
                            //index_candidates[tvk] = j;
                            tauVector[tvk] = alphasCAN[i];
                            tvk++;
                        }
                    }

                    // GCD
                    for (size_t j = 0; j < bucketSizeGCD_value; ++j) {
                        for (size_t i = 0; i < alphasGCD.size(); ++i) {
                            /*
                            if (j < 1) {
                                index_candidates[tvk] = -2;
                            } else {
                                index_candidates[tvk] = -3;
                            }
                            */
                            tauVector[tvk] = alphasGCD[i];
                            tvk++;
                        }
                    }

                    // fill the rest with dummy values from pre-sieving
                    for (; tvk < degree; ++tvk) {
                        tauVector[tvk] = alphasPS[0];
                    }
                }

                DBG("before eval poly");

                auto c_mpz_t = e.eval_poly(eit_poly, tauVector);

                DBG("after eval poly");

                // Need to compute presieve flags
                std::vector<mpz_class> c;
                std::tie(eitPS, flags_can, c) = computeEITAndFlags(c_mpz_t);

                
                eit_prunedPS = pruneAndReorderEIT(eitPS, flags_can, alphasPS.size());

                alphasPS_tick.insert(alphasPS_tick.begin(), mpz_class(4));
                for (int i = primesPS; i < primesPS + primesCAN; ++i) {
                    c_can[i - primesPS] = c[i];
                }

                for (int i = primesPS + primesCAN; i < primesPS + primesCAN + primesGCD; ++i) {
                    c_gcd[i - primesPS - primesCAN] = c[i];
                }

                std::vector<int> expectedFlagsCan;
                MessageType expectedFlagsType;
                ia >> expectedFlagsType;

                if (expectedFlagsType != MessageType::PS_SIEVING_FLAGS) {
                    std::cout << "Unexpected MessageType " << int(expectedFlagsType) << std::endl;
                    std::cout << "Expected: " << msgs[1 + int(MessageType::PS_SIEVING_FLAGS) - int(MessageType::ID_PARTY)] << std::endl;
                    std::cout << "Got: " << msgs[1 + int(expectedFlagsType) - int(MessageType::ID_PARTY)];
                    assert(false);
                }

                ia >> expectedFlagsCan;
                
                if (expectedFlagsCan.size() != flags_can.size()) {
                    std::cout << "Expected eit flags size " << expectedFlagsCan.size()
                              << " does not match compute flags size " << flags_can.size();
                    assert(false);
                }
                for (size_t i = 0; i < expectedFlagsCan.size(); i++) {
                    if (flags_can[i] != expectedFlagsCan[i]) {
                        std::cout << "Expected flag flags_can[" << i << "] = " << expectedFlagsCan[i]
                                  << " does not match computed flags[" << i << "] = " << flags_can[i];
                        assert(false);
                    }
                }
            }

            std::cout << "Pre-sieve: [SUCCESS]" << std::endl;
        }

        // Validate modulus candidate generation
        {
            auto bucketSize = lrint(floor(double(primesCAN) / double(alphasCAN.size())));
            if (bucketSize > eit_prunedPS.size()) {
                bucketSize = eit_prunedPS.size();
            }
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

            using PairOfVec = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>; 
            auto acc_pair = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>(std::vector<mpz_class>(primesCAN), std::vector<mpz_class>(primesCAN));
            auto ax_by_sum = validate<PairOfVec>(
                    MessageType::AX_BY_SHARES,
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
            auto ab = validate<std::vector<mpz_class>>(
                    MessageType::AXB_MINUS_BYA_SHARES,
                    [](const std::vector<mpz_class>& a, const std::vector<mpz_class>& b) {
                    assert(a.size() == b.size());

                    std::vector<mpz_class> result(a.size());

                    for(size_t i = 0; i < a.size(); ++i) {
                    result[i] = a[i] + b[i]; 
                    }
                    return result;
                    },
                    std::vector<mpz_class>(primesCAN)
                    );

            {
                int k = 0;
                for (int j = 0; j < bucketSize; ++j) {
                    for (int i = 0; i < alphasCAN.size(); ++i) {
                        ab[k] = math::mod(ab[k] + c_can[k] - ax_by_sum.first[k] * ax_by_sum.second[k], alphasCAN[i]);
                        ++k;
                    }
                }
            }

            // Validate MODULUS_CANDIDATES
            {

               

                std::vector<mpz_class> candidates_raw(bucketSize);

                int k = 0;
                auto pq_size = bucketSize;
                if (pq_size > eit_prunedPS.size()) {
                    pq_size = eit_prunedPS.size();
                }
                for (int i = 0; i < pq_size; ++i) {

                    std::vector<mpz_class> x(alphasCAN.size() + alphasPS_tick.size());

                    //DBG("copy from can");
                    for (int zz = 0; zz < alphasCAN.size(); ++zz) {
                        x[zz] = mpz_class(ab[i * alphasCAN.size() + zz]);
                    }
                    //DBG("copy from ps");
                    // copy from ps
                    for (int zz = 0; zz < alphasPS_tick.size(); ++zz) {
                        x[zz + alphasCAN.size()] = mpz_class(eit_prunedPS[i][zz]);
                    }

                    //DBG("before crt_reconstruct");
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
                            std::cout << "Candidate[" << k << "] = " << candidates_raw[k] << " is divisible by " << boost::math::prime(zz) << std::endl;
                            assert(false);
                        }
                        mpz_clear(result);
                    }
                    //DBG("after crt_reconstruct");
                    k++;
                }

                DBG("k = " << k);

                // Compute pq_size
                std::vector<mpz_class> candidates(k);
                for (int i = 0; i < k; ++i) {
                    candidates[i] = candidates_raw[i];
                }

                // Now compare with the expected one from file
                std::vector<mpz_class> expectedModulusCandidates;
                MessageType expectedModulusCandidatesType;
                ia >> expectedModulusCandidatesType;

                if (expectedModulusCandidatesType != MessageType::MODULUS_CANDIDATE) {
                    std::cout << "Unexpected MessageType " << int(expectedModulusCandidatesType) << std::endl;
                    std::cout << "Expected: " << msgs[1 + int(MessageType::MODULUS_CANDIDATE) - int(MessageType::ID_PARTY)] << std::endl;
                    std::cout << "Got: " << msgs[1 + int(expectedModulusCandidatesType) - int(MessageType::ID_PARTY)];
                    assert(false);
                }
                ia >> expectedModulusCandidates;

                if (expectedModulusCandidates.size() != candidates.size()) {
                    std::cout << "Size does not match" << std::endl;
                    std::cout << "Expected candidates size: " << expectedModulusCandidates.size() << std::endl;
                    std::cout << "Computed candidates size: " << candidates.size();
                    assert(false);
                }
                candidatesCAN = candidates;
            }



            // validate sync step
            validate<int>(
                    MessageType::SYNCHRONIZE_NOW,
                    [](const int& a, const int& b) {
                    return a+b;
                    },
                    0);

            std::cout << "Candidate generation: [SUCCESS]" << std::endl;
        }

        // Validate PostSieve
        {
            // read post sieve flags
            boost::dynamic_bitset<> postSieveFlags;
            MessageType postSieveType;

            ia >> postSieveType;
            if (postSieveType != MessageType::POST_SIEVE) {
                std::cout << "Unexpected MessageType " << int(postSieveType) << std::endl;
                std::cout << "Expected: " << msgs[1 + int(MessageType::POST_SIEVE) - int(MessageType::ID_PARTY)] << std::endl;
                std::cout << "Got: " << msgs[1 + int(postSieveType) - int(MessageType::ID_PARTY)];
                assert(false);
            }

            ia >> postSieveFlags;

            auto [Bs, Ms] = ligero::math::compute_m_b_vec(4096, 104729, 2);
            auto nb_threads = 90;
            auto postSieve = std::bind(ligero::math::test_factorizable_threaded, std::placeholders::_1, Bs, Ms, nb_threads);


            std::vector<mpz_class> nonBiPrimes;
            for (size_t i = 0; i < candidatesCAN.size(); i++) {
                if (postSieveFlags[i] == 1) {
                    nonBiPrimes.push_back(candidatesCAN[i]);
                }
            }

            auto discardFlagsToValidate = postSieve(std::vector<mpz_class>(nonBiPrimes.data(), nonBiPrimes.data() + nonBiPrimes.size()));

            // Check that all nonBiPrimes were eliminated
            if (discardFlagsToValidate.count() != nonBiPrimes.size()) {
                std::cout << "Post sieve validation failed" << std::endl;
                std::cout << "Total amount of non-biprimes which were not eliminated is: " << nonBiPrimes.size() - discardFlagsToValidate.count() << std::endl;
                std::cout << "The following non-biprimes were not eliminated:" << std::endl;
                for (size_t i = 0; i < nonBiPrimes.size(); i++) {
                    if (discardFlagsToValidate[i] == 0) {
                        std::cout << "index: " << i << " = " <<  nonBiPrimes[i] << std::endl;
                        assert(false);
                    }
                }
            }
            candidatesPostSieve = discardCandidates(candidatesCAN, postSieveFlags);
            std::cout << currentStep++ << ". " << "Post-Sieve: [SUCCESS]" << std::endl;
        }

        // Validate jacobi
        {
            auto candidates = candidatesPostSieve;
            validate<mpz_class>(
                    MessageType::GAMMA_RANDOM_SEED_SHARES,
                    [](const mpz_class& a,
                        const mpz_class& b) { return a ^ b; },
                    mpz_class(0)
                    );

            auto ggs = validate<std::vector<mpz_class>>(
                    MessageType::EXPONENTIATED_GAMMA_VALUE,
                    [](const std::vector<mpz_class>& a, const std::vector<mpz_class>& b) {

                    if (a.size() == 0) return b;

                    std::vector<mpz_class> c(b.size());
                    for (size_t i = 0; i < c.size(); ++i) {
                    c[i] = a[i] * b[i];
                    }
                    return c;
                    },
                    std::vector<mpz_class>()
                    );

            boost::dynamic_bitset<> discard (candidates.size());

            int discarded = 0;
            //DBG("my candidates = " << candidates);
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

            boost::dynamic_bitset<> discardFlagsJacobi;
            MessageType discardFlagsJacobiType;
            ia >> discardFlagsJacobiType;
            if (discardFlagsJacobiType != MessageType::DISCARD_FLAGS) {
                std::cout << "Unexpected MessageType " << int(discardFlagsJacobiType) << std::endl;
                std::cout << "Expected: " << msgs[1 + int(MessageType::DISCARD_FLAGS) - int(MessageType::ID_PARTY)] << std::endl;
                std::cout << "Got: " << msgs[1 + int(discardFlagsJacobiType) - int(MessageType::ID_PARTY)];
                assert(false);
            }
            ia >> discardFlagsJacobi;

            for (size_t i = 0; i < candidates.size(); i++) {
                if (discardFlagsJacobi[i] != discard[i]) {
                    std::cout << "Discard flags for jacobi do not match the computed ones" << std::endl;
                    std::cout << "Computed[" << i << "] = " << discard[i] << std::endl;
                    std::cout << "Expected[" << i << "] = " << discardFlagsJacobi[i] << std::endl;
                }
            }

            candidatesJacobi = discardCandidates(candidates, discardFlagsJacobi);
            std::cout << currentStep++ << ". " << "Jacobi flags validation passed" << std::endl;
            std::cout << "Jacobi: [SUCCESS]" << std::endl;
        }

        // Validate GCD and Jacobi
        {
            std::vector<mpz_class> ggsGCD;
            const int bucketSize = lrint(floor(double(primesGCD) / double(alphasGCD.size())));
            auto candidates = candidatesJacobi;
            candidates.resize(1);

            auto acc_pair = std::pair{std::pair{std::vector<mpz_class>(primesGCD), std::vector<mpz_class>(primesGCD)}, mpz_class(0)};

            auto result_value = validate<PairOfVecGCD> 
                (
                 MessageType::GCD_AX_BY_SHARES,
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

            auto [pair_axby, gammaSeed] = result_value;

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
            auto zcrt_value = validate<std::pair<std::vector<mpz_class>, std::vector<mpz_class>>>(
                    MessageType::AXB_MINUS_BYA_SHARES,
                    [](const std::pair<std::vector<mpz_class>, std::vector<mpz_class>>& pa,
                        const std::pair<std::vector<mpz_class>, std::vector<mpz_class>>& pb) {
                    auto [a, expGa] = pa;
                    auto [b, expGb] = pb;

                    if (a.size() == 0) return pb; 

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
                        std::pair{std::vector<mpz_class>(), std::vector<mpz_class>()}
            );

            std::vector<mpz_class> zCRTs;

            std::tie(zCRTs, ggsGCD) = zcrt_value;


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




            // TODO: check flags
            boost::dynamic_bitset<> discardFlagsGCD;
            MessageType discardFlagsGCDType;
            ia >> discardFlagsGCDType;

            if (discardFlagsGCDType != MessageType::DISCARD_FLAGS) {
                std::cout << "Unexpected MessageType " << int(discardFlagsGCDType) << std::endl;
                std::cout << "Expected: " << msgs[1 + int(MessageType::DISCARD_FLAGS) - int(MessageType::ID_PARTY)] << std::endl;
                std::cout << "Got: " << msgs[1 + int(discardFlagsGCDType) - int(MessageType::ID_PARTY)];
                assert(false);
            }

            ia >> discardFlagsGCD;

            for (size_t i = 0; i < candidates.size(); i++) {
                if (discardFlagsGCD[i] != discard[i]) {
                    std::cout << "Discard flags for GCD and jacobi do not match the computed ones" << std::endl;
                    std::cout << "Computed[" << i << "] = " << discard[i] << std::endl;
                    std::cout << "Expected[" << i << "] = " << discardFlagsGCD[i] << std::endl;
                }
            }

            std::cout << currentStep++ << ". " << "GCD Jacobi flags validation passed" << std::endl;
            std::cout << "GCD and Jacobi: [SUCCESS]" << std::endl;
        }

        std::cout << std::endl << ">>>>  VALIDATION PASSED <<<<" << std::endl;
    }
    ifs.close();
    return 0;
}
