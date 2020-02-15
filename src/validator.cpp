#include <sstream>
#include "Ligero.hpp"
#include "LatticeEncryption.hpp"
#include "Common.hpp"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/archive/binary_iarchive.hpp>


using namespace ligero;

constexpr auto degree = 1 << 16;
constexpr auto p = 9;
constexpr auto q = 21;

using Q = nfl::poly_p<uint64_t, degree, q>;

static boost::filesystem::ifstream ifs("script.data", std::ios::in | std::ios::binary);
static boost::archive::binary_iarchive ia(ifs, boost::archive::no_header);

static int numParties;
static int currentStep = 1;

template <typename T>
void validate(MessageType expectedMessageType,
        std::function<T(T, T)> op,
        T origAccumulator
        ) {
    //try {
        if (ifs.is_open()) {
            // Now read and validate each round
            MessageType t;
            ia >> t;

            int msgId = int(t);
            //std::cout << "t: " << msgId << std::endl;
            if (msgId > 0) msgId = 1 + msgId - int(MessageType::ID_PARTY);
            std::cout << currentStep++ << ". " << msgs[msgId] << std::endl;

            if (t != expectedMessageType) {
                LOG(FATAL) << "Unexpected MessageType. Expected message type: " << msgs[1 + int(expectedMessageType) - int(MessageType::ID_PARTY)]; 
            }
            T computedAccumulator = origAccumulator;
            bool initialized = false;
            for (size_t i = 0; i < numParties; ++i) {
                std::cout << "i = " << i << std::endl;
                if (ifs && ifs.peek() != EOF) {
                    std::pair<MessageType, T> zzz;
                    ia >> zzz;
                    auto [ti, x] = zzz;
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
            // TODO: Check that both are NOT zeros
            if (at != expectedMessageType || expectedAccumulator != computedAccumulator) {
                LOG(FATAL) << "Accumulated result did NOT match expected value";
            }
        } else {
            std::cout << "file is not good" << std::endl;
        }
        /*
    } catch (...) {
        std::cout << "caught exception" << std::endl;
        ifs.close();
    }*/

}


int main(int argc, char** argv)
{
    //try {
        if (ifs.is_open()) {

            // Read the number of parties
            ia >> numParties;
            std::cout << "numParties: " << numParties << std::endl;

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

                validate<Q>(
                        MessageType::PARTIAL_XY_MINUS_Z_SHARES,
                        std::plus<Q>(), 
                        Q{0}
                );

                std::cout << "Pre-sieve: [SUCCESS]" << std::endl;
            }

            // Validate modulus candidate generation
            {
                using PairOfVec = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>; 
                auto acc_pair = std::pair<std::vector<mpz_class>, std::vector<mpz_class>>(std::vector<mpz_class>(primesCAN), std::vector<mpz_class>(primesCAN));
                validate<PairOfVec>(
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

                validate<std::vector<mpz_class>>(
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

                // validate sync step
                validate<int>(
                        MessageType::MUTHU_ACK,
                        [](const int& a, const int& b) {
                        return a+b;
                        },
                        0);

                std::cout << "Candidate generation: [SUCCESS]" << std::endl;
            }

            // Validate jacobi
            {
                validate<mpz_class>(
                        MessageType::GAMMA_RANDOM_SEED_SHARES,
                        [](const mpz_class& a,
                            const mpz_class& b) { return a ^ b; },
                        mpz_class(0)
                );

                validate<std::vector<mpz_class>>(
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

                std::cout << "Jacobi: [SUCCESS]" << std::endl;
            }

            // Validate GCD and Jacobi
            {
                auto acc_pair = std::pair{std::pair{std::vector<mpz_class>(primesGCD), std::vector<mpz_class>(primesGCD)}, mpz_class(0)};

                validate<PairOfVecGCD> 
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


                validate<std::pair<std::vector<mpz_class>, std::vector<mpz_class>>>(
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
                std::cout << "GCD and Jacobi: [SUCCESS]" << std::endl;
            }

            std::cout << std::endl << ">>>>  VALIDATION PASSED <<<<" << std::endl;
            /*
               while (ifs && ifs.peek() != EOF) {
               MessageType t;
               ia >> t;
               int msgId = int(t);
               std::cout << "t: " << msgId << std::endl;
               if (msgId > 0) msgId = msgId - int(MessageType::ID_PARTY);
               std::cout << "header: " << msgs[msgId] << std::endl;

               }
               */
        }
        /*
    } catch (...) {
        std::cout << "caught exception" << std::endl;
        ifs.close();
    }
    */

    ifs.close();
    return 0;
}

