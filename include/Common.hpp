#pragma once

#include <boost/dynamic_bitset.hpp>
#include <boost/dynamic_bitset/serialization.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/serialization/array.hpp>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <sodium/crypto_generichash.h>

#include <iostream>
#include <gmpxx.h>

#include "nfl.hpp"
#include <sstream>

#define MAX_LINES_DISPLAY           10

template <class T, std::size_t N>
std::ostream& operator<< (std::ostream& o, const std::array<T, N>& arr)
{
    std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
    return o;
}

#include "easylogging++.h"

typedef boost::multiprecision::mpz_int MPInt;

typedef uint64_t RegularInteger;
typedef unsigned __int128 WideInteger;

enum class LinearConstraintFlag : uint64_t {
    Regular,
	SigmaProtocol
};

//==============================================================================
// Platform definitions and app configuration
INITIALIZE_EASYLOGGINGPP

#ifdef NDEBUG
  #define DBG(x)
#else
  #define DBG(x) do { LOG(DEBUG) << x << std::endl; } while (0)
#endif

#ifndef LIG_VECTOR_WIDTH
 #define LIG_VECTOR_WIDTH 1000
#endif

#define THROW_IF_NOT_INITIALIZED(X) initialized ? X : throw std::logic_error(kPublicKeyInvalidGetterInvocation);

// Timers Constraints
constexpr unsigned int comSizeBuffer        = 100;
constexpr unsigned int timeSizeBuffer       = 100;
constexpr unsigned int logLevels            = 4;

// Multi-Threading Physical Constraints
#define NB_MAX_THREADS                      90
#define MAX_SINGLE_THREADED_ITERATIONS      3

//==============================================================================
// Profiling    
typedef struct timepoint {
    size_t cummulative_communication;
    struct timespec timer;
} timepoint;

class timers {
protected:
    struct timespec anchor;
    std::unordered_map<std::string, timepoint> timer_start;
    bool active = false;

public:
    void initialize_timer() {
        this->active = true;
        clock_gettime(CLOCK_MONOTONIC, &this->anchor);
    }

    void begin(size_t level, std::string label, size_t cummulative_communication) {
        struct timespec time;
        
        if (!this->active) return;

        clock_gettime(CLOCK_MONOTONIC, &time);
 
        timepoint tmp;
        tmp.timer = time;
        tmp.cummulative_communication = cummulative_communication;
        timer_start.insert(std::pair<std::string, timepoint>(label, tmp));
    }

    void timespec_diff(struct timespec *start, struct timespec *stop, struct timespec *result)
    {
    if ((stop->tv_nsec - start->tv_nsec) < 0) {
            result->tv_sec = stop->tv_sec - start->tv_sec - 1;
            result->tv_nsec = stop->tv_nsec - start->tv_nsec + 1000000000L;
        } else {
            result->tv_sec = stop->tv_sec - start->tv_sec;
            result->tv_nsec = stop->tv_nsec - start->tv_nsec;
        }
    return;
    }

    // buf needs to store 30 characters
    int timespec2str(char *buf, uint len, struct timespec *ts) {
        uint ret;
        struct tm t;

        if ((ts->tv_nsec<0)||(ts->tv_sec<0)) {strcpy(buf,"00:00.000000"); return 0;}

        tzset();
        if (localtime_r(&(ts->tv_sec), &t) == NULL)
            return 1;

        ret = strftime(buf, len, "%M:%S", &t);
        if (ret == 0)
            return 2;
        len -= ret - 1;

        ret = snprintf(&buf[strlen(buf)], len, ".%06ld", ts->tv_nsec/1000L);

        if (ret >= len)
            return 3;

        return 0;
    }

    void end(size_t level, std::string label, size_t cummulative_communication) {
        struct timespec time;

        if (!this->active) return;

        clock_gettime(CLOCK_MONOTONIC, &time);
        
        std::unordered_map<std::string, timepoint>::iterator iterator = timer_start.find(label);
        if (iterator == timer_start.end()) {std::cout << "error, label " + label + " doesn't exist"; exit(0);}
        
        timepoint reference = iterator->second;
        timespec label_time; 
        timespec aggregate_time;        

        timespec_diff(&anchor, &time, &aggregate_time);
        timespec_diff(&reference.timer, &time, &label_time);

        char label_time_char[timeSizeBuffer];
        char aggregate_time_char[timeSizeBuffer];
        
        if (timespec2str(label_time_char, sizeof(label_time_char), &label_time) != 0) {
            strcpy(label_time_char,"invalid time");
        }

        if (timespec2str(aggregate_time_char, sizeof(aggregate_time_char), &aggregate_time) != 0) {
            strcpy(aggregate_time_char,"invalid time");
        }

        char cummulative_communication_char[comSizeBuffer];
        sprintf(cummulative_communication_char,"%0.02f MB", (double)cummulative_communication/1e6);

        char label_communication_char[comSizeBuffer];
        sprintf(label_communication_char,"%0.02f MB", (double)(cummulative_communication - reference.cummulative_communication)/1e6);

        timer_start.erase(iterator);

        std::string levels[logLevels];
        for (size_t i = 0; i < logLevels; i++) {
            levels[i].assign(" ");
        }

        levels[level-1].assign(label);
        std::string comma(",");

        std::string msg("‹");

        for (size_t i = 0; i < logLevels; i++) {
            msg += (levels[i] + comma);
        }

        msg += std::string(label_communication_char) 
            + comma 
            + std::string(cummulative_communication_char)
            + comma 
            + std::string(label_time_char)
            + comma 
            + std::string(aggregate_time_char)
            + comma;

        DBG(msg);
    }
};

//==============================================================================
/** All Public Data */

//==============================================================================

        /*
            All Public Data On Current Protocol
        */
        class PublicData {
            public:
                // Parameters
                size_t degree;
                size_t p;
                size_t q;
                size_t sigma;
                size_t lambda;
                size_t tau_limit_bit;
                size_t tau;
                size_t numCandidates;
                size_t modulusIdx;
                size_t ptNumberEvaluation;

                // Protocol
                // Alphas
                std::vector<uint64_t> ps;

                // Indices
                std::vector<uint64_t> indicesPS;
                std::vector<uint64_t> indicesCAN;
                std::vector<uint64_t> indicesGCD;

                // Round 3: Keygen
                std::vector<uint64_t> A;
                std::vector<uint64_t> bi;
                std::vector<uint64_t> b;

                // Round 4: Pre-sieving and triple generation
                std::vector<uint64_t> ci_1;
                std::vector<uint64_t> ci_2;

                // Round 5: Pre-sieving and triple generation
                std::vector<uint64_t> ci_1_prime;
                std::vector<uint64_t> ci_2_prime;
                std::vector<uint64_t> c_1;
                std::vector<uint64_t> c_2;

                // Round 6: Pre-sieving and triple generation
                bool special;
                std::vector<uint64_t> di;
                std::vector<uint64_t> c_1_prime;
                std::vector<uint64_t> c_2_prime;

                // Round 7: Candidates generation
                std::vector<uint64_t> ax_shares;
                std::vector<uint64_t> by_shares;
                std::vector<uint64_t> fixedCoefs;
                std::vector<uint64_t> coefsCAN;
                std::vector<uint64_t> cans;
                std::vector<uint64_t> prodcans;

                // Round 8: Candidates generation
                std::vector<uint64_t> ax;
                std::vector<uint64_t> by;
                std::vector<uint64_t> axby;

                // Round 10: Jacobi/Sigma Protocol
                //std::vector<uint64_t> sigmar;
                //std::vector<uint64_t> sigmax;
                //std::vector<uint64_t> sigmablind;
                //std::vector<uint64_t> sigmabeta;

                // Round 11: GCD
                std::vector<uint64_t> ax_shares_GCD;
                std::vector<uint64_t> by_shares_GCD;
                std::vector<uint64_t> coefsGCD;
                std::vector<uint64_t> gcds;
                std::vector<uint64_t> prodgcds;

                // Round 12: GCD result 
                std::vector<uint64_t> axGCD;
                std::vector<uint64_t> byGCD;
                std::vector<uint64_t> axbyGCD;

                // Roots Of Unity
                std::vector<uint64_t> roots;

                // Sigma protocol
                // std::vector<mpz_class> sigmaa;
                // std::vector<mpz_class> sigmae;
                // std::vector<mpz_class> sigmaz;
                std::vector<uint64_t> sigmaeGCD;
                std::vector<uint64_t> sigmazGCD;

                std::vector<uint64_t> log2Ds;
                std::vector<uint64_t> Cs;
                std::vector<uint64_t> Ds;

            private:
                friend class boost::serialization::access;

                // Implements seralization and deserialization of public data
                template <class Archive>
                void serialize(Archive& ar, const unsigned int version)
                {
                    ar & ps;

                    ar & degree;
                    ar & p;
                    ar & q;
                    ar & sigma;
                    ar & lambda;
                    ar & tau_limit_bit;
                    ar & tau;
                    ar & numCandidates;
                    ar & modulusIdx;
                    ar & roots;
                    ar & ptNumberEvaluation;

                    ar & special;
                    ar & A;
                    ar & bi;
                    ar & b;

                    ar & ci_1;
                    ar & ci_2;

                    ar & ci_1_prime;
                    ar & ci_2_prime;
                    ar & c_1;
                    ar & c_2;

                    ar & di;
                    ar & c_1_prime;
                    ar & c_2_prime;

                    ar & ax_shares;
                    ar & by_shares;
                    ar & fixedCoefs;
                    ar & coefsCAN;
                    ar & cans;
                    ar & prodcans;

                    ar & ax;
                    ar & by;
                    ar & axby;

                    ar & ax_shares_GCD;
                    ar & by_shares_GCD;
                    ar & gcds;
                    ar & coefsGCD;
                    ar & prodgcds;

                    ar & axGCD;
                    ar & byGCD;
                    ar & axbyGCD;

                    ar & roots;

                    // ar & sigmaa;
                    // ar & sigmae;
                    // ar & sigmaz;
                    ar & sigmaeGCD;
                    ar & sigmazGCD;
                    //ar & sigmar;
                    //ar & sigmax;
                    //ar & sigmablind;
                    //ar & sigmabeta;

                    ar & indicesPS;
                    ar & indicesCAN;
                    ar & indicesGCD;

                    ar & log2Ds;
                    ar & Cs;
                    ar & Ds;
                }
        };

        /*
            All Secret Data On Current Protocol
        */
        class SecretData {
            public:
                // Status
                bool isAvailable;

                // Protocol
                    // Round 3: Keygen
                    std::vector<uint64_t> siP;
                    std::vector<uint64_t> eiP;
                    std::vector<uint64_t> q_r3;

                    // Round 4: Pre-sieving and triple generation
                    std::vector<uint64_t> xi;
                    std::vector<uint64_t> ux;
                    std::vector<uint64_t> vx;
                    std::vector<uint64_t> wx;
                    std::vector<uint64_t> vxP;
                    std::vector<uint64_t> wxP;
                    std::vector<uint64_t> q_r4_1;
                    std::vector<uint64_t> q_r4_2;

                    // Round 5:
                    std::vector<uint64_t> zp; // Noise added to secure multiplication
                    std::vector<uint64_t> yi;
                    std::vector<uint64_t> zi;
                    std::vector<uint64_t> uz;
                    std::vector<uint64_t> vz;
                    std::vector<uint64_t> vzP;
                    std::vector<uint64_t> wz;
                    std::vector<uint64_t> wzP;
                    std::vector<uint64_t> q_r5_1;
                    std::vector<uint64_t> q_r5_2;

                    // Round 6:
                    std::vector<uint64_t> rNoise;
                    std::vector<uint64_t> q_r6;

                    // Round 7:
                    std::vector<uint64_t> x_sharesCAN;
                    std::vector<uint64_t> x_sharesPS;
                    std::vector<uint64_t> y_sharesCAN;
                    std::vector<uint64_t> y_sharesPS;
                    std::vector<uint64_t> z_sharesCAN;
                    std::vector<uint64_t> q_p_prod_r7;
                    std::vector<uint64_t> q_p_r7;
                    std::vector<uint64_t> q_q_prod_r7;
                    std::vector<uint64_t> q_q_r7;

		    // Round 8:
                    std::vector<uint64_t> q_r8;

		    // Round 11: 
                    std::vector<uint64_t> x_sharesGCD;
                    std::vector<uint64_t> y_sharesGCD;
                    std::vector<uint64_t> z_sharesGCD;
                    std::vector<uint64_t> q_p_prod_r11;
                    std::vector<uint64_t> q_q_prod_r11;
                    std::vector<uint64_t> q_pq_r11;
                    std::vector<uint64_t> r_CRTs;
                    std::vector<uint64_t> q_r_r11;

		    // Round 12:
                    std::vector<uint64_t> q_r12;

            // Sigma protocol:

            // Bit Decompositions
                    std::vector<std::vector<bool>> bitDecompositions;
                    std::vector<uint64_t> XplusCs;
                    std::vector<uint64_t> XplusDs;
                    std::vector<uint64_t> Xs;

		    // std::vector<mpz_class> sigmar;
                    std::vector<uint64_t> sigmarGCD;
                    std::vector<uint64_t> sigmaxGCD;
                    std::vector<uint64_t> sigmaqGCD;
                    std::vector<uint64_t> exp_q;
                    //std::vector<uint64_t> sigma_blind;

            private:
                friend class boost::serialization::access;

                // Implements seralization and deserialization of secret data
                template <class Archive>
                void serialize(Archive& ar, const unsigned int version)
                {
                    ar & siP;
                    ar & eiP;
                    ar & q_r3;

                    ar & xi;
                    ar & ux;
                    ar & vx;
                    ar & wx;
                    ar & vxP;
                    ar & wxP;
                    ar & q_r4_1;
                    ar & q_r4_2;

                    ar & zp;
                    ar & yi;
                    ar & zi;
                    ar & uz;
                    ar & vz;
                    ar & wz;
                    ar & vzP;
                    ar & wzP;
                    ar & q_r5_1;
                    ar & q_r5_2;

                    ar & rNoise;
                    ar & q_r6;

                    ar & x_sharesCAN;
                    ar & x_sharesPS;
                    ar & y_sharesCAN;
                    ar & y_sharesPS;
                    ar & y_sharesGCD;
                    ar & z_sharesCAN;
                    ar & q_p_prod_r7;
                    ar & q_p_r7;
                    ar & q_q_prod_r7;
                    ar & q_q_r7;

                    ar & q_r8;

                    ar & x_sharesGCD;
                    ar & y_sharesGCD;
                    ar & z_sharesGCD;
                    ar & q_p_prod_r11;
                    ar & q_q_prod_r11;
                    ar & q_pq_r11;
                    ar & q_r_r11;
                    ar & r_CRTs;

                    ar & q_r12;

                    // ar & sigmar;
                    ar & sigmarGCD;
                    ar & sigmaxGCD;
                    ar & sigmaqGCD;
                    ar & exp_q;
                    //ar & sigma_blind;

                    ar & isAvailable;
                    ar & bitDecompositions;
                    ar & Xs;
                    ar & XplusCs;
                    ar & XplusDs;

                }
        };

//==============================================================================
/** Misc Helper Functions */

//==============================================================================

void configureEasyLogging(std::string tag = "") {
    el::Configurations c;
    c.setToDefault();
    std::string config = std::string("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s.%g} %msg\n LOG_FLUSH_THRESHOLD = 1\n ENABLED = true\n FILENAME = \"test") + tag + std::string(".log\" \n TO_FILE = true \n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n LOG_FLUSH_THRESHOLD = 1\n *INFO:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false");
    c.parseFromText(config.c_str());
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
    el::Loggers::setDefaultConfigurations(c, true);
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
}

template<typename T>
void sendToFile(T& data, std::string filename) {
    std::stringstream ss;
    std::ofstream fileHandler(filename.c_str(), std::ios::binary);

    boost::archive::binary_oarchive oa(ss);
    oa << data;

    std::string msg = ss.str();

    fileHandler << msg;
    // Destructors will close both file handler and boost serializer
    system((std::string("touch ") + filename + std::string("_complete")).c_str());
}

template<typename T>
void fetchFromFile(T& data, std::string filename) {
    std::ifstream fileHandler(filename.c_str(), std::ios::binary);
    boost::archive::binary_iarchive ia(fileHandler);
    ia >> data;

    // Destructors will close both file handler and boost serializer 
}

inline bool doesFileExist(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

inline bool doesTouchFileExist(const std::string& name) {
  struct stat buffer;   
  return (stat ((name+std::string("_complete")).c_str(), &buffer) == 0); 
}

/* Data Structures */
class internalError : public std::runtime_error {using std::runtime_error::runtime_error;};
class failedTest : public std::runtime_error {using std::runtime_error::runtime_error;};

//==============================================================================
/** Stream overload for a std::array<T>. */

//==============================================================================

using namespace std::chrono_literals;

namespace ligero
{

    // Helper Class to Express Indices
    struct CandidateIndices {
        std::vector<size_t> ps;
        std::vector<size_t> can;
        std::vector<size_t> gcd;
    };

    // Constants
    constexpr auto kDefaultRegistrationWaitTime = 30s;
    constexpr int kTrialsThreshold = 10;
    constexpr int kMinAmountOfPartiesRequired = 2;
    constexpr int kDefaultZeroMQIOThreads = 15;
    constexpr int kMaxBitsRandomness = 128;
    constexpr double kDefaultChiStd = 3.2;
    const size_t kJacobiNumberOfRandomValues = 1000;

    constexpr auto primesPS = 43918; // 43200;
    constexpr auto primesCAN = 21600;
    constexpr auto primesGCD = 18; //180;
    constexpr auto primesTOTAL = primesPS + primesCAN + primesGCD;

    // Constant strings
    const std::string kPublicKeyInvalidGetterInvocation = "Getters should be used only after proper object instantiation.";

    // Application-wide typedefs
    typedef boost::multiprecision::mpz_int MPInt;
    typedef boost::multiprecision::cpp_int_backend<64, 64, boost::multiprecision::signed_magnitude, boost::multiprecision::checked, void> Int64_backend;
    typedef boost::multiprecision::number<Int64_backend> Int64;

    typedef boost::multiprecision::cpp_int_backend<512, 512, boost::multiprecision::signed_magnitude, boost::multiprecision::checked, void> Int512_backend;
    typedef boost::multiprecision::number<Int512_backend> Int512;

    typedef boost::uuids::uuid SocketId;

    typedef std::pair<std::pair<std::vector<mpz_class>, std::vector<mpz_class>>, mpz_class> PairOfVecGCD;


    // Protocol configuration values
    template <typename T>
    class ProtocolConfig {
    public:
        ProtocolConfig() : initialized(false) {};

        ProtocolConfig(int rows, int cols, double chiStd, T q, int numParties, int packingFactor)
            : _rows(rows), _cols(cols), _chiStd(chiStd), _q(q), _packingFactor(packingFactor), _numParties(numParties), initialized(true) {};

        ProtocolConfig(size_t numParties, size_t numPrimesP, size_t numPrimesQ, size_t degree,
                       size_t sigma, size_t lambda, size_t tauLimit, size_t pbs)
            : _numParties(numParties), _nb_primes_p(numPrimesP), _nb_primes_q(numPrimesQ),
              _degree(degree), _sigma(sigma), _lambda(lambda), _tau_bit_limit(tauLimit), _pbs(pbs), initialized(true)
            { }

        int rows() const { return THROW_IF_NOT_INITIALIZED(_rows); }
        int cols() const { return THROW_IF_NOT_INITIALIZED(_cols); }
        double chiStd() const { return THROW_IF_NOT_INITIALIZED(_chiStd); }
        T q() const { return THROW_IF_NOT_INITIALIZED(_q); }
        int packingFactor() const { return THROW_IF_NOT_INITIALIZED(_packingFactor); }

        // RingLWE getters
        size_t numParties() const { return THROW_IF_NOT_INITIALIZED(_numParties); }
        size_t& numParties() { return _numParties; }
        size_t numPrimesP() const { return _nb_primes_p; }
        size_t numPrimesQ() const { return _nb_primes_q; }
        size_t degree() const { return _degree; }
        size_t sigma() const { return _sigma; }
        size_t lambda() const { return _lambda; }
        size_t tauLimitBit() const { return _tau_bit_limit; }
        size_t pbs() const { return _pbs; }

    private:
        friend class boost::serialization::access;
        // Implements seralization and deserialization for ProtocolConfig
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & _rows;
            ar & _cols;
            ar & _chiStd;
            ar & _q;
            ar & _numParties;
            ar & _packingFactor;

            ar & _nb_primes_p;
            ar & _nb_primes_q;
            ar & _degree;
            ar & _sigma;
            ar & _lambda;
            ar & _tau_bit_limit;
            ar & _pbs;

            initialized = true;
        }

        int _rows;
        int _cols;
        double _chiStd;
        T _q;
        int _packingFactor;
        bool initialized;
        // RingLWE parameters
        size_t _numParties;
        size_t _nb_primes_p;
        size_t _nb_primes_q;
        size_t _degree;
        size_t _sigma;
        size_t _lambda;
        size_t _tau_bit_limit;
        size_t _pbs;
    };

    class PartyCommitment {
    public:

        PartyCommitment() {};

        PartyCommitment(
                std::string ipAddress,
                std::string ai,
                std::string seed1,
                std::string seed2
                )
            : _ipAddress(ipAddress), _aiHash(ai), _seed1hash(seed1), _seed2hash(seed2)
            { }

        std::string ipAddress() const { return _ipAddress; }
        std::string aiHash() const { return _aiHash; }
        std::string seed1hash() const { return _seed1hash; }
        std::string seed2hash() const { return _seed2hash; }

    private:

        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & _ipAddress;
            ar & _aiHash;
            ar & _seed1hash;
            ar & _seed2hash;
        }

        std::string _ipAddress;
        std::string _aiHash;
        std::string _seed1hash;
        std::string _seed2hash;
    };

    // Enumerating message types
    enum class MessageType : int {
        ABORT_PROTOCOL_VERSION_MISMATCH,
        ID_PARTY = PROTOCOL_VERSION_NUMBER,
        THROUGHPUT_TEST_START,
        THROUGHPUT_TEST_STOP,
        THROUGHPUT_TEST_SURVIVE,
        THROUGHPUT_TEST_KICKOUT,
        PROTOCOL_CONFIG,
        ENCRYPTION,
        ASSIGNMENT_P1,
        ASSIGNMENT_PN,
        ENCRYPTED_X_SHARES,
        ENCRYPTED_X_VALUE,
        ENCRYPTED_XY_SHARES,
        ENCRYPTED_XY_VALUE,
        ENCRYPTED_XY_PLUS_Z_SHARES,
        ENCRYPTED_XY_PLUS_Z_VALUE,
        PARTIAL_XY_MINUS_Z_SHARES,
        XY_MINUS_Z_VALUE,
        AX_BY_SHARES,
        GCD_AX_BY_SHARES,
        AX_BY_VALUE,
        AXB_MINUS_BYA_SHARES,
        AXB_MINUS_BYA_VALUE,
        AB_VALUE,
        ENCRYPTED_PS_A_SHARES,
        ENCRYPTED_PS_A_VALUES,
        ENCRYPTED_PS_E_SHARES,
        ENCRYPTED_PS_E_VALUES,
        PS_PARTIAL_E_SHARES,
        PS_SIEVING_FLAGS,
        REPEAT_PRESIEVING_SHARES,
        PERFORM_MODULUS_GENERATION,
        PRIME_SHARES,
        PRIME_CANDIDATE,
        MODULUS_SHARES,
        PARTIAL_MODULUS_SHARES,
        ENCRYPTED_MODULUS_CANDIDATE,
        MODULUS_CANDIDATE,
        POST_SIEVE,
        GAMMA_COMMITMENT,
        GAMMA_SHARES,
        GAMMA_RANDOM_SEED_SHARES,
        GAMMA_RANDOM_SEED_VALUE,
        GAMMA_VALUE,
        EXPONENTIATED_GAMMA_VALUE,
        DISCARD_FLAGS,
        GCD_RAND_SHARES,
        GCD_RAND_VALUE,
        GCD_Z_SHARES,
        GCD_Z_PARTIAL_SHARES,
        GCD_Z_VALUE,
        PUBLIC_KEY_A_SHARES,
        PUBLIC_KEY_A_VALUE,
        PUBLIC_KEY_B_SHARES,
        PUBLIC_KEY_B_VALUE,
        MESSAGE,
        P_CLEAR_DEBUG,
        Q_CLEAR_DEBUG,
        MUTHU_ACK,
        FOUND_MODULI,
        NO_MODULI,
        PROTOCOL_RESTART,
        END,
        GATHER_PROOFS,
        GATHER_PUBLIC_DATA
    };

std::vector<std::string> msgs = {
        "ABORT_PROTOCOL_VERSION_MISMATCH",
        "ID_PARTY = PROTOCOL_VERSION_NUMBER",
        "THROUGHPUT_TEST_START",
        "THROUGHPUT_TEST_STOP",
        "THROUGHPUT_TEST_SURVIVE",
        "THROUGHPUT_TEST_KICKOUT",
        "PROTOCOL_CONFIG",
        "ENCRYPTION",
        "ASSIGNMENT_P1",
        "ASSIGNMENT_PN",
        "ENCRYPTED_X_SHARES",
        "ENCRYPTED_X_VALUE",
        "ENCRYPTED_XY_SHARES",
        "ENCRYPTED_XY_VALUE",
        "ENCRYPTED_XY_PLUS_Z_SHARES",
        "ENCRYPTED_XY_PLUS_Z_VALUE",
        "PARTIAL_XY_MINUS_Z_SHARES",
        "XY_MINUS_Z_VALUE",
        "AX_BY_SHARES",
        "GCD_AX_BY_SHARES",
        "AX_BY_VALUE",
        "AXB_MINUS_BYA_SHARES",
        "AXB_MINUS_BYA_VALUE",
        "AB_VALUE",
        "ENCRYPTED_PS_A_SHARES",
        "ENCRYPTED_PS_A_VALUES",
        "ENCRYPTED_PS_E_SHARES",
        "ENCRYPTED_PS_E_VALUES",
        "PS_PARTIAL_E_SHARES",
        "PS_SIEVING_FLAGS",
        "REPEAT_PRESIEVING_SHARES",
        "PERFORM_MODULUS_GENERATION",
        "PRIME_SHARES",
        "PRIME_CANDIDATE",
        "MODULUS_SHARES",
        "PARTIAL_MODULUS_SHARES",
        "ENCRYPTED_MODULUS_CANDIDATE",
        "MODULUS_CANDIDATE",
        "POST_SIEVE",
        "GAMMA_COMMITMENT",
        "GAMMA_SHARES",
        "GAMMA_RANDOM_SEED_SHARES",
        "GAMMA_RANDOM_SEED_VALUE",
        "GAMMA_VALUE",
        "EXPONENTIATED_GAMMA_VALUE",
        "DISCARD_FLAGS",
        "GCD_RAND_SHARES",
        "GCD_RAND_VALUE",
        "GCD_Z_SHARES",
        "GCD_Z_PARTIAL_SHARES",
        "GCD_Z_VALUE",
        "PUBLIC_KEY_A_SHARES",
        "PUBLIC_KEY_A_VALUE",
        "PUBLIC_KEY_B_SHARES",
        "PUBLIC_KEY_B_VALUE",
        "MESSAGE",
        "P_CLEAR_DEBUG",
        "Q_CLEAR_DEBUG",
        "MUTHU_ACK",
        "FOUND_MODULI",
        "NO_MODULI",
        "PROTOCOL_RESTART",
        "END",
        "GATHER_PROOFS",
        "GATHER_PUBLIC_DATA"
    };

    /** Overload the addition operator for std::array types. */
    template <typename T, size_t N>
    std::array<T, N> operator+ (const std::array<T, N>& a, const std::array<T, N>& b) {
        std::array<T, N> result {};

        for (size_t i = 0; i < N; ++i)
            result[i] = a[i] + b[i];

        return result;
    }

    /** Overload the addition operator for std::array types with a scalar. */
    template <typename T, size_t N>
    std::array<T, N> operator+ (const std::array<T, N>& a, const T& b) {
        std::array<T, N> result {};

        for (size_t i = 0; i < N; ++i)
            result[i] = a[i] + b;

        return result;
    }

    /** Overload the multiplication operator for std::array types. */
    template <typename T, size_t N>
    std::array<T, N> operator* (const std::array<T, N>& a, const std::array<T, N>& b) {
        std::array<T, N> result {};

        for (size_t i = 0; i < N; ++i)
            result[i] = a[i] * b[i];

        return result;
    }

    /** Overload the multiplication operator for std::array types with a scalar. */
    template <typename T, size_t N>
    std::array<T, N> operator* (const std::array<T, N>& a, const T& b) {
        std::array<T, N> result {};

        for (size_t i = 0; i < N; ++i)
            result[i] = a[i] * b;

        return result;
    }

    /** Helper function for discarding values from a ColumnVector */
    template <typename NumberType>
    std::vector<NumberType> discardCandidates (
        const std::vector<NumberType>& candidates,
        const boost::dynamic_bitset<>& discardFlags
    ) {
        int numCandidates = candidates.size();
        int numSurvivors = numCandidates - discardFlags.count();

        assert(static_cast<int>(discardFlags.size()) == numCandidates);

        std::vector<NumberType> survivors (numSurvivors);

        // Discard any candidate that failed the test
        for (int j = 0, s = 0; j < candidates.size(); ++j) {
            if (discardFlags[j] == 0) {
                survivors[s++] = candidates[j];
            }
        }

        return survivors;
    }   

//==============================================================================
// Helper function for debug

        void displayVector(std::string s, std::vector<size_t> data) {

            DBG("----------------------------------------------------------");
            DBG(s << ":" << data.size() << " elements"); // Header

            if (data.size() > MAX_LINES_DISPLAY) {
                for (size_t i = 0; i < MAX_LINES_DISPLAY/2; i++) {
                    DBG(data[i]);
                }

                DBG("...");

                for (size_t i = MAX_LINES_DISPLAY/2; i>0 ; i--) {
                    DBG(data[data.size() - i]);
                }
            } else {
                for (size_t i = 0;i < data.size(); i++) {
                    DBG(data[i]);
                }
            }
            DBG("----------------------------------------------------------");
        }

        template <typename FieldT>
        void displayVector(std::string s, std::vector<FieldT> data) {

            DBG("----------------------------------------------------------");
            DBG(s << ":" << data.size() << " elements"); // Header

            if (data.size() > MAX_LINES_DISPLAY) {
                for (size_t i = 0; i < MAX_LINES_DISPLAY/2; i++) {
                    DBG(data[i].getValue());
                }

                DBG("...");

                for (size_t i = MAX_LINES_DISPLAY/2; i>0 ; i--) {
                    DBG(data[data.size() - i].getValue());
                }
            } else {
                for (size_t i = 0;i < data.size(); i++) {
                    DBG(data[i].getValue());
                }
            }
            DBG("----------------------------------------------------------");
        }

        template <typename FieldT>
        void displayArray(std::string s, FieldT *data, size_t size) {

            DBG("----------------------------------------------------------");
            DBG(s << ":" << size << " elements"); // Header

            if (size > MAX_LINES_DISPLAY) {
                for (size_t i = 0; i < MAX_LINES_DISPLAY/2; i++) {
                    DBG(data[i].getValue());
                }

                DBG("...");

                for (size_t i = MAX_LINES_DISPLAY/2; i>0 ; i--) {
                    DBG(data[size - i].getValue());
                }
            } else {
                for (size_t i = 0;i < size; i++) {
                    DBG(data[i].getValue());
                }
            }
            DBG("----------------------------------------------------------");
        }

        constexpr size_t upperBoundBits = 180ull;

        std::pair<mpz_class, uint64_t> roundToUpperPowerTwoLarge(mpz_class n) {
            size_t log2d = upperBoundBits;

            mpz_class powerOfN(1);
            mpz_class higherPower;

            for (uint64_t i = 0; i < upperBoundBits; i++) {
                if ((n & powerOfN) > 0) {
                    log2d = i;
                    higherPower = powerOfN;
                }
                powerOfN *= 2;
                }

            if (higherPower != n) {
                log2d++;
                higherPower *= 2;
            }

            return std::pair<mpz_class, uint64_t>(higherPower, log2d);
        }

        std::pair<uint64_t, uint64_t> roundToUpperPowerTwo(uint64_t n) {

            size_t log2d = 64;
            size_t higherPower;

            for (unsigned long long i = 0; i < 64; i++) {
                if ((n & (1ull << i)) > 0) {
                    log2d = i;
                    higherPower = (1ull << i);
                }
            }

            if (higherPower != n) {
                log2d++;
                higherPower *= 2;
            }

            return std::pair<uint64_t, uint64_t>(higherPower, log2d);
        }

        std::tuple<mpz_class, mpz_class, mpz_class, uint64_t> rephraseAs2Nct(mpz_class d) {

            // Rephrase the problem using the nearest 2^n 
            auto [ D, log2D ] = ligero::roundToUpperPowerTwoLarge(2*d);
            mpz_class C = mpz_class(D-d);

            // defining the constraint
            if (D == 0) {
                if (d == 0) {
                    DBG("Limit set to 0!");
                } else {
                    DBG("Overflow");
                }

                assert(d > 0);

            } else {
                DBG("log2D = " << log2D);
                DBG("d = " << d);
                DBG("C = " << C);
                DBG("D = " << D);
                DBG("d + C = D");
            }

            return std::tuple<mpz_class, mpz_class, mpz_class, uint64_t>(d,C,D,log2D);
        }

} // namespace ligero

//==============================================================================
namespace boost {
    namespace serialization {

        /** Implements serialization and deserialization for an MPInt type to
         *  an arbitrary Boost Archive.
         */
        template <class Archive>
        void serialize(Archive& ar, ligero::MPInt& x, const unsigned int version)
        {
            if (Archive::is_saving::value)
            {
                std::string s = mpz_get_str(nullptr, 32, x.backend().data());
                ar & s;
            }
            else
            {
                std::string s;

                ar & s;
                mpz_set_str(x.backend().data(), s.c_str(), 32);
            }
        }

        /** Implements serialization and deserialization for an NFL polynomial type to
         *  an arbitrary Boost Archive.
         */
        template <class Archive, typename T, size_t Degree, size_t NbModuli>
        void serialize(Archive& ar, nfl::poly<T, Degree, NbModuli>& x, const unsigned int version)
        {
            if (Archive::is_saving::value)
            {
                std::stringstream ss(std::stringstream::out | std::stringstream::binary);
                x.serialize_manually(ss);
                std::string s = ss.str();
                std::cout << "serialize: " << ss.str() << std::endl << std::endl;
                ar & s;
            }
            else
            {
                std::string s;
                ar & s;
                std::istringstream ssi(s);
                x.deserialize_manually(ssi);
            }
        }

        /** Implements serialization and deserialization for an NFL polynomial type to
         *  an arbitrary Boost Archive.
         */
        template <class Archive>
        void serialize(Archive& ar, mpz_class& x, const unsigned int version)
        {
            if (Archive::is_saving::value)
            {
                std::string s = x.get_str();
                ar & s;
            }
            else
            {
                std::string s;
                ar & s;
                x = mpz_class(s);
            }
        }
    } // namespace serialization
} // namespace boost
