#pragma once

// Disable Boost deprecated header message
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>

#include <boost/functional/hash.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/dynamic_bitset/serialization.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/unordered_map.hpp>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include <sodium/crypto_generichash.h>

#include <gmpxx.h>
#include <iostream>
#include <sstream>
#include <variant>
#include <zmq.hpp>

#include <cryptopp/files.h>
#include "cryptlib.h"
#include "eccrypto.h"
#include "ecp.h"
#include "hex.h"
#include "osrng.h"
#include "oids.h"

#include "easylogging++.h"
#include "nfl.hpp"
#include <sstream>

using namespace std::chrono_literals;
using namespace CryptoPP;

#define MAX_LINES_DISPLAY 10

typedef uint64_t RegularInteger;
typedef unsigned __int128 WideInteger;

/** Enum class for zero knowledge linear constraint **/
enum class LinearCombinationsStitching : uint64_t {
    Regular,
    StitchingNoScalar,      
    StitchingScalar         
};

/** Error handling */
enum class Error : size_t {
    TIMED_OUT,
    OUT_OF_SYNC,
    MODULUS_NOT_FOUND,
    DESERIALIZE_FAIL,
    TOO_FEW_PARTIES,
    TOO_MANY_RESTART,
    RESTART,
    KILLED_BY_COORDINATOR,
    UNKNOWN_ERROR
};

template <typename T, typename U = Error>
using expected = std::variant<T, U>;
using Unit = std::monostate;

/** Helper to conver Error to string error message */
std::string showError(Error e) {
    switch (e) {
        case Error::TIMED_OUT:
            return std::string("Receving timed out");
        case Error::OUT_OF_SYNC:
            return std::string("Out of synchronization");
        case Error::MODULUS_NOT_FOUND:
            return std::string("No modulus found");
        case Error::DESERIALIZE_FAIL:
            return std::string("Deserialization failed when receving message");
        case Error::TOO_FEW_PARTIES:
            return std::string("Too few parties left to continue the protocol");
        case Error::TOO_MANY_RESTART:
            return std::string("Too many restart");
        case Error::RESTART:
            return std::string("Restart the protocol");
        case Error::KILLED_BY_COORDINATOR:
            return std::string("Party killed by coordinator");
        case Error::UNKNOWN_ERROR:
        default:
            break;
    }
    return std::string("Unknown error: ") + std::to_string(size_t(e));
}

/** Check whether the result has error or a valid value */
template <typename T>
inline bool hasError(const expected<T>& result) {
    return std::holds_alternative<Error>(result);
}

/** Get error object from the result */
template <typename T>
inline Error getError(const expected<T>& result) {
    return std::get<Error>(result);
}

/** Get result from variable */
template <typename T>
inline auto getResult(const expected<T>& result) 
    -> decltype(std::get<T>(result)) 
{
    return std::get<T>(result);
}

/** Operator << overloading for output an array to stream
 * @param o Reference of output stream
 * @param arr Array to be copied
 * @return Reference of the output stream
 */
template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& o, const std::array<T, N>& arr) {
    std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
    return o;
}

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
/** Timer communication buffer size */
constexpr unsigned int comSizeBuffer = 100;
/** Timer time buffer size */
constexpr unsigned int timeSizeBuffer = 100;
/** Timer logging level */
constexpr unsigned int logLevels = 4;

// Multi-Threading Physical Constraints
#define NB_MAX_THREADS 90
#define MAX_SINGLE_THREADED_ITERATIONS 3

// Profiling

// Display
constexpr int restiveX = 1;
constexpr int restiveY = 25;

constexpr auto firstline = 5;
constexpr auto displayStride = 7;
constexpr auto displayOffset = 30;

//==============================================================================
// Profiling    
typedef struct timepoint {
    size_t cummulative_communication; /**< Cummulative communication cost */
    struct timespec timer;            /**< Structure holding an interval  */
} timepoint;

void printAtXY(int x,int y, std::string output)
{
    std::string display = std::string("%c[%d;%df") + output;
    printf(display.c_str(),0x1B,y,x);

    std::string cursor = std::string("%c[%d;%df\n");
    printf(cursor.c_str(),0x1B,restiveY,restiveX);
}

class timers {
protected:
    struct timespec anchor; /**< Anchor */
    std::unordered_map<std::string, timepoint>
    timer_start;     /**< Map of started timer */
    bool active = false; /**< Activation flag */

public:
    void initialize_timer() {
        this->active = true;
        clock_gettime(CLOCK_MONOTONIC, &this->anchor);
    }

   /** Start a new timer
    * @param level Timer indentation level in display
    * @param label Name of this timer
    * @param cummulative_communication Commulative communication cost
    */
    void begin(size_t level, std::string label, size_t cummulative_communication) {
        struct timespec time;

        if (!this->active) return;

        clock_gettime(CLOCK_MONOTONIC, &time);
        timepoint tmp;
        tmp.timer = time;
        tmp.cummulative_communication = cummulative_communication;
        timer_start.insert(std::pair<std::string, timepoint>(label, tmp));
    }

   /** Helper function for calculating the difference between two time points
    * @param[in] start The start time point
    * @param[in] stop The stop time point
    * @param[out] result Result we are writing to
    */
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

   /** Pretty printing a timespec structure to string
    * @param[out] buf Output buffer. Needs to store 30 characters
    * @param[in] len Length of the buffer
    * @param[in] ts Time interval
    * @return On success return 0
    */
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

   /** Stop a running timer. NOTE: The whole program will exit if wrong label
    * is given
    * @param level Timer indentation level in display
    * @param label Timer name, must be same with the running timer
    * @param cummulative_communication Cummulative communication
    */
    void end(size_t level, std::string label, size_t cummulative_communication, bool higherLogLevel = false) {
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

        levels[level - 1].assign(label);
        std::string comma(",");

        std::string msg("‹");

        for (size_t i = 0; i < logLevels; i++) {
            msg += (levels[i] + comma);
        }

        msg += std::string(label_communication_char) + comma +
               std::string(cummulative_communication_char) + comma +
               std::string(label_time_char) + comma +
               std::string(aggregate_time_char) + comma;

        if (higherLogLevel) LOG(INFO) << msg;
        else DBG(msg);
    }

    /** Simple helper to reset all timers */
    void reset() {
        LOG(INFO) << "Reset all timers as well";
        timer_start.clear();
        initialize_timer();
    }
};

        /*
            All Public Data On Current Protocol
        */
        class SigmaProtocolPublicData {
            public:
                std::vector<uint64_t> sigmaeGCD;
                std::vector<uint64_t> sigmazGCD;
                std::vector<uint64_t> sigmaaGCD;
                std::vector<uint64_t> sigmagGCD;

                // x, y, and z shares for failed candidates
                std::vector<uint64_t> xi;
                std::vector<uint64_t> yi;
                std::vector<uint64_t> zi;

            private:
                friend class boost::serialization::access;

                // Implements seralization and deserialization of public data
                template <class Archive>
                void serialize(Archive& ar, const unsigned int version)
                {
                    ar & sigmaeGCD;
                    ar & sigmazGCD;
                    ar & sigmaaGCD;
                    ar & sigmagGCD;

                    ar & xi;
                    ar & yi;
                    ar & zi;
                }
        };

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
                // size_t numCandidates;
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
                std::vector<uint64_t> finalModuli_GCD;

                // Roots Of Unity
                std::vector<uint64_t> roots;

                // Sigma protocol
                // std::vector<mpz_class> sigmaa;
                // std::vector<mpz_class> sigmae;
                // std::vector<mpz_class> sigmaz;
                std::vector<uint64_t> sigmaeGCD;
                std::vector<uint64_t> sigmazGCD;
                mpz_class gammaSeed;

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
                    ar & modulusIdx;
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
        	        ar & finalModuli_GCD;

                    ar & indicesPS;
                    ar & indicesCAN;
                    ar & indicesGCD;

                    ar & log2Ds;
                    ar & Cs;
                    ar & Ds;

                    ar & gammaSeed;
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
                    std::vector<uint64_t> ss_GCD;
                    std::vector<uint64_t> q_r12;

            // Cross-Moduli Integrity:
                    uint64_t              alpha;
                    uint64_t              W1;
                    uint64_t              W2;
                    std::vector<std::vector<std::vector<std::vector<uint64_t>>>> ai;
                    // ai are structured as follows:
                    // - one entry for the iteration idx
                    // - one entry per attached variable (ei,si,ux,uz)
                    // - one entry per block (degree/_blockSize)
                    // - _blockSize entries per block 

            // Bit Decompositions
                    std::vector<std::vector<bool>> bitDecompositions;
                    std::vector<uint64_t> XplusCs;
                    std::vector<uint64_t> XplusDs;
                    std::vector<uint64_t> Xs;

                    std::vector<uint64_t> sigmarGCD;
                    std::vector<uint64_t> sigmaxGCD;
                    std::vector<uint64_t> sigmaqGCD;
                    std::vector<uint64_t> exp_q;
                    std::vector<uint64_t> expqGCD;

            private:
                friend class boost::serialization::access;

                // Implements seralization and deserialization of secret data
                template <class Archive>
                void serialize(Archive& ar, const unsigned int version)
                {
                    ar & siP;
                    ar & eiP;
                    ar & q_r3;

                    ar & alpha;
                    ar & W1;
                    ar & W2;
                    ar & ai;

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

            		ar & ss_GCD;
                    ar & q_r12;

                    // ar & sigmar;
                    ar & sigmarGCD;
                    ar & sigmaxGCD;
                    ar & sigmaqGCD;
                    ar & exp_q;
                    ar & expqGCD;

                    ar & isAvailable;
                    ar & bitDecompositions;
                    ar & Xs;
                    ar & XplusCs;
                    ar & XplusDs;

    }
};

/** Misc Helper Functions */

/** Compute positive remainder for mpz_class type values */
mpz_class positiveRemainder(mpz_class & a, mpz_class & b) {
        return ((a % b + b) % b);
};

/** Extract data from the vector of mpz_class for a given set of indexes
 * into a vector of uint64_t values within a given ring defined by NbPrimesQ parameter */
template <size_t NbPrimesQ>
void dataExtractSubVector(std::vector<uint64_t> &assignee, std::vector<mpz_class> &target, std::vector<size_t> &idxs) {
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


/** Extract data from the vector of mpz_class values into a vector of 
 * uint64_t values within a given ring defined by NbPrimesQ parameter */
void dataExtractVector(std::vector<uint64_t> &assignee, std::vector<mpz_class> &target) {
        assignee.resize(target.size());
        for (size_t idx = 0; idx < target.size(); idx++) {
            assignee[idx] = mpz_get_ui(target[idx].get_mpz_t());
        }
    };


/** Extract data from the vector of mpz_class values in CRT representation 
 * into a vector of uint64_t values within a given ring defined by NbPrimesQ parameter */
template <size_t NbPrimesQ>
void dataExtractVectorCRT(std::vector<uint64_t> &assignee, std::vector<mpz_class> &target) {
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

/** Extract data from the value in Q into a vector of uint64_t values */
template <typename T, size_t Degree, size_t NbPrimesQ>
void dataExtractQ(std::vector<uint64_t> &assignee, nfl::poly_p<T, Degree, NbPrimesQ> & target) {
    assignee.assign(
            target.poly_obj().data(),
            target.poly_obj().data() + (size_t)target.nmoduli * (size_t)target.degree);
};

/** Serialization helper */
template<typename T>
std::string serialize(T obj) {
    std::stringstream ss;
    boost::archive::binary_oarchive oa(ss);

    oa << obj;

    return (ss.str());
};

/** Configure logging levels */
void configureEasyLogging(const char* fileChar = "default") {
    std::string file(fileChar);

    el::Configurations c;
    c.setToDefault();
    std::string config = std::string("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s.%g} %msg\n LOG_FLUSH_THRESHOLD = 1\n ENABLED = true\n FILENAME = \"") + file + std::string(".log\" \n TO_FILE = true \n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n LOG_FLUSH_THRESHOLD = 1\n *INFO:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false");
    c.parseFromText(config.c_str());
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
    el::Loggers::setDefaultConfigurations(c, true);
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
}

/** Helper function for serializing data in binary format and save to file
 * @param data Serializable data
 * @param filename File to save
 */
template <typename T>
void sendToFile(T& data, std::string filename) {
    std::stringstream ss;
    std::ofstream fileHandler(filename.c_str(), std::ios::binary);

    boost::archive::binary_oarchive oa(fileHandler);
    oa << data;
}

/** Helper function for serializing data in binary format and save to file
 * @param data string
 * @param filename File to save
 */
void sendToFileString(const std::string& data, std::string filename) {
    std::ofstream fileHandler(filename.c_str(), std::ios::binary);
    fileHandler << data;
}

/** Helper function for saving data to file
 * @param data message
 * @param filename File to save
 */
void sendToFileRaw(zmq::message_t* msg, std::string filename) {
    std::stringstream ss;
    std::ofstream fileHandler(filename.c_str(), std::ios::binary);

    ss.write(msg->data<char>(), msg->size());
    std::string msgOut = ss.str();
    fileHandler << msgOut;
}

/** Helper function for deserializing data from file
 * @param data[out] Destination to save fetched data
 * @param filename[in] File name to read from
 */
template <typename T>
void fetchFromFile(T& data, std::string filename) {
    std::ifstream fileHandler;
    fileHandler.open(filename.c_str(), std::ios::in | std::ios::binary);

    boost::archive::binary_iarchive ia(fileHandler);
    ia >> data;

    // Destructors will close both file handler and boost serializer
}

/** Helper function for checking if file exist
 * @param name File name to check
 * @return True if file exist
 */
inline bool doesFileExist(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

/** Helper function for checking if complete file exist
 * @param name File name to check
 * @return True if file exist
 */
inline bool doesTouchFileExist(const std::string& name) {
    struct stat buffer;
    return (stat((name + std::string("_complete")).c_str(), &buffer) == 0);
}

/** Runtime error */
class internalError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};
/** Runtime error */
class failedTest : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/** Stream overload for a std::array<T>. */

namespace ligero {
    // Constants
    constexpr auto kDefaultRegistrationWaitTime = 30s;
    constexpr int kTrialsThreshold = 10;
    constexpr int kMinAmountOfPartiesRequired = 2;
    constexpr int kDefaultZeroMQIOThreads = 15;
    constexpr int kMaxBitsRandomness = 128;
    constexpr double kDefaultChiStd = 3.2;
    const size_t kJacobiNumberOfRandomValues = 1000;

    constexpr auto primesPS = 43917; // 43200;
    constexpr auto primesCAN = 21600;
    constexpr auto primesGCD = 19;
    constexpr auto primesTOTAL = primesPS + primesCAN + primesGCD;

    // Constant strings
    const std::string kPublicKeyInvalidGetterInvocation = "Getters should be used only after proper object instantiation.";

    // Application-wide typedefs
    typedef std::pair<std::vector<mpz_class>, std::pair<std::vector<mpz_class>, std::vector<mpz_class>>> TripleVector;
    typedef boost::multiprecision::mpz_int MPInt;
    typedef boost::multiprecision::cpp_int_backend<64, 64, boost::multiprecision::signed_magnitude, boost::multiprecision::checked, void> Int64_backend;
    typedef boost::multiprecision::number<Int64_backend> Int64;

    typedef boost::multiprecision::cpp_int_backend<512, 512, boost::multiprecision::signed_magnitude, boost::multiprecision::checked, void> Int512_backend;
    typedef boost::multiprecision::number<Int512_backend> Int512;

    typedef boost::uuids::uuid SocketId;

    typedef std::pair<std::pair<std::vector<mpz_class>, std::vector<mpz_class>>, mpz_class> PairOfVecGCD;

    template <typename T, size_t Degree, size_t NbPrimesQ>
    class BadEventData {

        using Q = nfl::poly_p<T, Degree, NbPrimesQ>;
        public:
        BadEventData() : 
            jacobiSeed(mpz_class(0)),
            jacobiAndGCDSeed(mpz_class(0)),
            si(nfl::poly_p<T, Degree, NbPrimesQ>{0}),
            ei(nfl::poly_p<T, Degree, NbPrimesQ>{0}),
            bi(nfl::poly_p<T, Degree, NbPrimesQ>{0}),
            gcdRX(mpz_class(0)),
            gcdSS(std::vector<mpz_class>()) {}

        BadEventData(mpz_class _jacobiSeed, 
                mpz_class _jacobiAndGCDSeed,
                nfl::poly_p<T, Degree, NbPrimesQ> _si,
                nfl::poly_p<T, Degree, NbPrimesQ> _ei,
                nfl::poly_p<T, Degree, NbPrimesQ> _bi,
                mpz_class _gcdRX,
                std::vector<mpz_class> _gcdSS) 
            : jacobiSeed(_jacobiSeed),
            jacobiAndGCDSeed(_jacobiAndGCDSeed),
            si(_si),
            ei(_ei),
            bi(_bi),
            gcdRX(_gcdRX),
            gcdSS(_gcdSS)
            {}

        nfl::poly_p<T, Degree, NbPrimesQ> si, ei, bi;
        mpz_class jacobiSeed;
        mpz_class jacobiAndGCDSeed;
        mpz_class gcdRX;
        std::vector<mpz_class> gcdSS;

        private:
        friend class boost::serialization::access;
        // Implements seralization and deserialization of public data
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & jacobiSeed;
            ar & jacobiAndGCDSeed;
            ar & si;
            ar & ei;
            ar & bi;
            ar & gcdRX;
            ar & gcdSS;
        }

    };

    enum class ProtocolMode : int {
        NORMAL,
        REPLAY,
        RECORD
    };

    // Protocol configuration values
    template <typename T>
    class ProtocolConfig {
    public:
        ProtocolConfig() : initialized(false) {};

        ProtocolConfig(int rows, int cols, double chiStd, T q, int numParties, int packingFactor)
            : _rows(rows), _cols(cols), _chiStd(chiStd), _q(q), _packingFactor(packingFactor), _numParties(numParties), initialized(true) {};

        ProtocolConfig(size_t numParties, size_t numPrimesP, size_t numPrimesQ, size_t degree,
                       size_t sigma, size_t lambda, size_t tauLimit, size_t pbs, ProtocolMode protocolMode, ECDSA<ECP, SHA256>::PublicKey publicKey)
            : _numParties(numParties), _nb_primes_p(numPrimesP), _nb_primes_q(numPrimesQ),
              _degree(degree), _sigma(sigma), _lambda(lambda), _tau_bit_limit(tauLimit), _pbs(pbs), _protocolMode(protocolMode), _publicKey(publicKey), initialized(true)
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
        ProtocolMode protocolMode() const { return _protocolMode; }
        ECDSA<ECP, SHA256>::PublicKey publicKey() const { return _publicKey; }

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

            ar & _protocolMode;

            if (Archive::is_saving::value) {
                StringSink ss(_publicKeyString);
                _publicKey.Save(ss);
                ar & _publicKeyString;
            } else {
                ar & _publicKeyString;
                StringSource ss(_publicKeyString, true);
                _publicKey.Load(ss);
            }
            initialized = true;
        }

        // TODO: deprecate unused parameter in RingLWE
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
        ProtocolMode _protocolMode;
        std::string _publicKeyString;
        ECDSA<ECP, SHA256>::PublicKey _publicKey;
    };

    // All the data each party commits on registration
    class PartyCommitment {
    public:

        PartyCommitment() {};

        PartyCommitment(
                std::string ipAddress,
                std::vector<std::string> earlyW,
                std::string ai,
                std::string seed1,
                std::string seed2,
                ECDSA<ECP, SHA256>::PublicKey publicKey
                )
            : _ipAddress(ipAddress), _earlyW(earlyW), _aiHash(ai), _seed1hash(seed1), _seed2hash(seed2), _publicKey(publicKey)
            { }

        std::string ipAddress() const { return _ipAddress; }
        std::vector<std::string> earlyW() const { return _earlyW; }
        std::string aiHash() const { return _aiHash; }
        std::string seed1hash() const { return _seed1hash; }
        std::string seed2hash() const { return _seed2hash; }
        std::string publicKeyString() const { return _publicKeyString; }
        ECDSA<ECP, SHA256>::PublicKey publicKey() const { return _publicKey; }

    private:

        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & _ipAddress;
            ar & _earlyW;
            ar & _aiHash;
            ar & _seed1hash;
            ar & _seed2hash;

            if (Archive::is_saving::value) {
                StringSink ss(_publicKeyString);
                _publicKey.Save(ss);
                ar & _publicKeyString;
            } else {
                ar & _publicKeyString;
                StringSource ss(_publicKeyString, true);
                _publicKey.Load(ss);
            }
        }

        std::string                 _ipAddress;
        std::vector<std::string>    _earlyW;
        std::string                 _aiHash;
        std::string                 _seed1hash;
        std::string                 _seed2hash;
        std::string                 _publicKeyString;
        ECDSA<ECP, SHA256>::PublicKey _publicKey;
    };

    // Enumerating message types
    enum class MessageType : int {
        ABORT_PROTOCOL_VERSION_MISMATCH,
        ID_PARTY = PROTOCOL_VERSION_NUMBER,
        REPLAY_PROTOCOL_SHARES,
        REPLAY_PROTOCOL_RESPONSE,
        RECORD_PROTOCOL_SHARES,
        RECORD_PROTOCOL_RESPONSE,
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
        SYNCHRONIZE_NOW,
        FOUND_MODULI,
        NO_MODULI,
        NO_MODULI_VERIFICATION_SHARES,
        NO_MODULI_VERIFICATION_FINISHED,
        BAD_EVENT_DATA_SHARES,
        BAD_EVENT_DATA_RESPONSE,
        PROTOCOL_RESTART,
        PROTOCOL_KILLED,
        END,
        GATHER_PROOFS,
        GATHER_PUBLIC_DATA,
        GATHER_SIGMA_DATA,
        VERIFIER_READY,
        PUBLIC_DATA_TO_VERIFIER,
        PROOF_TO_VERIFIER,
        GATHER_REPORTS,
        VERIFIER_REPORT,
        RECYCLE_VERIFIER
    };

    // Stringified rounds names
    std::vector<std::string> msgs = {
        "ABORT_PROTOCOL_VERSION_MISMATCH",
        "ID_PARTY",
        "REPLAY_PROTOCOL_SHARES",
        "REPLAY_PROTOCOL_RESPONSE",
        "RECORD_PROTOCOL_SHARES",
        "RECORD_PROTOCOL_RESPONSE",
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
        "SYNCHRONIZE_NOW",
        "FOUND_MODULI",
        "NO_MODULI",
        "NO_MODULI_VERIFICATION_SHARES",
        "NO_MODULI_VERIFICATION_FINISHED",
        "BAD_EVENT_DATA_SHARES",
        "BAD_EVENT_DATA_RESPONSE",
        "PROTOCOL_RESTART",
        "PROTOCOL_KILLED",
        "END",
        "GATHER_PROOFS",
        "GATHER_PUBLIC_DATA",
        "GATHER_SIGMA_DATA",
        "VERIFIER_READY",
        "PUBLIC_DATA_TO_VERIFIER",
        "PROOF_TO_VERIFIER",
        "GATHER_REPORTS",
        "VERIFIER_REPORT",
        "RECYCLE_VERIFIER"
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

        // TODO: If a[i] is an EIgen::Matrix type, we need cwiseProduct here
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

    // Helper Class to Express Indices
    struct CandidateIndices {
        std::vector<size_t> ps;
        std::vector<size_t> can;
        std::vector<size_t> gcd;
    };


    CandidateIndices getCandidateIndices(size_t idx, const std::vector<size_t>& index_candidates) {

            std::vector<size_t> ps;
            size_t i;
            for (i = 0; i < primesPS; ++i) {
                if (index_candidates[i] == idx) {
                    ps.push_back(i);
                }
            }

            std::vector<size_t> can;
            for (; i < primesPS + primesCAN; ++i) {
                if (index_candidates[i] == idx) {
                    can.push_back(i);
                }
            }

            std::vector<size_t> gcd;
            for (; i < primesPS + primesCAN + primesGCD; ++i) {
                if (index_candidates[i] == -2) { // -2 is a special value used only for GCD. Only one candidate uses it
                    gcd.push_back(i);
                }
            }
            return CandidateIndices({ps, can, gcd});
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

        constexpr size_t upperBoundBits = 2080ull;

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
        
        std::string extractAll(std::vector<uint64_t>& v) {
            std::string tmp("");
            if (v.size()>0) {
                for (size_t idx = 0; idx < v.size()-1; idx++) {
                    tmp += std::to_string(v[idx]) + std::string(",");
                }
                tmp += std::to_string(v[v.size()-1]);
            } else {
                tmp.assign("-");
            }

            return tmp;
        }

        void dumpSigmaProtocolPublicData(SigmaProtocolPublicData& spd) {
            std::ofstream out("dump_sigmaProtocolPublicData");

            out << extractAll(spd.sigmaeGCD) << std::endl;
            out << extractAll(spd.sigmazGCD) << std::endl;

            out.close();
        }

        void dumpPublicData(PublicData& pd) {
            std::ofstream out("dump_publicData");

            out << pd.degree << std::endl;
            out << pd.p << std::endl;
            out << pd.q << std::endl;
            out << pd.sigma << std::endl;
            out << pd.lambda << std::endl;
            out << pd.tau_limit_bit << std::endl;
            out << pd.tau << std::endl;
            out << pd.modulusIdx << std::endl;
            out << pd.ptNumberEvaluation << std::endl;
            out << pd.special << std::endl;
            out << extractAll(pd.ps) << std::endl;
            out << extractAll(pd.indicesPS) << std::endl;
            out << extractAll(pd.indicesCAN) << std::endl;
            out << extractAll(pd.indicesGCD) << std::endl;
            out << extractAll(pd.A) << std::endl;
            out << extractAll(pd.bi) << std::endl;
            out << extractAll(pd.b) << std::endl;
            out << extractAll(pd.ci_1) << std::endl;
            out << extractAll(pd.ci_2) << std::endl;
            out << extractAll(pd.ci_1_prime) << std::endl;
            out << extractAll(pd.ci_2_prime) << std::endl;
            out << extractAll(pd.c_1) << std::endl;
            out << extractAll(pd.c_2) << std::endl;
            out << extractAll(pd.di) << std::endl;
            out << extractAll(pd.c_1_prime) << std::endl;
            out << extractAll(pd.c_2_prime) << std::endl;
            out << extractAll(pd.ax_shares) << std::endl;
            out << extractAll(pd.by_shares) << std::endl;
            out << extractAll(pd.coefsCAN) << std::endl;
            out << extractAll(pd.cans) << std::endl;
            out << extractAll(pd.prodcans) << std::endl;
            out << extractAll(pd.ax) << std::endl;
            out << extractAll(pd.by) << std::endl;
            out << extractAll(pd.axby) << std::endl;
            out << extractAll(pd.ax_shares_GCD) << std::endl;
            out << extractAll(pd.by_shares_GCD) << std::endl;
            out << extractAll(pd.coefsGCD) << std::endl;
            out << extractAll(pd.gcds) << std::endl;
            out << extractAll(pd.prodgcds) << std::endl;
            out << extractAll(pd.axGCD) << std::endl;
            out << extractAll(pd.byGCD) << std::endl;
            out << extractAll(pd.axbyGCD) << std::endl;
            out << extractAll(pd.finalModuli_GCD) << std::endl;
            out << extractAll(pd.roots) << std::endl;
            out << extractAll(pd.log2Ds) << std::endl;
            out << extractAll(pd.Cs) << std::endl;
            out << extractAll(pd.Ds) << std::endl;

            out.close();
        }

        struct coordinates {
            size_t bloc;
            size_t col;
            bool   isAvailable;
        };

} // namespace ligero

//==============================================================================
namespace boost {
    namespace serialization {

        template <class Archive>
        void serialize(Archive& ar, boost::uuids::uuid& id, const unsigned int version)
        {
            if (Archive::is_saving::value)
            {
                std::string s = boost::lexical_cast<std::string>(id);
                ar & s;
            }
            else
            {
                std::string s;
                ar & s;
                id = boost::lexical_cast<boost::uuids::uuid>(s);
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
