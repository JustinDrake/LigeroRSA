#pragma once

#include <functional>

#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <exrandom/discrete_normal_distribution.hpp>

#include <sodium/randombytes.h>
#include <gmpxx.h>
#include "Common.hpp"


/** Overriding abs for basic int type **/ 
namespace boost {
namespace multiprecision {
RegularInteger abs(RegularInteger in) {
    return in;
}
}
}

namespace ligero {
namespace math {


std::string computeHash(mpz_class val) {
    std::string s = val.get_str();
    std::string os = val.get_str();

    unsigned char hash[crypto_generichash_BYTES + 1];
    const unsigned char *message = reinterpret_cast<const unsigned char *> (s.c_str());
    crypto_generichash(hash, crypto_generichash_BYTES,
            message, s.length(),
            NULL, 0);
    hash[crypto_generichash_BYTES] = 0;
    DBG("For value os = " << os << " s = " << s << " generating hash = >>>" << hash << "<<<");
    return std::move(std::string(reinterpret_cast<char*>(hash)));
}

template <typename T, size_t Degree, size_t NbPrimesQ>
std::string computeHash(nfl::poly_p<T, Degree, NbPrimesQ> val) {
    std::stringstream ss(std::stringstream::out | std::stringstream::binary);
    val.serialize_manually(ss);
    std::string s = ss.str();

    unsigned char hash[crypto_generichash_BYTES + 1];
    const unsigned char *message = reinterpret_cast<const unsigned char *> (s.c_str());
    crypto_generichash(hash, crypto_generichash_BYTES,
            message, s.length(),
            NULL, 0);
    hash[crypto_generichash_BYTES] = 0;
    std::string result = std::string(reinterpret_cast<char*>(hash));
    return result;
}

std::string computeHash(std::pair<std::pair<std::vector<mpz_class>, std::vector<mpz_class>>, mpz_class> val) {
    mpz_class v = val.second;
    return computeHash(v);
}

template <typename T>
std::string computeHash(T val) {
    LOG(FATAL) << "computeHash for this type is not implemented yet";
}



/** Generate a vector of random numbers of 128-bit numbers **/
std::vector<mpz_class> generateRandomVector(mpz_class seed, size_t size, size_t numBits = 128) {
    std::vector<mpz_class> result(size);

    const int numBytes = lround(double(numBits) / 8);
    const size_t bufferSize = size * (numBytes + 1);

    char *buf = new char[bufferSize];

    size_t seedSize = (mpz_sizeinbase (seed.get_mpz_t(), 2) + CHAR_BIT-1) / CHAR_BIT;
    unsigned char seedArray[32U];
    std::vector<unsigned char> seedVector(seedSize);

    mpz_export(&seedVector[0], &seedSize, 1, 1, 0, 0, seed.get_mpz_t());
    for (size_t i = 0; i < seedSize && i < 32; ++i) {
        seedArray[i] = seedVector[i];
    }
    for (size_t i = seedSize; i < 32; ++i) {
        seedArray[i] = static_cast<unsigned char>(0);
    }

    //const unsigned char seed[32] = {0};
    randombytes_buf_deterministic((void *)buf, bufferSize, seedArray);

    size_t k = 0;
    for (size_t i = 0; i < result.size(); ++i) {
        assert(k < bufferSize);
        mpz_import(result[i].get_mpz_t(), numBytes, 1, 1, 0, 0, &buf[k]);
        k += numBytes;
        assert(result[i] > 0);
    }
    delete [] buf;

    return result;
}


mpz_class generateRandomValue(mpz_class seed, size_t numBits = 128) {
    mpz_class result;
    {
        std::vector t = generateRandomVector(seed, 1, numBits);
        result = t[0];
    }
    return result;
}

/** A quick wrapper around GMP's powm implementation for different types. */
mpz_class
powm (const mpz_class& b, const mpz_class& p, const mpz_class& m) {

    mpz_class x(b);
    mpz_powm(
        x.get_mpz_t(),
        x.get_mpz_t(),
        p.get_mpz_t(),
        m.get_mpz_t()
    );

    return x;
}

/** A quick wrapper around GMP's powm implementation for different types. */
template <typename NumberType>
NumberType
powm (const NumberType& b, const NumberType& p, const NumberType& m) {
    MPInt x (b);
    MPInt y (p);
    MPInt z (m);

    MPInt o;

    mpz_powm(
        o.backend().data(),
        x.backend().data(),
        y.backend().data(),
        z.backend().data()
    );

    return NumberType(o);
}

/** A quick wrapper around GMP's powm implementation for different types. */
RegularInteger
powm (const RegularInteger& b, const long p, const RegularInteger& m) {
    MPInt x (b);
    MPInt y (p);
    MPInt z (m);

//    MPInt o;

    mpz_powm(
        x.backend().data(),
        x.backend().data(),
        y.backend().data(),
        z.backend().data()
    );

    return RegularInteger(x);
}

Int64
divideAndRound(const Int64& a, const Int64& b) {
    Int64 quotient, remainder;
    boost::multiprecision::divide_qr(a, b, quotient, remainder);
    if (remainder * 2 >= b) {
        quotient++;
    }
    return quotient;
}

WideInteger pow(RegularInteger a, RegularInteger exp) {
    MPInt x(0);
    MPInt y(a);

    mpz_pow_ui(
        x.backend().data(), 
        y.backend().data(), 
        exp);

    return WideInteger(x);
}

WideInteger pow(RegularInteger a, int exp) {
    MPInt x(0);
    MPInt y(a);

    mpz_pow_ui(
        x.backend().data(), 
        y.backend().data(), 
        exp);

    return WideInteger(x);
}

template<typename NumberType>
NumberType pow(NumberType a, uint64_t exp) {
    MPInt x(0);
    MPInt y(a);

    mpz_pow_ui(
        x.backend().data(), 
        y.backend().data(), 
        exp);

    return NumberType(x);
}

MPInt
divideAndRound(const MPInt& a, const MPInt& b) {
    MPInt quotient, remainder;
    boost::multiprecision::divide_qr(a, b, quotient, remainder);
    if (remainder * 2 >= b) {
        quotient++;
    }
    return quotient;
}

int64_t
divideAndRound(int64_t a, int64_t b) {
    int64_t remainder = a % b;
    int64_t result = static_cast<int64_t>(a / b);
    if (remainder * 2 >= b) {
        result++;
    }
    return result;
}

mpz_class
mod(const mpz_class& a, const mpz_class& b) {
    return powm(a, 1, b);
}

MPInt
mod(MPInt a, MPInt b) {
    return boost::multiprecision::powm(a, 1, b);
}

RegularInteger
mod(RegularInteger a, RegularInteger b) {
// return static_cast<RegularInteger>(boost::multiprecision::powm(MPInt(a), 1, MPInt(b)));
   return static_cast<RegularInteger>((static_cast<WideInteger>(a % b) + static_cast<WideInteger>(b)) % static_cast<WideInteger>(b));
}

RegularInteger
mod(WideInteger a, WideInteger b) {
// return static_cast<RegularInteger>(boost::multiprecision::powm(MPInt(a), 1, MPInt(b)));
   return (((a % b) + b) % b);
}

Int64
mod(Int64 a, Int64 b) {
    return boost::multiprecision::powm(a, 1, b);
}

int64_t
mod(int64_t a, int64_t b) {
    return ((a % b) + b) % b;
}

template <typename NumberType>
struct ChiDist {

    // This constant defines distributions standard deviation

    std::mt19937 g;
    double chiStd;

    ChiDist(int seed, double _chiStd) : g(seed), chiStd(_chiStd){ }
    ChiDist(int seed) : g(seed), chiStd(ligero::kDefaultChiStd){ }

    ChiDist() : g(std::random_device()()), chiStd(ligero::kDefaultChiStd){}

    NumberType operator()() {
        exrandom::discrete_normal_distribution N(0, 1, lround(chiStd * 10), 10);
        return boost::multiprecision::abs(NumberType(N(g)));
    }
};

template <typename FieldT>
std::vector<FieldT>
sampleRandomVectorOverField(size_t m, FieldT q, double chiStd) {
    //std::function<typename FieldT::underlyingType> 
    auto dist = ligero::math::ChiDist<typename FieldT::underlyingType>(std::random_device()(), chiStd);
    std::vector<FieldT> e(m);   
    
    for (int i = 0; i < static_cast<int>(m); i++) {
        typename FieldT::underlyingType tmp;
        do {
            tmp = static_cast<typename FieldT::underlyingType>(dist());
        } while (tmp>q.getValue());

        e[i] = FieldT(tmp);
    }

    return e;
}

/** A quick wrapper around GMP's modular inverse implementation for different types.
 *
 *  Note that the implementation here promotes arbitrary number types to our MPInt
 *  type so that we can call into `mpz_invert` directly with mpz types. This allows
 *  us to leverage GMP's existing implementation, but obviously has overhead. May
 *  be an optimization opportunity later if we want to consider writing our own modular
 *  inverse implementation.
 */
mpz_class
mod_inverse (const mpz_class& a, const mpz_class& b) {
    mpz_t z;
    mpz_init(z);

    int r = mpz_invert(z, a.get_mpz_t(), b.get_mpz_t());

    if (r == 0)
        throw std::runtime_error("Modular inverse does not exist for the given operands.");

    auto result = mpz_class(z);
    mpz_clear(z);

    return result;
}


/** Deconstruct input `v` by Chinese Remainder Theorem given a vector of moduli
 *  in the input `alphas`.
 */
std::vector<mpz_class> crt_deconstruct (mpz_class v, const std::vector<mpz_class>& alphas)
{
    std::vector<mpz_class> ds (alphas.size());

    for (int i = 0; i < alphas.size(); ++i) {
        ds[i] = mod(v, alphas[i]);
    }

    return ds;
}

/** Reconstruct a value `v` by Chinese Remainder Theorem given a vector of its
 *  components and a vector of moduli.
 */
mpz_class crt_reconstruct (const std::vector<mpz_class>& ds, std::vector<mpz_class>& coeffs, const std::vector<mpz_class>& alphas)
{
    if (ds.size() != alphas.size())
        throw std::runtime_error("Reconstruction vector lengths don't match.");

    coeffs.resize(ds.size());

    //mpz_class p = alphas.prod();
    mpz_class p = 1;
    for (size_t i = 0; i < alphas.size(); ++i) {
        p *= alphas[i];
    }

    for (size_t i = 0; i < ds.size(); ++i) {
        mpz_class pa = p / alphas[i];
        mpz_class x = mod(pa, alphas[i]);
        coeffs[i] = pa * mod_inverse(x, alphas[i]);
    }

    mpz_class dotProduct = 0;
    for (size_t i = 0; i < ds.size(); ++i) {
        dotProduct += ds[i] * coeffs[i];
    }

    return dotProduct % p;
}

/** A quick helper function for calculating the product of a vector of numbers. */
mpz_class vectorProduct(const std::vector<mpz_class>& xs) {
    mpz_class product = 1;

    for (auto& x : xs)
        product = product * x;

    return product;
}
/* fixed version */
std::pair<std::vector<mpz_class>, std::vector<size_t>>
fixed_bucket_n_primes (int productBitThreshold, size_t degree, int tauLimitBit)
{
    int number_of_buckets = std::ceil(double(productBitThreshold)/double(tauLimitBit)); // This rounds down. 
    std::vector<mpz_class> alphas(number_of_buckets);
    std::vector<size_t> bucketSize(number_of_buckets);


    gmp_randclass grc(gmp_randinit_default);
    grc.seed(0);

    for (int i = 0; i < number_of_buckets; ++i) {
        mpz_class r = grc.get_z_bits(tauLimitBit);
        while (!boost::multiprecision::miller_rabin_test(MPInt(r.get_mpz_t()), 64) || mpz_sizeinbase(r.get_mpz_t(),2) < tauLimitBit) {
            r = grc.get_z_bits(tauLimitBit);
        }
        alphas[i] = r;
        // DBG("alphas[" << i << "] = " << r);
    }


    for (size_t i = 0; i < number_of_buckets; ++i) {
        bucketSize[i] = lrint(floor(double(degree) / double(number_of_buckets)));
        DBG("bucketSize[ " << i << "] = " << bucketSize[i]);
    }

    return std::pair{alphas, bucketSize};
}


/** Implements a prime bucketing algorithm, returning a vector of numbers such that
 *  each number is a product of subsequent primes that is at least 18bits, and such
 *  that the product of numbers in the vector is at least `productBitThreshold` bits.
 */
std::pair<std::vector<mpz_class>, std::vector<size_t>>
balanced_bucket_n_primes (int productBitThreshold, size_t degree, int tauLimitBit, int startJ = 1)
{
    int number_of_buckets = std::ceil(double(productBitThreshold)/double(tauLimitBit)); // This rounds down. 
    std::vector<mpz_class> alphas(number_of_buckets);
    std::vector<mpz_class> primes;
    mpz_class alpha = 1;
    int j = startJ;
    while(1)
    {
        alpha*= boost::math::prime(j);
        primes.push_back(boost::math::prime(j));
        if(mpz_sizeinbase(alpha.get_mpz_t(), 2) >= static_cast<unsigned int>(productBitThreshold)) break;
        j++;
    }
    DBG("Maximum size of prime = " << mpz_sizeinbase((primes[primes.size()-1]).get_mpz_t(),2));

    std::vector<double> weights(number_of_buckets);
    for ( int i = 0; i < number_of_buckets; i++)
    {
      alphas[i] = 1;
      weights[i] = 1;
    }
    double sum = 0.0;
    for ( int i = primes.size()-1;i >= 0;  i--)
    //for ( int i = 0; i < primes.size(); i++)
    {
      double maxweight = weights[0]; 
      int maxpos = 0;
      for (int j = 1; j < number_of_buckets; j++)
      {
	      mpz_class temp = alphas[j];
        if(weights[j] > maxweight && mpz_sizeinbase(temp.get_mpz_t(),2) <= tauLimitBit)
        {
          maxweight = weights[j];
          maxpos = j;
        }
      }
      alphas[maxpos] *= primes[i];
      weights[maxpos] *= (double(1.0) - (double(1)/primes[i].get_d())); 
    }
    for ( int i = 0; i < number_of_buckets; i++)
    {
      weights[i] = double(1)/weights[i];
      sum += weights[i];
    }

/*    while(!primes.empty())
    {
        double weight=1;
        auto it = primes.begin();
        alpha = *it;
        weight *= (double(1.0) - (double(1)/alpha.get_d()));
        it = primes.erase(it);
        while(!primes.empty() && mpz_sizeinbase(alpha.get_mpz_t(), 2) < tauLimitBit )
        {
            weight *= (double(1.0) - (double(1)/(primes.back()).get_d()));
            alpha *= primes.back();
            primes.pop_back();
        }
        alphas.push_back(alpha);
        weight = double(1)/weight;
        sum += weight;
        weights.push_back(weight);

    }*/

    std::vector<size_t> bucketSize(weights.size());
    size_t total = 0;
    for (size_t i = 0; i < weights.size() - 1; ++i) {
        bucketSize[i] = lround(double(weights[i] * double(degree)) / double(sum));
        total += bucketSize[i];
        DBG("bucketSize[ " << i << "] = " << bucketSize[i]);
    }
    bucketSize[bucketSize.size() - 1] = degree - total;
    DBG("bucketSize[ " << bucketSize.size() - 1 << "] = " << bucketSize[bucketSize.size() - 1]);

    return std::pair{alphas, bucketSize};
    }

}//namespace math
}//namespace ligero
