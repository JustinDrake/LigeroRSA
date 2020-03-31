/**
 * Finite field F[p], for prime p of NumberType type.
 *
 * A proth prime is of the form: p = a * 2^n + 1
 */

#pragma once

#include "Common.hpp"
#include "Math.hpp"

#include "nfl/params.hpp"

using namespace nfl;

namespace ligero {

/*
  Modular Multiplication for 64 bit integers
  */

struct mulmod {
uint64_t operator()(uint64_t x, uint64_t y, size_t cm) const
{
    auto const p = params<uint64_t>::P[cm];
    auto const shift = params<uint64_t>::kModulusRepresentationBitsize;
    ASSERT_STRICTMOD((x<p) && (y<p));
    __uint128_t  res = (__uint128_t )x * y;
    __uint128_t  q = ((__uint128_t )params<uint64_t>::Pn[cm] * (res >> shift)) + (res<<2) ;
    uint64_t r  = res - (q>>shift) * p;
    if (r >= p) { r -= p; }
    ASSERT_STRICTMOD(r == ((__uint128_t )(x) * y) % p);
    return r;
}
};

// ================================================================================================
// Class Structure

// Assumptions on NumberType:
// - A << overload must be available for this type.

template<typename NumberType, const NumberType& modulus>
class Fp {
    protected:
        static NumberType modulus_;
        NumberType value;

    public:
        using underlyingType = NumberType;

        Fp () {
        }

        Fp(NumberType p) { 
            this->value = ligero::math::mod(p, modulus);
        }

        bool operator==(const Fp& other) const;
        bool operator!=(const Fp& other) const;

        Fp& operator=(const NumberType& other);
        Fp& operator=(const Fp& other);
        Fp& operator+=(const Fp& other);
        Fp& operator-=(const Fp& other);
        Fp& operator*=(const Fp& other);

        Fp operator+(const Fp& other) const;
        Fp operator-(const Fp& other) const;
        Fp operator*(const Fp& other) const;
        Fp operator-() const;

        Fp operator^(const unsigned long pow) const;

        Fp inverse() const;
        void print() const;

        static Fp one() {return Fp(NumberType(1));}
        NumberType getValue() const {return this->value;}
        NumberType& getValueRef() const {return this->value;}
        static NumberType getModulus() {return modulus_;}
};

template<typename NumberType, const NumberType& modulus>
class Fpp : public Fp<NumberType, modulus> {
    public:
        static NumberType a;
        static NumberType n;
        static NumberType omega;
        static std::vector<NumberType> rootsOfUnity;
    public:
        Fpp(Fp<NumberType, modulus> p) {
            this->value = ligero::math::mod(p.getValue(), modulus);
        }

        Fpp() {
        }

        Fpp(NumberType p) {
            this->value = ligero::math::mod(NumberType(p), modulus);
        }

        Fpp& operator=(const Fpp<NumberType, modulus>& other);
        Fpp& operator=(const NumberType& other);

        void populate(NumberType value) {this->value = value;};
        static void preprocessing();                        // Initialize directly through Prime

        // Accessors
        static NumberType getOmega(size_t domain_size);     // Domain size must divide composite domain size 
        static NumberType getOmega();                       // Domain size must divide composite domain size 
        static size_t getCompositeDomainSize();
        static NumberType getN();
        static Fpp one() {return Fpp(NumberType(1));}

        static Fpp<NumberType, modulus> randomVector();
        static void randomVector(Fpp<NumberType, modulus> *,size_t, bool non_zero = false);
        static void randomPartialMatrix(Fpp<NumberType, modulus> *,size_t, size_t, size_t, bool non_zero = false);
        static void randomPartialMatrixColumns(Fpp<NumberType, modulus> *,size_t, size_t, size_t, bool non_zero = false);
        static void zeroSumNonZRandomVector(Fpp<NumberType, modulus> *vector, size_t colsRandomness);

        // Boost Serialization
        friend class boost::serialization::access;

        private:
            // Implements seralization and deserialization
            template <class Archive>
            void serialize(Archive& ar, const unsigned int version)
            {
                ar & this->value;
            }
};

template<const RegularInteger& modulus>
class Fpp_Fixed : public Fpp<RegularInteger, modulus> {

    public:
        static size_t modulusIdx_;

        Fpp_Fixed(RegularInteger p) { 
            this->value = ligero::math::mod(p, modulus);
        }

        Fpp_Fixed() {
        }

        static Fpp_Fixed one() {return Fpp_Fixed(RegularInteger(1));}

        Fpp_Fixed& operator+=(const Fpp_Fixed& other);
        Fpp_Fixed& operator-=(const Fpp_Fixed& other);
        Fpp_Fixed& operator*=(const Fpp_Fixed& other);

        Fpp_Fixed operator+(const Fpp_Fixed& other) const;
        Fpp_Fixed operator-(const Fpp_Fixed& other) const;
        Fpp_Fixed operator*(const Fpp_Fixed& other) const;
        Fpp_Fixed operator-() const;

        Fpp_Fixed operator^(const unsigned long pow) const;

        Fpp_Fixed inverse() const;

        static void preprocessing(size_t overrideN = 0);                        // Initialize directly through Prime

};

// ===============================================================================================
// Operators 

template<typename NumberType, const NumberType& modulus>
bool Fp<NumberType, modulus>::operator==(const Fp& other) const
{
 return (this->value == other.value);
}

template<typename NumberType, const NumberType& modulus>
bool Fp<NumberType, modulus>::operator!=(const Fp& other) const
{
 return (this->value != other.value);
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus>& Fp<NumberType, modulus>::operator=(const Fp& other)
{
 this->value = ligero::math::mod(other.value, modulus);
 return *this;
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus>& Fp<NumberType, modulus>::operator+=(const Fp& other)
{
 this->value = ligero::math::mod(this->value + other.value, modulus);
 return *this;
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus>& Fp<NumberType, modulus>::operator-=(const Fp& other)
{
 this->value = ligero::math::mod(this->value - other.value + modulus, modulus);
 return *this;
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus>& Fp<NumberType, modulus>::operator*=(const Fp& other)
{
 this->value = ligero::math::mod(this->value * other.value, modulus);
 return *this;
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus> Fp<NumberType, modulus>::operator+(const Fp& other) const
{
 return (ligero::math::mod(this->value + other.value, modulus));
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus> Fp<NumberType, modulus>::operator-(const Fp<NumberType, modulus>& other) const
{
 return (ligero::math::mod(this->value - other.value + modulus, modulus));
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus> Fp<NumberType, modulus>::operator*(const Fp& other) const
{
 return (ligero::math::mod(this->value * other.value, modulus));
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus> Fp<NumberType, modulus>::operator-() const
{
 return (ligero::math::mod(-this->value, modulus));
}

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus> Fp<NumberType, modulus>::operator^(const unsigned long pow) const
{
 return (ligero::math::powm<NumberType>(this->value,pow,modulus));
}


// Optimized for 64-bits numbers
// =============================
template<const RegularInteger& modulus>
Fpp_Fixed<modulus>& Fpp_Fixed<modulus>::operator+=(const Fpp_Fixed& other)
{
 this->value = static_cast<RegularInteger>(ligero::math::mod(static_cast<WideInteger>(this->value) + static_cast<WideInteger>(other.value), static_cast<WideInteger>(modulus)));
 return *this;
}

template<const RegularInteger& modulus>
Fpp_Fixed<modulus>& Fpp_Fixed<modulus>::operator-=(const Fpp_Fixed& other)
{
 this->value = static_cast<RegularInteger>(ligero::math::mod(static_cast<WideInteger>(this->value) - static_cast<WideInteger>(other.value) + static_cast<WideInteger>(modulus), static_cast<WideInteger>(modulus)));
 return *this;
}

template<const RegularInteger& modulus>
Fpp_Fixed<modulus>& Fpp_Fixed<modulus>::operator*=(const Fpp_Fixed& other)
{
 this->value = (mulmod{}(this->value, other.value, this->modulusIdx_));
 return *this;
}

template<const RegularInteger& modulus>
Fpp_Fixed<modulus> Fpp_Fixed<modulus>::operator+(const Fpp_Fixed& other) const
{
 return (static_cast<RegularInteger>(ligero::math::mod(static_cast<WideInteger>(this->value) + static_cast<WideInteger>(other.value), static_cast<WideInteger>(modulus))));
}

template<const RegularInteger& modulus>
Fpp_Fixed<modulus> Fpp_Fixed<modulus>::operator-(const Fpp_Fixed& other) const
{
 return (static_cast<RegularInteger>(ligero::math::mod(static_cast<WideInteger>(this->value) - static_cast<WideInteger>(other.value) + static_cast<WideInteger>(modulus), static_cast<WideInteger>(modulus))));
}

template<const RegularInteger& modulus>
Fpp_Fixed<modulus> Fpp_Fixed<modulus>::operator*(const Fpp_Fixed& other) const
{
 return (mulmod{}(this->value, other.value, this->modulusIdx_));
}

template<const RegularInteger& modulus>
Fpp_Fixed<modulus> Fpp_Fixed<modulus>::operator^(const unsigned long pow) const
{
 return (ligero::math::powm<RegularInteger>(static_cast<RegularInteger>(this->value),static_cast<RegularInteger>(pow),static_cast<RegularInteger>(modulus)));
}

template<const RegularInteger& modulus>
void Fpp_Fixed<modulus>::preprocessing(size_t overrideN) {

    DBG("------------------------------------------------------------------");

    // Infer factorization, naive algorithm
    uint64_t nLocal = 0, nFinal = 0;
    RegularInteger aFinal;

    if (overrideN >0) {

        nFinal = overrideN;
        aFinal = (modulus-1)/(pow(RegularInteger(2),nFinal));

        DBG("preprocessing " << aFinal << " " << nFinal << " " << modulus << std::endl);

        assert((pow(RegularInteger(2),nFinal)*aFinal + 1 == modulus));

    } else {

        RegularInteger val = ligero::math::pow(RegularInteger(2),static_cast<uint64_t>(nLocal));

        while (val < modulus) {
            if ((ligero::math::mod(RegularInteger(modulus-1), val)) == RegularInteger(0)) nFinal = nLocal;
            
            nLocal++;
            val = ligero::math::pow(RegularInteger(2),static_cast<uint64_t>(nLocal));
            if (ligero::math::pow(MPInt(2),static_cast<uint64_t>(nLocal))>MPInt(modulus)) break; 
        }

        aFinal = (modulus-1)/(pow(RegularInteger(2),nFinal));
        while ((nLocal>0)&&(pow(RegularInteger(2),nFinal)*aFinal + 1 != modulus)) {
            nLocal--;
            aFinal = (modulus-1)/pow(RegularInteger(2),nLocal);
        }

        if (nFinal==0) throw std::runtime_error("Misformed Proth Prime Field, cannot find a p = a.2^n + 1 decomposition.\n");
    }

    Fpp_Fixed::modulus_ = modulus;
    Fpp_Fixed::n = nFinal;
    Fpp_Fixed::a = aFinal;

    DBG("Proth Identified p = " << Fpp_Fixed::a << " x 2^" << Fpp_Fixed::n << " + 1");

    // find a primitive root of unity 
    Fpp_Fixed elt(RegularInteger(1));
    Fpp_Fixed pow_elt(RegularInteger(1));
    RegularInteger size = Fpp_Fixed<modulus>::getCompositeDomainSize();
    DBG("Composite Domain Size = " << size);
    DBG("Searching for a primitive root of unity:");

    bool foundPrimitiveRoot = false;
    Fpp_Fixed<modulus> minusOne(modulus-RegularInteger(1));
    DBG("-1 = " << minusOne.getValue() << " mod " << modulus);
    RegularInteger val = ligero::math::pow(RegularInteger(2),static_cast<uint64_t>(nFinal));

    do {
        if (ligero::math::mod(RegularInteger(modulus-1), val) == RegularInteger(0)) nFinal = nLocal;
        elt += Fpp_Fixed(RegularInteger(1));
        pow_elt = ligero::math::powm(elt.getValue(),RegularInteger(aFinal*size/2),modulus);

        if (pow_elt == minusOne) {
            foundPrimitiveRoot = true;
            break;
        }
    } while ((!foundPrimitiveRoot) && (elt != minusOne) );
    
    if (!foundPrimitiveRoot) throw std::runtime_error("Could not find a primitive root of unity.\n");

    Fpp_Fixed::omega = ligero::math::powm<RegularInteger>(elt.getValue(),Fpp_Fixed::a,modulus);
    DBG("ω = " << Fpp_Fixed::omega);
    DBG("");
    DBG("Pre-processing completed - Proth Prime Field Constructed.");
    // 
    DBG("------------------------------------------------------------------");

};

// ===============================================================================================
// Public Methods 

template<typename NumberType, const NumberType& modulus>
Fp<NumberType, modulus> Fp<NumberType, modulus>::inverse() const {
    return ligero::math::powm<NumberType>(static_cast<NumberType>(this->value),static_cast<NumberType>(-1),static_cast<NumberType>(modulus));
};

template<const RegularInteger& modulus>
Fpp_Fixed<modulus> Fpp_Fixed<modulus>::inverse() const {
    return Fpp_Fixed<modulus>(ligero::math::powm(static_cast<RegularInteger>(this->value),int(-1),modulus));
};

// ===============================================================================================
// Static Members

template<typename NumberType, const NumberType& modulus> 
NumberType Fp<NumberType, modulus>::modulus_;

/**
 * Finite field F[Pp], for Proth prime p of MPInt type.
 *
 * A proth prime is of the form: p = a * 2^n + 1
 */


// ================================================================================================
// Class Structure


// ===============================================================================================
// Static Members 

template<typename NumberType, const NumberType& modulus> 
NumberType Fpp<NumberType, modulus>::a;

template<typename NumberType, const NumberType& modulus> 
NumberType Fpp<NumberType, modulus>::n;

template<typename NumberType, const NumberType& modulus> 
NumberType Fpp<NumberType, modulus>::omega;

template<typename NumberType, const NumberType& modulus> 
std::vector<NumberType> Fpp<NumberType, modulus>::rootsOfUnity;

template<const RegularInteger& modulus> 
size_t Fpp_Fixed<modulus>::modulusIdx_;

template<typename NumberType, const NumberType& modulus>
Fpp<NumberType, modulus>& Fpp<NumberType, modulus>::operator=(const Fpp& other)
{
 this->value = ligero::math::mod(other.value, modulus);
 return *this;
}

template<typename NumberType, const NumberType& modulus>
Fpp<NumberType, modulus>& Fpp<NumberType, modulus>::operator=(const NumberType& other)
{
 this->value = ligero::math::mod(other, modulus);
 return *this;
}

// ===============================================================================================
// Static Members s

template<typename NumberType, const NumberType& modulus> 
void Fpp<NumberType, modulus>::preprocessing() {

    // 
    DBG("------------------------------------------------------------------");

    // Infer factorization, naive algorithm
    uint64_t nLocal = 0, nFinal = 0;

    NumberType val = ligero::math::pow(NumberType(2),static_cast<uint64_t>(nLocal));

    while (val < modulus) {
        if ((ligero::math::mod(NumberType(modulus-1), val)) == NumberType(0)) nFinal = nLocal;
        
        nLocal++;
        val = ligero::math::pow(NumberType(2),static_cast<uint64_t>(nLocal));
        if (ligero::math::pow(MPInt(2),static_cast<uint64_t>(nLocal))>MPInt(modulus)) break; 
    }

    NumberType aFinal = (modulus-1)/(pow(NumberType(2),nFinal));
    while ((nLocal>0)&&(pow(NumberType(2),nFinal)*aFinal + 1 != modulus)) {
        nLocal--;
        aFinal = (modulus-1)/pow(NumberType(2),nLocal);
    }

    if (nFinal==0) throw std::runtime_error("Misformed Proth Prime Field, cannot find a p = a.2^n + 1 decomposition.\n");

    Fpp::modulus_ = modulus;
    Fpp::n = nFinal;
    Fpp::a = aFinal;

    DBG("Proth Identified p = " << a << " x 2^" << n << " + 1");

    // find a primitive root of unity 
    Fpp elt(NumberType(1));
    Fpp pow_elt(NumberType(1));
    NumberType size = Fpp<NumberType, modulus>::getCompositeDomainSize();
    DBG("Composite Domain Size = " << size);
    DBG("Searching for a primitive root of unity:");

    bool foundPrimitiveRoot = false;
    Fpp<NumberType, modulus> minusOne(modulus-NumberType(1));
    DBG("-1 = " << minusOne.getValue() << " mod " << modulus);
    val = ligero::math::pow(NumberType(2),static_cast<uint64_t>(nFinal));

    do {
        if (ligero::math::mod(NumberType(modulus-1), val) == NumberType(0)) nFinal = nLocal;
        elt += Fpp(NumberType(1));
        pow_elt = ligero::math::powm(elt.getValue(),NumberType(aFinal*size/2),modulus);

        if (pow_elt == minusOne) {
            foundPrimitiveRoot = true;
            break;
        }
    } while ((!foundPrimitiveRoot) && (elt != minusOne) );
    
    if (!foundPrimitiveRoot) throw std::runtime_error("Could not find a primitive root of unity.\n");

    Fpp::omega = ligero::math::powm<NumberType>(elt.getValue(),Fpp::a,modulus);
    DBG("ω = " << Fpp::omega);
    DBG("");
    DBG("Pre-processing completed - Proth Prime Field Constructed.");

    Fpp::rootsOfUnity.reserve(static_cast<size_t>(size));
    Fpp::rootsOfUnity.push_back(1);

    pow_elt = Fpp(1);
    
    Fpp p_omega(omega);
    for (size_t idx = 1; idx< static_cast<size_t>(size) ;idx++) {
        pow_elt = pow_elt * p_omega;
        Fpp::rootsOfUnity.push_back(pow_elt.getValue());
    }

    assert(Fpp::rootsOfUnity[0] == 1);
    assert(Fpp::rootsOfUnity[1] == (p_omega^1).getValue());
    assert(Fpp::rootsOfUnity[2] == (p_omega^2).getValue());

    assert(Fpp::rootsOfUnity[size/2] == (p_omega^(size/2)).getValue());
    assert(Fpp::rootsOfUnity[size/4] == (p_omega^(size/4)).getValue());

    // 
    DBG("------------------------------------------------------------------");

};

template<typename NumberType, const NumberType& modulus> 
NumberType Fpp<NumberType, modulus>::getOmega(size_t sub_domain) {

            NumberType ret = ligero::math::powm<NumberType>(Fpp<NumberType, modulus>::omega,NumberType(Fpp::getCompositeDomainSize()/sub_domain),modulus);

            return ret;
}

template<typename NumberType, const NumberType& modulus> 
NumberType Fpp<NumberType, modulus>::getOmega() {

            return Fpp<NumberType, modulus>::omega;
}

template<typename NumberType, const NumberType& modulus> 
size_t Fpp<NumberType, modulus>::getCompositeDomainSize() {
            return static_cast<size_t>(ligero::math::pow(NumberType(2),static_cast<uint64_t>(Fpp<NumberType, modulus>::n)));
}

template<typename NumberType, const NumberType& modulus> 
NumberType Fpp<NumberType, modulus>::getN() {
            return Fpp<NumberType, modulus>::n;
}

// Computes a non-null random vector
template<typename NumberType, const NumberType& modulus>
void Fpp<NumberType, modulus>::randomVector(Fpp<NumberType, modulus> *array, size_t elements, bool non_zero) {
    randomPartialMatrix(array, 1, elements, elements, non_zero);
}

// Computes a non-null random matrix
template<typename NumberType, const NumberType& modulus>
void Fpp<NumberType, modulus>::randomPartialMatrix(Fpp<NumberType, modulus> *array, size_t rows, size_t cols, size_t colsRandomness, bool non_zero) {
    boost::random::independent_bits_engine<boost::random::mt19937, 64, uint64_t> engine;
    boost::random::random_device rd;
    boost::random::seed_seq seq = { rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
    engine.seed(seq);
    
    // Populate the Matrix
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j< colsRandomness; j++) {
            array[i*cols + j].populate(NumberType(ligero::math::mod(engine(), static_cast<NumberType>(modulus-1)) + (non_zero ? 1 : 0)));
        }
    }

    return;
};

// Computes a random vector
template<typename NumberType, const NumberType& modulus>
void Fpp<NumberType, modulus>::zeroSumNonZRandomVector(Fpp<NumberType, modulus> *vector, size_t colsRandomness) {
    boost::random::independent_bits_engine<boost::random::mt19937, 64, uint64_t> engine;
    boost::random::random_device rd;
    boost::random::seed_seq seq = { rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
    engine.seed(seq);
    
    // Populate the Vector
    for (size_t j = 0; j< colsRandomness; j++) {
        NumberType pick(engine());
        while ((pick == 0)||(pick >= modulus)) {
            pick = engine();
        }
        
        vector[j].populate(pick);
    }

    return;
};

// Computes a non-null random element
	template<typename NumberType, const NumberType& modulus>
	void Fpp<NumberType, modulus>::randomPartialMatrixColumns(Fpp<NumberType, modulus> *array, size_t rows, size_t cols, size_t colsRandomness, bool non_zero) {
	    boost::random::independent_bits_engine<boost::random::mt19937, 64, uint64_t> engine;
	    boost::random::random_device rd;
	    boost::random::seed_seq seq = { rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
	    engine.seed(seq);
	
	    // Populate the Matrix
	    for (size_t i = 0; i < rows; i++) {
	        for (size_t j = 0; j< colsRandomness; j++) {
	            array[j*rows + i].populate(NumberType(ligero::math::mod(engine(), static_cast<NumberType>(modulus-1)) + (non_zero ? 1 : 0)));
	        }
	    }
	
	    return;
	};
	
template<typename NumberType, const NumberType& modulus>
void Fp<NumberType, modulus>::print() const {
    std::cout << this->value;
};

}
