#pragma once

#include <fstream>

#include "Common.hpp"
#include "Math.hpp"

/* Naive Evalution for Testing Purpose */
template<typename FieldT>
void naiveEvaluation(FieldT *evals, FieldT *coeffs, size_t k, std::string label = "")
{
    FieldT p = FieldT::one();
    size_t idx = 0;

    // this should be equal to -1
    FieldT pow_elt(ligero::math::powm(p.getValue(),ligero::math::powm(FieldT(2).getValue(), FieldT::n, FieldT::getModulus()),FieldT::getModulus()));

    DBG("Naive Evaluation");
    DBG("p:" << p.getValue());
    DBG("n:" << FieldT::n);
    DBG("2^n:" << ligero::math::powm(FieldT(2).getValue(), FieldT::n, FieldT::getModulus()));

    DBG("omega:" << FieldT::omega);
    DBG("stride:" << FieldT::getOmega(k));
    DBG("pow:" << pow_elt.getValue());

    std::ofstream check(label);
    do 
    {
        FieldT v(0);
        
        for (size_t i = k; i--;)
        {
            v *= p;
            v += coeffs[i];
        }

        evals[idx++] = v;

        p *= FieldT::getOmega(k);
        check << static_cast<uint64_t>(p.getValue()) << std::endl;
        if (p == FieldT::one()) break;

    } while (true);

    check.close();
}

template<typename FieldT>
void radix2_FFT(FieldT *a, size_t a_size, bool isDirectReverse, std::vector<FieldT>& roots);

/** Reverse the order of all bits **/
size_t bitreverse(size_t n, const size_t l) {
    size_t r = 0;
    for (size_t k = 0; k < l; ++k)
    {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    return r;
}

template<typename FieldT>
void FFT(FieldT *a, size_t a_size, FieldT shift, std::vector<FieldT>& roots)
{
    // std::cout << "reg FFT:" << std::endl;    
    for (size_t i = 0; i < a_size; ++i)
    {
        a[i] *= (shift^i);
    }

    radix2_FFT<FieldT>(a, a_size, true, roots);
}

template<typename FieldT>
void FFT(FieldT *a, size_t a_size, std::vector<FieldT>& roots)
{
    // std::cout << "reg FFT:" << std::endl;    
    radix2_FFT<FieldT>(a, a_size, true, roots);
}

template<typename FieldT>
void iFFT(FieldT *a, size_t a_size, FieldT shift, std::vector<FieldT>& roots)
{
    // std::cout << "inverse FFT:" << std::endl;
    radix2_FFT<FieldT>(a, a_size, false, roots);

    const FieldT invshift = shift.inverse();
    const FieldT sconst = FieldT(a_size).inverse();

    for (size_t i = 0; i < a_size; ++i)
    {
        a[i] *= (invshift^i) * sconst;
    }
}

template<typename FieldT>
void iFFT(FieldT *a, size_t a_size, std::vector<FieldT>& roots)
{
    // std::cout << "inverse FFT:" << std::endl;
    radix2_FFT<FieldT>(a, a_size, false, roots);

    const FieldT sconst = FieldT(a_size).inverse();

    for (size_t i = 0; i < a_size; ++i)
    {
        a[i] *= sconst;
    }
}

/** Performs a standard DIT, radix-2 computation **/
template<typename FieldT>
void radix2_FFT(FieldT *a, size_t a_size, bool directReverseSwitch, std::vector<FieldT>& roots)
{
    FieldT omega;
    long long int stride;
    long long int sizeLookupTable = FieldT::getCompositeDomainSize();

    if (directReverseSwitch) {
        omega = FieldT(FieldT::getOmega(a_size));
        stride = (FieldT::getCompositeDomainSize()/a_size);
    } else {
        omega = FieldT(FieldT::getOmega(a_size)).inverse();
        stride = -(FieldT::getCompositeDomainSize()/a_size);
    }

    const size_t n = a_size, logn = log2(n);
    if (n != (1u << logn)) throw internalError(std::string("FFT error: n should be a power of 2: n = ") + std::to_string(n));

    /* swapping in place (from Storer's book) */
    for (size_t k = 0; k < n; ++k)
    {
        const size_t rk = bitreverse(k, logn);
        if (k < rk)
            std::swap(a[k], a[rk]);
    }

    size_t m = 1; // invariant: m = 2^{s-1}
    size_t size = 2ull << FieldT::n;
    for (size_t s = 1; s <= logn; ++s)
    {
        asm volatile  ("/* pre-inner */");
        for (size_t k = 0; k < n; k += 2*m)
        {
            long long int lookupIdx = 0;
            for (size_t j = 0; j < m; ++j)
            {
                FieldT t = roots[lookupIdx] * a[k+j+m];
                a[k+j+m] = a[k+j] - t;
                a[k+j] = a[k+j] + t;

                lookupIdx = (lookupIdx + stride * (n/(2*m)) + sizeLookupTable) % sizeLookupTable;
            }
        }
        asm volatile ("/* post-inner */");

        m *= 2;
    }
}
//        stride = (FieldT::getCompositeDomainSize()/a_size);
//    long long int size = FieldT::getCompositeDomainSize();

/** Performs a standard DIT, radix-2 computation **/
void radix2_FFT_64(uint64_t *a, size_t size, long long int stride, bool directReverseSwitch, size_t prime)
{
    uint64_t omega;
    WideInteger wprime = static_cast<WideInteger>(prime);

    auto mod = [&](RegularInteger&& b) -> RegularInteger {
        return (((b % prime) + prime) % prime);
    };

    auto modW = [&](WideInteger&& b) -> RegularInteger {
        return static_cast<RegularInteger>(((b % wprime) + wprime) % wprime);
    };

    if (directReverseSwitch) omega = a[size];
     else omega = a[2*size-1];

    const size_t n = size, logn = log2(n);
    if (n != (1u << logn)) throw internalError(std::string("FFT error: n should be a power of 2: n = ") + std::to_string(n));

    /* swapping in place (from Storer's book) */
    for (size_t k = 0; k < n; ++k)
    {
        const size_t rk = bitreverse(k, logn);
        if (k < rk)
            std::swap(a[k], a[rk]);
    }

    size_t m = 1; // invariant: m = 2^{s-1}
    for (size_t s = 1; s <= logn; ++s)
    {
        asm volatile  ("/* pre-inner */");
        for (size_t k = 0; k < n; k += 2*m)
        {
            long long int lookupIdx = 0;
            for (size_t j = 0; j < m; ++j)
            {
                uint64_t t = a[size+lookupIdx]*a[k+j+m];
                a[k+j+m] = a[k+j] - t;
                a[k+j] = a[k+j] + t;

                lookupIdx = (lookupIdx + stride * (n/(2*m)) + size) % size;
            }
        }
        asm volatile ("/* post-inner */");

        m *= 2;
    }
}

/** Performs a standard DIT, radix-2 computation **/
template<typename FieldT>
void radix2_FFT_bench(uint64_t *a, size_t a_size, bool directReverseSwitch, std::vector<uint64_t>& roots)
{
    FieldT omega;
    long long int stride;
    long long int sizeLookupTable = FieldT::getCompositeDomainSize();

    if (directReverseSwitch) {
        omega = FieldT(FieldT::getOmega(a_size));
        stride = (FieldT::getCompositeDomainSize()/a_size);
    } else {
        omega = FieldT(FieldT::getOmega(a_size)).inverse();
        stride = -(FieldT::getCompositeDomainSize()/a_size);
    }

    const size_t n = a_size, logn = log2(n);
    if (n != (1u << logn)) throw internalError(std::string("FFT error: n should be a power of 2: n = ") + std::to_string(n));

    /* swapping in place (from Storer's book) */
    for (size_t k = 0; k < n; ++k)
    {
        const size_t rk = bitreverse(k, logn);
        if (k < rk)
            std::swap(a[k], a[rk]);
    }

    size_t m = 1; // invariant: m = 2^{s-1}
    size_t size = 2ull << FieldT::n;
    for (size_t s = 1; s <= logn; ++s)
    {
        asm volatile  ("/* pre-inner */");
        for (size_t k = 0; k < n; k += 2*m)
        {
            long long int lookupIdx = 0;
            for (size_t j = 0; j < m; ++j)
            {
                uint64_t t = roots[lookupIdx] * a[k+j+m];
                a[k+j+m] = a[k+j] - t;
                a[k+j] = a[k+j] + t;

                lookupIdx = (lookupIdx + stride * (n/(2*m)) + sizeLookupTable) % sizeLookupTable;
            }
        }
        asm volatile ("/* post-inner */");

        m *= 2;
    }
}

