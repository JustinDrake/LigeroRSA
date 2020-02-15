#pragma once

#include <fstream>

#include "Common.hpp"
#include "Math.hpp"

#include "nfl/params.hpp"

using namespace nfl;

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

/*
  Number Theoretic Transform
  */

class NTT {
  protected:
    uint64_t
        *phis,
        *shoupphis,
        *invpoly_times_invphis,
        *shoupinvpoly_times_invphis,
        *omegas,
        *shoupomegas,
        *invomegas,
        *shoupinvomegas,
        invpolyDegree;
    
    size_t 
        _J,
        _indexModulus,
        _degree;

    void (*_compute) (uint64_t*,uint64_t const*);

    struct ntt_loop_body
    {
      public:
        uint64_t _p;

        ntt_loop_body(uint64_t const p) {
            _p = p;
        }

        inline void operator()(uint64_t* x0, uint64_t* x1, uint64_t const* winvtab, uint64_t const* wtab) const {
            uint64_t u0 = *x0;
            uint64_t u1 = *x1;

            uint64_t t0 = u0 + u1;
            t0 -= ((t0 >= 2*_p) ? (2*_p) : 0);

            uint64_t t1 = u0 - u1 + 2*_p;

            uint64_t q = ((__uint128_t) t1 * (*winvtab)) >> params<uint64_t>::kModulusRepresentationBitsize;
            uint64_t t2 = t1 * (*wtab) - q * _p;

            *x0 = t0;
            *x1 = t2;
        }
    };

    static constexpr uint64_t getModulus(size_t n) { return params<uint64_t>::P[n]; }

    size_t run(uint64_t* x, const uint64_t* &wtab, const uint64_t* &winvtab, const uint64_t p) {
        ntt_loop_body body(p);

        for (size_t w = 0; w < this->_J; w++) {
        const size_t M = 1 << w;
        const size_t N = _degree >> w;

        for (size_t r = 0; r < M; r++) {
            for (size_t i = 0; i < N/2; i += 2) {
            body(&x[N * r + i + 0], &x[N * r + i + 0 + N/2], &winvtab[i + 0], &wtab[i + 0]);
            body(&x[N * r + i + 1], &x[N * r + i + 1 + N/2], &winvtab[i + 1], &wtab[i + 1]);
            }
        }
        wtab += N / 2;
        winvtab += N / 2;
        }

        return 1<<this->_J;
    }


  bool ntt(uint64_t* x, const uint64_t* wtab, const uint64_t* winvtab, uint64_t const p)
  {
  #ifdef CHECK_STRICTMOD
    for (size_t i = 0 ; i < degree ; i++)
    {
      ASSERT_STRICTMOD(x[i] < p);
    }
  #endif

  #ifdef NTT_STRICTMOD
    uint64_t* x_orig = x;
  #endif

    if (_degree == 1)
      return true;

    // special case
    if (_degree == 2)
    {
      uint64_t u0 = x[0];
      uint64_t u1 = x[1];
      uint64_t t0 = u0 + u1;
      uint64_t t1 = u0 - u1;
      t0 -= (t0 >= 2*p) ? (2*p) : 0;
      t1 += ((typename std::make_signed<uint64_t>::type) t1 < 0) ? (2*p) : 0;
      x[0] = t0;
      x[1] = t1;
      return true;
    }
    const size_t M = run(x, wtab, winvtab, p);

    typedef typename std::make_signed<uint64_t>::type signed_uint64_t;
    // last two layers
    for (size_t r = 0; r < M; r++, x += 4)
    {
      uint64_t u0 = x[0];
      uint64_t u1 = x[1];
      uint64_t u2 = x[2];
      uint64_t u3 = x[3];

      uint64_t v0 = u0 + u2;
      v0 -= (v0 >= 2*p) ? (2*p) : 0;
      uint64_t v2 = u0 - u2;
      v2 += ((signed_uint64_t) v2 < 0) ? (2*p) : 0;

      uint64_t v1 = u1 + u3;
      v1 -= (v1 >= 2*p) ? (2*p) : 0;
      uint64_t t = u1 - u3 + 2*p;

      uint64_t q = ((__uint128_t) t * winvtab[1]) >> params<uint64_t>::kModulusRepresentationBitsize;
      uint64_t v3 = t * wtab[1] - q * p;

      uint64_t z0 = v0 + v1;
      z0 -= (z0 >= 2*p) ? (2*p) : 0;
      uint64_t z1 = v0 - v1;
      z1 += ((signed_uint64_t) z1 < 0) ? (2*p) : 0;

      uint64_t z2 = v2 + v3;
      z2 -= (z2 >= 2*p) ? (2*p) : 0;
      uint64_t z3 = v2 - v3;
      z3 += ((signed_uint64_t) z3 < 0) ? (2*p) : 0;

      x[0] = z0;
      x[1] = z1;
      x[2] = z2;
      x[3] = z3;
    }

  #ifdef NTT_STRICTMOD
    for (size_t i = 0; i < _degree; i++)
    {
      x_orig[i]-= ((x_orig[i]>=p)? p : 0);
      ASSERT_STRICTMOD(x_orig[i] < p);
    }
  #endif

  return true;
}

// Inverse NTT: replaces NTT values representation by the classic
// coefficient representation, parameters are the same as for ntt() with
// inv_wtab, inv_winvtab tables for the inverse operation and invK the
// inverse of the polynomialDegree
bool inv_ntt(uint64_t * x, const uint64_t* const inv_wtab, const uint64_t* const inv_winvtab, const uint64_t invK, uint64_t const p)
{
  alignas(32) uint64_t y[_degree+1];

  if (_degree == 1)
    return true;

  // bit-reverse
  _compute(y, x);

  ntt(y, inv_wtab, inv_winvtab, p);

  // bit-reverse again
  _compute(x, y);
  return true;
}

// This function is private but it is kept here to group all the functions
// by David Harvey
// This function just computes the powers of the root of unity and its inverse
// It is used by initialize() and results are used by ntt and inv_ntt
inline void prep_wtab(uint64_t* wtab, uint64_t* wtabshoup, uint64_t w, size_t cm)
{
  auto const p = getModulus(cm);
  unsigned K = _degree;

  while (K >= 2)
  {
    __uint128_t wi = 1;     // w^i
    for (size_t i = 0; i < K/2; i++)
    {
      *wtab++ = wi;
      *wtabshoup++ = ((__uint128_t) wi << params<uint64_t>::kModulusRepresentationBitsize) / p;
      wi = mulmod{}(static_cast<uint64_t>(wi), w, cm);
    }
    w = mulmod{}(w, w, cm);
    K /= 2;
  }

}

public:
    ~NTT() {
	    delete[] phis;
	    delete[] shoupphis;
	    delete[] invpoly_times_invphis;
	    delete[] shoupinvpoly_times_invphis;
	    delete[] omegas;
	    delete[] invomegas;
    }
    NTT(size_t indexModulus, size_t degree, void (*compute) (uint64_t*,uint64_t const*)) : _compute(compute) {

        // allocate all 
        phis                        = new uint64_t[degree],
        shoupphis                   = new uint64_t[degree],
        invpoly_times_invphis       = new uint64_t[degree],
        shoupinvpoly_times_invphis  = new uint64_t[degree],
        omegas                      = new uint64_t[degree * 2],
        invomegas                   = new uint64_t[degree * 2],
        shoupomegas;

        // Take in inputs
        _indexModulus   = indexModulus;
        _degree         = degree;
        _J              = log2(degree)-2;

        // initialize all 
        uint64_t phi, invphi, omega, invomega, temp;
        shoupomegas = omegas + degree;
        shoupinvomegas = invomegas + degree;

        // We start by computing phi
        // The roots in the array are primitve
        // X=2*params<T>::kMaxPolyDegree-th roots
        // Squared log2(X)-log2(degree) times they become
        // degree-th roots as required by the NTT
        // But first we get phi = sqrt(omega) squaring them
        // log2(X/2)-log2(degree) times
        phi = params<uint64_t>::primitive_roots[indexModulus];
        for (unsigned int i = 0 ; i < static_log2<params<uint64_t>::kMaxPolyDegree>::value - log2(degree) ; i++) {
            phi = mulmod{}(phi, phi, indexModulus);
          }

        // Now that temp = phi we initialize the array of phi**i values
        // Initialized to phi**0
        temp = 1;
        for (unsigned int i = 0 ; i < degree ; i++) {
            phis[i] = temp;
            shoupphis[i] = ((__uint128_t) temp << params<uint64_t>::kModulusRepresentationBitsize) / getModulus(indexModulus);
            // phi**(i+1)
            temp = mulmod{}(temp, phi, indexModulus);
          }

        // At the end of the loop temp = phi**degree

        // Computation of invphi
        // phi**(2*polydegree)=1 -> temp*phi**(degree-1) = phi**(-1)
        invphi = mulmod{}(temp, phis[degree-1], indexModulus);

        // Computation of the inverse of degree using the inverse of kMaxPolyDegree
        invpolyDegree = mulmod{}(
                  params<uint64_t>::invkMaxPolyDegree[indexModulus],
                  static_cast<uint64_t>(params<uint64_t>::kMaxPolyDegree/degree), 
                  indexModulus);

        // Now we can compute the table invpoly_times_invphis
        temp = invpolyDegree;
        for (unsigned int i = 0 ; i < degree ; i++) {
              invpoly_times_invphis[i] = temp;
              shoupinvpoly_times_invphis[i] = ((__uint128_t) temp << params<uint64_t>::kModulusRepresentationBitsize)
                  / getModulus(indexModulus);
              // This is invpolyDegree*invphi**(i+1)
              temp = mulmod{}(temp, invphi, indexModulus);
            }

        // For the omegas it is easy, we just use the function of David Harvey modified for our needs
        omega = mulmod{}(phi, phi, indexModulus);
        prep_wtab(omegas, shoupomegas, omega, indexModulus);

        // And again for the invomegas
        invomega = mulmod{}(invphi, invphi, indexModulus);
        prep_wtab(invomegas, shoupinvomegas, invomega, indexModulus);
    }

    inline void inv_ntt(uint64_t *data)
    {
        inv_ntt(data, invomegas, shoupinvomegas, invpolyDegree, getModulus(_indexModulus));

        for (size_t idx = 0; idx < _degree; idx++) {
            data[idx] = mulmod{}(data[idx], invpoly_times_invphis[idx], _indexModulus);
        }
    }

    inline void ntt(uint64_t *data) {
        
        for (size_t idx = 0; idx < _degree; idx++) {
            data[idx] = mulmod{}(data[idx], phis[idx], _indexModulus);
        }

        ntt(data, omegas, shoupomegas, getModulus(_indexModulus));
    }

};
