#include <thread>
#include <future>

#include <iostream>
#include <tuple>
#include <array>

#include <gmpxx.h>
#include <nfl.hpp>

#include "Common.hpp"
#include "FiniteFields.hpp"
#include "FastFourierTransform.hpp"
using namespace ligero;

typedef RegularInteger NumberType;
NumberType p1(4611686018326724609ULL);
typedef Fpp_Fixed<p1> FieldT;

typedef unsigned int uint128_t __attribute__((mode(TI)));

__attribute__((noinline)) void ntt_new(uint64_t* x, uint64_t* wtab, uint64_t* winvtab,
    unsigned k, uint64_t p)
{
  if (k == 1)
    return;

  // special case
  if (k == 2)
  {
    uint64_t u0 = x[0];
    uint64_t u1 = x[1];
    uint64_t t0 = u0 + u1;
    uint64_t t1 = u0 - u1;
    t0 -= (t0 >= 2*p) ? (2*p) : 0;
    t1 += ((int64_t) t1 < 0) ? (2*p) : 0;
    x[0] = t0;
    x[1] = t1;
    return;
  }

  size_t N = k;                       // size of block
  size_t M = 1;                       // number of blocks

  for (; N > 4; N /= 2, M *= 2)
  {
    uint64_t* x0 = x;
    uint64_t* x1 = x + N/2;
    for (size_t r = 0; r < M; r++, x0 += N, x1 += N)
    {
      ptrdiff_t i = N/2 - 2;
      do
      {
        {
          uint64_t u0 = x0[i+1];
          uint64_t u1 = x1[i+1];

          uint64_t t0 = u0 + u1;
          t0 -= ((t0 >= 2*p) ? (2*p) : 0);

          uint64_t t1 = u0 - u1 + 2*p;

          uint64_t q = ((uint128_t) t1 * winvtab[i+1]) >> 64;
          uint64_t t2 = t1 * wtab[i+1] - q * p;

          x0[i+1] = t0;
          x1[i+1] = t2;
        }
        {
          uint64_t u0 = x0[i];
          uint64_t u1 = x1[i];

          uint64_t t0 = u0 + u1;
          t0 -= ((t0 >= 2*p) ? (2*p) : 0);

          uint64_t t1 = u0 - u1 + 2*p;

          uint64_t q = ((uint128_t) t1 * winvtab[i]) >> 64;
          uint64_t t2 = t1 * wtab[i] - q * p;

          x0[i] = t0;
          x1[i] = t2;
        }

        i -= 2;
      }
      while (i >= 0);
    }
    wtab += N/2;
    winvtab += N/2;
  }

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
    v2 += ((int64_t) v2 < 0) ? (2*p) : 0;

    uint64_t v1 = u1 + u3;
    v1 -= (v1 >= 2*p) ? (2*p) : 0;
    uint64_t t = u1 - u3 + 2*p;

    uint64_t q = ((uint128_t) t * winvtab[1]) >> 64;
    uint64_t v3 = t * wtab[1] - q * p;

    uint64_t z0 = v0 + v1;
    z0 -= (z0 >= 2*p) ? (2*p) : 0;
    uint64_t z1 = v0 - v1;
    z1 += ((int64_t) z1 < 0) ? (2*p) : 0;

    uint64_t z2 = v2 + v3;
    z2 -= (z2 >= 2*p) ? (2*p) : 0;
    uint64_t z3 = v2 - v3;
    z3 += ((int64_t) z3 < 0) ? (2*p) : 0;

    x[0] = z0;
    x[1] = z1;
    x[2] = z2;
    x[3] = z3;
  }

}

namespace nfl { namespace tests {

template <class P>
class poly_tests_proxy
{
  using value_type = typename P::value_type;

  public:
  static inline bool ntt(value_type* x, const value_type* wtab, const value_type* winvtab, value_type const p)
  {
    return P::core::ntt(x, wtab, winvtab, p);
  }

  static inline value_type* get_omegas(P& p) { return &p.base.omegas[0][0]; }
  static inline value_type* get_shoupomegas(P& p) { return &p.base.shoupomegas[0][0]; }
};

} // tests

} // nfl

//==============================================================================
int main(int argc, char** argv)
{

    // Configure easylogging
    el::Configurations c;
    c.setToDefault();
    std::string config = std::string("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s.%g} %msg\n LOG_FLUSH_THRESHOLD = 1\n ENABLED = true\n FILENAME = \"test") + std::string(".log\" \n TO_FILE = true \n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n LOG_FLUSH_THRESHOLD = 1\n *INFO:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false");
    c.parseFromText(config.c_str());
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
    el::Loggers::setDefaultConfigurations(c, true);
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);

    timers myTimers;
    myTimers.initialize_timer();

    myTimers.begin(1,"pop1",0);

    nfl::poly_p<uint64_t, 2048, 1> referencePoly(nfl::uniform{});

    for (size_t i=0; i<2048; i++) {
      referencePoly.poly_obj().data()[i] = 0;
    }

    referencePoly.poly_obj().data()[1] = 1;
    referencePoly.poly_obj().ntt_pow_phi();
  
    for (size_t i=0; i<10; i++) {
      std::cout << referencePoly.poly_obj().data()[i] << std::endl;
    }

    myTimers.end(1,"pop1",0);

    myTimers.begin(1,"NTT",0);

    referencePoly.ntt_pow_phi();

    myTimers.end(1,"NTT",0);

    myTimers.begin(1,"pop2",0);

    p1=nfl::params<uint64_t>::P[0];

    FieldT::preprocessing(16);

    FieldT data[65536], copy[65536];
    std::vector<FieldT> roots(65536, FieldT(0));

    // Regular FieldT computation
    uint64_t datafast[2*65536];
    std::vector<uint64_t> rootsfast(65536, uint64_t(0));

    for (size_t i = 0; i< 65536; i++) {data[i] = referencePoly.poly_obj().data()[i];}
    for (size_t i = 0; i< 65536; i++) {copy[i] = referencePoly.poly_obj().data()[i];}
    for (size_t j = 0; j< 65536; j++) {roots[j]= FieldT(FieldT::rootsOfUnity[j]);}

    myTimers.end(1,"pop2",0);

    myTimers.begin(1,"FFT",0);

	FFT<FieldT>(&data[0], 65536, roots);

    myTimers.end(1,"FFT",0);

    // 64-bit idealized computation

    for (size_t i = 0; i< 65536; i++) {datafast[i] = referencePoly.poly_obj().data()[i];}
    for (size_t j = 0; j< 65536; j++) {rootsfast[j]= uint64_t(FieldT::rootsOfUnity[j]);}

    myTimers.begin(1,"FFT64",0);

    radix2_FFT_bench<FieldT>(&datafast[0], 65536, true, rootsfast);

    myTimers.end(1,"FFT64",0);

    // Compact Data Structure
    for (size_t i = 0; i< 65536; i++) {datafast[i] = referencePoly.poly_obj().data()[i];}
    for (size_t j = 0; j< 65536; j++) {datafast[65536 + j]= uint64_t(FieldT::rootsOfUnity[j]);}

    myTimers.begin(1,"FFT64_CC",0);

    radix2_FFT_64(&datafast[0], 65536, 1, true, 4611686018326724609ULL);

    myTimers.end(1,"FFT64_CC",0);

    myTimers.begin(1,"NTT1",0);

    using poly_t = nfl::poly_from_modulus<uint64_t, 65536, 62>;
    using poly_proxy_t = nfl::tests::poly_tests_proxy<poly_t>;

    // REC FFTs

    poly_t p;

    std::cout << "prime:" << p.get_modulus(0) << std::endl;

    for (size_t i = 0; i< 65536; i++) {p(0,i) = 0;}
    p(0,1) = 1;
    ntt_new(&p(0, 0), poly_proxy_t::get_omegas(p), poly_proxy_t::get_shoupomegas(p), 65536, p.get_modulus(0));
    // poly_proxy_t::ntt(&p(0,0), poly_proxy_t::get_omegas(p), poly_proxy_t::get_shoupomegas(p), p.get_modulus(0));

    for (size_t i =0; i< 30;i++) {std::cout << "eval:" << p(0,i) << std::endl;}
    for (size_t i =0; i< 30;i++) {std::cout << "roots:" << poly_proxy_t::get_omegas(p)[i] << std::endl;}
    for (size_t i =0; i< 30;i++) {std::cout << "roots (mycalc):" << FieldT::rootsOfUnity[i] << std::endl;}

    for (size_t i =0; i< 30;i++) {
        bool found = false;
        size_t jj;
        RegularInteger *roots = &FieldT::rootsOfUnity[0];
        for (size_t j =0; j< 65536;j++) {if ((uint64_t)FieldT::rootsOfUnity[j] == (uint64_t)p(0,i)) {found = true;jj =j; break;}} 
        assert(found);
        if (found) std::cout << "eval_idx:" << i << " = " << jj << std::endl;
        else std::cout << p(0,i) << " not found!" << std::endl;
        }
        std::cout << "root1= " << FieldT::rootsOfUnity[32768] << std::endl;
        std::cout << "root2= " << p(0,2) << std::endl;

    for (size_t i = 0; i< 65536; i++) {std::cout << i << std::endl; assert(p(0,i) == data[i].getValue());}

    std::cout << "prime:" << p.get_modulus(0) << std::endl;
    std::cout << "omega:" << poly_proxy_t::get_omegas(p)[10001] << std::endl;
    for (size_t idx =0; idx < 65536; idx++) {if (FieldT::rootsOfUnity[idx] == poly_proxy_t::get_shoupomegas(p)[0]) std::cout << "yep:" << idx << std::endl;}
    std::cout << "omega:" << FieldT::rootsOfUnity[10001] << std::endl;
    myTimers.end(1,"NTT1",0);

}
