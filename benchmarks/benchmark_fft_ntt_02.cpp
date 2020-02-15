#include <thread>
#include <future>

#include <iostream>
#include <tuple>
#include <array>

#include <gmpxx.h>
#include <nfl.hpp>
#include "nfl/ops.hpp"

#include <type_traits>
#include <vector>
#include <numeric>
#include <algorithm>

#include "Common.hpp"
#include "FiniteFields.hpp"
#include "FastFourierTransform.hpp"
#include "NumberTheoreticTransform.hpp"
using namespace ligero;

constexpr auto degree = (1ull << 11);
constexpr auto q = 21;

using namespace nfl;

//==============================================================================
int main(int argc, char** argv)
{
    timers myTimers;
    myTimers.initialize_timer();

    // nfl::poly_p<uint64_t, degree, 1> referencePoly(nfl::uniform{});

    // myTimers.begin(1,"NTT_64_full",0);
    std::ofstream roots("roots");

    // Reconciling vs naive evaluation

    // Build Evaluations
    std::vector<uint64_t> dataDegOne(degree,0);
    std::vector<uint64_t> dataRandom(degree,0);
    std::vector<uint64_t> dataLarger(2*degree,0);
    std::vector<uint64_t> flatOmegaTable(2*degree,0);

    dataDegOne[1]=1;
    dataRandom[2]=564;
    dataRandom[10]=4632;
    dataRandom[26]=35;
    dataRandom[175]=654;

    std::vector<uint64_t> dataPoly(dataRandom);

    NTT<degree> degOne(0, &dataDegOne[0]);
    NTT<degree> random(0, &dataRandom[0]);
    NTT<2*degree> larger(0, &dataLarger[0]);

    // Evaluating f(x) = x
    degOne.ntt();

    // Build a flat table with roots of unity
    // Build the map to the roots
    uint64_t root = 1;
    for (size_t i = 0; i < degree; i++) {
      flatOmegaTable.emplace_back(root);
      root = mulmod{}(root, larger.omegas[1],0);
    }

    // Compute the values of random

    // Compare with NTT
    random.ntt();

    size_t match = 0;
    for (size_t i = 0; i < degree; i++) {
        uint64_t eval = 0;
        for (size_t m = 0; m < degree; m++) {
          uint64_t term = 1;
          if (dataPoly[m] > 0) {
            for (size_t k = 0; k < m; k++) {
              term = mulmod{}(term, degOne.data[i],0);
            }
            eval = static_cast<uint64_t>((static_cast<__uint128_t>(eval) + static_cast<__uint128_t>(term) * static_cast<__uint128_t>(dataPoly[m])) % static_cast<__uint128_t>(params<uint64_t>::P[0]));
          }
        }

        assert(eval == random.data[i]);
        if (eval == random.data[i]) match++;
    }

    std::cout << "evaluation test successful: # matches = " << match << std::endl;

  //   // Checking whether the roots are distinct
  //  for (size_t i = 10; i < 11; i++) {
      constexpr size_t i = 14;

      constexpr size_t smallDegree = 1ull << i;
      constexpr size_t largeDegree = 1ull << (i+1);

      std::vector<uint64_t> dataSmall(smallDegree,0);
      std::vector<uint64_t> dataLarge(largeDegree,0);

      dataSmall[1]=1;
      dataLarge[1]=1;

      NTT<smallDegree> polySmall(0, &dataSmall[0]);
      NTT<largeDegree> polyLarge(0, &dataLarge[0]);

      polySmall.ntt();
      polyLarge.ntt();

      for (size_t k = 0; k < smallDegree; k++) {
        for (size_t j = 0; j < largeDegree; j++) {
          if (polySmall.data[k] == polyLarge.data[j]) {
            roots.close();
            std::cout << "domain test failed: overlapping evaluation domains" << std::endl;
            return 0;
          }
        }
      }

      std::cout << "domain test successful: domains of size " << smallDegree << " and " << largeDegree << " are disjoint." << std::endl;

    // }

    roots.close();
    return 0;
}


