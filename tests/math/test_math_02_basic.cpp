#include <iostream>
#include <sodium/randombytes.h>
#include "gtest/gtest.h"

#include "Ligero.hpp"


using namespace ligero;
using namespace math;

using Q = nfl::poly_p<uint64_t, 1<<16, 21>;

TEST(BasicMath, SameSeedSameNumber) {

    std::vector<mpz_class> v = generateRandomVector(mpz_class(1), 1);
    std::vector<mpz_class> reference = generateRandomVector(mpz_class(1), 1);
    ASSERT_EQ(reference[0], v[0]);
}

TEST(BasicMath, DifferentSeedDifferentNumber) {

    std::vector<mpz_class> v = generateRandomVector(mpz_class(1), 1);
    std::vector<mpz_class> reference = generateRandomVector(mpz_class(2), 1);
    ASSERT_NE(reference[0], v[0]);
}

TEST(BasicMath, NumberOfBitsShouldMatch) {

    mpz_class v = generateRandomVector(mpz_class(1), 1, 175)[0];
    unsigned int numBits = static_cast<unsigned int>(mpz_sizeinbase(v.get_mpz_t(), 2));
    DBG("size in bits = " << numBits);
    gmp_randclass grc(gmp_randinit_default);
    grc.seed(std::random_device()());
    mpz_class x = grc.get_z_bits(175);
    DBG("grc size in bits = " << mpz_sizeinbase(x.get_mpz_t(), 2));
    ASSERT_LT(numBits, 175);
}

TEST(BasicMath, SameValueSameHashMPZ) {
    ASSERT_EQ(computeHash(mpz_class(42)), computeHash(mpz_class(42)));
}

TEST(BasicMath, SameValueSameHashQ) {
    auto v = Q(nfl::uniform{});
    ASSERT_EQ(computeHash(v), computeHash(v));
}

TEST(BasicMath, DifferentValueDifferenHashMPZ) {
    mpz_class v1 = mpz_class(42);
    mpz_class v2 = mpz_class(10);
    ASSERT_NE(v1, v2);
    ASSERT_NE(computeHash(v1), computeHash(v2));
}

TEST(BasicMath, DifferentValueDifferentHashQ) {
    auto v1 = Q(nfl::uniform{});
    auto v2 = Q(nfl::uniform{});
    ASSERT_NE(v1, v2);
    ASSERT_NE(computeHash(v1), computeHash(v2));
}
int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
