#include <iostream>
#include "gtest/gtest.h"
#include "Ligero.hpp"



TEST(CRTTest, TestBucketing) {
    size_t degree = 1 << 15;
    auto [alphas, bucketSize] = ligero::math::balanced_bucket_n_primes(1000, degree);
    ASSERT_EQ(alphas.size(), 8);
    ASSERT_EQ(bucketSize.size(), 8);

    size_t sum = 0;
    for (size_t i = 0; i < bucketSize.size(); ++i) {
        sum += bucketSize[i];
    }
    ASSERT_EQ(sum, degree);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
