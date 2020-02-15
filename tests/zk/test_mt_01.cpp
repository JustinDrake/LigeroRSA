/* Third Party */
#include <sodium/randombytes.h>
#include "gtest/gtest.h"

/* Ligero */
#include "Common.hpp"
#include "FiniteFields.hpp"
#include "SecretSharing.hpp"

using namespace ligero;

typedef MPInt NumberType2;
const NumberType2 p2(449);
typedef Fpp<NumberType2, p2> FieldT2;

typedef RegularInteger NumberType1;
const NumberType1 p1(449);
typedef Fpp_Fixed<p1> FieldT1;

TEST(TestMT, libsodiumRandomGenerator) {

    char *buf = new char[1ULL << 30ULL];

    /* With a seed */
    const unsigned char seed[32] = {0};
      randombytes_buf_deterministic((void *)buf, 1ULL << 30ULL, seed);
    // for (size_t i = 1; i< (1ULL << 30ULL); i++) {
    //   randombytes_buf_deterministic((void *)buf, 1, seed);
    // }    
    /* Without one */
    randombytes_buf((void *)buf, 1ULL << 30ULL);
}


int main(int argc, char** argv)
{

  // Configure easylogging
  el::Configurations c;
  c.setToDefault();
  c.parseFromText("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s.%g} %msg\n FILENAME = \"test_mt_01.log\" \n TO_FILE = true \n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n*INFO:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = false");
  el::Loggers::setDefaultConfigurations(c, true);

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
