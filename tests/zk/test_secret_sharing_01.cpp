#include "Common.hpp"
#include "FiniteFields.hpp"
#include "SecretSharing.hpp"
#include "gtest/gtest.h"

#define OPTIMIZED

using namespace ligero;

// Define the precision
#ifdef OPTIMIZED
typedef RegularInteger NumberType;
// const NumberType p(97ULL);
const NumberType p(10514644122014433281ULL);
//const NumberType p(10514644122014433281ULL);
typedef Fpp_Fixed<p> FieldT;
#else
typedef MPInt NumberType;   
const NumberType p(10514644122014433281ULL);
typedef Fpp<NumberType, p> FieldT;
#endif

TEST(TestSecretSharing, integrityTest) {

    try {
        // Parameters
        size_t l = 1 << 11;               // Secret Size
        size_t k = 1 << 12;               // Degree
        size_t n = 1 << 13;               // Number of shares
        size_t c = 1;                     // Number of rows

        FieldT* secret = new FieldT[n*c]; 
        FieldT* copy = new FieldT[n*c]; 

        // Initialize Fpp
        FieldT::preprocessing(log2(n)+1);

        // Initialize Secret Sharing Interface
        ligero::SecretSharing<FieldT> localSecretSharing((size_t)l, (size_t)k, (size_t)n);

        // for (size_t i = 0 ; i < n*c; i++) { 
        //     secret[i].populate(NumberType(1));
        // }

        DBG("generating c random secrets of size " << l << ":");
        FieldT::randomPartialMatrix(secret, c, n, l, true);

        // Share the secrets
        // for (size_t i = 0 ; i < c; i++)
        //     print("Original Secret", &secret[i*n], n);

        // Copy the secret
        for (size_t i = 0 ; i < n*c; i++) {
            copy[i] = secret[i];
        }

        DBG("done with generation");

        localSecretSharing.shareMany(secret, c);

        for (size_t i = 0 ; i < 10; i++)
            DBG("Shared Secret:" << copy[i].getValue());

        // Reconstruct the secret
        localSecretSharing.reconstructMany(secret, c);

        for (size_t i = 0 ; i < 10; i++ )
            DBG("Reconstructed Secret:" << secret[i].getValue());

        DBG("done reconstructing");

        int success = 1;

        for (size_t i = 0 ; i < l; i++) {
            for (size_t j = 0 ; j < 1; j++) {
                if (copy[i+j*n] != secret[i+j*n]) {
                    success = 0;
                    // std::cout << "failed at " << i << "," << j << std::endl;
                    break;
                }
            }
        }

        ASSERT_EQ(success,1);

    } catch(std::exception& e) {
        DBG("Aborting: " << e.what());
    }

}

int main(int argc, char** argv)
{

  // Configure easylogging
  el::Configurations c;
  c.setToDefault();
  c.parseFromText("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s.%g} %msg\n FILENAME = \"test_secret_sharing_01.log\" \n TO_FILE = true \n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n*INFO:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = false");
  el::Loggers::setDefaultConfigurations(c, true);

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}