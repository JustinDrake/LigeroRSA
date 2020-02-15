#include "Common.hpp"
#include <stdio.h>
#include <stdlib.h>
#include "FiniteFields.hpp"
#include "SecretSharingNTT.hpp"
#include "gtest/gtest.h"

#define OPTIMIZED

using namespace ligero;

constexpr size_t l = 1ull << 12;               // Secret Size
constexpr size_t k = 1ull << 13;               // Degree
constexpr size_t n = 1ull << 16;               // Number of shares
constexpr size_t c = 3000;                     // Number of rows

TEST(TestSecretSharing, integrityTest) {

    try {
        uint64_t *secret = static_cast<uint64_t*>(calloc(n*c, sizeof(uint64_t)));
        uint64_t *copy = static_cast<uint64_t*>(calloc(n*c, sizeof(uint64_t))); 

        // Initialize Secret Sharing Interface
        ligero::SecretSharingNTT localSecretSharing((size_t)l, (size_t)k, (size_t)n, permut<k>::compute, permut<n>::compute);

        DBG("generating c random secrets of size " << l << ":");
        secret[3] = 585457478;

        // Copy the secret
        for (size_t i = 0 ; i < n*c; i++) {
            copy[i] = secret[i];
        }

        DBG("done with generation");

        localSecretSharing.shareMany(secret, c);

        for (size_t i = 0 ; i < 10; i++)
            DBG("Shared Secret:" << copy[i]);

        // Reconstruct the secret
        localSecretSharing.reconstructMany(secret, c);

        for (size_t i = 0 ; i < 10; i++ )
            DBG("Reconstructed Secret:" << secret[i]);

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