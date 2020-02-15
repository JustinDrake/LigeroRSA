#include "gtest/gtest.h"

#include "Common.hpp"
#include "FiniteFields.hpp"
#include "SecretSharing.hpp"

using namespace ligero;

typedef MPInt NumberType2;
const NumberType2 p2(4611686018326724609ULL);
typedef Fpp<NumberType2, p2> FieldT2;

typedef RegularInteger NumberType1;
const NumberType1 p1(4611686018326724609ULL);
typedef Fpp_Fixed<p1> FieldT1;

TEST(TestFFT, FFT_MPInt_64bits_Poly_Order_1) {

    // Initialize Fpp
    FieldT1::preprocessing();
    FieldT2::preprocessing();

    // Size of the Domain
    size_t dim1 = FieldT1::getCompositeDomainSize();
    size_t dim2 = FieldT2::getCompositeDomainSize();

    // Target Evaluation Domain
    size_t evaluationDim = 1ull << 15;

    // Build Random Vector of Polynomial Coefficients
    FieldT1 coefs1[evaluationDim]; 
    FieldT2 coefs2[evaluationDim]; 

    for (size_t i = 0 ; i < evaluationDim; i++) {
        coefs1[i].populate(NumberType1(0));
        coefs2[i].populate(NumberType2(0));
    }

    // Generating Sample Polynomials
    DBG("Generating 1 polynomial of degree " << 1 << ".");
    coefs1[0] = NumberType1(0);
    coefs1[1] = NumberType1(1);

    DBG("Generating 1 polynomial of degree " << 1 << ".");
    coefs2[0] = NumberType2(0);
    coefs2[1] = NumberType2(1);

    FieldT1 evals1[evaluationDim]; 
    FieldT2 evals2[evaluationDim]; 

    DBG("Minus1 : " << FieldT1::getModulus()-NumberType1(1));

   naiveEvaluation(evals1, coefs1, evaluationDim, std::string("64"));
   FFT<FieldT1>(coefs1, evaluationDim, FieldT1(1));
    
   naiveEvaluation(evals2, coefs2, evaluationDim, std::string("MPInt"));
   FFT<FieldT2>(coefs2, evaluationDim, FieldT2(1));

    int test = 0;

    for (size_t idx = 0 ; idx < evaluationDim ; idx++) {
        if (!(
            (MPInt(coefs1[idx].getValue()) == MPInt(evals1[idx].getValue())) 
            &&(MPInt(coefs2[idx].getValue() == MPInt(evals2[idx].getValue()))) 
            &&(MPInt(coefs1[idx].getValue() == MPInt(coefs2[idx].getValue()))))
        )
            {test = -1;DBG("error at idx : " << idx);}
    }

    ASSERT_EQ(test,0);
}

TEST(TestFFT, iFFT_MPInt_64bits_Poly_Order_1) {
    // Initialize Fpp
    FieldT1::preprocessing();
    FieldT2::preprocessing();

    // Size of the Domain
    size_t dim1 = FieldT1::getCompositeDomainSize();
    size_t dim2 = FieldT2::getCompositeDomainSize();

    // Target Evaluation Domain
    size_t evaluationDim = 1ull << 15;

    // Build Random Vector of Polynomial Coefficients
    FieldT1 evals1[evaluationDim]; 
    FieldT2 evals2[evaluationDim]; 

    for (size_t i = 0 ; i < evaluationDim; i++) {
        evals1[i].populate(NumberType1(0));
        evals2[i].populate(NumberType2(0));
    }

    // Generating Sample Polynomials
    DBG("Generating 1 eval of dimension " << 2 << ".");
    evals1[0] = NumberType1(0);
    evals1[1] = NumberType1(1);

    DBG("Generating 1 eval of dimension " << 2 << ".");
    evals2[0] = NumberType2(0);
    evals2[1] = NumberType2(1);

    iFFT<FieldT1>(evals1, evaluationDim, FieldT1(1));
    iFFT<FieldT2>(evals2, evaluationDim, FieldT2(1));

    int test = 0;

    for (size_t idx = 0 ; idx < evaluationDim ; idx++) {
        if (!(
            (MPInt(evals1[idx].getValue()) == MPInt(evals2[idx].getValue())) 
        ))
            {test = -1;DBG("error at idx : " << idx);}
    }

    ASSERT_EQ(test,0);
}

TEST(TestFFT, FFT_MPInt_64bits_Poly_Order_Dom_Size) {
    // Initialize Fpp
    FieldT1::preprocessing();
    FieldT2::preprocessing();

    // Size of the Domain
    size_t dim1 = FieldT1::getCompositeDomainSize();
    size_t dim2 = FieldT2::getCompositeDomainSize();

    // Target Evaluation Domain
    size_t evaluationDim = 1ull << 15;

    // Build Random Vector of Polynomial Coefficients
    FieldT1 coefs1[evaluationDim]; 
    FieldT2 coefs2[evaluationDim]; 

    for (size_t i = 0 ; i < evaluationDim; i++) {
        coefs1[i].populate(NumberType1(0));
        coefs2[i].populate(NumberType2(0));
    }

    // Generating Sample Polynomials
    DBG("Generating 1 polynomial of degree " << evaluationDim << ".");
    FieldT1::randomPartialMatrix(coefs1, 1, evaluationDim, evaluationDim, true);

    DBG("Generating 1 polynomial of degree " << evaluationDim << ".");
    for (size_t i = 0 ; i < evaluationDim; i++) {
        coefs2[i] = FieldT2(coefs1[i].getValue());
    }    

    FieldT1 evals1[evaluationDim]; 
    FieldT2 evals2[evaluationDim]; 

   naiveEvaluation(evals1, coefs1, evaluationDim);
   FFT<FieldT1>(coefs1, evaluationDim, FieldT1(1));

   naiveEvaluation(evals2, coefs2, evaluationDim);
   FFT<FieldT2>(coefs2, evaluationDim, FieldT2(1));

    int test = 0;

    for (size_t idx = 0 ; idx < evaluationDim ; idx++) {
        if (!(
            (MPInt(coefs1[idx].getValue()) == MPInt(evals1[idx].getValue())) 
            &&(MPInt(coefs2[idx].getValue() == MPInt(evals2[idx].getValue()))) 
            &&(MPInt(coefs1[idx].getValue() == MPInt(coefs2[idx].getValue()))))
        )
            {test = -1;DBG("error at idx : " << idx);}
    }

    ASSERT_EQ(test,0);
}

int main(int argc, char** argv)
{

  // Configure easylogging
  el::Configurations c;
  c.setToDefault();
  c.parseFromText("*GLOBAL:\n FORMAT = %datetime{%H:%m:%s.%g} %msg\n FILENAME = \"test_fft_01.log\" \n TO_FILE = true \n TO_STANDARD_OUTPUT = false\n*DEBUG:\n \n TO_FILE = true\n TO_STANDARD_OUTPUT = false\n*INFO:\n \n TO_FILE = false\n TO_STANDARD_OUTPUT = false");
  el::Loggers::setDefaultConfigurations(c, true);

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
