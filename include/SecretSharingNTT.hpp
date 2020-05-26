#pragma once

#include <math.h>
#include <time.h>
#include <boost/multiprecision/gmp.hpp>

#include "FiniteFields.hpp"
#include "NumberTheoreticTransform.hpp"

namespace ligero {

/**
 *  Secret Sharing Interface
 *   Handles the sharing and reconstruction of the private witness, as well as
 * support for multiple tests for inter-block, intra-block and quadratic
 * constraints
 */
class SecretSharingNTT {
 public:
  /* Constructors allowing some flexibility in the choice of domains */
  ~SecretSharingNTT() {
    delete _smallDomain;
    delete _largeDomain;
  }
  SecretSharingNTT(size_t modulusIdx, size_t lSecretLength, size_t kDegree,
                   size_t nShares,
                   void (*computeSmall)(uint64_t *, uint64_t const *),
                   void (*computeLarge)(uint64_t *, uint64_t const *))
      : modulusIdx_(modulusIdx),
        lSecretLength_(lSecretLength),
        kDegree_(kDegree),
        nNumberShares_(nShares) {
    this->dCompositeDomainSize_ = nShares;

    _smallDomain = new NTT(modulusIdx, this->kDegree_, computeSmall);
    _largeDomain = new NTT(modulusIdx, this->nNumberShares_, computeLarge);
  };

  void share(uint64_t *secret);
  void reconstruct(uint64_t *eval, bool expanded = false);

  void shareMany(uint64_t *secrets, size_t lrows);
  void reconstructMany(uint64_t *secrets, size_t lrows);

  bool degreeTest(uint64_t *eval);
  bool zeroTest(uint64_t *eval, bool expanded = false);
  bool zeroSumTest(uint64_t *eval);

 protected:
  void padPolynominal(uint64_t *secret);
  void padIntrablocRandomness(uint64_t *secret);

  size_t modulusIdx_;
  size_t lSecretLength_;
  size_t kDegree_;
  size_t nNumberShares_;
  size_t dCompositeDomainSize_;

  NTT *_smallDomain;
  NTT *_largeDomain;
};

template <typename FieldT>
class SecretSharingInterface : SecretSharingNTT {
 public:
  /* Constructors allowing some flexibility in the choice of domains */
  SecretSharingInterface(size_t modulusIdx, size_t lSecretLength,
                         size_t kDegree, size_t nShares,
                         void (*computeSmall)(uint64_t *, uint64_t const *),
                         void (*computeLarge)(uint64_t *, uint64_t const *))
      : SecretSharingNTT(modulusIdx, lSecretLength, kDegree, nShares,
                         computeSmall, computeLarge) {
    this->dCompositeDomainSize_ = nShares;
  };

  using SecretSharingNTT::degreeTest;
  using SecretSharingNTT::reconstruct;
  using SecretSharingNTT::share;
  using SecretSharingNTT::zeroSumTest;
  using SecretSharingNTT::zeroTest;

  void share(FieldT *secret);
  void reconstruct(FieldT *eval, bool expanded = false);

  void padMany(FieldT *secrets, size_t lrows);
  void padIntrablocRandomness(FieldT *secret);
  void shareMany(FieldT *secrets, size_t lrows);
  void reconstructMany(FieldT *secrets, size_t lrows);

  bool degreeTest(FieldT *eval);
  bool zeroTest(FieldT *eval, bool expanded = false);
  bool zeroSumTest(FieldT *eval);
  FieldT sumReconstruct(FieldT *eval);

 protected:
  void padSecret(FieldT *secret);
  void padPolynominal(FieldT *secret);

  // adapter
  std::vector<uint64_t> demote(FieldT *eval, size_t length);
  void promote(uint64_t *eval, size_t length, FieldT *dest);
};

/** Implementation
/*
=============================================================================================
 */

/** Pads a block to k degree with numbers generated at random from the
 * underlying field
 * @param secret witness block to be padded
 */
template <typename FieldT>
inline void SecretSharingInterface<FieldT>::padSecret(FieldT *secret) {
  FieldT::randomVector(secret + this->lSecretLength_,
                       this->kDegree_ - this->lSecretLength_, true);
}

/** Pads randomness for intrabloc constraints up to degree k with zeroes
 * @param secret randomness block to be padded
 */
template <typename FieldT>
inline void SecretSharingInterface<FieldT>::padIntrablocRandomness(
    FieldT *secret) {
  for (size_t i = this->lSecretLength_; i < kDegree_; i++) {
    secret[i] = FieldT(0);
  }
}

/** Pads polynomial coefficients from k degree to n number of shares
 * @param secret array of polynomial coefficients to be padded
 */
template <typename FieldT>
inline void SecretSharingInterface<FieldT>::padPolynominal(FieldT *secret) {
  for (size_t i = this->kDegree_; i < nNumberShares_; i++) {
    secret[i] = FieldT(0);
  }
}

/** Pads multiple blocks
 * @param secret witness blocks to be padded
 * @param lrows number of blocks
 */
template <typename FieldT>
void SecretSharingInterface<FieldT>::padMany(FieldT *secret, size_t lrows) {
  /* Add randomness for padding */
  for (size_t i = 0; i < lrows; i++) {
    this->padSecret(secret + i * this->nNumberShares_);
  }
}

/** Secret-sharing interface for a block in the witness: this demotes field
 * elements into 64-bit numbers and calls the specialized 64-bit secret-sharing
 * method
 * @param secret witness block to be padded
 */
template <typename FieldT>
void SecretSharingInterface<FieldT>::share(FieldT *secret) {
  std::vector<uint64_t> secret64 = demote(secret, this->nNumberShares_);
  this->share(&secret64[0]);
  promote(&secret64[0], this->nNumberShares_, secret);
}

/** Reconstructs the secret corresponding to a given shared block
 * @param eval array of field elements representing a shared block
 * @param expandedDegree flag indicating that the degree of the polynomial for
 * this bloc is more than k
 */
template <typename FieldT>
void SecretSharingInterface<FieldT>::reconstruct(FieldT *eval,
                                                 bool expandedDegree) {
  std::vector<uint64_t> eval64 = demote(eval, this->nNumberShares_);
  this->reconstruct(&eval64[0], expandedDegree);
  promote(&eval64[0], this->kDegree_, eval);
}

/** Demotes an array of field elements to 64-bit numbers
 * @param eval array of field elements
 * @param length length of the array to be demoted
 * @return array of demoted 64-bit numbers
 */
template <typename FieldT>
std::vector<uint64_t> SecretSharingInterface<FieldT>::demote(FieldT *eval,
                                                             size_t length) {
  std::vector<uint64_t> temp(length);
  for (size_t idx = 0; idx < length; idx++) {
    temp[idx] = static_cast<uint64_t>(eval[idx].getValue());
  }

  return temp;
}

/** promotes an array of 64-bit numbers to field elements
 * @param eval array of 64-bit numbers
 * @param length length of the array to be promoted
 * @param dest location of the promoted content
 */
template <typename FieldT>
void SecretSharingInterface<FieldT>::promote(uint64_t *eval, size_t length,
                                             FieldT *dest) {
  for (size_t idx = 0; idx < length; idx++) {
    dest[idx] = FieldT(eval[idx]);
  }
}

/** performs low degree testing on the shared block provided
 * @param eval pointer to the shared block of field elements
 * @return boolean indicating whether the test succeeded or failed
 */
template <typename FieldT>
bool SecretSharingInterface<FieldT>::degreeTest(FieldT *eval) {
  std::vector<uint64_t> eval64 = demote(eval, this->nNumberShares_);
  return this->degreeTest(&eval64[0]);
}

/** tests whether the polynomial corresponding to the shared block provided
 * evaluates to 0 on each point of the secret domain
 * @param eval pointer to the shared block of field elements
 * @param largerDegree boolean indicates whether the underlying polynomial is of
 * degree > k
 * @return boolean indicating whether the test succeeded or failed
 */
template <typename FieldT>
bool SecretSharingInterface<FieldT>::zeroTest(FieldT *eval, bool largerDegree) {
  std::vector<FieldT> localCopy(eval, eval + this->nNumberShares_);
  reconstruct(&localCopy[0], largerDegree);

  for (size_t i = 0; i < this->lSecretLength_; i++) {
    if (!(localCopy[i] == FieldT(0))) return false;
  }

  return true;
}

/** tests whether the sum of the evaluations for the polynomial corresponding to
 * the shared block over all points of the secret domain is 0
 * @param eval pointer to the shared block of field elements
 * @return boolean indicating whether the test succeeded or failed
 */
template <typename FieldT>
bool SecretSharingInterface<FieldT>::zeroSumTest(FieldT *eval) {
  std::vector<FieldT> localCopy(eval, eval + this->nNumberShares_);
  reconstruct(&localCopy[0], true);

  FieldT sum = FieldT(0);
  for (size_t i = 0; i < this->lSecretLength_; i++) {
    sum += localCopy[i];
  }

  if (sum == FieldT(0))
    return true;
  else
    return false;
}

/** computes and returns the sum of the evaluations for the polynomial
 * corresponding to the shared block over all points of the secret domain
 * @param eval pointer to the shared block of field elements
 * @return sum of evaluations over the secret domain
 */
template <typename FieldT>
FieldT SecretSharingInterface<FieldT>::sumReconstruct(FieldT *eval) {
  std::vector<FieldT> localCopy(eval, eval + this->nNumberShares_);
  reconstruct(&localCopy[0], true);

  FieldT sum = FieldT(0);
  for (size_t i = 0; i < this->lSecretLength_; i++) {
    sum += localCopy[i];
  }

  return sum;
}

/** performs secret-sharing for multiple blocks at once. The sharing is
 * performed in-place, i.e. the shared block will be replacing the secret in
 * memory
 * @param secret pointer to a set of private blocks of field elements
 * @param lrows total number of blocks to share
 */
template <typename FieldT>
void SecretSharingInterface<FieldT>::shareMany(FieldT *secret, size_t lrows) {
  // std::cout << "Performing FFTs: " << lrows << std::endl;
  for (size_t i = 0; i < lrows; i++) {
    share(secret + i * this->nNumberShares_);
  }
}

/** reconstructs multiple private blocs based on the corresponding shared
 * elements at once. The sharing is performed in-place, i.e. the secrets will be
 * replacing the shared blocks in memory
 * @param eval pointer to a set of shared blocks of field elements
 * @param lrows total number of blocks to share
 */
template <typename FieldT>
void SecretSharingInterface<FieldT>::reconstructMany(FieldT *eval,
                                                     size_t lrows) {
  for (size_t i = 0; i < lrows; i++) {
    reconstruct(eval + i * this->nNumberShares_);
  }
}

/** pads intrabloc randomness with zeroes from the length of a block to the
 * polynomial degree
 * @param randomness pointer to a private block of field elements
 */
void SecretSharingNTT::padIntrablocRandomness(uint64_t *randomness) {
  for (size_t i = this->lSecretLength_; i < kDegree_; i++) {
    randomness[i] = uint64_t(0);
  }
}

/** pads polynomial coefficients with zeroes from degree k to number of shares n
 * @param coefs pointer to the coefficients
 */
void SecretSharingNTT::padPolynominal(uint64_t *coefs) {
  for (size_t i = this->kDegree_; i < nNumberShares_; i++) {
    coefs[i] = uint64_t(0);
  }
}

/** secret-sharing of a demoted vector of 64-bit secrets. The secret sharing is
 * performed in place so secrets are replaced by the shared block in memory
 * @param data pointer to the demoted vector of secrets
 */
void SecretSharingNTT::share(uint64_t *data) {
  // Interpolate Poly from Secret
  _smallDomain->inv_ntt(data);

  this->padPolynominal(data);

  // Evaluate Publicly Shared Data from Poly
  _largeDomain->ntt(data);
}

/** reconstructing a secret block based on the input shared block
 * so secrets are replaced by the shared block in memory
 * @param eval point to the shared block
 * @param exapandedDegree boolean indicating whether the polynomial linked to
 * this block has a degree > k
 */
void SecretSharingNTT::reconstruct(uint64_t *eval, bool expandedDegree) {
  size_t degree = kDegree_;

  // Interpolate Poly from Publicly Shared Data
  _largeDomain->inv_ntt(eval);

  // if (expandedDegree) degree = kDegree_*2;
  // Because we evaluate on roots of unity, we can use coefs[i] - coefs[i+k]
  // on k coefficients instead of a 2k evaluation
  if (expandedDegree) {
    for (size_t i = 0; i < kDegree_; i++) {
      if (eval[i + kDegree_] > eval[i]) {
        eval[i] = eval[i] - eval[i + kDegree_] +
                  params<uint64_t>::P[this->modulusIdx_];
      } else {
        eval[i] = eval[i] - eval[i + kDegree_];
      }
    }
  }

  // Evaluate Poly Into Secret
  _smallDomain->ntt(eval);
}

/** performs low degree testing on the demoted shared block provided
 * @param eval pointer to the shared block of demoted 64-bit numbers
 * @return boolean indicating whether the test succeeded or failed
 */
bool SecretSharingNTT::degreeTest(uint64_t *eval) {
  // Interpolate Poly from Publicly Shared Data
  std::vector<uint64_t> localCopy(eval, eval + this->nNumberShares_);
  _largeDomain->inv_ntt(&localCopy[0]);

  for (size_t i = kDegree_; i < nNumberShares_; i++) {
    if (!(localCopy[i] == uint64_t(0))) {
      DBG("val:" << localCopy[i] << ",degree>" << i);
      return false;
    }
  }

  return true;
}

/** tests whether the polynomial corresponding to the demoted shared block
 * provided evaluates to 0 on each point of the secret domain
 * @param eval pointer to the demoted shared block of 64-bit numbers
 * @param largerDegree boolean indicates whether the underlying polynomial is of
 * degree > k
 * @return boolean indicating whether the test succeeded or failed
 */
bool SecretSharingNTT::zeroTest(uint64_t *eval, bool largerDegree) {
  std::vector<uint64_t> localCopy(eval, eval + this->nNumberShares_);
  reconstruct(&localCopy[0], largerDegree);

  for (size_t i = 0; i < this->lSecretLength_; i++) {
    if (!(localCopy[i] == uint64_t(0))) return false;
  }

  return true;
}

/** tests whether the sum of the evaluations for the polynomial corresponding to
 * the demoted shared block over all points of the secret domain is 0
 * @param eval pointer to the demoted shared block of 64-bit numbers
 * @return boolean indicating whether the test succeeded or failed
 */
bool SecretSharingNTT::zeroSumTest(uint64_t *eval) {
  std::vector<uint64_t> localCopy(eval, eval + this->nNumberShares_);
  reconstruct(&localCopy[0], true);

  uint64_t sum = uint64_t(0);
  for (size_t i = 0; i < this->lSecretLength_; i++) {
    sum += localCopy[i];
  }

  if (sum == uint64_t(0))
    return true;
  else
    return false;
}

/** performs secret-sharing for multiple demoted blocks at once. The sharing is
 * performed in-place, i.e. the shared block will be replacing the secret in
 * memory
 * @param secret pointer to a set of demoted private blocks of 64-bit numbers
 * @param lrows total number of blocks to share
 */
void SecretSharingNTT::shareMany(uint64_t *secret, size_t lrows) {
  // std::cout << "Performing FFTs: " << lrows << std::endl;
  for (size_t i = 0; i < lrows; i++) {
    share(secret + i * this->nNumberShares_);
  }
}

/** reconstructs multiple private blocs based on the corresponding demoted
 * shared elements at once. The sharing is performed in-place, i.e. the secrets
 * will be replacing the shared blocks in memory
 * @param eval pointer to a set of demoted shared blocks of 64-bit numbers
 * @param lrows total number of blocks to share
 */
void SecretSharingNTT::reconstructMany(uint64_t *secret, size_t lrows) {
  for (size_t i = 0; i < lrows; i++) {
    reconstruct(secret + i * this->nNumberShares_);
  }
}

}  // namespace ligero
