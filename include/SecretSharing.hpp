#pragma once

#include <boost/multiprecision/gmp.hpp>
#include <math.h>
#include <time.h>

#include "FiniteFields.hpp"
#include "FastFourierTransform.hpp"

namespace ligero {
template <typename FieldT>
class 		SecretSharing {
	public:
		/* Constructors allowing some flexibility in the choice of domains */
		SecretSharing(size_t lSecretLength, size_t kDegree, size_t nShares) : lSecretLength_(lSecretLength), kDegree_(kDegree), nNumberShares_(nShares) {
			// In this implementation, we use an evaluation domain that is twice the size of the number of shares, and we discard
			this->dCompositeDomainSize_ = nShares;

			// // embed roots of unity
			// rootsOfUnity_.resize(FieldT::rootsOfUnity.size());
			// for (size_t idx = 0; idx < FieldT::rootsOfUnity.size(); idx++) {
			// 	rootsOfUnity_[idx]= FieldT(FieldT::rootsOfUnity[idx]);
			// }
		};

		int share(FieldT *secret);
		int reconstruct(FieldT *eval, bool expanded = false);

		int padMany(FieldT *secrets, size_t lrows);
		int padIntrablocRandomnessMany(FieldT *secret, size_t lrows);
		int shareMany(FieldT *secrets, size_t lrows);
		int reconstructMany(FieldT *secrets, size_t lrows);

		bool degreeTest(FieldT *eval);
		bool zeroTest(FieldT *eval, bool expanded = false);
		bool zeroSumTest(FieldT *eval);

	protected:
		int padSecret(FieldT *secret);
		int padPolynominal(FieldT *secret);
		int padShares(FieldT *secret);
		int padIntrablocRandomness(FieldT *secret);

		size_t	lSecretLength_;
		size_t	kDegree_;
		size_t	nNumberShares_;
		size_t	dCompositeDomainSize_;
		std::vector<FieldT> rootsOfUnity_;
};

template <typename FieldT>
int SecretSharing<FieldT>::padSecret(FieldT *secret) {
	FieldT::randomVector(secret + this->lSecretLength_, this->kDegree_ - this->lSecretLength_, true);
	return 0;
}

template <typename FieldT>
int SecretSharing<FieldT>::padIntrablocRandomness(FieldT *secret) {
	for (size_t i = this->lSecretLength_; i < kDegree_ ;i++) {
		secret[i] = FieldT(0); 
	}

	return 0;
}

template <typename FieldT>
int SecretSharing<FieldT>::padPolynominal(FieldT *secret) {
	for (size_t i = this->kDegree_; i < nNumberShares_;i++) {
		secret[i] = FieldT(0);
	}

	return 0;
}

template <typename FieldT>
int SecretSharing<FieldT>::padShares(FieldT *secret) {
	for (size_t i = this->nNumberShares_; i < dCompositeDomainSize_ ;i++) {
		secret[i] = FieldT(0);
	}

	return 0;
}

template <typename FieldT>
int SecretSharing<FieldT>::padMany(FieldT *secret, size_t lrows) {

	// Add randomness for padding
	for (size_t i = 0; i < lrows; i++) {
		this->padSecret(secret + i * this->nNumberShares_);
	}

	return 0;
}

template <typename FieldT>
int SecretSharing<FieldT>::padIntrablocRandomnessMany(FieldT *secret, size_t lrows) {

	// Add randomness for padding
	for (size_t i = 0; i < lrows; i++) {
		this->padIntrablocRandomness(secret + i * this->nNumberShares_);
	}

	return 0;
}

template <typename FieldT>
int SecretSharing<FieldT>::share(FieldT *secret) {

	// Interpolate a polynomial of degree k for the secret
    iFFT<FieldT>(secret, kDegree_, rootsOfUnity_);

	// Evaluate the polynomial on a coset within the composite domain (size numberShares)
	this->padPolynominal(secret);
	FFT<FieldT>(secret, nNumberShares_, FieldT(FieldT::getOmega(nNumberShares_)), rootsOfUnity_);

	return 0;
}

template <typename FieldT>
inline int SecretSharing<FieldT>::reconstruct(FieldT *eval, bool expandedDegree) {
	size_t degree = kDegree_;
	iFFT<FieldT>(eval, nNumberShares_, FieldT(FieldT::getOmega(nNumberShares_)), rootsOfUnity_);

	// if (expandedDegree) degree = kDegree_*2; 
	// Because we evaluate on roots of unity, we can use coefs[i] + coefs[i+k]  
	// on k coefficients instead of a 2k evaluation
	if (expandedDegree) {
		for (size_t i = 0; i < kDegree_; i++) {eval[i] = eval[i] + eval[i + kDegree_];}
	}

    FFT<FieldT>(eval, degree, rootsOfUnity_);

	// rebuild accounting for the stride implied by the larger degree
	// size_t stride = degree/this->kDegree_;
	// for (size_t col = 0; col < this->lSecretLength_; col++) {
	// 	eval[col] = eval[col*stride];
	// }

	return 0;
}

template <typename FieldT>
bool SecretSharing<FieldT>::degreeTest(FieldT *eval) {

	std::vector<FieldT> localCopy(eval, eval + this->nNumberShares_);
	iFFT<FieldT>(&localCopy[0], nNumberShares_, FieldT::getOmega(nNumberShares_), rootsOfUnity_);

	for (size_t i = kDegree_; i < nNumberShares_; i++) {
		if (!(localCopy[i] == FieldT(0))) {
			DBG("val:"<< localCopy[i].getValue() << ",degree>" << i);
			return false;
		}
	}

	return true;
}

template <typename FieldT>
bool SecretSharing<FieldT>::zeroTest(FieldT *eval, bool largerDegree) {

	std::vector<FieldT> localCopy(eval, eval + this->nNumberShares_);
	reconstruct(&localCopy[0], largerDegree);

	for (size_t i = 0; i < this->lSecretLength_; i++) {
		if (!(localCopy[i] == FieldT(0))) return false;
	}

	return true;
}

template <typename FieldT>
bool SecretSharing<FieldT>::zeroSumTest(FieldT *eval) {
	
	std::vector<FieldT> localCopy(eval, eval+ this->nNumberShares_);
	reconstruct(&localCopy[0], true);

	FieldT sum = FieldT(0);
	for (size_t i = 0; i < this->lSecretLength_; i++) {
		sum += localCopy[i];
	}

	if (sum == FieldT(0)) return true;
	else return false;
}

template <typename FieldT>
int SecretSharing<FieldT>::shareMany(FieldT *secret, size_t lrows) {

	// std::cout << "Performing FFTs: " << lrows << std::endl; 
	for (size_t i = 0; i < lrows; i++) {
		share(secret + i * this->nNumberShares_);
	}

	return 0;
}

template <typename FieldT>
int SecretSharing<FieldT>::reconstructMany(FieldT *secret, size_t lrows) {

	for (size_t i = 0; i < lrows; i++) {
		reconstruct(secret + i * this->nNumberShares_);
	}

	return 0;
}


}
