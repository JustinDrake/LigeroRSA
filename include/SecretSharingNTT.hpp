#pragma once

#include <boost/multiprecision/gmp.hpp>
#include <math.h>
#include <time.h>

#include "FiniteFields.hpp"
#include "FastFourierTransform.hpp"
#include "NumberTheoreticTransform.hpp"

namespace ligero {

class SecretSharingNTT {
	public:
		/* Constructors allowing some flexibility in the choice of domains */
		~SecretSharingNTT() {
			delete _smallDomain;
			delete _largeDomain;
		}
		SecretSharingNTT(size_t modulusIdx, size_t lSecretLength, size_t kDegree, size_t nShares, void (*computeSmall) (uint64_t*,uint64_t const*), void (*computeLarge) (uint64_t*,uint64_t const*)) : modulusIdx_(modulusIdx), lSecretLength_(lSecretLength), kDegree_(kDegree), nNumberShares_(nShares) {
			this->dCompositeDomainSize_ = nShares;

		    _smallDomain = new NTT(modulusIdx, this->kDegree_, computeSmall);
		    _largeDomain = new NTT(modulusIdx, this->nNumberShares_, computeLarge);
		};

		int share(uint64_t *secret);
		int reconstruct(uint64_t *eval, bool exapnded = false);

		int padMany(uint64_t *secrets, size_t lrows);
		int shareMany(uint64_t *secrets, size_t lrows);
		int reconstructMany(uint64_t *secrets, size_t lrows);

		bool degreeTest(uint64_t *eval);
		bool zeroTest(uint64_t *eval, bool expanded = false);
		bool zeroSumTest(uint64_t *eval);
		bool zeroSumTestSigmaProtocol(uint64_t *eval, std::vector<uint64_t> beta);

	protected:
		int padSecret(uint64_t *secret);
		int padPolynominal(uint64_t *secret);
		int padShares(uint64_t *secret);
		int padIntrablocRandomness(uint64_t *secret);

		size_t  modulusIdx_;
		size_t	lSecretLength_;
		size_t	kDegree_;
		size_t	nNumberShares_;
		size_t	dCompositeDomainSize_;

		NTT 	*_smallDomain;
		NTT		*_largeDomain;
};

template<typename FieldT>
class SecretSharingInterface: SecretSharingNTT {
	public:
		/* Constructors allowing some flexibility in the choice of domains */
		SecretSharingInterface(size_t modulusIdx, size_t lSecretLength, size_t kDegree, size_t nShares, void (*computeSmall) (uint64_t*,uint64_t const*), void (*computeLarge) (uint64_t*,uint64_t const*)) : SecretSharingNTT(modulusIdx, lSecretLength, kDegree, nShares, computeSmall, computeLarge) {
			this->dCompositeDomainSize_ = nShares;
		};

		using SecretSharingNTT::degreeTest;
		using SecretSharingNTT::zeroTest;
		using SecretSharingNTT::zeroSumTest;
		using SecretSharingNTT::share;
		using SecretSharingNTT::reconstruct;

		int share(FieldT *secret);
		int reconstruct(FieldT *eval, bool expanded = false);

		int padMany(FieldT *secrets, size_t lrows);
		int padIntrablocRandomness(FieldT *secret);
		int shareMany(FieldT *secrets, size_t lrows);
		int reconstructMany(FieldT *secrets, size_t lrows);

		bool degreeTest(FieldT *eval);
		bool zeroTest(FieldT *eval, bool expanded = false);
		bool zeroSumTestSigmaProtocol(FieldT *eval, std::vector<uint64_t> beta);
		bool zeroSumTest(FieldT *eval);
		FieldT sumReconstruct(FieldT *eval);

	protected:
		int padSecret(FieldT *secret);
		int padPolynominal(FieldT *secret);
		int padShares(FieldT *secret);
		
		// adapter
		std::vector<uint64_t> demote(FieldT *eval, size_t length);
		void promote(uint64_t *eval, size_t length, FieldT *dest);
};

// Implementation
// =============================================================================================

// Secret Sharing Interface
template <typename FieldT>
inline int SecretSharingInterface<FieldT>::padSecret(FieldT *secret) {
	FieldT::randomVector(secret + this->lSecretLength_, this->kDegree_ - this->lSecretLength_, true);
	return 0;
}

template <typename FieldT>
inline int SecretSharingInterface<FieldT>::padIntrablocRandomness(FieldT *secret) {
	for (size_t i = this->lSecretLength_; i < kDegree_ ;i++) {
		secret[i] = FieldT(0); 
	}

	return 0;
}

template <typename FieldT>
inline int SecretSharingInterface<FieldT>::padPolynominal(FieldT *secret) {
	for (size_t i = this->kDegree_; i < nNumberShares_;i++) {
		secret[i] = FieldT(0);
	}

	return 0;
}

template <typename FieldT>
inline int SecretSharingInterface<FieldT>::padShares(FieldT *secret) {
	for (size_t i = this->nNumberShares_; i < dCompositeDomainSize_ ;i++) {
		secret[i] = FieldT(0);
	}

	return 0;
}

template <typename FieldT>
int SecretSharingInterface<FieldT>::padMany(FieldT *secret, size_t lrows) {

	// Add randomness for padding
	for (size_t i = 0; i < lrows; i++) {
		this->padSecret(secret + i * this->nNumberShares_);
	}

	return 0;
}

template <typename FieldT>
int SecretSharingInterface<FieldT>::share(FieldT *secret) {

	std::vector<uint64_t> secret64 = demote(secret, this->nNumberShares_); 
	this->share(&secret64[0]);
	promote(&secret64[0], this->nNumberShares_, secret); 

	return 0;
}

// this method is used when the resulting degree of the polynomial is larger than the domain space of the secret
template <typename FieldT>
int SecretSharingInterface<FieldT>::reconstruct(FieldT *eval, bool expandedDegree) {
	std::vector<uint64_t> eval64 = demote(eval, this->nNumberShares_); 
	this->reconstruct(&eval64[0], expandedDegree);
	promote(&eval64[0], this->kDegree_, eval); 

	return 0;
}

template <typename FieldT>
std::vector<uint64_t> SecretSharingInterface<FieldT>::demote(FieldT *eval, size_t length) {
	std::vector<uint64_t> temp(length);
	for (size_t idx = 0; idx < length; idx++) {
		temp[idx] = static_cast<uint64_t>(eval[idx].getValue());
	}

	return temp;
}

template <typename FieldT>
void SecretSharingInterface<FieldT>::promote(uint64_t *eval, size_t length, FieldT *dest) {
	for (size_t idx = 0; idx < length; idx++) {
		dest[idx] = FieldT(eval[idx]);
	}
}

template <typename FieldT>
bool SecretSharingInterface<FieldT>::degreeTest(FieldT *eval) {

	std::vector<uint64_t> eval64 = demote(eval,this->nNumberShares_); 
	return this->degreeTest(&eval64[0]);
}

template <typename FieldT>
bool SecretSharingInterface<FieldT>::zeroTest(FieldT *eval, bool largerDegree) {

	std::vector<FieldT> localCopy(eval, eval + this->nNumberShares_);
	reconstruct(&localCopy[0], largerDegree);

	for (size_t i = 0; i < this->lSecretLength_; i++) {
		if (!(localCopy[i] == FieldT(0))) return false;
	}

	return true;
}

template <typename FieldT>
bool SecretSharingInterface<FieldT>::zeroSumTest(FieldT *eval) {
	
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
FieldT SecretSharingInterface<FieldT>::sumReconstruct(FieldT *eval) {
	
	std::vector<FieldT> localCopy(eval, eval+ this->nNumberShares_);
	reconstruct(&localCopy[0], true);

	FieldT sum = FieldT(0);
	for (size_t i = 0; i < this->lSecretLength_; i++) {
		sum += localCopy[i];
	}

	return sum;
}

template <typename FieldT>
int SecretSharingInterface<FieldT>::shareMany(FieldT *secret, size_t lrows) {

	// std::cout << "Performing FFTs: " << lrows << std::endl; 
	for (size_t i = 0; i < lrows; i++) {
		share(secret + i * this->nNumberShares_);
	}

	return 0;
}

template <typename FieldT>
int SecretSharingInterface<FieldT>::reconstructMany(FieldT *secret, size_t lrows) {

	for (size_t i = 0; i < lrows; i++) {
		reconstruct(secret + i * this->nNumberShares_);
	}

	return 0;
}

template <typename FieldT>
bool SecretSharingInterface<FieldT>::zeroSumTestSigmaProtocol(FieldT *eval, std::vector<uint64_t> beta) {
	std::vector<uint64_t> eval64 = demote(eval, this->nNumberShares_); 
	
	return true;
}


// Secret Sharing NTT

int SecretSharingNTT::padSecret(uint64_t *secret) {

	// uint64_t::randomVector(secret + this->lSecretLength_, this->kDegree_ - this->lSecretLength_, true);
	return 0;
}

int SecretSharingNTT::padIntrablocRandomness(uint64_t *secret) {

	for (size_t i = this->lSecretLength_; i < kDegree_ ;i++) {
		secret[i] = uint64_t(0); 
	}

	return 0;
}

int SecretSharingNTT::padPolynominal(uint64_t *secret) {

	for (size_t i = this->kDegree_; i < nNumberShares_;i++) {
		secret[i] = uint64_t(0);
	}

	return 0;
}

int SecretSharingNTT::padShares(uint64_t *secret) {

	for (size_t i = this->nNumberShares_; i < dCompositeDomainSize_ ;i++) {
		secret[i] = uint64_t(0);
	}

	return 0;
}

int SecretSharingNTT::padMany(uint64_t *secret, size_t lrows) {

	// Add randomness for padding
	for (size_t i = 0; i < lrows; i++) {
		this->padSecret(secret + i * this->nNumberShares_);
	}

	return 0;
}

int SecretSharingNTT::share(uint64_t *data) {

	// Interpolate Poly from Secret
	_smallDomain->inv_ntt(data);

	this->padPolynominal(data);

	// Evaluate Publicly Shared Data from Poly
	_largeDomain->ntt(data);

	return 0;
}

// this method is used when the resulting degree of the polynomial is larger than the domain space of the secret
int SecretSharingNTT::reconstruct(uint64_t *eval, bool expandedDegree) {
	size_t degree = kDegree_;

	// Interpolate Poly from Publicly Shared Data
	_largeDomain->inv_ntt(eval);

	// if (expandedDegree) degree = kDegree_*2; 
	// Because we evaluate on roots of unity, we can use coefs[i] + coefs[i+k]  
	// on k coefficients instead of a 2k evaluation
	if (expandedDegree) {
		for (size_t i = 0; i < kDegree_; i++) {
				if (eval[i+kDegree_] > eval[i]) {
					eval[i] = eval[i] - eval[i + kDegree_] + params<uint64_t>::P[this->modulusIdx_];
				} else {
					eval[i] = eval[i] - eval[i + kDegree_];
				}
			}
	}

	// Evaluate Poly Into Secret
	_smallDomain->ntt(eval);

	// rebuild accounting for the stride implied by the larger degree
	size_t stride = degree/this->kDegree_;
	for (size_t col = 0; col < this->lSecretLength_; col++) {
		eval[col] = eval[col*stride];
	}

	return 0;
}

bool SecretSharingNTT::degreeTest(uint64_t *eval) {

	// Interpolate Poly from Publicly Shared Data
	std::vector<uint64_t> localCopy(eval, eval + this->nNumberShares_);
	_largeDomain->inv_ntt(&localCopy[0]);

	for (size_t i = kDegree_; i < nNumberShares_; i++) {
		if (!(localCopy[i] == uint64_t(0))) {
			DBG("val:"<< localCopy[i] << ",degree>" << i);
			return false;
		}
	}

	return true;
}

bool SecretSharingNTT::zeroTest(uint64_t *eval, bool largerDegree) {

	std::vector<uint64_t> localCopy(eval, eval + this->nNumberShares_);
	reconstruct(&localCopy[0], largerDegree);

	for (size_t i = 0; i < this->lSecretLength_; i++) {
		if (!(localCopy[i] == uint64_t(0))) return false;
	}

	return true;
}

bool SecretSharingNTT::zeroSumTest(uint64_t *eval) {
	
	std::vector<uint64_t> localCopy(eval, eval+ this->nNumberShares_);
	reconstruct(&localCopy[0], true);

	uint64_t sum = uint64_t(0);
	for (size_t i = 0; i < this->lSecretLength_; i++) {
		sum += localCopy[i];
	}

	if (sum == uint64_t(0)) return true;
	else return false;
}

bool SecretSharingNTT::zeroSumTestSigmaProtocol(uint64_t *eval, std::vector<uint64_t> beta) {

	std::vector<uint64_t> localCopy(eval, eval+ this->nNumberShares_);
	reconstruct(&localCopy[0], true);

	for (size_t i = 0; i < this->lSecretLength_; i++) {
		if (localCopy[i] != beta[i]) return false;
	}

	return true;
}

int SecretSharingNTT::shareMany(uint64_t *secret, size_t lrows) {

	// std::cout << "Performing FFTs: " << lrows << std::endl; 
	for (size_t i = 0; i < lrows; i++) {
		share(secret + i * this->nNumberShares_);
	}

	return 0;
}

int SecretSharingNTT::reconstructMany(uint64_t *secret, size_t lrows) {

	for (size_t i = 0; i < lrows; i++) {
		reconstruct(secret + i * this->nNumberShares_);
	}

	return 0;
}

}
