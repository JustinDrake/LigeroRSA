#pragma once

/** STL */
#include <algorithm>
#include <vector>

/** Third Party */
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

/** Ligero */
#include "Common.hpp"
#include "Hash.hpp"
#include "Math.hpp"
#include "MerkleTree.hpp"
#include "SecretSharingNTT.hpp"
#include "Transport.hpp"

#include "protocol/Builder.hpp"
#include "protocol/Constraint.hpp"
#include "protocol/ConstraintSystem.hpp"
#include "protocol/ExpressNPStatement.hpp"

/**
 * ==============================================================================
 */
/** ZKSnark Definitions */
#define OPEN_COLUMNS 2

#define STATE_UPDATE_LINEAR_COMBOS 1
#define STATE_UPDATE_OPEN_COLS 2

#define TEST_INTERLEAVED 0
#define TEST_LINEAR_INTERBLOC_CONSTRAINTS 1
#define TEST_LINEAR_INTRABLOC_CONSTRAINTS 2
#define TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS 3
#define TEST_QUADRATIC_CONSTRAINTS 4
#define OPEN_COLUMNS 5

std::vector<std::string> testDescription = {
    "TEST_INTERLEAVED",
    "TEST_LINEAR_INTERBLOC_CONSTRAINTS",
    "TEST_LINEAR_INTRABLOC_CONSTRAINTS",
    "TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS",
    "TEST_QUADRATIC_CONSTRAINTS",
    "OPEN_COLUMNS"};

constexpr size_t numberOfRounds = 1;
constexpr size_t numberOfIterations = 3;
constexpr uint8_t nbBlindingRows = 6;

/** These are the offset for the location various blinding rows from the # of
 * rows of the main witness */
std::vector<int8_t> offsetIntrablocBlinding = {-1, -2, -3};

std::vector<int8_t> offsetStitchingBlinding = {-4, -5, -6};

std::vector<int8_t> offsetLDTBlinding = {-6, -7, -8};

/**
 * ==============================================================================
 */
/** Structured Objects */
namespace ligero {
namespace zksnark {

template <typename FieldT>
std::vector<FieldT> pick(std::vector<FieldT> data,
                         std::vector<size_t> indices) {
  std::vector<FieldT> tmp(indices.size());
  transform(indices.begin(), indices.end(), tmp.begin(),
            [&](size_t index) { return data[index]; });

  return tmp;
}

/** Data Structures */

/** Fiat-Shamir Transformation: We generate "random" choices by the Verifier
 *  based on the values of the partial transcript
 */
template <typename FieldT>
class VerifierSimulation {
  uint64_t *&_randomness;
  size_t _modulusIdx;

 public:
  VerifierSimulation(uint64_t *&randomness, size_t modulusIdx,
                     const hash::HashIntegerSeedGeneratorF<FieldT> &generator,
                     const hash::HashStringSeedGeneratorF &generatorString)
      : _randomness(randomness),
        _modulusIdx(modulusIdx),
        _hashIntegerSeedGenerator(generator),
        _hashStringSeedGenerator(generatorString){};

  const hash::HashIntegerSeedGeneratorF<FieldT> &_hashIntegerSeedGenerator;
  const hash::HashStringSeedGeneratorF &_hashStringSeedGenerator;

  /** Consume preprocessed randomness to yield random field elements
   * @param randomness reference to the randomness matrix to populate
   * @param ncols number of columns in the matrix
   * @param nrows number of rows in the matrix
   * @return 0 if successful
   */
  int pickFieldTupleIntra(std::vector<FieldT> &randomness, size_t ncols,
                          size_t nrows) {
    size_t idx = 0;
    randomness.resize(nrows * ncols, FieldT(0));
    for (size_t row = 0; row < nrows; row++) {
      for (size_t col = 0; col < ncols; col++) {
        randomness[row * ncols + col] = FieldT(
            ligero::math::scaleToPrime(*this->_randomness++, _modulusIdx));
      }
    }

    return 0;
  }

  /** Consume preprocessed randomness to yield random field elements
   * @param randomness reference to the randomness vector to populate
   * @param n number of elements to generate
   * @return 0 if successful
   */
  int pickFieldTuple(std::vector<FieldT> &randomness, const size_t n) {
    randomness.resize(n);
    for (size_t idx = 0; idx < n; idx++) {
      randomness[idx] =
          FieldT(ligero::math::scaleToPrime(*this->_randomness++, _modulusIdx));
    }

    return 0;
  }

  /** Consume preprocessed randomness to yield random field elements
   * @param randomness reference to the randomness vector to populate
   * @return 0 if successful
   */
  int pickFieldTuple(std::vector<FieldT> &randomness) {
    randomness.push_back(
        FieldT(ligero::math::scaleToPrime(*this->_randomness++, _modulusIdx)));
    return 0;
  }

  /** Consume preprocessed randomness to yield random 64-bit unsigned integers
   * @param randomness reference to the randomness vector to populate
   * @param n number of elements to populate
   * @param upperBound integers will be in the interval [0,upperBound)
   * @return 0 if successful
   */
  int pickUintTuple(std::vector<size_t> &randomness, const size_t n,
                    size_t upperBound) {
    randomness.resize(n);
    for (size_t idx = 0; idx < n; idx++) {
      randomness[idx] =
          ligero::math::scaleToUpperBound(*this->_randomness++, upperBound);
    }

    return 0;
  }
};

/**
 *  Parameters for Interleaved Reed-Solomon Code:
 *  - each row is a linear code over F consisting of n-tuples (p(η1), ... ,
 *  p(ηn)) given η ∈ Fn and p is a polynomial of degree less than k.
 *  - m rows are interleaved.
 */
struct InterleavedRSCode {
  size_t n;
  size_t k;
  size_t m; /** Private witness */
};

/**
 *  Full transcript for a single round
 */
template <typename FieldT>
class roundFSTranscript {
 public:
  hash::digest proverCommitment;
  hash::digest proverCommitment_early;
  std::vector<std::unordered_map<size_t, std::vector<FieldT>>>
      proverResponsesLinearCombinations;
  std::vector<size_t> verifierColumnQueries;
  hash::sparseMultipleDecommitment<FieldT> proverDecommitment;
  hash::sparseMultipleDecommitment<FieldT> proverDecommitment_early;

 private:
  friend class boost::serialization::access;

  /** Implements seralization and deserialization for a round of transcript */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar & this->proverCommitment;
    ar & this->proverCommitment_early;
    ar & this->proverResponsesLinearCombinations;
    ar & this->proverDecommitment;
    ar & this->proverDecommitment_early;
    _initialized = true;
  }

  bool _initialized;
};

/** This structure holds on to the indices of the blocs where variables can be
 * found across multiple witnesses */
struct multiWitnessesReferences {
  std::vector<size_t> si_early_BlocIdxs;
  std::vector<size_t> si_BlocIdxs;

  std::vector<size_t> ei_early_BlocIdxs;
  std::vector<size_t> ei_BlocIdxs;

  std::vector<size_t> ri_early_BlocIdxs;
  std::vector<size_t> ri_BlocIdxs;

  size_t size(PublicData &p) {
    return (this->si_BlocIdxs.size() + this->ei_BlocIdxs.size() +
            ((p.modulusIdx < p.p) ? 1 : 0));
  }
};

template <typename FieldT>
using FullFSTranscript = std::vector<roundFSTranscript<FieldT>>;

/**
 * ==============================================================================
 */
/** Prover Column Commitment Scheme */
template <typename FieldT>
class ProverColumnCS : public ligero::MT_CommitmentScheme<FieldT> {
 public:
  ProverColumnCS(ligero::MTParameters<FieldT> pcs, FieldT *U)
      : ligero::MT_CommitmentScheme<FieldT>(
            pcs.hashInnerNodeF, pcs.hashLeafContentF, pcs.hashZKRandomnessF,
            pcs.digestLength, pcs.zkByteSize, pcs.leavesNumber,
            pcs.leavesEntropy),
        _U(U) {}

  hash::digest performInitialCommit(FieldT *U) { return (this->commit(U)); }

  hash::sparseMultipleDecommitment<FieldT> performProverDecommit(
      std::vector<size_t> &q) {
    return (this->decommit(q));
  } /** This includes the values picked as well as the MT authentication path,
       which are both added to the transcript */

  FieldT *_U;
};

template <typename FieldT>
class VerifierColumnCS : public ligero::MT_CommitmentScheme<FieldT> {
 public:
  VerifierColumnCS(ligero::MTParameters<FieldT> pcs)
      : ligero::MT_CommitmentScheme<FieldT>(
            pcs.hashInnerNodeF, pcs.hashLeafContentF, pcs.hashZKRandomnessF,
            pcs.digestLength, pcs.zkByteSize, pcs.leavesNumber,
            pcs.leavesEntropy) {}
};

/**
 * ==============================================================================
 */
/** Generic Non Interactive Protocol */
template <typename FieldT>
class NonInteractiveArgument {
 public:
  Statements _statement;

  /** Protocol Parameters */
  InterleavedRSCode &_rs;
  InterleavedRSCode _rs_early;
  size_t _l;
  size_t _t;
  hash::HashIntegerSeedGeneratorF<FieldT> _gen;
  hash::HashStringSeedGeneratorF _genS;
  ligero::MTParameters<FieldT> _mtParams;
  ligero::MTParameters<FieldT> _mtParams_early;

  /** Functionality */
  FullFSTranscript<FieldT> _transcript;
  VerifierSimulation<FieldT> _vs;

  /** Constraint System */
  ConstraintSystem<FieldT> _constraints_early;
  ConstraintSystem<FieldT> _constraints;

  /** Current Round */
  size_t _round;

  /** Timing Information */
  timers _myTimer;

  /** Inputs */
  PublicData &_publicData;
  SigmaProtocolPublicData &_sigmaProtocolPublicData;
  SecretData &_secretData;

  /** Randomness */
  uint64_t *_randomness;
  uint64_t *_alloc;

  /** Verifying Integrity Multi-Witness */
  multiWitnessesReferences _splitWitness;

  /** Size of the Constraint System Randomness */
  size_t _constraintSystemRandomnessSize;

  /** Cross-Moduli References */
  std::vector<size_t> _uxBlocks;
  std::vector<size_t> _uzBlocks;

  std::vector<size_t> _eiBlocks;
  std::vector<size_t> _siBlocks;

  std::vector<std::vector<size_t>> _varBlocks;
  uint64_t _stitching;

  /** Bit Reverse Templates */
  void (*_computeSmall)(uint64_t *, uint64_t const *);
  void (*_computeLarge)(uint64_t *, uint64_t const *);

  /** Maintaining a Hash */
  std::string _hashedPublicData;

 public:
  ~NonInteractiveArgument() {
    cleanup();
    cleanupEarly();
  }

  /** Constructor for the Generic Non Interactive Argument */
  NonInteractiveArgument(
      Statements statement, InterleavedRSCode &irs, size_t t, size_t l,
      const hash::HashIntegerSeedGeneratorF<FieldT> &generator,
      const hash::HashStringSeedGeneratorF &generatorString,
      ligero::MTParameters<FieldT> &mt, PublicData &pdata,
      SigmaProtocolPublicData &spdata, SecretData &sdata,
      void (*computeSmall)(uint64_t *, uint64_t const *),
      void (*computeLarge)(uint64_t *, uint64_t const *),
      std::string &hashedPublicData)
      : _statement(statement),
        _rs(irs),
        _t(t),
        _l(l),
        _gen(generator),
        _genS(generatorString),
        _mtParams(mt),
        _mtParams_early(mt),
        _vs(this->_randomness, pdata.modulusIdx, generator, generatorString),
        _publicData(pdata),
        _sigmaProtocolPublicData(spdata),
        _secretData(sdata),
        _computeSmall(computeSmall),
        _computeLarge(computeLarge),
        _hashedPublicData(hashedPublicData){

            /** We need to increase the leaves entropy by one to account for the
                     blinding row */
            /** _mtParams.leavesEntropy++; */
        };

  void cleanup() {
    for (constraint *p : this->_constraints.constraints) {
      delete p;
    }
  }

  void cleanupEarly() {
    for (constraint *p : this->_constraints_early.constraints) {
      delete p;
    }
  }

  /** Here we compute the product of the transpose of the randomness vector r,
   * namely the row vector (rT), with the linear constraint matrix A.
   * The linear constraint matrix A is defined as : Ax = b, where x ∈ Fml, b ∈
   * Fml and A ∈ Fml x ml
   * In this implementation, for a sparse representation of the matrix A, we
   * use the transformaton_constraint
   * class, each object representing the relationship between two elements of
   * the extended witness.
   * The resulting matrix (rT) x A will then be multiplied pointwise by our
   * secret matrix U,
   * subsequently shared, with a limited number of columns being opened based
   * on security preferences.
   * @param out reference to the (rT) x A matrix, stored in vector form
   * @param constraint reference to the transformation constraint processed
   * @param rml dense matrix of randomness
   * @param targetBlocIdx block index target for the generative constraint
   * @param processed reference to the set of blocks connected to transformation
   * constraints
   */
  void calculateRTAmatrix(std::vector<FieldT> &out,
                          transformation_constraint<FieldT> &constraint,
                          std::vector<FieldT> &rml, size_t targetBlocIdx,
                          std::unordered_set<size_t> &processed) {
    if (processed.find(targetBlocIdx) == processed.end()) {
      processed.insert(targetBlocIdx);
    }

    /** Browse through each component */
    for (size_t idx = 0; idx < constraint.scalar.size(); idx++) {
      size_t sourceWidx = (constraint.block[idx] * this->_rs.n) +
                          constraint.source_position[idx];
      size_t targetWidx =
          (targetBlocIdx * this->_rs.n) + constraint.target_position[idx];

      out[sourceWidx] +=
          constraint.scalar[idx] *
          rml[(targetBlocIdx * this->_l) + constraint.target_position[idx]];
      out[targetWidx] +=
          FieldT(FieldT::getModulus() - 1) *
          rml[(targetBlocIdx * this->_l) + constraint.target_position[idx]];

      if (processed.find(constraint.block[idx]) == processed.end()) {
        processed.insert(constraint.block[idx]);
      }
    }
  }

  /** Similar functionality, except on the (rT) x A' matrix here to generate
   * linear combinations used for the stitching of proofs accross moduli
   * @param out reference to the (rT) x A' matrix, stored in vector form
   * @param constraint reference to the transformation constraint processed
   * @param rml dense matrix of randomness
   * @param targetBlocIdx block index target for the generative constraint
   * @param processed reference to the set of blocks connected to transformation
   * @param iteration the index of the iteration in the proof
   * @param populate flag indicating whether we are storing randomness or
   * reusing what is already stored
   */
  void calculateRTAStitchingMatrix(
      std::vector<FieldT> &out, transformation_constraint<FieldT> &constraint,
      std::vector<FieldT> &rml, size_t targetBlocIdx,
      std::unordered_set<size_t> &processed, size_t iteration,
      bool populate = false) {
    if (constraint.stitching ==
        LinearCombinationsStitching::StitchingNoScalar) {
      /** Browse through each component */
      for (size_t idx = 0; idx < constraint.scalar.size(); idx++) {
        size_t sourceWidx = (constraint.block[idx] * this->_rs.n) +
                            constraint.source_position[idx];

        out[sourceWidx] += constraint.scalar[idx];

        if (processed.find(constraint.block[idx]) == processed.end()) {
          processed.insert(constraint.block[idx]);
        }
      }
    } else {
      auto findVarIdx = [&](size_t index) -> std::pair<size_t, size_t> {
        std::pair<size_t, size_t> out;
        for (size_t i = 0; i < this->_varBlocks.size(); i++) {
          auto it = find(this->_varBlocks[i].begin(), this->_varBlocks[i].end(),
                         index);

          if (it != this->_varBlocks[i].end()) {
            out = {i, it - this->_varBlocks[i].begin()};
            break;
          }
        }

        return out;
      };

      auto [varIdx, blockIdx] = findVarIdx(constraint.block[0]);

      if (this->_publicData.modulusIdx == 0) {
        /** Record randomness */
        for (size_t idx = 0; idx < constraint.scalar.size(); idx++) {
          size_t sourceWidx = (constraint.block[idx] * this->_rs.n) +
                              constraint.source_position[idx];
          FieldT scalar(rml[(targetBlocIdx * this->_l) +
                            constraint.target_position[idx]]);

          out[sourceWidx] += scalar;
          this->_secretData.ai[iteration][varIdx][blockIdx]
                              [constraint.source_position[idx]] =
              scalar.getValue();

          if (processed.find(constraint.block[idx]) == processed.end()) {
            processed.insert(constraint.block[idx]);
          }
        }
      } else {
        /** Use recorded randomness */
        for (size_t idx = 0; idx < constraint.scalar.size(); idx++) {
          size_t sourceWidx = (constraint.block[idx] * this->_rs.n) +
                              constraint.source_position[idx];

          out[sourceWidx] +=
              FieldT(this->_secretData.ai[iteration][varIdx][blockIdx]
                                         [constraint.source_position[idx]]);

          if (processed.find(constraint.block[idx]) == processed.end()) {
            processed.insert(constraint.block[idx]);
          }
        }
      }
    }
  }

  /** Generate Seed from Public Data and Transcript
   * @param transcript the full transcript up to this point
   * @return a hash of the combined information
   */
  std::vector<unsigned char> generateFSSeed(
      FullFSTranscript<FieldT> &transcript) {
    std::string msg;

    msg += this->_hashedPublicData;
    msg += std::to_string(this->_publicData.modulusIdx);
    msg += serialize(transcript[0].proverCommitment);
    msg += serialize(transcript[0].proverDecommitment);

    for (size_t idx = 0;
         idx < transcript[0].proverResponsesLinearCombinations.size(); idx++) {
      msg += serialize(
          transcript[0]
              .proverResponsesLinearCombinations[idx][TEST_INTERLEAVED]);
      msg += serialize(transcript[0].proverResponsesLinearCombinations
                           [idx][TEST_LINEAR_INTERBLOC_CONSTRAINTS]);
      msg += serialize(transcript[0].proverResponsesLinearCombinations
                           [idx][TEST_LINEAR_INTRABLOC_CONSTRAINTS]);
      msg += serialize(transcript[0].proverResponsesLinearCombinations
                           [idx][TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS]);
      msg += serialize(
          transcript[0]
              .proverResponsesLinearCombinations[idx]
                                                [TEST_QUADRATIC_CONSTRAINTS]);
    }

    return (this->_vs._hashStringSeedGenerator(msg.c_str(), msg.size()));
  }

  /** Bulk generation of randomness based on the state of the transcript up to
   * this point
   * @param sizeRandomness number of random field elements to generate at once
   * @param transcript the full transcript up to this point
   */
  void pickFiatShamirRandomness(size_t sizeRandomness,
                                FullFSTranscript<FieldT> &transcript) {
    std::vector<unsigned char> seed = this->generateFSSeed(transcript);
    this->_randomness =
        static_cast<uint64_t *>(malloc(sizeRandomness * sizeof(uint64_t)));
    randombytes_buf_deterministic((void *)this->_randomness,
                                  sizeRandomness * sizeof(uint64_t), &seed[0]);
    this->_alloc = this->_randomness;
  }

  /** These will be implemented by specialized classes, as the logic will differ
   * here between the prover and the verifier */
  virtual void runRound(size_t, bool) = 0;

  /** Running the protocol
   * @param nbRounds number of rounds included in the protocol (1 for the RSA
   * Ceremony)
   * @param earlyOnly if true, all that will be done is to generate the early
   * witness, no proof generated at this point
   */
  void runProtocol(size_t nbRounds, bool earlyOnly = false) {
    for (size_t round = 0; round < nbRounds; round++)
      runRound(round, earlyOnly);
  }

  /** Pick challenges for the simulated verifier based on Fiat-Shamir
   * transformation
   * @param round index of the current protocol round
   * @param phase the type of test we are generating challenges for
   * @param r reference to the vector of challenges
   * @param size target size for the challenge vector
   */
  void pickVerifierChallenges(size_t round, unsigned int phase,
                              vector<FieldT> &r, size_t size = 0) {
    /** Pick randomness */
    switch (phase) {
      case TEST_INTERLEAVED:
        this->_vs.pickFieldTuple(r, size); /** r ∈ Fm */
        break;

      case TEST_QUADRATIC_CONSTRAINTS:
      case TEST_LINEAR_INTERBLOC_CONSTRAINTS:
        this->_vs.pickFieldTuple(r);
        break;

      case TEST_LINEAR_INTRABLOC_CONSTRAINTS:
        this->_vs.pickFieldTupleIntra(r, this->_l, size); /** r ∈ Fml */
        break;

      default:
        throw internalError("incorrect phase for pickVerifierChallenges: " +
                            std::to_string(phase));
        break;
    }
  }

  /** Pick challenges for the simulated verifier for opening columns based on
   * Fiat-Shamir transformation
   * @param round index of the current protocol round
   * @param phase the type of test we are generating challenges for (here
   * opening columns)
   * @param r reference to the vector of challenges
   */
  void pickVerifierChallenges(size_t round, unsigned int phase,
                              std::vector<size_t> &r) {
    /** Pick randomness */
    switch (phase) {
      case OPEN_COLUMNS:
        this->_vs.pickUintTuple(r, this->_t, this->_rs.n);

        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());

        break;

      default:
        throw internalError("incorrect phase for pickVerifierChallenges: " +
                            std::to_string(phase));
        break;
    }
  }

  /** Interleaved: prover helper function to compute linear combinations over
   * the full witness as the prover response to the verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param interleaved_r verifiers challenges
   * @param U pointer to the witness
   */
  void interbloc(std::vector<FieldT> &res,
                 const std::vector<FieldT> &interleaved_r, const FieldT *U) {
    if (res.size() == 0) {
      res.resize(this->_rs.n, FieldT(0));
    }

    for (size_t col = 0; col < this->_rs.n; col++) {
      for (size_t row = 0; row < interleaved_r.size(); row++) {
        res[col] += U[col + (row * this->_rs.n)] * interleaved_r[row];
      }
    }
  }

  /** Interleaved: verifier helper function to compute linear combinations over
   * the opened columns as the prover response to the verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param interleaved_r verifiers challenges
   * @param Uq matrix of opened columns
   */
  void interbloc(std::vector<FieldT> &res,
                 const std::vector<FieldT> &interleaved_r,
                 const std::vector<std::vector<FieldT>> &Uq) {
    if (res.size() == 0) {
      res.resize(Uq.size(), FieldT(0));
    }

    for (size_t col = 0; col < Uq.size(); col++) {
      for (size_t row = 0; row < interleaved_r.size(); row++) {
        res[col] += Uq[col][row] * interleaved_r[row];
      }
    }
  }

  /** Interbloc Linear Constraint: prover helper function to compute linear
   * combinations over the full witness as the prover response to the
   * verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param block indices of the input blocks involved in the linear constraint
   * @param scalar scalars for each of the blocks
   * @param target index of the target block for the linear constraint (output)
   * @param r verifiers challenges
   * @param U pointer to the witness
   */
  void interbloc(std::vector<FieldT> &res, const std::vector<size_t> &block,
                 const std::vector<FieldT> &scalar, const size_t target,
                 const FieldT &r, const FieldT *U) {
    if (res.size() == 0) {
      res.resize(this->_rs.n, FieldT(0));
    }

    for (size_t col = 0; col < this->_rs.n; col++) {
      for (size_t idx = 0; idx < block.size(); idx++) {
        res[col] += U[col + (block[idx] * this->_rs.n)] * scalar[idx] * r;
      }
      res[col] = res[col] - U[col + (target * this->_rs.n)] * r;
    }
  }

  /** Interbloc Linear Constraint: verifier helper function to compute linear
   * combinations over opened columns as the prover response to the
   * verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param block indices of the input blocks involved in the linear constraint
   * @param scalar scalars for each of the blocks
   * @param target index of the target block for the linear constraint (output)
   * @param r verifiers challenges
   * @param Uq matrix of opened columns
   */
  void interbloc(std::vector<FieldT> &res, const std::vector<size_t> &block,
                 const std::vector<FieldT> &scalar, const size_t target,
                 const FieldT &r, const std::vector<std::vector<FieldT>> &Uq) {
    if (res.size() == 0) {
      res.resize(Uq.size(), FieldT(0));
    }

    for (size_t col = 0; col < Uq.size(); col++) {
      for (size_t idx = 0; idx < block.size(); idx++) {
        res[col] += Uq[col][block[idx]] * scalar[idx] * r;
      }

      res[col] = res[col] - Uq[col][target] * r;
    }
  }

  /** Interbloc Linear Constraint Across Witnesses : prover helper
   * function to compute linear combinations across witnesses to prove
   * consistency
   * @param res reference to the vector containing the linear combinations
   * @param block_early indices of the target blocks in the early witness
   * @param block indices of the corresponding target blocks in the main witness
   * @param r verifiers challenges
   * @param U_early pointer to the early witness
   * @param U pointer to the witness
   */
  void splitWitnessInterbloc(std::vector<FieldT> &res,
                             const size_t &block_early, const size_t &block,
                             const FieldT &r, FieldT *U_early, FieldT *U) {
    if (res.size() == 0) {
      res.resize(this->_rs.n, FieldT(0));
    }

    for (size_t col = 0; col < this->_rs.n; col++) {
      res[col] += U_early[col + (block_early * this->_rs.n)] * r;
      res[col] -= U[col + (block * this->_rs.n)] * r;
    }
  }

  /** Interbloc Linear Constraint Across Witnesses : verifier helper
   * function to compute linear combinations over opened columns across
   * witnesses to prove consistency
   * @param res reference to the vector containing the linear combinations
   * @param block_early indices of the target blocks in the early witness
   * @param block indices of the corresponding target blocks in the main witness
   * @param r verifiers challenges
   * @param Uq_early matrix of opened columns in the early witness
   * @param Uq matrix of opened columns in the main witness
   */
  void splitWitnessInterbloc(std::vector<FieldT> &res,
                             const size_t &block_early, const size_t &block,
                             const FieldT &r,
                             const std::vector<std::vector<FieldT>> &Uq_early,
                             const std::vector<std::vector<FieldT>> &Uq) {
    if (res.size() == 0) {
      res.resize(Uq.size(), FieldT(0));
    }

    for (size_t col = 0; col < Uq.size(); col++) {
      res[col] += Uq_early[col][block_early] * r;
      res[col] -= Uq[col][block] * r;
    }
  }

  /** Quadratic Constraint: prover helper function to compute linear
   * combinations over the full witness as the prover response to the
   * verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param block_pl indices of the blocks containing the left-handside operand
   * of the quadratic constraint
   * @param block_pr indices of the blocks containing the right-handside operand
   * of the quadratic constraint
   * @param scalar scalars for each of the blocks
   * @param target index of the target block for the quadratic constraint
   * (output)
   * @param r verifier's challenges
   * @param U pointer to the witness
   */
  void interbloc(std::vector<FieldT> &res, const std::vector<size_t> &block_pl,
                 const std::vector<size_t> &block_pr,
                 const std::vector<FieldT> &scalar, const size_t target,
                 const FieldT &r, const FieldT *U) {
    if (res.size() == 0) {
      res.resize(this->_rs.n, FieldT(0));
    }

    for (size_t col = 0; col < this->_rs.n; col++) {
      for (size_t idx = 0; idx < block_pl.size(); idx++) {
        res[col] += U[col + (block_pl[idx] * this->_rs.n)] *
                    U[col + (block_pr[idx] * this->_rs.n)] * scalar[idx] * r;
      }

      res[col] = res[col] - U[col + (target * this->_rs.n)] * r;
    }
  }

  /** Quadratic Constraint: verifier helper function to compute linear
   * combinations over opened columns as the prover response to the
   * verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param block_pl indices of the blocks containing the left-handside operand
   * of the quadratic constraint
   * @param block_pr indices of the blocks containing the right-handside operand
   * of the quadratic constraint
   * @param scalar scalars for each of the blocks
   * @param target index of the target block for the quadratic constraint
   * (output)
   * @param r verifier's challenges
   * @param Uq matrix of opened columns
   */
  void interbloc(std::vector<FieldT> &res, const std::vector<size_t> &block_pl,
                 const std::vector<size_t> &block_pr,
                 const std::vector<FieldT> &scalar, const size_t target,
                 const FieldT &r, const std::vector<std::vector<FieldT>> &Uq) {
    if (res.size() == 0) {
      res.resize(Uq.size(), FieldT(0));
    }

    for (size_t col = 0; col < Uq.size(); col++) {
      for (size_t idx = 0; idx < block_pl.size(); idx++) {
        res[col] +=
            Uq[col][block_pl[idx]] * Uq[col][block_pr[idx]] * scalar[idx] * r;
      }

      res[col] = res[col] - Uq[col][target] * r;
    }
  }

  /** Transformation Constraint: verifier helper function to compute linear
   * combinations over opened columns as the prover response to the
   * verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param r verifiers challenges
   * @param Uq matrix of opened columns in the main witness
   * @param Q positions of the opened columns
   * @param processed reference to the set of blocks involved in this linear
   * constraint
   * @param iteration index of the current iteration within the current round
   */
  void intrabloc(std::vector<FieldT> &res, const std::vector<FieldT> &r,
                 const std::vector<std::vector<FieldT>> &Uq,
                 const vector<size_t> &Q, std::unordered_set<size_t> &processed,
                 uint8_t iteration) {
    if (res.size() == 0) {
      res.resize(Uq.size(), FieldT(0));
    }

    for (size_t col = 0; col < Uq.size(); col++) {
      for (size_t row : processed) {
        res[col] += Uq[col][row] * r[Q[col] + (row * this->_rs.n)];
      }

      res[col] += Uq[col][this->_rs.m + offsetIntrablocBlinding[iteration]];
    }
  }

  /** Stitching: verifier helper function to compute linear
   * combinations over opened columns as the prover response to the
   * verifier's challenges
   * @param res reference to the vector containing the linear combinations
   * @param r verifiers challenges
   * @param Uq matrix of opened columns in the main witness
   * @param Q positions of the opened columns
   * @param processed reference to the set of blocks involved in this linear
   * constraint
   * @param iteration index of the current iteration within the current round
   */
  void intrablocStitching(std::vector<FieldT> &res,
                          const std::vector<FieldT> &r,
                          const std::vector<std::vector<FieldT>> &Uq,
                          const vector<size_t> &Q,
                          std::unordered_set<size_t> &processed,
                          uint8_t iteration) {
    if (res.size() == 0) {
      res.resize(Uq.size(), FieldT(0));
    }

    for (size_t col = 0; col < Uq.size(); col++) {
      for (size_t row : processed) {
        res[col] += Uq[col][row] * r[Q[col] + (row * this->_rs.n)];
      }

      res[col] += Uq[col][this->_rs.m + offsetStitchingBlinding[iteration]];
    }
  }
};

template <typename FieldT>
class Prover : public NonInteractiveArgument<FieldT> {
  /** Extended Witness */
  FieldT *_w_early;
  FieldT *_w;
  ligero::ZeroMQClientTransport &_transport;

 public:
  /** Constructor for the Prover */
  Prover(Statements statement, InterleavedRSCode &irs, size_t t,
         size_t blocSize,
         const hash::HashIntegerSeedGeneratorF<FieldT> &generator,
         const hash::HashStringSeedGeneratorF &generatorString,
         ligero::MTParameters<FieldT> &mt, PublicData &pdata,
         SigmaProtocolPublicData &spdata, SecretData &sdata,
         void (*computeSmall)(uint64_t *, uint64_t const *),
         void (*computeLarge)(uint64_t *, uint64_t const *),
         ligero::ZeroMQClientTransport &transport,
         std::string &hashedPublicData)
      : NonInteractiveArgument<FieldT>(
            statement, irs, t, blocSize, generator, generatorString, mt, pdata,
            spdata, sdata, computeSmall, computeLarge, hashedPublicData),
        _transport(transport) {}

  /** Updating an aggregate linear combination vector with the prover's response
   * to the simulated verifier's challenges for a specific constraint
   * @param round current round
   * @param target target block for the constraint handled
   * @param ct pointer to the constraint being handled
   * @param type type of underlying constraint
   * @param r verifier challenges
   * @param rU reference to an aggregate linear combination vector
   */
  inline void obtainProverResponses(size_t round, size_t target, constraint *ct,
                                    int type, const vector<FieldT> &r,
                                    vector<FieldT> &rU) {
    ligero::SecretSharingInterface<FieldT> localSecretSharing(
        this->_publicData.modulusIdx, (size_t)this->_l, (size_t)this->_rs.k,
        (size_t)this->_rs.n, this->_computeSmall, this->_computeLarge);

    switch (type) {
      case TEST_LINEAR_INTERBLOC_CONSTRAINTS: {
        linear_constraint<FieldT> constraint =
            *static_cast<const linear_constraint<FieldT> *>(ct);
        this->interbloc(rU, constraint.block, constraint.scalar, target,
                        r.back(), this->_w);

#ifndef NDEBUG

        {
          std::string report(std::to_string(target) +
                             std::string(",test,interbloc,blocs,"));
          for (size_t idx = 0; idx < constraint.block.size() - 1; idx++) {
            report += std::to_string(constraint.block[idx]) + std::string(",");
          }
          report += std::to_string(constraint.block.back());

          DBG(report);
        }

        {
          std::string report(std::to_string(target) +
                             std::string(",test,interbloc,coefs,"));
          for (size_t idx = 0; idx < constraint.scalar.size() - 1; idx++) {
            report += std::to_string(constraint.scalar[idx].getValue()) +
                      std::string(",");
          }
          report += std::to_string(constraint.scalar.back().getValue());

          DBG(report);
        }

#endif

      } break;

      case TEST_LINEAR_INTRABLOC_CONSTRAINTS:
        break;

      case TEST_QUADRATIC_CONSTRAINTS: {
        quadratic_constraint<FieldT> constraint =
            *static_cast<const quadratic_constraint<FieldT> *>(ct);
        this->interbloc(rU, constraint.left_block, constraint.right_block,
                        constraint.scalar, target, r.back(), this->_w);

#ifndef NDEBUG

        {
          std::string report(std::to_string(target) +
                             std::string(",quadratic,left_bloc,"));
          for (size_t idx = 0; idx < constraint.left_block.size() - 1; idx++) {
            report +=
                std::to_string(constraint.left_block[idx]) + std::string(",");
          }
          report += std::to_string(constraint.left_block.back());

          DBG(report);
        }

        {
          std::string report(std::to_string(target) +
                             std::string(",quadratic,right_bloc,"));
          for (size_t idx = 0; idx < constraint.right_block.size() - 1; idx++) {
            report +=
                std::to_string(constraint.right_block[idx]) + std::string(",");
          }
          report += std::to_string(constraint.right_block.back());

          DBG(report);
        }

        {
          std::string report(std::to_string(target) +
                             std::string(",quadratic,coefs,"));
          for (size_t idx = 0; idx < constraint.scalar.size() - 1; idx++) {
            report += std::to_string(constraint.scalar[idx].getValue()) +
                      std::string(",");
          }
          report += std::to_string(constraint.scalar.back().getValue());

          DBG(report);
        }

#endif
      } break;
      default:
        throw std::runtime_error("incorrect test type for prover responses:" +
                                 std::to_string(type));
        break;
    }
  }

  /** Builds the constraint system and witness for the statement currently
   * stored, then generates a ZKP and sends it to the coordinator
   * @param nbRounds total number of rounds (1 for the RSA Ceremony)
   * @param earlyCommitmentOnly if true, all we will do is build and commit the
   * merkle tree for the early witness
   */
  void produceArgumentOfKnowledge(size_t nbRounds, bool earlyCommitmentOnly) {
    this->_myTimer.initialize_timer();
    this->_myTimer.begin(1, "Producing Proof", 0);

    /** Initialize Fpp */
    /** =========================== */

    auto assessBlocIdx = [&](auto &builder,
                             std::vector<size_t> idxs) -> std::vector<size_t> {
      std::vector<size_t> ret;
      for (auto idx : idxs) {
        ret.emplace_back(builder.proof_values_location[idx]);
      }

      return ret;
    };

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      /** Build a first witness/constraint system for Early Commitment */
      /** ================================================================== */

      this->_rs_early = this->_rs;
      this->_myTimer.begin(
          2, "Build Extended Witness & Constraint System For Early Commitment",
          0);
      builder<FieldT> helper(this->_l);
      vector<vector<FieldT>> earlyCommitments2DVector;

      expressNPStatement<FieldT> expressStatementEarlyCommitment(
          this->_publicData, this->_sigmaProtocolPublicData, this->_secretData,
          helper, &earlyCommitments2DVector, this->_l, this->_randomness);
      expressStatementEarlyCommitment.buildConstraintSet(
          Statements::RSACeremony_EarlyCommitments);

      this->_myTimer.begin(3, "compile", 0);
      this->_constraints_early = helper.compile(earlyCommitments2DVector);
      this->_myTimer.end(3, "compile", 0);

      /** Dimensioning the circuit accordingly */
      /** We reserve one additional row for blinding */
      this->_mtParams_early.leavesEntropy = this->_rs_early.m =
          earlyCommitments2DVector.size() +
          this->_constraints_early.constant_blocks + 1;

      /** Record Position of Variables in the Witness */
      this->_splitWitness.si_early_BlocIdxs = assessBlocIdx(
          helper, expressStatementEarlyCommitment._siCoefsBlocIds);
      this->_splitWitness.ei_early_BlocIdxs = assessBlocIdx(
          helper, expressStatementEarlyCommitment._eiCoefsBlocIds);
      if (this->_publicData.modulusIdx < this->_publicData.p)
        this->_splitWitness.ri_early_BlocIdxs = assessBlocIdx(
            helper, {expressStatementEarlyCommitment._riCoefsBlocId});

      /** Embedding the early commit witness */
      /** ================================== */

      this->_myTimer.begin(3, "Embedding the early extended witness", 0);

      /** Migrate the extended witness to an expandable matrix */
      FieldT *earlyWitnessCArray =
          new FieldT[this->_rs_early.n * this->_rs_early.m]{FieldT(0)};

      {
        /** Populate the private witness */
        size_t idx = 0;
        for (auto &row : earlyCommitments2DVector) {
          memcpy(static_cast<void *>(
                     &earlyWitnessCArray[(idx++) * this->_rs_early.n]),
                 static_cast<void *>(&row[0]), row.size() * sizeof(FieldT));
          row.clear();
        }

        /** Add public constants */
        size_t topConstantRow =
            this->_rs_early.m - this->_constraints.constant_blocks -
            1; /** Accounting for constants and blinding row */
        for (size_t row = 0; row < this->_constraints_early.constant_blocks;
             row++) {
          for (size_t col = 0; col < this->_l; col++) {
            earlyWitnessCArray[(topConstantRow + row) * this->_rs_early.n +
                               col] =
                this->_constraints_early
                    .constants[col * this->_constraints.constant_blocks + row];
          }
        }
      }

      this->_w_early = earlyWitnessCArray;

      this->_myTimer.end(3, "Embedding the early extended witness", 0);
      this->_myTimer.end(
          2, "Build Extended Witness & Constraint System For Early Commitment",
          0);
    }

    if (earlyCommitmentOnly) {
      /** Dimension Randomness for the Constraint System */
      this->_constraintSystemRandomnessSize = determineCSrandomnessSize(
          this->_statement, this->_publicData.ptNumberEvaluation);

      this->_myTimer.begin(2, "Generating the Proof", 0);
      this->runProtocol(nbRounds, earlyCommitmentOnly);
      this->_myTimer.end(2, "Generating the Proof", 0);

      return;
    }

    /** Building the constraint system and populating the extended witness */
    /** ================================================================== */

    this->_myTimer.begin(2, "Build Extended Witness & Constraint System", 0);
    vector<vector<FieldT>> extendedWitness2DVector;

    {
      builder<FieldT> bob(this->_l);

      this->_myTimer.begin(3, "Express NP statement", 0);
      expressNPStatement<FieldT> expressStatement(
          this->_publicData, this->_sigmaProtocolPublicData, this->_secretData,
          bob, &extendedWitness2DVector, this->_l, this->_randomness);
      expressStatement.buildConstraintSet(this->_statement, false);
      this->_myTimer.end(3, "Express NP statement", 0);

      this->_myTimer.begin(3, "compile", 0);
      this->_constraints = bob.compile(
          extendedWitness2DVector); /** Main functionality to retain */
      this->_myTimer.end(3, "compile", 0);

      this->_myTimer.end(2, "Build Extended Witness & Constraint System", 0);
      DBG(this->_constraints.constraints.size());

      /** Details of the witness constructed */
      /** We need to add */
      DBG("witness rows:" << extendedWitness2DVector.size());
      DBG("witness cols:" << extendedWitness2DVector[0].size());
      this->_mtParams.leavesEntropy = this->_rs.m =
          extendedWitness2DVector.size() + this->_constraints.constant_blocks +
          nbBlindingRows;

      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        /** Record Position of Variables in the Witness */
        this->_splitWitness.si_BlocIdxs =
            assessBlocIdx(bob, expressStatement._siBlocs);
        this->_splitWitness.ei_BlocIdxs =
            assessBlocIdx(bob, expressStatement._eiBlocs);
        this->_splitWitness.ri_BlocIdxs =
            assessBlocIdx(bob, {expressStatement._riCoefsBlocId});
      }

      /** Here checking that we got the referencing right */
      this->_uxBlocks = assessBlocIdx(bob, expressStatement._uxBlocs);
      this->_uzBlocks = assessBlocIdx(bob, expressStatement._uzBlocs);

      this->_siBlocks = assessBlocIdx(bob, expressStatement._siBlocs);
      this->_eiBlocks = assessBlocIdx(bob, expressStatement._eiBlocs);

      /** Storing for Cross-Moduli Integrity Checking */
      this->_varBlocks.emplace_back(this->_siBlocks);
      this->_varBlocks.emplace_back(this->_eiBlocks);
      this->_varBlocks.emplace_back(this->_uxBlocks);
      this->_varBlocks.emplace_back(this->_uzBlocks);
    }

    /** Embedding the extended witness */
    /** ============================== */

    this->_myTimer.begin(2, "Embedding the extended witness", 0);

    /** Migrate the extended witness to an expandable matrix */
    FieldT *extendedWitnessCArray =
        new FieldT[this->_rs.n * this->_rs.m]{FieldT(0)};

    {
      /** Populate the private witness */
      size_t idx = 0;
      for (auto &row : extendedWitness2DVector) {
        memcpy(
            static_cast<void *>(&extendedWitnessCArray[(idx++) * this->_rs.n]),
            static_cast<void *>(&row[0]), row.size() * sizeof(FieldT));
        row.clear();
      }

      /** Add public constants */
      size_t topConstantRow =
          this->_rs.m - this->_constraints.constant_blocks -
          nbBlindingRows; /** Accounting for constants and blinding row */
      for (size_t row = 0; row < this->_constraints.constant_blocks; row++) {
        for (size_t col = 0; col < this->_l; col++) {
          extendedWitnessCArray[(topConstantRow + row) * this->_rs.n + col] =
              this->_constraints
                  .constants[col * this->_constraints.constant_blocks + row];
        }
      }
    }

    this->_w = extendedWitnessCArray;

    this->_myTimer.end(2, "Embedding the extended witness", 0);

    /** Debug: dump of the extended witness and constraints [first 10 elements
     * of each row] */
    /**
     * ===================================================================================
     */
    std::stringstream ss;
    ss << boost::uuids::random_generator()();
    std::string id(ss.str());
    constexpr auto nDisplayCol = 20;

    for (size_t row = 0; row < this->_rs.m; row++) {
      std::string line = id + std::string(",dump_witness,");

      line += std::to_string(row) + std::string(",");
      for (size_t col = 0; col < nDisplayCol; col++) {
        char buf[100];
        sprintf(buf, "%lu",
                (uint64_t)this->_w[col + row * this->_rs.n].getValue());
        line += std::string(buf) + std::string(",");
      }

      DBG(line);
    }

    /** Testing constraints */
    size_t target = 0;

#ifndef NDEBUG
    /** Dump of the constraint system for debug */
    for (auto constraint : this->_constraints.constraints) {
      switch (constraint->type) {
        case linear: {
          linear_constraint<FieldT> ct =
              *static_cast<const linear_constraint<FieldT> *>(constraint);

          std::string s =
              id + std::string(",dump_constraints,linear_constraint,block,");
          for (size_t i = 0; i < std::min(nDisplayCol, (int)ct.block.size());
               i++) {
            s += std::to_string(ct.block[i]) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,scalar,");
          for (size_t i = 0; i < std::min(nDisplayCol, (int)ct.scalar.size());
               i++) {
            char buf[100];
            sprintf(buf, "%lu", ct.scalar[i].getValue());
            s += std::string(buf) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,target_block,");
          s += std::to_string(this->_constraints.targets[target]) +
               std::string(",");
          DBG(s);
        } break;

        case transformation: {
          transformation_constraint<FieldT> ct =
              *static_cast<const transformation_constraint<FieldT> *>(
                  constraint);

          std::string s =
              id +
              std::string(",dump_constraints,transformation,source_block,");
          for (size_t i = 0; i < std::min(nDisplayCol, (int)ct.block.size());
               i++) {
            s += std::to_string(ct.block[i]) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,source_pos,");
          for (size_t i = 0;
               i < std::min(nDisplayCol, (int)ct.source_position.size()); i++) {
            char buf[100];
            sprintf(buf, "%lu", ct.source_position[i]);
            s += std::string(buf) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,target_pos,");
          for (size_t i = 0;
               i < std::min(nDisplayCol, (int)ct.target_position.size()); i++) {
            char buf[100];
            sprintf(buf, "%lu", ct.target_position[i]);
            s += std::string(buf) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,target_block,");
          s += std::to_string(this->_constraints.targets[target]) +
               std::string(",");
          DBG(s);
        } break;

        case quadratic: {
          quadratic_constraint<FieldT> ct =
              *static_cast<const quadratic_constraint<FieldT> *>(constraint);

          std::string s =
              id +
              std::string(",dump_constraints,quadratic_constraint,left_block,");
          for (size_t i = 0;
               i < std::min(nDisplayCol, (int)ct.left_block.size()); i++) {
            s += std::to_string(ct.left_block[i]) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,right_bloc,");
          for (size_t i = 0;
               i < std::min(nDisplayCol, (int)ct.right_block.size()); i++) {
            s += std::to_string(ct.right_block[i]) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,scalar,");
          for (size_t i = 0; i < std::min(nDisplayCol, (int)ct.scalar.size());
               i++) {
            char buf[100];
            sprintf(buf, "%lu", ct.scalar[i].getValue());
            s += std::string(buf) + std::string(",");
          }
          DBG(s);

          s = id + std::string(",dump_constraints,-,target_block,");
          s += std::to_string(this->_constraints.targets[target]) +
               std::string(",");
          DBG(s);
        } break;

        default:
          throw std::runtime_error("incorrect constraint type for test:" +
                                   std::to_string(constraint->type));
          break;
      }

      target++;
    }
#endif

    /** Generating and sending the proof */
    /** ================================ */

    /** Dimension Randomness for the Constraint System */
    this->_constraintSystemRandomnessSize = determineCSrandomnessSize(
        this->_statement, this->_publicData.ptNumberEvaluation);

    this->_myTimer.begin(2, "Generating the Proof", 0);
    this->runProtocol(nbRounds, earlyCommitmentOnly);
    this->_myTimer.end(2, "Generating the Proof", 0);

    this->_myTimer.begin(2, "Sending the Proof", 0);
    this->sendArgumentOfKnowledge();
    this->_myTimer.end(2, "Sending the Proof", 0);

    /** Clean Up */
    delete[] extendedWitnessCArray;

    this->_myTimer.end(1, "Producing Proof", 0);
  }

  /** Once the coordinator announces it has provided verifiers with
   * the data linked to the sigma protocol then we can start sending our proofs
   */
  void sendArgumentOfKnowledge() {
#ifndef NDEBUG
    DBG("Linear Combinations Sizes Sent By Prover");
    DBG("========================================");
    DBG(this->_transcript[0]
            .proverResponsesLinearCombinations[0][TEST_INTERLEAVED]
            .size());
    DBG(this->_transcript[0]
            .proverResponsesLinearCombinations
                [0][TEST_LINEAR_INTERBLOC_CONSTRAINTS]
            .size());
    DBG(this->_transcript[0]
            .proverResponsesLinearCombinations
                [0][TEST_LINEAR_INTRABLOC_CONSTRAINTS]
            .size());
    DBG(this->_transcript[0]
            .proverResponsesLinearCombinations
                [0][TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS]
            .size());
    DBG(this->_transcript[0]
            .proverResponsesLinearCombinations[0][TEST_QUADRATIC_CONSTRAINTS]
            .size());
#endif

    /** Send Data and the argument to the Coordinator */

    if (this->_publicData.modulusIdx == 0) {
      auto msg = _transport.awaitReply();
      if (hasError(msg)) {
        throw std::runtime_error(showError(getError(msg)));
      }
      assert(getResult(msg) == MessageType::GATHER_PROOFS);
    }

    _transport.send(MessageType::GATHER_PROOFS, this->_transcript, false);
  }

  /** This section will add multiple blinding rows:
   * two for LDT testing
   * two for regular intra-bloc constraints
   * two for stitching proofs together */
  void addBlindingRows() {
    for (uint8_t blindingRowIdx = 0; blindingRowIdx < nbBlindingRows;
         blindingRowIdx++) {
      std::vector<FieldT> blindingRow(this->_l);

      FieldT::randomVector(&blindingRow[0], this->_l - 1, true);
      FieldT sum = FieldT(0);
      for (size_t idx = 0; idx < this->_l - 1; idx++) {
        sum += blindingRow[idx];
      }
      blindingRow[this->_l - 1] = FieldT(0) - sum;

      memcpy(this->_w + ((this->_rs.m - blindingRowIdx - 1) * this->_rs.n),
             (const void *)&blindingRow[0], this->_l * sizeof(FieldT));
    }
  }

  /** Prover runs a round of the protocol
   * (there is only one round in total for the RSA Ceremony)
   *
   * @param round current round
   * @param earlyOnly if true, only populate and generate the commitment for the
   * early witness
   */
  void runRound(size_t round, bool earlyOnly) {
    /** Expand the transcript */
    this->_transcript.emplace_back(roundFSTranscript<FieldT>());

    /** Initialize secret sharing interface */
    ligero::SecretSharingInterface<FieldT> localSecretSharing(
        this->_publicData.modulusIdx, (size_t)this->_l, (size_t)this->_rs.k,
        (size_t)this->_rs.n, this->_computeSmall, this->_computeLarge);

    /** Have the "Verifier" pick randomness for all tests according to the
     * Fiat-Shamir transform
     */
    std::vector<FieldT> interleaved_r, interleaved_r_early,
        linearConstraintsInterBloc_r, quadraticConstraints_r,
        linearConstraintsIntraBloc_ri, linearConstraintsSigma_r,
        linearConstraintsStitchingIntraBloc_ri;

    /** Early witness
     * ========================================================================================
     */

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      this->_myTimer.begin(3, "Sharing early witness", 0);

      /** Pad the early witness with randomness */
      localSecretSharing.padMany(this->_w_early, this->_rs_early.m);

      /** Share the early witness */
      localSecretSharing.shareMany(this->_w_early, this->_rs_early.m);

      this->_myTimer.end(3, "Sharing early witness", 0);
    }

    /** Commit to the early witness */
    ProverColumnCS<FieldT> proverCCS_early(this->_mtParams_early,
                                           this->_w_early);

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      /** Merkel Tree Commitment */
      this->_myTimer.begin(3, "Commit Early Witness MT", 0);

      this->_transcript[round].proverCommitment_early =
          proverCCS_early.performInitialCommit(this->_w_early);

      this->_myTimer.end(3, "Commit Early Witness MT", 0);
    }

    /** End here if all we wanted was the early commitment */
    if (earlyOnly) return;

    /** Moving on to the main witness */
    /* =========================================================================
     */

    /** Adding blinding */
    addBlindingRows();

    this->_myTimer.begin(3, "Sharing witness", 0);

    /** Pad the extended witness with randomness */
    localSecretSharing.padMany(this->_w, this->_rs.m);

    /** Share the extended witness */
    localSecretSharing.shareMany(this->_w, this->_rs.m);

    this->_myTimer.end(3, "Sharing witness", 0);

    /** Merkel Tree Commitment */
    this->_myTimer.begin(3, "Commit MT", 0);

    ProverColumnCS<FieldT> proverCCS(this->_mtParams, this->_w);
    this->_transcript[round].proverCommitment =
        proverCCS.performInitialCommit(this->_w);

    this->_myTimer.end(3, "Commit MT", 0);

    this->_myTimer.begin(3, "Estimate Randomness", 0);
    bool hasTransformation = false;

    /** Estimate All Randomness */
    size_t sizeRandomness = this->_rs.m; /** LDT testing */

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      sizeRandomness +=
          this->_rs_early.m + this->_splitWitness.size(this->_publicData);
    }

    for (auto constraint : this->_constraints.constraints) {
      switch (constraint->type) {
        case linear:
          sizeRandomness++;
          break;

        case transformation:
          if (!hasTransformation) {
            hasTransformation = true;
            sizeRandomness += this->_rs.m * this->_l;
          }
          break;

        case quadratic:
          sizeRandomness++;
          break;

        default:
          throw std::runtime_error(
              "incorrect constraint type to obtain prover responses:" +
              std::to_string(constraint->type));
          break;
      }
    }

    this->_myTimer.end(3, "Estimate Randomness", 0);

    for (size_t iteration = 0; iteration < numberOfIterations; iteration++) {
      std::vector<FieldT> interleaved_rU, linearConstraintsInterBloc_rU,
          linearConstraintsSigma_rU, linearConstraintsIntraBloc_riU,
          linearConstraintsStitchingIntraBloc_riU, quadraticConstraints_rU;
      std::unordered_map<size_t, std::vector<FieldT>> linearCombinations;

      /** Update the state (= establish a new seed for Fiat-Shamir randomness)
       * and pick all randomness in a bulk */
      sizeRandomness += this->_constraintSystemRandomnessSize;

      this->_myTimer.begin(3, "Pick Randomness", 0);
      this->pickFiatShamirRandomness(sizeRandomness, this->_transcript);
      this->_myTimer.end(3, "Pick Randomness", 0);

      /** Rebuild the constraint system based on this round's randomness
       * We are locally scoping objects so as to invoke their destructor
       * immediately; they won't be needed here
       */
      {
        this->cleanup();

        builder<FieldT> bob(this->_l);
        vector<vector<FieldT>> extendedWitness2DVector;

        this->_myTimer.begin(3, "Express NP statement", 0);
        expressNPStatement<FieldT> expressStatement(
            this->_publicData, this->_sigmaProtocolPublicData,
            this->_secretData, bob, &extendedWitness2DVector, this->_l,
            this->_randomness);
        expressStatement.buildConstraintSet(this->_statement);
        this->_myTimer.end(3, "Express NP statement", 0);

        this->_myTimer.begin(3, "compile", 0);
        this->_constraints = bob.compile(
            extendedWitness2DVector, false); /** Main functionality to retain */
        this->_myTimer.end(3, "compile", 0);
      }

      /** LDT testing */
      this->_myTimer.begin(3, "LDT Testing Core Witness", 0);
      this->pickVerifierChallenges(round, TEST_INTERLEAVED, interleaved_r,
                                   this->_rs.m);
      this->interbloc(interleaved_rU, interleaved_r, this->_w);
      this->_myTimer.end(3, "LDT Testing Core Witness", 0);

      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        this->_myTimer.begin(3, "LDT Testing Early Commit", 0);
        this->pickVerifierChallenges(round, TEST_INTERLEAVED,
                                     interleaved_r_early, this->_rs_early.m);
        this->interbloc(interleaved_rU, interleaved_r_early, this->_w_early);
        this->_myTimer.end(3, "LDT Testing Early Commit", 0);
      }

      /** transformation constraints, sparse randomness */
      std::unordered_set<size_t> processed, processedStitched;
      std::vector<FieldT> localRandomness;

      /** Testing constraints */
      size_t target = 0;

      /** Intrabloc constraints only when there are transformations */
      this->_myTimer.begin(3, "Pick Randomness", 0);
      this->pickVerifierChallenges(round, TEST_LINEAR_INTRABLOC_CONSTRAINTS,
                                   localRandomness, this->_rs.m);
      this->_myTimer.end(3, "Pick Randomness", 0);

      if (hasTransformation) {
        linearConstraintsIntraBloc_ri.assign(this->_rs.n * this->_rs.m,
                                             FieldT(0));
      }

      /** Stitching constraints for all moduli */
      linearConstraintsStitchingIntraBloc_ri.assign(this->_rs.n * this->_rs.m,
                                                    FieldT(0));

      this->_myTimer.begin(3, "Linear Combinations for Constraints Testing", 0);

      for (auto constraint : this->_constraints.constraints) {
        switch (constraint->type) {
          case linear: {
            this->_myTimer.begin(3, "Linear Combination", 0);

            this->pickVerifierChallenges(
                round, TEST_LINEAR_INTERBLOC_CONSTRAINTS,
                linearConstraintsInterBloc_r, this->_rs.m);
            this->obtainProverResponses(
                round, this->_constraints.targets[target], constraint,
                TEST_LINEAR_INTERBLOC_CONSTRAINTS, linearConstraintsInterBloc_r,
                linearConstraintsInterBloc_rU);

            this->_myTimer.end(3, "Linear Combination", 0);
          } break;

          case transformation: {
            transformation_constraint<FieldT> transfCt =
                *static_cast<const transformation_constraint<FieldT> *>(
                    constraint);

            this->_myTimer.begin(3, "Calc Matrix", 0);
            if (transfCt.stitching == LinearCombinationsStitching::Regular) {
              this->calculateRTAmatrix(
                  linearConstraintsIntraBloc_ri, transfCt, localRandomness,
                  this->_constraints.targets[target], processed);
            } else {
              this->calculateRTAStitchingMatrix(
                  linearConstraintsStitchingIntraBloc_ri, transfCt,
                  localRandomness, this->_constraints.targets[target],
                  processedStitched, iteration);
            }
            this->_myTimer.end(3, "Calc Matrix", 0);
          } break;

          case quadratic: {
            this->_myTimer.begin(3, "Quadratic", 0);
            this->pickVerifierChallenges(round, TEST_QUADRATIC_CONSTRAINTS,
                                         quadraticConstraints_r, this->_rs.m);
            this->obtainProverResponses(
                round, this->_constraints.targets[target], constraint,
                TEST_QUADRATIC_CONSTRAINTS, quadraticConstraints_r,
                quadraticConstraints_rU);
            this->_myTimer.end(3, "Quadratic", 0);
          } break;

          default:
            throw std::runtime_error(
                "incorrect constraint type to obtain prover responses:" +
                std::to_string(constraint->type));
            break;
        }

        target++;
      }

      auto vectorInterBloc =
          [&](std::vector<FieldT> &linearConstraintsInterBloc_rU,
              const std::vector<size_t> &v_early,
              const std::vector<size_t> &v) -> void {
        for (size_t idx = 0; idx < v.size(); idx++) {
          this->splitWitnessInterbloc(
              linearConstraintsInterBloc_rU, v_early[idx], v[idx],
              FieldT(ligero::math::scaleToPrime(*this->_randomness++,
                                                this->_publicData.modulusIdx)),
              this->_w_early, this->_w);
        }
      };

      /** In the specific case of the full ceremony,
       * add inter-bloc constraints tying the variables
       * across the split witness
       */
      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        vectorInterBloc(linearConstraintsInterBloc_rU,
                        this->_splitWitness.si_early_BlocIdxs,
                        this->_splitWitness.si_BlocIdxs);
        vectorInterBloc(linearConstraintsInterBloc_rU,
                        this->_splitWitness.ei_early_BlocIdxs,
                        this->_splitWitness.ei_BlocIdxs);
        if (this->_publicData.modulusIdx < this->_publicData.p)
          vectorInterBloc(linearConstraintsInterBloc_rU,
                          this->_splitWitness.ri_early_BlocIdxs,
                          this->_splitWitness.ri_BlocIdxs);
      }

      this->_myTimer.end(3, "Linear Combinations for Constraints Testing", 0);

      /** Sharing & Committing */
      /* ================================================================== */
      this->_myTimer.begin(
          3, "Sharing + linear combination intrabloc randomness", 0);

      if (hasTransformation) {
        /** Perform linear combinations */
        linearConstraintsIntraBloc_riU.resize(this->_rs.n, FieldT(0));

        for (size_t row : processed) {
          /** Pad randomness */
          localSecretSharing.padIntrablocRandomness(
              &linearConstraintsIntraBloc_ri[this->_rs.n * row]);

          /** Secret Share Randomness First */
          localSecretSharing.share(
              &linearConstraintsIntraBloc_ri[this->_rs.n * row]);

          for (size_t col = 0; col < this->_rs.n; col++) {
            linearConstraintsIntraBloc_riU[col] +=
                this->_w[col + (row * this->_rs.n)] *
                linearConstraintsIntraBloc_ri[col + (row * this->_rs.n)];
          }
        }

        /** Adding Blinding */
        for (size_t col = 0; col < this->_rs.n; col++) {
          linearConstraintsIntraBloc_riU[col] +=
              this->_w[col +
                       ((this->_rs.m + offsetIntrablocBlinding[iteration]) *
                        this->_rs.n)];
        }

#ifndef NDEBUG
        /** Sanity check */
        std::vector<FieldT> check;
        for (size_t col = 0; col < this->_rs.n; col++) {
          check.push_back(this->_w[col + ((this->_rs.m +
                                           offsetIntrablocBlinding[iteration]) *
                                          this->_rs.n)]);
        }

        localSecretSharing.reconstruct(&check[0], false);

        FieldT sanityCheck(0);
        for (size_t col = 0; col < this->_l; col++) {
          sanityCheck += check[col];
        }

        assert(sanityCheck == FieldT(0));
#endif
      }

      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        /** Perform linear combinations */
        linearConstraintsStitchingIntraBloc_riU.resize(this->_rs.n, FieldT(0));

        for (size_t row : processedStitched) {
          /** Pad randomness */
          localSecretSharing.padIntrablocRandomness(
              &linearConstraintsStitchingIntraBloc_ri[this->_rs.n * row]);

          /** Secret Share Randomness First */
          localSecretSharing.share(
              &linearConstraintsStitchingIntraBloc_ri[this->_rs.n * row]);

          for (size_t col = 0; col < this->_rs.n; col++) {
            linearConstraintsStitchingIntraBloc_riU[col] +=
                this->_w[col + (row * this->_rs.n)] *
                linearConstraintsStitchingIntraBloc_ri[col +
                                                       (row * this->_rs.n)];
          }
        }

        /** Adding Blinding */
        for (size_t col = 0; col < this->_rs.n; col++) {
          linearConstraintsStitchingIntraBloc_riU[col] +=
              this->_w[col +
                       ((this->_rs.m + offsetStitchingBlinding[iteration]) *
                        this->_rs.n)];
        }
      }

      this->_myTimer.end(3, "Sharing + linear combination intrabloc randomness",
                         0);

      /** Prepare for a new protocol round */
      this->_round = round;

      /** Populate the transcript */
      linearCombinations.insert(std::pair<int, std::vector<FieldT>>(
          TEST_INTERLEAVED, interleaved_rU));
      linearCombinations.insert(std::pair<int, std::vector<FieldT>>(
          TEST_LINEAR_INTERBLOC_CONSTRAINTS, linearConstraintsInterBloc_rU));
      linearCombinations.insert(std::pair<int, std::vector<FieldT>>(
          TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS,
          linearConstraintsStitchingIntraBloc_riU));
      linearCombinations.insert(std::pair<int, std::vector<FieldT>>(
          TEST_LINEAR_INTRABLOC_CONSTRAINTS, linearConstraintsIntraBloc_riU));
      linearCombinations.insert(std::pair<int, std::vector<FieldT>>(
          TEST_QUADRATIC_CONSTRAINTS, quadraticConstraints_rU));
      this->_transcript[round].proverResponsesLinearCombinations.emplace_back(
          linearCombinations);

#ifndef NDEBUG
      /** Display Randomness */
      displayVector(",Prover,local, Verifier LDT Combinations Challenge:",
                    interleaved_r);
      displayVector(
          ",Prover,local, Verifier Linear Interbloc Constraints Combinations "
          "Challenge:",
          linearConstraintsInterBloc_r);
      if (hasTransformation)
        displayVector(
            ",Prover,local, Verifier Linear Intrabloc Constraints Combinations "
            "Challenge:",
            linearConstraintsIntraBloc_ri);
      displayVector(
          ",Prover,local, Verifier Quadratic Constraints Combinations "
          "Challenge:",
          quadraticConstraints_r);

      /** Display Linear Combinations */
      displayVector(",Prover,transcript, Prover LDT Response:", interleaved_rU);
      displayVector(
          ",Prover,transcript, Prover Interbloc Constraints Response:",
          linearConstraintsInterBloc_rU);
      displayVector(
          ",Prover,transcript, Prover Intrabloc Stitching Constraints "
          "Response:",
          linearConstraintsStitchingIntraBloc_riU);
      displayVector(
          ",Prover,transcript, Prover Sigma Protocol Constraints Response:",
          linearConstraintsSigma_rU);
      displayVector(
          ",Prover,transcript, Prover Intrabloc Constraints Response:",
          linearConstraintsIntraBloc_riU);
      displayVector(
          ",Prover,transcript, Prover Quadratic Constraints Response:",
          quadraticConstraints_rU);
#endif

      /** Free Randomness */
      free(this->_alloc);
    }

    /** Verifier Checking Certain Columns */
    /* ==================================================================== */

    /** Source Randomness */
    this->_myTimer.begin(3, "Pick Randomness to Open Up Columns", 0);
    this->pickFiatShamirRandomness(this->_t, this->_transcript);
    this->_myTimer.end(3, "Pick Randomness to Open Up Columns", 0);

    this->_myTimer.begin(3, "Decommit Columns", 0);

    /** Open columns */
    this->pickVerifierChallenges(
        round, OPEN_COLUMNS, this->_transcript[round].verifierColumnQueries);
    displayVector(",Prover,local, Opening Columns:",
                  this->_transcript[round].verifierColumnQueries);

    /** Free Randomness */
    free(this->_alloc);

    /** Test out the columns */
    for (size_t idx = 0;
         idx < this->_transcript[round].verifierColumnQueries.size(); idx++) {
      assert(this->_transcript[round].verifierColumnQueries[idx] <
             (size_t)this->_rs.n);
    }

    /** Add to transcript round */
    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      this->_transcript[round].proverDecommitment_early =
          proverCCS_early.performProverDecommit(
              this->_transcript[round].verifierColumnQueries);
    }

    this->_transcript[round].proverDecommitment =
        proverCCS.performProverDecommit(
            this->_transcript[round].verifierColumnQueries);

#ifndef NDEBUG
    displayVector(",Prover,transcript, Column 1 Content:",
                  this->_transcript[round].proverDecommitment.contents.front());
    displayVector(",Prover,transcript, Column q Content:",
                  this->_transcript[round].proverDecommitment.contents.back());

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      displayVector(
          ",Prover,transcript, Column 1 Content:",
          this->_transcript[round].proverDecommitment_early.contents.front());
      displayVector(
          ",Prover,transcript, Column q Content:",
          this->_transcript[round].proverDecommitment_early.contents.back());
    }
#endif

    this->_myTimer.end(3, "Decommit Columns", 0);
  }
};

template <typename FieldT>
class Verifier : public NonInteractiveArgument<FieldT> {
  const size_t _l;

 public:
  /** Constructor for the verifier */
  Verifier(Statements statement, InterleavedRSCode &irs, size_t t,
           size_t blockSize,
           const hash::HashIntegerSeedGeneratorF<FieldT> &generator,
           const hash::HashStringSeedGeneratorF &generatorString,
           ligero::MTParameters<FieldT> &mt, PublicData &pdata,
           SigmaProtocolPublicData &spdata, SecretData &sdata,
           void (*computeSmall)(uint64_t *, uint64_t const *),
           void (*computeLarge)(uint64_t *, uint64_t const *),
           std::string &hashedPublicData)
      : _l(blockSize),
        NonInteractiveArgument<FieldT>(
            statement, irs, t, blockSize, generator, generatorString, mt, pdata,
            spdata, sdata, computeSmall, computeLarge, hashedPublicData) {
    /** By definition, the verifier does not have access to the secret data */
    this->_secretData.isAvailable = false;
  };

  /** Run tests */
  void test(int type, vector<FieldT> &rU,
            ligero::SecretSharingInterface<FieldT> &secretSharing) {
    bool success = true;

    switch (type) {
      case TEST_INTERLEAVED:
        success = secretSharing.degreeTest(&rU[0]);
        break;
      case TEST_LINEAR_INTERBLOC_CONSTRAINTS:
        success = secretSharing.zeroTest(&rU[0]);
        break;
      case TEST_LINEAR_INTRABLOC_CONSTRAINTS:
        success = secretSharing.zeroSumTest(&rU[0]);
        break;
      case TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS:
        this->_stitching = secretSharing.sumReconstruct(&rU[0]).getValue();
        break;
      case TEST_QUADRATIC_CONSTRAINTS:
        success = secretSharing.zeroTest(&rU[0], true);
        break;
      default:
        throw std::runtime_error("incorrect test type:" + std::to_string(type));
        break;
    }

    if (success)
      DBG("Success:" << testDescription[type]);
    else
      throw failedTest(
          (std::string("test failed:") + testDescription[type]).c_str());
  }

  /** Compares linear combinations calculated by the verifier for the opened
   * columns with those provided by the prover for the entire witness
   * @param rU reference to linear combinations for the entire witness provided
   * by the prover
   * @param Q indices of the opened columns
   * @param c_rUq reference to the linear combinations for the opened columns
   */
  bool compare(const vector<FieldT> &rU, const vector<size_t> &Q,
               const vector<FieldT> &c_rUq) {
    bool success = true;
    for (size_t i = 0; i < Q.size(); i++) {
      if (!(rU[Q[i]] == c_rUq[i])) {
        return false;
      }
    }

    return success;
  }

  /** Checks consistency of the linear combinations provided by the prover with
   * the calculations of the verifier for the opened columns
   * @param r simulated verifier's challenges for the main witness
   * @param r_early simulated verifier's challenges for the early witness
   * @param rU reference to linear combinations for the entire witness provided
   * by the prover
   * @param Q indices of the opened columns
   * @param Uq reference to the matrix of opened columns for the main witness
   * @param Uq_early reference to the matrix of opened columns for the early
   * witness
   */
  void checkConsistency(const vector<FieldT> &r, const vector<FieldT> &r_early,
                        const vector<FieldT> &rU, const vector<size_t> &Q,
                        const std::vector<std::vector<FieldT>> &Uq,
                        const std::vector<std::vector<FieldT>> &Uq_early) {
    std::vector<FieldT> c_rUq;

    this->interbloc(c_rUq, r, Uq);
    this->interbloc(c_rUq, r_early, Uq_early);

    if (compare(rU, Q, c_rUq))
      DBG("Success: consistency " << testDescription[TEST_INTERLEAVED]);
    else
      throw failedTest((std::string("consistency check failed:") +
                        testDescription[TEST_INTERLEAVED])
                           .c_str());
  }

  /** verifier helper function updating the prover's response to the simulated
   * verifier's challenges for the opened columns
   * @param type the type of underlying constraint
   * @param c_rUq reference to the prover's response (linear combinations for
   * the opened columns as a vector)
   * @param r simulated verifier's challenges for the main witness
   * @param Q indices of the opened columns
   * @param Uq reference to the matrix of opened columns for the main witness
   * @param processed set of blocks involved in the linear constraint
   * @param iteration current iteration within the current round
   * witness
   */
  void linearCombinationPerConstraint(
      int type, std::vector<FieldT> &c_rUq, std::vector<FieldT> &r,
      const vector<size_t> &Q, const std::vector<std::vector<FieldT>> &Uq,
      std::unordered_set<size_t> &processed, uint8_t iteration) {
    /** Have The Secret Sharing Interface Ready */
    ligero::SecretSharingInterface<FieldT> localSecretSharing(
        this->_publicData.modulusIdx, (size_t)this->_l, (size_t)this->_rs.k,
        (size_t)this->_rs.n, this->_computeSmall, this->_computeLarge);

    for (auto idx : processed) {
      localSecretSharing.padIntrablocRandomness(&r[0] + (idx * this->_rs.n));
      localSecretSharing.share(&r[0] + (idx * this->_rs.n));
    }

    switch (type) {
      case TEST_LINEAR_INTRABLOC_CONSTRAINTS:
        this->intrabloc(c_rUq, r, Uq, Q, processed, iteration);
        break;
      case TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS:
        this->intrablocStitching(c_rUq, r, Uq, Q, processed, iteration);
        break;
      default:
        throw internalError("incorrect consistency check type:" +
                            testDescription[type]);
        break;
    }
  }

  /** verifier helper function updating the prover's response to the simulated
   * verifier's challenges for the opened columns for a specific constraint
   * @param type the type of underlying constraint
   * @param c_rUq reference to the prover's response (linear combinations for
   * the opened columns as a vector)
   * @param ct pointer to the contraint being handled
   * @param r simulated verifier's challenges for the main witness
   * @param Q indices of the opened columns
   * @param Uq reference to the matrix of opened columns for the main witness
   * @param processed set of blocks involved in the linear constraint
   * @param iteration current iteration within the current round
   * witness
   */
  void linearCombinationPerConstraint(
      int type, std::vector<FieldT> &c_rUq, constraint *ct, size_t target,
      const FieldT &r, const vector<size_t> &Q,
      const std::vector<std::vector<FieldT>> &Uq) {
    switch (type) {
      case TEST_QUADRATIC_CONSTRAINTS: {
        quadratic_constraint<FieldT> constraint =
            *static_cast<const quadratic_constraint<FieldT> *>(ct);
        this->interbloc(c_rUq, constraint.left_block, constraint.right_block,
                        constraint.scalar, target, r, Uq);
      } break;
      case TEST_LINEAR_INTERBLOC_CONSTRAINTS: {
        linear_constraint<FieldT> constraint =
            *static_cast<const linear_constraint<FieldT> *>(ct);
        this->interbloc(c_rUq, constraint.block, constraint.scalar, target, r,
                        Uq);
      } break;
      default:
        throw internalError("incorrect consistency check type:" +
                            testDescription[type]);
        break;
    }
  }

  /** Run a round of the protocol for the verifier
   * @param round the current round
   * @param earlyOnly if true, only generate and commit to the merkle tree for
   * the early witness
   */
  void runRound(size_t round, bool earlyOnly) {
    /** Have The Secret Sharing Interface Ready */
    ligero::SecretSharingInterface<FieldT> localSecretSharing(
        this->_publicData.modulusIdx, (size_t)this->_l, (size_t)this->_rs.k,
        (size_t)this->_rs.n, this->_computeSmall, this->_computeLarge);

    /** Prepare this transcript round */
    FullFSTranscript<FieldT> replicatingTranscript;

    if (round > 0) {
      replicatingTranscript.resize(round - 1);

      for (size_t i = 0; i < round; i++) {
        replicatingTranscript[i].proverCommitment.swap(
            this->_transcript[i].proverCommitment);
        replicatingTranscript[i].proverResponsesLinearCombinations.swap(
            this->_transcript[i].proverResponsesLinearCombinations);
        replicatingTranscript[i].verifierColumnQueries.swap(
            this->_transcript[i].verifierColumnQueries);

        replicatingTranscript[i].proverDecommitment.randomnessHashes.swap(
            this->_transcript[i].proverDecommitment.randomnessHashes);
        replicatingTranscript[i].proverDecommitment.auxiliaryHashes.swap(
            this->_transcript[i].proverDecommitment.auxiliaryHashes);
        replicatingTranscript[i].proverDecommitment.contentHashes.swap(
            this->_transcript[i].proverDecommitment.contentHashes);
        replicatingTranscript[i].proverDecommitment.contents.swap(
            this->_transcript[i].proverDecommitment.contents);
        replicatingTranscript[i].proverDecommitment.positions.swap(
            this->_transcript[i].proverDecommitment.positions);

        if ((this->_statement == Statements::RSACeremony) ||
            (this->_statement == Statements::RSACeremony_Rounds3to6)) {
          replicatingTranscript[i].proverCommitment_early.swap(
              this->_transcript[i].proverCommitment_early);

          replicatingTranscript[i]
              .proverDecommitment_early.randomnessHashes.swap(
                  this->_transcript[i]
                      .proverDecommitment_early.randomnessHashes);
          replicatingTranscript[i]
              .proverDecommitment_early.auxiliaryHashes.swap(
                  this->_transcript[i]
                      .proverDecommitment_early.auxiliaryHashes);
          replicatingTranscript[i].proverDecommitment_early.contentHashes.swap(
              this->_transcript[i].proverDecommitment_early.contentHashes);
          replicatingTranscript[i].proverDecommitment_early.contents.swap(
              this->_transcript[i].proverDecommitment_early.contents);
          replicatingTranscript[i].proverDecommitment_early.positions.swap(
              this->_transcript[i].proverDecommitment_early.positions);
        }
      }
    }

    roundFSTranscript<FieldT> current;
    replicatingTranscript.emplace_back(current);

    roundFSTranscript<FieldT> &transcriptRound = this->_transcript[round];

    /** Instantiate Tests */
    VerifierColumnCS<FieldT> verifierCCS_early(this->_mtParams_early);
    VerifierColumnCS<FieldT> verifierCCS(this->_mtParams);

    this->_myTimer.begin(3, "Estimate Randomness", 0);
    bool hasTransformation = false;

    /** Estimate All Randomness */
    size_t sizeRandomness = this->_rs.m; /** LDT testing */
    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      sizeRandomness +=
          this->_rs_early.m + this->_splitWitness.size(this->_publicData);
    }

    for (auto constraint : this->_constraints.constraints) {
      switch (constraint->type) {
        case linear:
          sizeRandomness++;
          break;

        case transformation:
          if (!hasTransformation) {
            hasTransformation = true;
            sizeRandomness += this->_rs.m * this->_l;
          }
          break;

        case quadratic:
          sizeRandomness++;
          break;

        default:
          throw std::runtime_error(
              "incorrect constraint type to obtain prover responses:" +
              std::to_string(constraint->type));
          break;
      }
    }

    this->_myTimer.end(3, "Estimate Randomness", 0);

    /** Open columns */
    /* =========================================================================
     */

    /** Replicate the state of the prover transcript */
    replicatingTranscript.back().proverCommitment.swap(
        transcriptRound.proverCommitment);
    replicatingTranscript.back().proverResponsesLinearCombinations.swap(
        transcriptRound.proverResponsesLinearCombinations);

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      replicatingTranscript.back().proverCommitment_early.swap(
          transcriptRound.proverCommitment_early);
    }

    /** Source Randomness */
    this->_myTimer.begin(3, "Pick Randomness to Open Up Columns", 0);
    this->pickFiatShamirRandomness(this->_t, replicatingTranscript);
    this->_myTimer.end(3, "Pick Randomness to Open Up Columns", 0);

    /** Pick columns */
    this->pickVerifierChallenges(
        round, OPEN_COLUMNS,
        replicatingTranscript.back().verifierColumnQueries);

    vector<size_t> &Q = replicatingTranscript.back().verifierColumnQueries;
    displayVector(",Verifier,local, Opening Columns:", Q);

    /** Free Randomness */
    free(this->_alloc);

    DBG("Verifying decommitment root for main witness.");

    /** Now verify the decommitment vs. initial commitment */
    if (!verifierCCS.verify(
            replicatingTranscript.back().proverCommitment,
            transcriptRound.proverDecommitment,
            replicatingTranscript.back().verifierColumnQueries)) {
      throw failedTest("openColumns|contents");
    }

    DBG("Verification successful for decommitment root for main witness.");
    DBG("Verifying decommitment root for early witness.");

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      if (!verifierCCS_early.verify(
              replicatingTranscript.back().proverCommitment_early,
              transcriptRound.proverDecommitment_early,
              replicatingTranscript.back().verifierColumnQueries)) {
        throw failedTest("openColumns|contents");
      }
    }

    DBG("Verification successful for decommitment root for early witness.");

    /** The transcript is already populated, see what we read in */
    std::vector<std::vector<FieldT>> &Uq =
        transcriptRound.proverDecommitment.contents;
    std::vector<std::vector<FieldT>> &Uq_early =
        transcriptRound.proverDecommitment_early.contents;

    /** Debug: Opened Column Content */
    displayVector(",Verifier,transcript, First Opened Column:", Uq[0]);

    /** Unwind the linear combinations */
    replicatingTranscript.back().proverResponsesLinearCombinations.swap(
        transcriptRound.proverResponsesLinearCombinations);

    /** Verifying Linear Combinations */
    /* =========================================================================
     */

    for (size_t iteration = 0; iteration < numberOfIterations; iteration++) {
      /** Replicate the state of the prover transcript */
      if (iteration > 0) {
        std::unordered_map<size_t, std::vector<FieldT>> a;

        replicatingTranscript.back()
            .proverResponsesLinearCombinations.push_back(a);
        replicatingTranscript.back()
            .proverResponsesLinearCombinations.back()
            .swap(this->_transcript[round]
                      .proverResponsesLinearCombinations[iteration - 1]);
      }

      /** Pick Randomness for this iteration */
      sizeRandomness += this->_constraintSystemRandomnessSize;

      this->_myTimer.begin(3, "Pick Randomness", 0);
      this->pickFiatShamirRandomness(sizeRandomness, replicatingTranscript);
      this->_myTimer.end(3, "Pick Randomness", 0);

      /** Rebuild the constraint system based on this round's randomness
       * We are locally scoping objects so as to invoke their destructor
       * immediately; they won't be needed here
       */
      {
        this->cleanup();

        builder<FieldT> bob(this->_l);
        vector<vector<FieldT>> extendedWitness2DVector;

        this->_myTimer.begin(3, "Express NP statement", 0);
        expressNPStatement<FieldT> expressStatement(
            this->_publicData, this->_sigmaProtocolPublicData,
            this->_secretData, bob, &extendedWitness2DVector, this->_l,
            this->_randomness);

        expressStatement.buildConstraintSet(this->_statement);
        this->_myTimer.end(3, "Express NP statement", 0);

        this->_myTimer.begin(3, "compile", 0);
        this->_constraints = bob.compile(
            extendedWitness2DVector, false); /** Main functionality to retain */
        this->_myTimer.end(3, "compile", 0);
      }

      /** Have the "Verifier" pick randomness for all tests according to the
       * Fiat-Shamir transform
       */
      std::vector<FieldT> interleaved_r, interleaved_r_early,
          linearConstraintsInterBloc_r, quadraticConstraints_r,
          linearConstraintsSigma_r, linearConstraintsIntraBloc_ri,
          linearConstraintsStitchingIntraBloc_ri;

      /** Read provers' responses from the transcript */
      std::vector<FieldT> &interleaved_rU =
          transcriptRound
              .proverResponsesLinearCombinations[iteration][TEST_INTERLEAVED];
      std::vector<FieldT> &linearConstraintsInterBloc_rU =
          transcriptRound.proverResponsesLinearCombinations
              [iteration][TEST_LINEAR_INTERBLOC_CONSTRAINTS];
      std::vector<FieldT> &linearConstraintsStitchingIntraBloc_riU =
          transcriptRound.proverResponsesLinearCombinations
              [iteration][TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS];
      std::vector<FieldT> &linearConstraintsIntraBloc_riU =
          transcriptRound.proverResponsesLinearCombinations
              [iteration][TEST_LINEAR_INTRABLOC_CONSTRAINTS];
      std::vector<FieldT> &quadraticConstraints_rU =
          transcriptRound
              .proverResponsesLinearCombinations[iteration]
                                                [TEST_QUADRATIC_CONSTRAINTS];

      /** Basic Sanity Checks */
      assert((interleaved_rU.size() == this->_rs.n) ||
             (interleaved_rU.size() == 0));
      assert((linearConstraintsInterBloc_rU.size() == this->_rs.n) ||
             (linearConstraintsInterBloc_rU.size() == 0));
      assert((linearConstraintsIntraBloc_riU.size() == this->_rs.n) ||
             (linearConstraintsIntraBloc_riU.size() == 0));
      assert((linearConstraintsStitchingIntraBloc_riU.size() == this->_rs.n) ||
             (linearConstraintsStitchingIntraBloc_riU.size() == 0));
      assert((quadraticConstraints_rU.size() == this->_rs.n) ||
             (quadraticConstraints_rU.size() == 0));

      /** Pick Randomness for Low Degree Testing */
      this->pickVerifierChallenges(round, TEST_INTERLEAVED, interleaved_r,
                                   this->_rs.m);

      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        this->pickVerifierChallenges(round, TEST_INTERLEAVED,
                                     interleaved_r_early, this->_rs_early.m);
      }

      /** Low Degree Testing Debug */
      displayVector(",Verifier,local, Verifier LDT Combinations Challenge:",
                    interleaved_r);
      displayVector(",Verifier,transcript, Prover LDT Response:",
                    interleaved_rU);

      /** Low Degree Test */
      test(TEST_INTERLEAVED, interleaved_rU, localSecretSharing);

      /** Testing constraints */
      size_t target = 0;
      std::unordered_set<size_t> processed, processedStitched;
      std::vector<FieldT> localRandomness;
      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        linearConstraintsStitchingIntraBloc_ri.assign(this->_rs.n * this->_rs.m,
                                                      FieldT(0));
      }

      /** Intrabloc constraints */
      this->_myTimer.begin(3, "Pick Randomness", 0);
      this->pickVerifierChallenges(round, TEST_LINEAR_INTRABLOC_CONSTRAINTS,
                                   localRandomness, this->_rs.m);
      this->_myTimer.end(3, "Pick Randomness", 0);

      if (hasTransformation) {
        linearConstraintsIntraBloc_ri.assign(this->_rs.n * this->_rs.m,
                                             FieldT(0));
      }

      /** Generating randomness for testing constraints */
      for (auto constraint : this->_constraints.constraints) {
        switch (constraint->type) {
          case linear:
            this->pickVerifierChallenges(
                round, TEST_LINEAR_INTERBLOC_CONSTRAINTS,
                linearConstraintsInterBloc_r, this->_rs.m);
            break;

          case transformation: {
            transformation_constraint<FieldT> transfCt =
                *static_cast<const transformation_constraint<FieldT> *>(
                    constraint);
            if (transfCt.stitching == LinearCombinationsStitching::Regular) {
              this->calculateRTAmatrix(
                  linearConstraintsIntraBloc_ri, transfCt, localRandomness,
                  this->_constraints.targets[target], processed);
            } else {
              this->calculateRTAStitchingMatrix(
                  linearConstraintsStitchingIntraBloc_ri, transfCt,
                  localRandomness, this->_constraints.targets[target],
                  processedStitched, iteration, true);
            }

          } break;

          case quadratic:
            this->pickVerifierChallenges(round, TEST_QUADRATIC_CONSTRAINTS,
                                         quadraticConstraints_r, this->_rs.m);
            break;

          default:
            throw std::runtime_error("incorrect constraint type for test:" +
                                     std::to_string(constraint->type));
            break;
        }

        target++;
      }

      /** Testing */

      /** Aggregate test for each constraint */
      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        test(TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS,
             linearConstraintsStitchingIntraBloc_riU, localSecretSharing);
      }

      if (linearConstraintsInterBloc_r.size() > 0) {
        displayVector(
            ",Verifier,local, Verifier Linear Interbloc Constraints "
            "Combinations Challenge:",
            linearConstraintsInterBloc_r);
        displayVector(
            ",Verifier,transcript, Prover Interbloc Constraints Response:",
            linearConstraintsInterBloc_rU);
        test(TEST_LINEAR_INTERBLOC_CONSTRAINTS, linearConstraintsInterBloc_rU,
             localSecretSharing);
      }

      if (hasTransformation) {
        displayVector(
            ",Verifier,local, Verifier Linear Interbloc Constraints "
            "Combinations Challenge:",
            linearConstraintsIntraBloc_ri);
        displayVector(
            ",Verifier,transcript, Prover Interbloc Constraints Response:",
            linearConstraintsIntraBloc_riU);
        test(TEST_LINEAR_INTRABLOC_CONSTRAINTS, linearConstraintsIntraBloc_riU,
             localSecretSharing);
      }

      if (quadraticConstraints_r.size() > 0) {
        displayVector(
            ",Verifier,local, Verifier Quadratic Constraints Combinations "
            "Challenge:",
            quadraticConstraints_r);
        displayVector(
            ",Verifier,transcript, Prover Quadratic Constraints Response:",
            quadraticConstraints_rU);
        test(TEST_QUADRATIC_CONSTRAINTS, quadraticConstraints_rU,
             localSecretSharing);
      }

      /** Checking consistency with decommitment */
      checkConsistency(interleaved_r, interleaved_r_early, interleaved_rU, Q,
                       Uq, Uq_early);

      /** Read provers' responses from the transcript */
      std::vector<FieldT> linearConstraintsInterBloc_c_rUq;
      std::vector<FieldT> linearConstraintsInterBlocSigmaProtocol_c_rUq;
      std::vector<FieldT> linearConstraintsIntraBloc_c_riUq;
      std::vector<FieldT> linearConstraintsStitchingIntraBloc_c_riUq;
      std::vector<FieldT> quadraticConstraints_c_rUq;

      size_t idxLinear = 0, idxTransformation = 0, idxQuadratic = 0;
      target = 0;

      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        linearCombinationPerConstraint(
            TEST_LINEAR_INTRABLOC_STITCHING_CONSTRAINTS,
            linearConstraintsStitchingIntraBloc_c_riUq,
            linearConstraintsStitchingIntraBloc_ri, Q, Uq, processedStitched,
            iteration);
      }

      if (hasTransformation) {
        linearCombinationPerConstraint(TEST_LINEAR_INTRABLOC_CONSTRAINTS,
                                       linearConstraintsIntraBloc_c_riUq,
                                       linearConstraintsIntraBloc_ri, Q, Uq,
                                       processed, iteration);
      }

      for (auto &constraint : this->_constraints.constraints) {
        switch (constraint->type) {
          case linear:
            linearCombinationPerConstraint(
                TEST_LINEAR_INTERBLOC_CONSTRAINTS,
                linearConstraintsInterBloc_c_rUq, constraint,
                this->_constraints.targets[target],
                linearConstraintsInterBloc_r[idxLinear++], Q, Uq);
            break;

          case transformation:
            break;

          case quadratic:
            linearCombinationPerConstraint(
                TEST_QUADRATIC_CONSTRAINTS, quadraticConstraints_c_rUq,
                constraint, this->_constraints.targets[target],
                quadraticConstraints_r[idxQuadratic++], Q, Uq);
            break;

          default:
            throw std::runtime_error(
                "incorrect constraint type for checking consistency:" +
                std::to_string(constraint->type));
            break;
        }

        target++;
      }

      /** Split Witness Handling */
      auto vectorInterBloc =
          [&](std::vector<FieldT> &linearConstraintsInterBloc_c_rUq,
              const std::vector<size_t> &v_early,
              const std::vector<size_t> &v) -> void {
        for (size_t idx = 0; idx < v.size(); idx++) {
          this->splitWitnessInterbloc(
              linearConstraintsInterBloc_c_rUq, v_early[idx], v[idx],
              FieldT(ligero::math::scaleToPrime(*this->_randomness++,
                                                this->_publicData.modulusIdx)),
              Uq_early, Uq);
        }
      };

      /** In the specific case of the full ceremony,
       * add inter-bloc constraints tying the variables
       * across the split witness
       */
      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        if (linearConstraintsInterBloc_c_rUq.size() == 0) {
          linearConstraintsInterBloc_c_rUq.resize(Uq.size(), FieldT(0));
        }

        vectorInterBloc(linearConstraintsInterBloc_c_rUq,
                        this->_splitWitness.si_early_BlocIdxs,
                        this->_splitWitness.si_BlocIdxs);
        vectorInterBloc(linearConstraintsInterBloc_c_rUq,
                        this->_splitWitness.ei_early_BlocIdxs,
                        this->_splitWitness.ei_BlocIdxs);
        if (this->_publicData.modulusIdx < this->_publicData.p)
          vectorInterBloc(linearConstraintsInterBloc_c_rUq,
                          this->_splitWitness.ri_early_BlocIdxs,
                          this->_splitWitness.ri_BlocIdxs);
      }

#ifndef NDEBUG
      displayVector(
          ",Verifier,transcript, Verifier Linear Interbloc Constraints "
          "Combinations Challenge Open Col:",
          linearConstraintsInterBloc_c_rUq);
      displayVector(
          ",Verifier,transcript, Verifier Linear Sigma Protocol Constraints "
          "Combinations Challenge Open Col:",
          linearConstraintsInterBlocSigmaProtocol_c_rUq);
      displayVector(
          ",Verifier,transcript, Verifier Linear Intrabloc Constraints "
          "Combinations Challenge Open Col:",
          linearConstraintsIntraBloc_c_riUq);
      displayVector(
          ",Verifier,transcript, Verifier Linear Intrabloc Constraints "
          "Combinations Challenge Open Col:",
          linearConstraintsStitchingIntraBloc_c_riUq);
      displayVector(
          ",Verifier,transcript, Verifier Quadratic Constraints Combinations "
          "Challenge Open Col:",
          quadraticConstraints_c_rUq);
#endif

      /** Aggregate test for each constraint */
      if (linearConstraintsInterBloc_r.size() > 0) {
        if (compare(linearConstraintsInterBloc_rU, Q,
                    linearConstraintsInterBloc_c_rUq))
          DBG("Success: consistency linearConstraintsInterBloc");
        else
          throw failedTest(
              (std::string(
                   "consistency check failed: linearConstraintsInterBloc"))
                  .c_str());
      }

      if (linearConstraintsIntraBloc_ri.size() > 0) {
        if (compare(linearConstraintsIntraBloc_riU, Q,
                    linearConstraintsIntraBloc_c_riUq))
          DBG("Success: consistency linearConstraintsIntraBloc");
        else
          throw failedTest(
              (std::string(
                   "consistency check failed: linearConstraintsIntraBloc"))
                  .c_str());
      }

      if (quadraticConstraints_r.size() > 0) {
        if (compare(quadraticConstraints_rU, Q, quadraticConstraints_c_rUq))
          DBG("Success: consistency quadraticConstraints");
        else
          throw failedTest(
              (std::string("consistency check failed: quadraticConstraints"))
                  .c_str());
      }

      if ((this->_statement == Statements::RSACeremony) ||
          (this->_statement == Statements::RSACeremony_Rounds3to6)) {
        if (compare(linearConstraintsStitchingIntraBloc_riU, Q,
                    linearConstraintsStitchingIntraBloc_c_riUq))
          DBG("Success: consistency linearConstraintsStitchingIntraBloc");
        else
          throw failedTest((std::string("consistency check failed: "
                                        "linearConstraintsStitchingIntraBloc"))
                               .c_str());
      }

      /** Free Randomness */
      free(this->_alloc);
    }
  }

  /** Verify A Zero-Knowledge Proof
   * @param t full transcript provided by the prover
   */
  void verifyProof(FullFSTranscript<FieldT> &t) {
    this->_myTimer.initialize_timer();
    this->_transcript = t;

    /** Initialize Fpp */
    /** =========================== */

    auto assessBlocIdx = [&](auto &builder,
                             std::vector<size_t> idxs) -> std::vector<size_t> {
      std::vector<size_t> ret;
      for (auto idx : idxs) {
        ret.emplace_back(builder.proof_values_location[idx]);
      }

      return ret;
    };

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      /** Build a first witness/constraint system for Early Commitment */
      /** ================================================================== */

      this->_rs_early = this->_rs;
      this->_myTimer.begin(
          2, "Build Extended Witness & Constraint System For Early Commitment",
          0);

      builder<FieldT> helper(this->_l);
      vector<vector<FieldT>> earlyCommitments2DVector;

      expressNPStatement<FieldT> expressStatementEarlyCommitment(
          this->_publicData, this->_sigmaProtocolPublicData, this->_secretData,
          helper, &earlyCommitments2DVector, this->_l, this->_randomness);
      expressStatementEarlyCommitment.buildConstraintSet(
          Statements::RSACeremony_EarlyCommitments);

      this->_myTimer.begin(3, "compile", 0);
      this->_constraints_early = helper.compile(earlyCommitments2DVector);
      this->_myTimer.end(3, "compile", 0);

      this->_myTimer.end(
          2, "Build Extended Witness & Constraint System For Early Commitment",
          0);

      this->_splitWitness.si_early_BlocIdxs = assessBlocIdx(
          helper, expressStatementEarlyCommitment._siCoefsBlocIds);
      this->_splitWitness.ei_early_BlocIdxs = assessBlocIdx(
          helper, expressStatementEarlyCommitment._eiCoefsBlocIds);
      if (this->_publicData.modulusIdx < this->_publicData.p)
        this->_splitWitness.ri_early_BlocIdxs = assessBlocIdx(
            helper, {expressStatementEarlyCommitment._riCoefsBlocId});

      this->_mtParams_early.leavesEntropy = this->_rs_early.m =
          earlyCommitments2DVector.size() +
          this->_constraints_early.constant_blocks + 1;
    }

    /** Building the constraint system and populating the extended witness */
    /** ================================================================== */

    builder<FieldT> bob(this->_l);
    vector<vector<FieldT>> extendedWitness2DVector;

    expressNPStatement<FieldT> expressStatement(
        this->_publicData, this->_sigmaProtocolPublicData, this->_secretData,
        bob, &extendedWitness2DVector, this->_l, this->_randomness);
    expressStatement.buildConstraintSet(this->_statement, false);

    this->_constraints = bob.compile(
        extendedWitness2DVector); /** Main functionality to retain */
    DBG(this->_constraints.constraints.size());

    /** Dimensioning the circuit accordingly */
    this->_mtParams.leavesEntropy = this->_rs.m =
        extendedWitness2DVector.size() + this->_constraints.constant_blocks +
        nbBlindingRows;

    if ((this->_statement == Statements::RSACeremony) ||
        (this->_statement == Statements::RSACeremony_Rounds3to6)) {
      this->_splitWitness.si_BlocIdxs =
          assessBlocIdx(bob, expressStatement._siBlocs);
      this->_splitWitness.ei_BlocIdxs =
          assessBlocIdx(bob, expressStatement._eiBlocs);
      this->_splitWitness.ri_BlocIdxs =
          assessBlocIdx(bob, {expressStatement._riCoefsBlocId});
    }

    /** Here checking that we got the referencing right */
    this->_uxBlocks = assessBlocIdx(bob, expressStatement._uxBlocs);
    this->_uzBlocks = assessBlocIdx(bob, expressStatement._uzBlocs);

    this->_siBlocks = assessBlocIdx(bob, expressStatement._siBlocs);
    this->_eiBlocks = assessBlocIdx(bob, expressStatement._eiBlocs);

    /** Storing for Cross-Moduli Integrity Checking */
    this->_varBlocks.emplace_back(this->_siBlocks);
    this->_varBlocks.emplace_back(this->_eiBlocks);
    this->_varBlocks.emplace_back(this->_uxBlocks);
    this->_varBlocks.emplace_back(this->_uzBlocks);

    /** Generating and sending the proof */
    /** ================================ */

    /** Dimension Randomness for the Constraint System */
    this->_constraintSystemRandomnessSize = determineCSrandomnessSize(
        this->_statement, this->_publicData.ptNumberEvaluation);

    DBG("Processing New Protocol Round.");

    this->_myTimer.begin(2, "Verifying the Proof", 0);
    this->runProtocol(numberOfRounds);
    this->_myTimer.end(2, "Verifying the Proof", 0);
  }
};

}  // namespace zksnark
}  // namespace ligero
