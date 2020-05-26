#pragma once
#include <vector>
#include "FiniteFields.hpp"
#include "protocol/Constraint.hpp"

/** Defines the constraint system for the Ligero-style ZKP
 */
template <typename FieldT>
class ConstraintSystem {
 public:
  ConstraintSystem(const size_t bs, const size_t pb);
  ConstraintSystem();
  std::vector<constraint *> constraints;
  std::vector<FieldT> constants;
  std::vector<FieldT> constants_shares;
  std::vector<size_t> targets;
  std::vector<size_t> special;
  size_t block_size;
  size_t proof_blocks;
  size_t constant_blocks;
};

/** Defines the constraint system for the Ligero-style ZKP
 * @param bs block size for the proof system (this is the size of the private
 * blocks)
 * @param pb number of blocks in the witness
 */
template <typename FieldT>
ConstraintSystem<FieldT>::ConstraintSystem(const size_t bs, const size_t pb)
    : block_size(bs), proof_blocks(pb) {}

/** Default destructor for the constraint system
 */
template <typename FieldT>
ConstraintSystem<FieldT>::ConstraintSystem() {}