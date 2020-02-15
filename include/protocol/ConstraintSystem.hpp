#pragma once
#include <vector>
#include "protocol/Constraint.hpp"
#include "FiniteFields.hpp"

template<typename FieldT>
class ConstraintSystem {
public:
	ConstraintSystem(const size_t bs, const size_t pb);
	ConstraintSystem();
	std::vector<constraint*> constraints;
	std::vector<FieldT> constants;
	std::vector<FieldT> constants_shares;
	std::vector<size_t> targets;
	std::vector<size_t> special;
	size_t block_size;
	size_t proof_blocks;
	size_t constant_blocks;
};

template<typename FieldT>
ConstraintSystem<FieldT>::ConstraintSystem(const size_t bs, const size_t pb):
	block_size(bs),
	proof_blocks(pb) { 
}

template<typename FieldT>
ConstraintSystem<FieldT>::ConstraintSystem() {}
