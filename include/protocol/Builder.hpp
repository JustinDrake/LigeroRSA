#pragma once
#include <vector>

#include "protocol/ConstraintSystem.hpp"
#include "protocol/WeightedUnionFind.hpp"

#include "FiniteFields.hpp"

template <typename FieldT>
class builder {
 public:
  std::vector<size_t> proof_values_location;
  std::vector<constraint*> constraints;
  size_t QUAD_COUNTER = 0;
  const size_t block_size;

  builder(size_t block_sz);

  size_t add(constraint* c, bool special = false);
  void ensure_equality(size_t a, size_t b);

  /** destructively compile the circuit, reducing blocks */
  ConstraintSystem<FieldT> compile();
  ConstraintSystem<FieldT> compile(std::vector<std::vector<FieldT>>& values,
                                   bool populate = true);
  ConstraintSystem<FieldT> build_system(weightedUnionFind& representative,
                                        const std::vector<bool>& removable);

 private:
  std::vector<size_t> special_blocks;
  std::vector<std::vector<size_t>> equality_graph;
  void determine_removable(std::vector<bool>& removable) const;
  void rewire_nonremovable(const std::vector<bool>& removable);
  void rewire_inverse(const size_t idx, const std::vector<bool>& removable);
  void rewire_linear(const size_t idx, const std::vector<bool>& removable);
  void rewire_quadratic(const size_t idx, const std::vector<bool>& removable);
  void rewire_transformation(const size_t idx,
                             const std::vector<bool>& removable);
  std::vector<bool> find_and_rewire();
  weightedUnionFind calc_equivalence_classes();

  std::vector<std::vector<FieldT>> compute_proof_values(
      const std::vector<std::vector<FieldT>>& variable_values,
      const std::vector<bool>& removable, weightedUnionFind& representative);
};

using std::cout;
using std::endl;

using std::array;
using std::make_pair;
using std::pair;
using std::tuple;
using std::vector;

template <typename FieldT>
builder<FieldT>::builder(size_t block_sz) : block_size(block_sz) {}

template <typename FieldT>
size_t builder<FieldT>::add(constraint* c, bool special) {
  constraints.emplace_back(c);
  equality_graph.push_back(vector<size_t>());

  if (special) {
    special_blocks.push_back(constraints.size() - 1);
  }
  return constraints.size() - 1;
}

/**
 * ==============================================================================================
 */
/** classify variable 'type' for use with ensure_equality */
/** see which takes precedence for staying */
/** var = variable */
/** res = reserved */
/** con = constant */
/** otr = other */
/** NEW = save both, check later that they subtract to 0 */
/** LOW = choose one with lowest index */
/** */
/**       var res con otr */
/**      +--------------- */
/** 0 var|NEW NEW con var */
/** 1 res|NEW NEW NEW res */
/** 2 con|con NEW NEW con */
/** 3 otr|var res con LOW */
/**
 * ==============================================================================================
 */

/** populates the equality graph to reflect that blocks indexed by a and b are
 * equal
 * @param a first block
 * @param b second block
 */
template <typename FieldT>
void builder<FieldT>::ensure_equality(size_t a, size_t b) {
  assert(a < constraints.size());
  assert(b < constraints.size());
  equality_graph[a].push_back(b);
  equality_graph[b].push_back(a);
}

/** checks whether elem is in a sorted vector v
 * @param elem field element
 * @param v sorted vector of field elements
 * @return true if the elem is in the vector
 */
template <class T>
inline bool in_sorted_vector(const vector<T>& v, T elem) {
  auto it = lower_bound(v.begin(), v.end(), elem);
  if (it == v.end()) return false;
  return *it == elem;
}

/**
 * ==============================================================================================
 */
/** >>> STAGE 1 : */
/** Sweep through blocks and determine which constraints can be deleted */
/** Can be deleted iff */
/** 1) not in a nontrivial equality class */
/** 2) is nonterminal */
/** 3) meets the below table. format is TYPE : (followed by only) */
/** LINEAR : INVERSE, LINEAR, QUADRATIC, TRANSFORM */
/** TRANSFORM : LINEAR, TRANSFORM <-- note that a following linear is 'promoted'
 */
/** to transform We figure this out, and 'skip' unnecessary constraints by */
/** changing their surroundings */
/** >>> STAGE 2 : actually rewire these constraints */
/** >>> STAGE 3 : Sweep through blocks and calculate values for only blocks
 * we're */
/** keepin */
/** >>> STAGE 4 : Determine which block(s) of an equality class will be kept, */
/** 				 add new constraints if necessary, and change */
/** constraints which refer to such blocks */
/**
 * ==============================================================================================
 */

/** +----------------------+ */
/** | BEGINNING OF STAGE 1 | */
/** +----------------------+ */

/** determines which blocks are removable and which are not
 * @param removable vector of booleans to be populated according to whether
 * blocks are removable
 */
template <typename FieldT>
void builder<FieldT>::determine_removable(vector<bool>& removable) const {
  removable = vector<bool>(constraints.size(), true);

  /** is the block produced by this constraint relied upon afterwards? */
  vector<bool> terminal(constraints.size(), true);

  const size_t csz = constraints.size();
  for (size_t cons_idx = 0; cons_idx < csz; ++cons_idx) {
    /** can't remove things in a nontrivial equivalence class */
    /** or if they are special */
    if (!equality_graph[cons_idx].empty() ||
        !in_sorted_vector(special_blocks, cons_idx)) {
      removable[cons_idx] = false;
    }

    removable[cons_idx] = true;
    const constraint* cons = constraints[cons_idx];

    switch (cons->type) {
      /** inverse constraint makes transformation blocks unremovable */
      case ConstraintType::inverse: {
        removable[cons_idx] = false;
        const inverse_constraint<FieldT>* i_cons =
            static_cast<const inverse_constraint<FieldT>*>(cons);
        for (size_t idx : i_cons->block) {
          if (constraints[idx]->type == ConstraintType::transformation) {
            removable[idx] = false;
          }
          terminal[idx] = false;
        }
        break;
      }

      /** linear constraint just marks blocks as nonterminal */
      case ConstraintType::linear: {
        const linear_constraint<FieldT>* l_cons =
            static_cast<const linear_constraint<FieldT>*>(cons);
        for (const size_t idx : l_cons->block) {
          terminal[idx] = false;
        }
        break;
      }

      /** quadratic constraint makes transformation blocks unremovable */
      case ConstraintType::quadratic: {
        removable[cons_idx] = false;
        const quadratic_constraint<FieldT>* q_cons =
            static_cast<const quadratic_constraint<FieldT>*>(cons);
        for (size_t idx : q_cons->left_block) {
          if (constraints[idx]->type == ConstraintType::transformation) {
            removable[idx] = false;
          }
          terminal[idx] = false;
        }
        for (size_t idx : q_cons->right_block) {
          if (constraints[idx]->type == ConstraintType::transformation)
            removable[idx] = false;
          terminal[idx] = false;
        }
        break;
      }

      /** transformation constraint just marks blocks as nonterminal */
      case ConstraintType::transformation: {
        const transformation_constraint<FieldT>* t_cons =
            static_cast<const transformation_constraint<FieldT>*>(cons);
        for (size_t idx : t_cons->block) {
          terminal[idx] = false;
        }
        break;
      }

      /** variables and constants are not removable */
      default: {
        removable[cons_idx] = false;
        break;
      }
    }
  }

  /** final sweep to set terminal blocks as non-removable */
  for (size_t cons_idx = 0; cons_idx < constraints.size(); ++cons_idx)
    removable[cons_idx] = removable[cons_idx] & !terminal[cons_idx];
}

/** +----------------------+ */
/** | BEGINNING OF STAGE 2 | */
/** +----------------------+ */

/** rewire a inverse constraint depending on what is removable
 * @param idx index of the constraint
 * @param removable reference to the vector of booleans indicating removability
 * of the blocks
 */
template <typename FieldT>
void builder<FieldT>::rewire_inverse(size_t idx,
                                     const vector<bool>& removable) {
  typedef pair<size_t, FieldT> data;
  vector<data> new_data;
  inverse_constraint<FieldT>* const cons =
      static_cast<inverse_constraint<FieldT>*>(constraints[idx]);

  for (size_t pos = 0; pos < cons->block.size(); ++pos) {
    /** remove linear constraints */
    const size_t old_block = cons->block[pos];
    const FieldT old_scalar = cons->scalar[pos];
    const constraint* const sub_cons = constraints[old_block];

    if (removable[old_block] && sub_cons->type == ConstraintType::linear) {
      const linear_constraint<FieldT>* const sub_lcons =
          static_cast<const linear_constraint<FieldT>*>(sub_cons);

      /** forward data from linear constraint, multiplied by scalar */
      for (size_t sub_pos = 0; sub_pos < sub_lcons->block.size(); ++sub_pos) {
        const size_t new_block = sub_lcons->block[sub_pos];
        const FieldT new_scalar = sub_lcons->scalar[sub_pos];

        new_data.emplace_back(new_block, new_scalar * old_scalar);
      }

    } else {
      new_data.emplace_back(old_block, old_scalar);
    }
  }

  /** sort by block so we can reduce redundant entries when combining */
  std::sort(
      new_data.begin(), new_data.end(),
      [](const data& a, const data& b) -> bool { return a.first < b.first; });

  /** push new data back into constraint */
  cons->block.clear();
  cons->scalar.clear();

  for (const data& entry : new_data) {
    if (cons->block.empty() || cons->block.back() != entry.first) {
      cons->block.emplace_back(entry.first);
      cons->scalar.emplace_back(entry.second);
    } else {
      /** combine scalar with old data */
      cons->scalar.back() += entry.second;
    }
    /** remove zero scalars */
    if (!cons->block.empty() && cons->scalar.back() == FieldT(0)) {
      cons->block.pop_back();
      cons->scalar.pop_back();
    }
  }
}

/** Converts an interblock linear constaint into a transformation
 * @param cons interblock constraint to be converted
 * @param block_size size of blocks in the witness
 * @return transformation constraint generated
 */
template <typename FieldT>
transformation_constraint<FieldT>* convert_linear_to_transform(
    const linear_constraint<FieldT>* const cons, const size_t block_size) {
  transformation_constraint<FieldT>* const t_cons =
      new transformation_constraint<FieldT>();

  for (size_t pos = 0; pos < cons->block.size(); ++pos) {
    for (size_t loc = 0; loc < block_size; ++loc) {
      /** doesn't convert scalars */
      t_cons->emplace(cons->block[pos], loc, loc, cons->scalar[pos]);
    }
  }
  delete cons;
  return t_cons;
}

/** rewires a quadratic constraint depending on what is removable
 * @param idx index for the constraint
 * @param removable reference to the removability vector
 */
template <typename FieldT>
void builder<FieldT>::rewire_quadratic(size_t idx,
                                       const vector<bool>& removable) {
  typedef tuple<size_t, size_t, FieldT> data;
  vector<data> new_data;

  quadratic_constraint<FieldT>* const cons =
      static_cast<quadratic_constraint<FieldT>*>(constraints[idx]);

  for (size_t pos = 0; pos < cons->left_block.size(); ++pos) {
    /** remove linear constraints */
    size_t old_left_block = cons->left_block[pos];
    size_t old_right_block = cons->right_block[pos];
    const FieldT old_scalar = cons->scalar[pos];

    const constraint* l_cons = constraints[old_left_block];
    const constraint* r_cons = constraints[old_right_block];

    bool remove_left =
        removable[old_left_block] && l_cons->type == ConstraintType::linear;
    bool remove_right =
        removable[old_right_block] && r_cons->type == ConstraintType::linear;

    if (remove_left && remove_right) {
      /** both are removable linear */
      const linear_constraint<FieldT>* const l_lcons =
          static_cast<const linear_constraint<FieldT>*>(l_cons);
      const linear_constraint<FieldT>* const r_lcons =
          static_cast<const linear_constraint<FieldT>*>(r_cons);

      for (size_t lpos = 0; lpos < l_lcons->block.size(); ++lpos) {
        const FieldT lscalar = l_lcons->scalar[lpos];

        for (size_t rpos = 0; rpos < r_lcons->block.size(); ++rpos) {
          const FieldT rscalar = r_lcons->scalar[rpos];
          new_data.emplace_back(l_lcons->block[lpos], r_lcons->block[rpos],
                                (lscalar * rscalar) * old_scalar);
        }
      }
    } else if (remove_left || remove_right) {
      /** we will assume the left is to be removed, so swap if this is false */
      if (remove_right) {
        std::swap(l_cons, r_cons);
        std::swap(remove_left, remove_right);
        std::swap(old_left_block, old_right_block);
      }
      const linear_constraint<FieldT>* const l_lcons =
          static_cast<const linear_constraint<FieldT>*>(l_cons);
      for (size_t lpos = 0; lpos < l_lcons->block.size(); ++lpos) {
        const FieldT lscalar = l_lcons->scalar[lpos];
        new_data.emplace_back(l_lcons->block[lpos], old_right_block,
                              lscalar * old_scalar);
      }
    } else {
      new_data.emplace_back(old_left_block, old_right_block, old_scalar);
    }
  }

  /** flip left and right in new_data so that left < right */
  /** this is done to further reduce duplicates */
  for (data& entry : new_data) {
    if (std::get<0>(entry) > std::get<1>(entry)) {
      std::swap(std::get<0>(entry), std::get<1>(entry));
    }
  }

  /** sort by block so we can reduce redundant entries when combining */
  std::sort(new_data.begin(), new_data.end(),
            [](const data& a, const data& b) -> bool {
              if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) < std::get<0>(b);
              return std::get<1>(a) < std::get<1>(b);
            });

  /** push new data back into constraint */
  cons->left_block.clear();
  cons->right_block.clear();
  cons->scalar.clear();

  for (const data& entry : new_data) {
    if (cons->left_block.empty() ||
        cons->left_block.back() != std::get<0>(entry) ||
        cons->right_block.back() != std::get<1>(entry)) {
      cons->left_block.emplace_back(std::get<0>(entry));
      cons->right_block.emplace_back(std::get<1>(entry));
      cons->scalar.emplace_back(std::get<2>(entry));
    } else {
      /** combine scalar with old data */
      cons->scalar.back() += std::get<2>(entry);
    }
    /** remove zero scalars */
    if (!cons->left_block.empty() && cons->scalar.back() == FieldT(0)) {
      cons->left_block.pop_back();
      cons->right_block.pop_back();
      cons->scalar.pop_back();
    }
  }
}

/** rewires a transformation constraint depending on what is removable
 * @param idx index for the constraint
 * @param removable reference to the removability vector
 */
template <typename FieldT>
void builder<FieldT>::rewire_transformation(const size_t idx,
                                            const vector<bool>& removable) {
  transformation_constraint<FieldT>* const cons =
      static_cast<transformation_constraint<FieldT>*>(constraints[idx]);
  typedef tuple<size_t, size_t, size_t, FieldT> data;
  vector<data> new_data;

  /** we want to be able to look up all (source, pos, scalar) by [target, */
  /** target_position], where 'target' is anything we map FROM this will be */
  /** useful for combining */
  std::map<pair<size_t, size_t>, vector<tuple<size_t, size_t, FieldT>>> lookup;

  /** contains all blocks that 'lookup' contains data about */
  /** NOT nessecarily all blocks that the transformation depends on! */
  std::set<size_t> seen_blocks;
  /** fill the lookup */
  for (size_t pos = 0; pos < cons->block.size(); ++pos) {
    const size_t old_block = cons->block[pos];
    const constraint* const sub_cons = constraints[old_block];
    if (removable[old_block] && sub_cons->type == ConstraintType::linear) {
      if (seen_blocks.count(old_block)) continue;
      seen_blocks.insert(old_block);
      /** add each individual mapping implied by a linear constraint */
      const linear_constraint<FieldT>* const sub_lcons =
          static_cast<const linear_constraint<FieldT>*>(sub_cons);
      for (size_t sub_pos = 0; sub_pos < sub_lcons->block.size(); ++sub_pos) {
        const size_t sub_block = sub_lcons->block[sub_pos];
        const FieldT sub_scalar = sub_lcons->scalar[sub_pos];
        /** location in the block */
        for (size_t loc = 0; loc < block_size; ++loc) {
          lookup[make_pair(old_block, loc)].emplace_back(sub_block, loc,
                                                         sub_scalar);
        }
      }
    } else if (removable[old_block] &&
               sub_cons->type == ConstraintType::transformation) {
      if (seen_blocks.count(old_block)) continue;
      seen_blocks.insert(old_block);
      /** add each individual mapping implied by a transformation constraint */
      const transformation_constraint<FieldT>* const sub_tcons =
          static_cast<const transformation_constraint<FieldT>*>(sub_cons);
      for (size_t sub_pos = 0; sub_pos < sub_tcons->block.size(); ++sub_pos) {
        const size_t sub_block = sub_tcons->block[sub_pos];
        const size_t sub_source_position = sub_tcons->source_position[sub_pos];
        const size_t sub_target_position = sub_tcons->target_position[sub_pos];
        const FieldT& sub_scalar = sub_tcons->scalar[sub_pos];

        lookup[make_pair(old_block, sub_target_position)].emplace_back(
            sub_block, sub_source_position, sub_scalar);
      }
    }
  }

  /** populate new data and perform the 'matrix multiplication' */
  for (size_t pos = 0; pos < cons->block.size(); ++pos) {
    const size_t old_block = cons->block[pos];
    const size_t old_source_position = cons->source_position[pos];
    const size_t old_target_position = cons->target_position[pos];
    const FieldT old_scalar = cons->scalar[pos];

    /** only seen */
    if (seen_blocks.count(old_block)) {
      auto it = lookup.find(make_pair(old_block, old_source_position));
      if (it == lookup.end()) continue;
      /** iterate through list of transformations ending up in slot */
      /** (old_block, old_source_position) and combine */
      for (const auto& entry : it->second) {
        new_data.emplace_back(std::get<0>(entry), std::get<1>(entry),
                              old_target_position,
                              old_scalar * std::get<2>(entry));
      }
    } else {
      new_data.emplace_back(old_block, old_source_position, old_target_position,
                            old_scalar);
    }
  }

  /** sort so we can reduce redundant entries when combining */
  /** IMPORTANT : notice that we intentionally sort by source_position first */
  /** this is is so that when we later apply randomness, if a transformation */
  /** pulls from ITSELF (as it might!) it will only add scalars to locations */
  /** which have already adjusted all scalars of elements which map into it
   * under */
  /** the transformation */
  std::sort(new_data.begin(), new_data.end(),
            [](const data& a, const data& b) -> bool {
              if (std::get<1>(a) != std::get<1>(b))
                return std::get<1>(a) < std::get<1>(b);
              if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) < std::get<0>(b);
              return std::get<2>(a) < std::get<2>(b);
            });

  /** push new data back into constraint */
  cons->block.clear();
  cons->source_position.clear();
  cons->target_position.clear();
  cons->scalar.clear();

  /** remove duplicates and zeroes */
  for (const data& entry : new_data) {
    if (cons->block.empty() || cons->block.back() != std::get<0>(entry) ||
        cons->source_position.back() != std::get<1>(entry) ||
        cons->target_position.back() != std::get<2>(entry)) {
      cons->emplace(std::get<0>(entry), std::get<1>(entry), std::get<2>(entry),
                    std::get<3>(entry));
    } else {
      /** combine scalar with old data */
      cons->scalar.back() += std::get<3>(entry);
    }
    /** remove zero scalars */
    if (!cons->block.empty() && cons->scalar.back() == FieldT(0)) {
      cons->block.pop_back();
      cons->source_position.pop_back();
      cons->target_position.pop_back();
      cons->scalar.pop_back();
    }
  }
}

/** rewire a linear constraint depending on what is removable
 * we take an index instead of a linear constraint because the type may change
 * @param idx index for the constraint
 * @param removable reference to the removability vector
 */
template <typename FieldT>
void builder<FieldT>::rewire_linear(const size_t idx,
                                    const vector<bool>& removable) {
  linear_constraint<FieldT>* const cons =
      static_cast<linear_constraint<FieldT>*>(constraints[idx]);

  /** before anything, we check to see if it depends on any removable */
  /** transformation constraints if so, we simply convert it to a transform and
   */
  /** call that rewire function */

  for (size_t pos = 0; pos < cons->block.size(); ++pos) {
    const size_t old_block = cons->block[pos];
    const constraint* sub_cons = constraints[old_block];
    if (removable[old_block] &&
        sub_cons->type == ConstraintType::transformation) {
      constraints[idx] = convert_linear_to_transform(cons, block_size);
      rewire_transformation(idx, removable);
      return;
    }
  }

  /** we have only linear dependencies to worry about */
  /** this code is EXACTLY the same as for the inverse case */
  typedef pair<size_t, FieldT> data;
  vector<data> new_data;
  for (size_t pos = 0; pos < cons->block.size(); ++pos) {
    /** remove linear constraints */
    const size_t old_block = cons->block[pos];
    const FieldT old_scalar = cons->scalar[pos];

    const constraint* sub_cons = constraints[old_block];
    if (removable[old_block] && sub_cons->type == ConstraintType::linear) {
      const linear_constraint<FieldT>* sub_lcons =
          static_cast<const linear_constraint<FieldT>*>(sub_cons);

      /** forward data from linear constraint, multiplied by scalar */
      for (size_t sub_pos = 0; sub_pos < sub_lcons->block.size(); ++sub_pos) {
        const size_t new_block = sub_lcons->block[sub_pos];
        const FieldT new_scalar = sub_lcons->scalar[sub_pos];

        new_data.emplace_back(new_block, new_scalar * old_scalar);
      }
    } else {
      new_data.emplace_back(old_block, old_scalar);
    }
  }

  /** sort by block so we can reduce redundant entries when combining */
  std::sort(
      new_data.begin(), new_data.end(),
      [](const data& a, const data& b) -> bool { return a.first < b.first; });

  /** push new data back into constraint */
  cons->block.clear();
  cons->scalar.clear();
  for (const data& entry : new_data) {
    if (cons->block.empty() || cons->block.back() != entry.first) {
      cons->block.emplace_back(entry.first);
      cons->scalar.emplace_back(entry.second);
    } else {
      /** combine scalar with old data */
      cons->scalar.back() += entry.second;
    }
    /** remove zero scalars */
    if (!cons->block.empty() && cons->scalar.back() == FieldT(0)) {
      cons->block.pop_back();
      cons->scalar.pop_back();
    }
  }
}

/** rewire constraints that are non removable
 * @param removable reference to the removability vector
 */
template <typename FieldT>
void builder<FieldT>::rewire_nonremovable(const vector<bool>& removable) {
  for (size_t cons_idx = 0; cons_idx < constraints.size(); ++cons_idx) {
    const constraint* const cons = constraints[cons_idx];
    /* no variable, constant, or reserved will ever have to be rewired */
    if (cons->type == ConstraintType::variable ||
        cons->type == ConstraintType::constant ||
        cons->type == ConstraintType::reserved) {
      continue;
    }
    switch (cons->type) {
      // reduce linear constraints
      case ConstraintType::inverse:
        rewire_inverse(cons_idx, removable);
        break;
      case ConstraintType::linear:
        rewire_linear(cons_idx, removable);
        break;
      case ConstraintType::quadratic:
        rewire_quadratic(cons_idx, removable);
        break;
      case ConstraintType::transformation:
        rewire_transformation(cons_idx, removable);
        break;
      default:
        assert(false);
    }
  }
}

/** identify removable constraints and rewire nonremovable ones
 * @param return removability vector
 */
template <typename FieldT>
std::vector<bool> builder<FieldT>::find_and_rewire() {
  vector<bool> removable;
  // determine which constraints are removable;
  determine_removable(removable);
  // rewire constraints to remove removable ones
  // rewire_nonremovable(removable);
  return removable;
}

/** +----------------------+ */
/** | BEGINNING OF STAGE 3 | */
/** +----------------------+ */

/** compute proof values
 * @param variable_values matrix reference to populated variable values
 * @param removable removability vector
 * @param representative weighted union find structure
 * @return matrix of proof values
 */
template <typename FieldT>
vector<vector<FieldT>> builder<FieldT>::compute_proof_values(
    const vector<vector<FieldT>>& variable_values,
    const vector<bool>& removable, weightedUnionFind& representative) {
  vector<vector<FieldT>> proof_values;
  /** on the initial pass, we just determine the size of proof_values */
  size_t proof_values_sz = 0;
  for (size_t idx = 0; idx < constraints.size(); ++idx) {
    /** as removable constraints have been removed, we can simply skip them */
    /** we can also skip constant constraints */
    if (!removable[idx] && constraints[idx]->type != ConstraintType::constant)
      ++proof_values_sz;
  }

  proof_values.reserve(proof_values_sz);
  proof_values_location.resize((size_t)constraints.size(), (size_t)-1);

  // get the block of values
  auto get_block_vals = [&](size_t idx) -> const vector<FieldT>& {
    idx = representative.id(idx);
    if (constraints[idx]->type == ConstraintType::constant) {
      return static_cast<const constant_constraint<FieldT>*>(constraints[idx])
          ->value;
    } else {
      return proof_values[proof_values_location[idx]];
    }
  };

  size_t vars_index = 0;
  /** now compute the block values */
  for (size_t idx = 0; idx < constraints.size(); ++idx) {
    /** as removable constraints have been removed, we can simply skip them */
    /** we can also skip constant constraints */
    if (removable[idx] || constraints[idx]->type == ConstraintType::constant) {
      continue;
    }

    /** if you're not your own representative we don't need to compute you. */
    /** However, if you are a variable we do need to increment vars_index */
    if (representative.id(idx) != idx) {
      proof_values_location[idx] =
          proof_values_location[representative.id(idx)];
      if (constraints[idx]->type == ConstraintType::variable) {
        ++vars_index;
      }
      continue;
    }

    switch (constraints[idx]->type) {
      /** compute and inverse block */
      case ConstraintType::inverse: {
        const inverse_constraint<FieldT>* const cons =
            static_cast<const inverse_constraint<FieldT>*>(constraints[idx]);

        vector<FieldT> ans(block_size, FieldT(0));
        for (size_t idx = 0; idx < cons->block.size(); ++idx) {
          const auto& block = get_block_vals(cons->block[idx]);

#pragma omp parallel for
          for (size_t loc = 0; loc < block_size; ++loc) {
            assert(loc < block.size());
            assert(idx < cons->scalar.size());
            ans[loc] += (block[loc] * cons->scalar[idx]);
          }
        }

#pragma omp parallel for
        for (size_t loc = 0; loc < block_size; ++loc) {
          assert(ans[loc] != FieldT(0));
          ans[loc] = ans[loc].inverse();
        }

        proof_values.push_back(ans);
        proof_values_location[idx] = proof_values.size() - 1;
        break;
      }
      /** compute a linear block */
      case ConstraintType::linear: {
        const linear_constraint<FieldT>* const cons =
            static_cast<const linear_constraint<FieldT>*>(constraints[idx]);

        vector<FieldT> ans(block_size, FieldT(0));
        for (size_t idx = 0; idx < cons->block.size(); ++idx) {
          const auto& block = get_block_vals(cons->block[idx]);

#pragma omp parallel for
          for (size_t loc = 0; loc < block_size; ++loc) {
            ans[loc] += (block[loc] * cons->scalar[idx]);
          }
        }

        proof_values.push_back(ans);
        proof_values_location[idx] = proof_values.size() - 1;
        break;
      }
      case ConstraintType::quadratic: {
        const quadratic_constraint<FieldT>* const cons =
            static_cast<const quadratic_constraint<FieldT>*>(constraints[idx]);

        vector<FieldT> ans(block_size, FieldT(0));
        for (size_t idx = 0; idx < cons->left_block.size(); ++idx) {
          const auto& lblock = get_block_vals(cons->left_block[idx]);
          const auto& rblock = get_block_vals(cons->right_block[idx]);

#pragma omp parallel for
          for (size_t loc = 0; loc < block_size; ++loc) {
            ans[loc] = (lblock[loc] * rblock[loc]) * cons->scalar[idx];
          }
        }

        proof_values.push_back(ans);
        proof_values_location[idx] = proof_values.size() - 1;
        break;
      }
      case ConstraintType::variable:
        proof_values.push_back(variable_values[vars_index++]);
        proof_values_location[idx] = proof_values.size() - 1;
        break;
      case ConstraintType::transformation: {
        const transformation_constraint<FieldT>* const cons =
            static_cast<const transformation_constraint<FieldT>*>(
                constraints[idx]);

        vector<FieldT> ans(block_size, FieldT(0));
        for (size_t j = 0; j < cons->block.size(); ++j) {
          ans[cons->target_position[j]] +=
              (get_block_vals(cons->block[j])[cons->source_position[j]] *
               cons->scalar[j]);
        }

        proof_values.push_back(ans);
        proof_values_location[idx] = proof_values.size() - 1;
        break;
      }
      default:
        break;
    }
  }

  return proof_values;
}

/** generates the weighted union find structure used to determine the correct
 * representative. Constant takes precedence, else first thing
 * @return weighted union find structure
 */
template <typename FieldT>
weightedUnionFind builder<FieldT>::calc_equivalence_classes() {
  weightedUnionFind uf(constraints.size());
  for (size_t i = 0; i < constraints.size(); ++i) {
    if (constraints[i]->type == ConstraintType::constant) {
      uf.change_weight(i, 1);
    }
  }

  for (size_t i = 0; i < constraints.size(); ++i) {
    if (equality_graph[i].empty()) continue;
    for (size_t j : equality_graph[i]) {
      uf.merge(i, j);
    }
  }
  return uf;
}

/** +----------------------+ */
/** | BEGINNING OF STAGE 4 | */
/** +----------------------+ */

/** Building the constraint system: just filter through all the constraints we
 * kept and replace blocks with their reps
 * @param representative reference to the structure to determine representatives
 * @param removable reference to the removability vector
 * @return fully formed constraint system
 */
template <typename FieldT>
ConstraintSystem<FieldT> builder<FieldT>::build_system(
    weightedUnionFind& representative, const vector<bool>& removable) {
  /** figure out how many blocks will be in our final proof */
  size_t proof_blocks = 0;

  for (size_t idx = 0; idx < constraints.size(); ++idx) {
    if (!removable[idx] && constraints[idx]->type != ConstraintType::constant &&
        representative.id(idx) == idx)
      ++proof_blocks;
  }

  /** the constraint system we're compiling */
  ConstraintSystem<FieldT> cs(block_size, proof_blocks);

  /** where does each block end up, or INF if not determined yet */
  const size_t INF = size_t(-1);
  vector<size_t> eventual_location(constraints.size(), INF);

  /** we need to figure out where constants will be first, since in */
  /** some equivalence classes, the lowest block is not the representative */
  /** if it's a constant instead */
  vector<vector<FieldT>> saved_constants;

  for (size_t idx = 0; idx < constraints.size(); ++idx) {
    if (constraints[idx]->type == ConstraintType::constant) {
      if (representative.id(idx) == idx) {
        /** reserve a spot AFTER all the proof blocks */
        const constant_constraint<FieldT>* cons =
            static_cast<const constant_constraint<FieldT>*>(constraints[idx]);
        eventual_location[idx] = saved_constants.size() + proof_blocks;
        saved_constants.push_back(cons->value);
      } else {
        /** we must have already sent the rep somewhere */
        eventual_location[idx] = eventual_location[representative.id(idx)];
      }
    }
  }

  /** put in all the saved constants */
  cs.constant_blocks = saved_constants.size();
  cs.constants.resize(cs.constant_blocks * block_size);

  for (size_t loc = 0; loc < block_size; ++loc) {
    for (size_t block = 0; block < cs.constant_blocks; ++block) {
      cs.constants[loc * cs.constant_blocks + block] =
          saved_constants[block][loc];
    }
  }

  std::function<void(vector<size_t>&)> relabel = [&](vector<size_t>& v) {
    for (size_t& idx : v) {
      idx = eventual_location[representative.id(idx)];
    }
  };

  size_t assigned = 0;
  /** process all the non-constant blocks */
  for (size_t idx = 0; idx < constraints.size(); ++idx) {
    constraint* cons = constraints[idx];
    /** delete constraints which won't be moved to the system */
    if (removable[idx]) {
      delete cons;
      continue;
    }

    if (eventual_location[representative.id(idx)] == INF) {
      eventual_location[representative.id(idx)] = assigned++;
    }
    eventual_location[idx] = eventual_location[representative.id(idx)];

    switch (cons->type) {
      /** convert inverse into quadratic */
      case ConstraintType::inverse: {
        inverse_constraint<FieldT>* inv_cons =
            static_cast<inverse_constraint<FieldT>*>(cons);
        relabel(inv_cons->block);
        quadratic_constraint<FieldT>* q = new quadratic_constraint<FieldT>(
            {inv_cons->block,
             vector<size_t>(inv_cons->block.size(), eventual_location[idx])});

        cs.constraints.push_back(q);
        cs.targets.push_back(eventual_location[1]);
        delete cons;
        break;
      }

      case ConstraintType::linear:
        relabel(static_cast<linear_constraint<FieldT>*>(cons)->block);
        cs.constraints.push_back(cons);
        cs.targets.push_back(eventual_location[idx]);
        break;

      case ConstraintType::quadratic: {
        quadratic_constraint<FieldT>* c =
            static_cast<quadratic_constraint<FieldT>*>(cons);
        relabel(c->left_block);
        relabel(c->right_block);

        cs.constraints.push_back(cons);
        cs.targets.push_back(eventual_location[idx]);
      } break;

      case ConstraintType::transformation:
        relabel(static_cast<transformation_constraint<FieldT>*>(cons)->block);
        cs.constraints.push_back(cons);
        cs.targets.push_back(eventual_location[idx]);
        break;
        /** we are done with these descriptors now that blocks and constraints
         */
        /** are not one-and-one */
        [[fallthrough]];
      case ConstraintType::reserved:
      case ConstraintType::constant:
      case ConstraintType::variable:
        delete cons;
        break;
    }
  }

  /** add in special blocks */
  cs.special.resize(special_blocks.size());
  for (size_t idx = 0; idx < special_blocks.size(); ++idx) {
    cs.special[idx] = eventual_location[special_blocks[idx]];
  }

  return cs;
}

/** +-----------------------+ */
/** | Tying it all together | */
/** +-----------------------+ */

/** Perform all steps above to yield an optimized constraint system
 * @param values reference to the matrix of values in the circuit
 * @param populate flag specifying whether to populate values in the circuit
 * @return fully formed constraint system
 */
template <typename FieldT>
ConstraintSystem<FieldT> builder<FieldT>::compile(
    vector<vector<FieldT>>& values, bool populate) {
  vector<bool> removable = find_and_rewire();
  weightedUnionFind representative = calc_equivalence_classes();

  if (populate) {
    values = compute_proof_values(values, removable, representative);
  }

  DBG("Extended Witness");
  DBG("================");
  for (size_t i = 0; i < values.size(); i++) {
    std::string line = std::to_string(i) + std::string(":");

    for (size_t j = 0; j < std::min((int)values[i].size(), 10); j++) {
      line += std::to_string(values[i][j].getValue()) + std::string(" ");
    }

    DBG(line);
  }

  auto x = build_system(representative, removable);

  return x;
}

/** similar functionality */
template <typename FieldT>
ConstraintSystem<FieldT> builder<FieldT>::compile() {
  vector<bool> removable = find_and_rewire();
  weightedUnionFind representative = calc_equivalence_classes();
  return build_system(representative, removable);
}

template <typename FieldT>
struct vars {
  std::map<size_t, size_t> mapping;
  vector<vector<FieldT>>* data;
  const size_t block_size;

  vars(vector<vector<FieldT>>* data, const size_t block_size)
      : data(data), block_size(block_size) {}

  void add(size_t block) {
    mapping[block] = data->size();
    data->push_back(vector<FieldT>(block_size, FieldT(0)));
  }

  void add(size_t block, const vector<FieldT>& row) {
    mapping[block] = data->size();
    data->push_back(row);
  }

  vector<FieldT>& operator[](size_t block) {
    assert(mapping.count(block));
    return (*data)[mapping[block]];
  }
};