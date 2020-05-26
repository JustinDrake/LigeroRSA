#pragma once
#include <bits/stdint-uintn.h>
#include <cassert>
#include <vector>

using std::vector;

enum ConstraintType {
  inverse,
  linear,
  quadratic,
  constant,
  variable,
  transformation,
  reserved
};

std::unordered_map<ConstraintType, std::string> constraintMap = {
    {inverse, "inverse"},     {linear, "linear"},
    {quadratic, "quadratic"}, {constant, "constant"},
    {variable, "variable"},   {transformation, "transformation"},
    {reserved, "reserved"},
};

class constraint {
 public:
  ConstraintType type;
  virtual size_t size() const { return 0; }
  virtual ~constraint() {}
};

/** X = (s1*A1 + ... + sn*An)^(-1) */
template <typename FieldT>
class inverse_constraint : public constraint {
 public:
  inverse_constraint(const std::vector<size_t> &b);
  inverse_constraint() { type = ConstraintType::inverse; }
  std::vector<size_t> block;
  std::vector<FieldT> scalar;
  size_t size() const { return block.size(); }
  void emplace(size_t b, FieldT s = 1);
};

/** X = s*A + t*B */
template <typename FieldT>
class linear_constraint : public constraint {
 public:
  linear_constraint(const std::vector<size_t> &b);
  linear_constraint(const std::vector<size_t> &b, const std::vector<FieldT> &s);
  linear_constraint() { type = ConstraintType::linear; }

  std::vector<size_t> block;
  std::vector<FieldT> scalar;
  size_t size() const { return block.size(); }
  void emplace(size_t b, FieldT s = 1);
};

/** X = s*AB + t*CD */
template <typename FieldT>
class quadratic_constraint : public constraint {
 public:
  quadratic_constraint(const std::vector<size_t> &l,
                       const std::vector<size_t> &r);
  quadratic_constraint(const std::vector<size_t> &l,
                       const std::vector<size_t> &r,
                       const std::vector<FieldT> &s);
  quadratic_constraint() { type = ConstraintType::quadratic; }
  std::vector<size_t> left_block;
  std::vector<size_t> right_block;
  std::vector<FieldT> scalar;
  size_t size() const { return left_block.size(); }
  void emplace(size_t a, size_t b, FieldT s = 1);
};

/** X = [fixed block] */
template <typename FieldT>
class constant_constraint : public constraint {
 public:
  constant_constraint(const std::vector<FieldT> &f);
  constant_constraint() { type = ConstraintType::constant; }
  std::vector<FieldT> value;
  size_t size() const { return value.size(); }
  void emplace(FieldT s);
};

/** X = [indeterminate] */
template <typename FieldT>
class variable_constraint : public constraint {
 public:
  variable_constraint() { type = ConstraintType::variable; }
  size_t size() const { return 0; }
};

/** X = linear_transform(some blocks) */
template <typename FieldT>
class transformation_constraint : public constraint {
 public:
  transformation_constraint(const std::vector<size_t> &b,
                            const std::vector<size_t> &s,
                            const std::vector<size_t> &t);

  transformation_constraint() {
    type = ConstraintType::transformation;
    stitching = LinearCombinationsStitching::Regular;
  }

  std::vector<size_t> block;
  std::vector<size_t> source_position;
  std::vector<size_t> target_position;
  std::vector<FieldT> scalar;
  LinearCombinationsStitching stitching;
  size_t size() const { return block.size(); }
  void emplace(size_t b, size_t src, size_t tgt, FieldT s = 1);
};

/** A block that doesn't have a value yet but is not a variable
 * this will be computed by the builder
 */
template <typename FieldT>
class reserved_constraint : public constraint {
 public:
  reserved_constraint() { type = ConstraintType::reserved; }
  reserved_constraint(const std::vector<size_t> &d);
  /** generic data to decode */
  std::vector<size_t> data;
  size_t size() const { return data.size(); }
  void emplace(size_t size);
};

/** Constructor for inverse constraint: instanciates for a set of private blocks
 * b
 * @param b vector of indices of blocks for which to calculate inverses, with a
 * default scalar of 1
 */
template <typename FieldT>
inverse_constraint<FieldT>::inverse_constraint(const vector<size_t> &b)
    : block(b), scalar(block.size(), FieldT(1)) {
  type = ConstraintType::inverse;
}

/** Add an incremental block to an inverse constraint
 * @param b index of the block to be added
 * @param s scalar for this block
 */
template <typename FieldT>
void inverse_constraint<FieldT>::emplace(size_t b, FieldT s) {
  block.emplace_back(b);
  scalar.emplace_back(s);
}

/** Constructor for linear interblock constraint: instanciates for a set of
 * private blocks b with an assumed scalar of 1 for each block
 * @param b vector of indices of blocks, with a default scalar of 1
 */
template <typename FieldT>
linear_constraint<FieldT>::linear_constraint(const vector<size_t> &b)
    : block(b), scalar(b.size(), FieldT(1)) {
  type = ConstraintType::linear;
}

/** Alternate constructor for linear interblock constraint: instanciates for a
 * set of private blocks b and related scalars s
 * @param b vector of indices of blocks
 * @param s vector of scalars for each block
 */
template <typename FieldT>
linear_constraint<FieldT>::linear_constraint(const vector<size_t> &b,
                                             const vector<FieldT> &s)
    : block(b), scalar(s) {
  type = ConstraintType::linear;
}

/** Add an incremental block to an interblock constraint
 * @param b index of the block to be added
 * @param s scalar for this block
 */
template <typename FieldT>
void linear_constraint<FieldT>::emplace(size_t b, FieldT s) {
  block.emplace_back(b);
  scalar.emplace_back(s);
}

/** Constructor for quadratic constraint: instanciates for a set of
 * private blocks l and r (resp. left handside and right handside of the
 * multiplication) assuming scalars of 1 for each multiplication
 * @param l vector of left handside blocks
 * @param r vector of right handside blocks
 */
template <typename FieldT>
quadratic_constraint<FieldT>::quadratic_constraint(const vector<size_t> &l,
                                                   const vector<size_t> &r)
    : left_block(l), right_block(r), scalar(l.size(), FieldT(1)) {
  type = ConstraintType::quadratic;
}

/** Alternate constructor for quadratic constraint: instanciates for a set of
 * private blocks l and r (resp. left handside and right handside of the
 * multiplication) with scalars s attached to each multiplication
 * @param l vector of left handside blocks
 * @param r vector of right handside blocks
 * @param s vector of scalars applied to each multiplication
 */
template <typename FieldT>
quadratic_constraint<FieldT>::quadratic_constraint(const vector<size_t> &l,
                                                   const vector<size_t> &r,
                                                   const vector<FieldT> &s)
    : left_block(l), right_block(r), scalar(s) {
  type = ConstraintType::quadratic;
}

/** Add an incremental block to a quadratic constraint
 * @param a index of the left handside block
 * @param b index of the right handside block
 * @param s scalar for this multiplication
 */
template <typename FieldT>
void quadratic_constraint<FieldT>::emplace(size_t a, size_t b, FieldT s) {
  left_block.emplace_back(a);
  right_block.emplace_back(b);
  scalar.emplace_back(s);
}

/** Constructor for constant constraint: this adds constants to the constraint
 * system that are publicly known
 * @param f vector of field elements to be added as constant constraints
 */
template <typename FieldT>
constant_constraint<FieldT>::constant_constraint(const vector<FieldT> &f)
    : value(f) {
  type = ConstraintType::constant;
}

/** Add an incremental block to a constant constraint
 * @param s Field element to be added as a constant
 */
template <typename FieldT>
void constant_constraint<FieldT>::emplace(FieldT s) {
  value.emplace_back(s);
}

/** Constructor for a transformation constraint: this type of constraint ties in
 * two elements of the witness that are not necessarily located in the same
 * column: y = Σ x * s
 * The constructor assumes s = 1
 * @param b vector of block indices containing x
 * @param s vector of columns containing x
 * @param t vector of columns containing y
 *
 */
template <typename FieldT>
transformation_constraint<FieldT>::transformation_constraint(
    const vector<size_t> &b, const vector<size_t> &s, const vector<size_t> &t)
    : block(b), source_position(s), target_position(t), scalar(b.size(), 1) {
  assert(b.size() == s.size());
  assert(b.size() == t.size());
  type = ConstraintType::transformation;
}

/** Add an incremental block to a transformation constraint
 * @param b source block
 * @param src source column
 * @param tgt target column
 * @param s scalar to be applied to the source
 */
template <typename FieldT>
void transformation_constraint<FieldT>::emplace(size_t b, size_t src,
                                                size_t tgt, FieldT s) {
  block.emplace_back(b);
  source_position.emplace_back(src);
  target_position.emplace_back(tgt);
  scalar.emplace_back(s);
}

/** Constructor for reserved constraint
 * @param d vector of indices of blocks being reserved
 */
template <typename FieldT>
reserved_constraint<FieldT>::reserved_constraint(const std::vector<size_t> &d) {
  type = ConstraintType::reserved;
  data = d;
}

/** adds a block to a reserve contraint
 * @param x index of the block to be added
 */
template <typename FieldT>
void reserved_constraint<FieldT>::emplace(size_t x) {
  data.emplace_back(x);
}