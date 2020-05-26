#pragma once
#include <cstdint>
#include <vector>

/** union finds where the id-representative is of the highest weight
 * only guaranteed correct if changing weights for top-level nodes
 */
class weightedUnionFind {
 public:
  weightedUnionFind(size_t size);
  size_t id(size_t a);
  void change_weight(size_t a, uint8_t weight);
  void merge(size_t a, size_t b);

 private:
  std::vector<size_t> nodes;
  std::vector<uint8_t> weights;
};

/** Constructor: intializes the nodes structure
 * @param size size of the structure to initialize
 */
weightedUnionFind::weightedUnionFind(size_t size) {
  nodes.resize(size);
  for (size_t a = 0; a < size; ++a) nodes[a] = a;
  weights.resize(size);
}

/** changes the weight for a specific node
 * @param a node index
 * @param weight new weight
 */
void weightedUnionFind::change_weight(size_t a, uint8_t weight) {
  assert(nodes[a] == a);
  weights[a] = weight;
}

/** returns node identifier
 * @param a node index
 * @return node identifier
 */
size_t weightedUnionFind::id(size_t a) {
  if (nodes[a] == a) return a;
  return nodes[a] = id(nodes[a]);
}

/** merging two nodes
 * @param a first node index
 * @param b second node index
 */
void weightedUnionFind::merge(size_t a, size_t b) {
  a = id(a);
  b = id(b);
  if (weights[a] == weights[b]) {
    if (a < b)
      nodes[b] = a;
    else
      nodes[a] = b;
  } else if (weights[a] > weights[b]) {
    nodes[b] = a;
  } else {
    nodes[a] = b;
  }
}
