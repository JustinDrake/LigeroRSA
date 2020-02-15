#pragma once
#include <vector>
#include <cstdint>

// union find where the id-representative is of the highest weight
// only guaranteed correct if changing weights for top-level nodes
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

weightedUnionFind::weightedUnionFind(size_t size) {
	nodes.resize(size);
	for(size_t a = 0; a < size; ++a)
		nodes[a] = a;
	weights.resize(size);
}

void weightedUnionFind::change_weight(size_t a, uint8_t weight) {
	assert(nodes[a] == a);
	weights[a] = weight;
}

size_t weightedUnionFind::id(size_t a) {
	if(nodes[a] == a) return a;
	return nodes[a] = id(nodes[a]);
}

void weightedUnionFind::merge(size_t a, size_t b) {
	a = id(a);
	b = id(b);
	if(weights[a] == weights[b]) {
		if(a < b) nodes[b] = a;
		else nodes[a] = b;
	} else if(weights[a] > weights[b]) {
		nodes[b] = a;
	 } else {
		nodes[a] = b;
	 }
}
