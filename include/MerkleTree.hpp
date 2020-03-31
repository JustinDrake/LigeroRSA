/* 
    Ligero Merkle Tree Interface 

    streamlined interface, enables commitment (constructing a tree), decommitment (providing an authentication path) 
    and verification (recomputing the root based on the authentication path provided, and verify against the original root commitment)
    */

#pragma once

/* sortedHashTableTL */
#include <string>
#include <vector>
#include <functional>

/* Ligero */
#include "Common.hpp"
#include "Math.hpp"
#include "Hash.hpp"

namespace ligero {

    template<typename FieldT> 
    struct MTParameters {
         hash::HashInnerNodeF hashInnerNodeF;
         hash::HashLeafContentF<FieldT> hashLeafContentF;
         hash::HashZKRandomnessF hashZKRandomnessF;
         size_t digestLength;
         size_t zkByteSize;
         size_t leavesNumber;
         size_t leavesEntropy;
    };

/*
    Merkle Tree Commitment Scheme Class Definition 
    Inner Nodes Structure Representation:
    
    [ Tree structure [leaves = 32]:               |0*|**|****|********|****************|********************************]
    [ Layer idx                                     1  2   3      4           5                        6                ]
    [ Layer length                                  1  2   4      8           16                       32               ]
    
    The top layer of the tree has dimension one, as it only contains the root. The total inner node length (ommitting the top layer) is 
    _leavesNumber - 1, as per the diagram above.  
    */

template<typename FieldT> 
class MT_CommitmentScheme {

protected:
    size_t                              _digestLength;
    size_t                              _leavesNumber;
    FieldT*                             _leavesContent;
    size_t                              _leavesEntropy;
    
    hash::digest*                       _nodes;
    hash::HashLeafContentF<FieldT>      _hashLeavesContentF;
    hash::HashInnerNodeF                _hashInnerNodeF;
    hash::HashZKRandomnessF             _hashZKRandomnessF;

    std::vector<std::vector<uint8_t>>   _zkRandomness;
    size_t                              _zkByteSize;
    bool                                _nodeAllocation;

public:
    MT_CommitmentScheme(
                const hash::HashInnerNodeF& hashInnerNodeF,
                const hash::HashLeafContentF<FieldT>& hashLeafContentF,
                const hash::HashZKRandomnessF& hashZKRandomnessF,
                const size_t digestLength,
                const size_t zkByteSize,
                const size_t leavesNumber,      // number of servers
                const size_t leavesEntropy      // number of blocks
                );

    ~MT_CommitmentScheme() {
        if (this->_nodeAllocation) {delete[] this->_nodes;}
    };

    hash::digest commit(
                const FieldT* leavesContent);

    hash::sparseMultipleDecommitment<FieldT> decommit(const std::vector<size_t> & decommitPositions);
    bool verify(const hash::digest &root, const hash::sparseMultipleDecommitment<FieldT>& decommitment, const std::vector<size_t> & decommitPositions);
};

/* 
    Methods Implementation 
    */

template<typename FieldT>
MT_CommitmentScheme<FieldT>::MT_CommitmentScheme(
               const hash::HashInnerNodeF& hashInnerNodeF,
               const hash::HashLeafContentF<FieldT>& hashLeafContentF,
               const hash::HashZKRandomnessF& hashZKRandomnessF,
               const size_t digestLength,               
               const size_t zkByteSize, 
               const size_t leavesNumber,
               const size_t leavesEntropy
               ):   _digestLength(digestLength),
                    _hashLeavesContentF(hashLeafContentF),
                    _hashInnerNodeF(hashInnerNodeF),
                    _hashZKRandomnessF(hashZKRandomnessF),
                    _zkByteSize(zkByteSize),
                    _leavesNumber(leavesNumber),
                    _leavesEntropy(leavesEntropy)
    {
        this->_nodeAllocation = false;
        this->_leavesNumber   = leavesNumber;
        this->_leavesEntropy  = leavesEntropy;
    }

template<typename FieldT>
hash::digest MT_CommitmentScheme<FieldT>::commit(
               const FieldT* leavesContent) {
        /*
            Initialize Variables for this commitment
            */

        this->_leavesContent    = (FieldT*)leavesContent;
        size_t& leavesNumber    = this->_leavesNumber;
        size_t& leavesEntropy   = this->_leavesEntropy;

        /* 
            Each layer i of the MT occupies m(i) memory such that m(i) = m(i-1)/2, with layer 0 being the input (leaves)  
            Hence the MT requires 2 * n - 1 (as the element left of the root can be discarded) 
            */

        this->_nodes = new hash::digest[leavesNumber * 2 - 1];
        this->_nodeAllocation = true;

        /* We generate randomness on the heap, incurring a minor performance hit */
        uint8_t *tmp = new uint8_t[this->_zkByteSize * this->_leavesNumber];
        randombytes_buf(tmp, this->_leavesNumber * this->_zkByteSize);

        std::vector<FieldT> leaf(leavesEntropy);

        /* First hash all leaves */
        for (size_t idx = 0 ; idx < leavesNumber ; idx++) {

            /* Additional Randomness for ZK */
            this->_zkRandomness.emplace_back(tmp + idx * this->_zkByteSize, tmp + (idx + 1) * this->_zkByteSize);
            hash::digest randDigest = this->_hashZKRandomnessF((const uint8_t *)&this->_zkRandomness.back()[0], this->_zkByteSize, this->_digestLength);

            /* Leaf */
            for (size_t row = 0; row < leavesEntropy; row++) {leaf[row]=*(_leavesContent+idx+(row*this->_leavesNumber));}
            hash::digest leafDigest = this->_hashLeavesContentF(&leaf[0], leavesEntropy, this->_digestLength);

            /* Final digest */
            _nodes[leavesNumber + idx - 1] = this->_hashInnerNodeF(leafDigest, randDigest, this->_digestLength);
        }

        /* Clean up */
        delete[] tmp;

        /* Then build the inner leaves */
        size_t nbLeavesThisLayer = leavesNumber;

        /* While we have not reached the Root */
        while (nbLeavesThisLayer > 1) {

            /* ZK: produce randomness for each node */
            /* Produce a hash for each node         */
            for (size_t idx = nbLeavesThisLayer ; idx < 2 * nbLeavesThisLayer; idx+=2) {
                _nodes[idx/2-1] = this->_hashInnerNodeF(_nodes[idx-1],_nodes[idx], this->_digestLength);
            }            

            /* sortedHashTabletep forward to the next layer */
            nbLeavesThisLayer = nbLeavesThisLayer >> 1;
        }
    
        return (this->_nodes[0]);
    }

/*
    Convention on the decommitment: positions are [ 0 , nbLeaves - 1 ] 
 */
template<typename FieldT>
hash::sparseMultipleDecommitment<FieldT> MT_CommitmentScheme<FieldT>::decommit(const std::vector<size_t> &positions)
{
    hash::sparseMultipleDecommitment<FieldT> result;

    /* 
        Prepare a sorted set of positions 
        and remove potential duplicates 
        */
    result.positions = positions;
    std::vector<size_t> inputSortedPositions(positions); 

    std::sort(inputSortedPositions.begin(), inputSortedPositions.end());
    inputSortedPositions.erase(std::unique(inputSortedPositions.begin(), inputSortedPositions.end()), inputSortedPositions.end()); 

    /*  
        Now generate the digest for combination of original content and randomness
        */
    std::vector<size_t> sortedPositions; 
    sortedPositions.reserve(inputSortedPositions.size());

    for (auto pos : inputSortedPositions)
    {
        /* Add Content with ZK Randomness */
        const hash::digest randomDigest = this->_hashZKRandomnessF((const uint8_t *)&this->_zkRandomness[pos][0], this->_zkByteSize, this->_digestLength);
        result.randomnessHashes.emplace_back(randomDigest);

        /* One time adustment for initial tree leaves */
        sortedPositions.emplace_back(pos + this->_leavesNumber - 1);

        /* Leaf */
        std::vector<FieldT> leaf(this->_leavesEntropy);
        for (size_t row = 0; row < this->_leavesEntropy; row++) {leaf[row]=*(_leavesContent+pos+(row*this->_leavesNumber));}
        const hash::digest contentDigest = this->_hashLeavesContentF(&leaf[0], _leavesEntropy, this->_digestLength);
        result.contents.emplace_back(leaf);
        result.contentHashes.emplace_back(contentDigest);

        assert(this->_hashInnerNodeF(result.contentHashes.back(), result.randomnessHashes.back(), this->_digestLength) == this->_nodes[pos + this->_leavesNumber - 1]);
    };

    /* 
        Add auxiliary hashes for each query and skip overlaps
        */
    while (true) 
    {
        auto it = sortedPositions.begin();

        /* have we arrived at the root yet? */
        if (*it == 0 && it == --sortedPositions.end()) break;

        std::vector<size_t> newsortedPositions;
        while (it != sortedPositions.end())
        {
            const size_t posIt = *it;
            auto nextIt = ++it;

            /* Always process parent.                                        */
            /* sortedPositionshould be positioned in ( (posIt + 3)/2 - 1 )   */
            newsortedPositions.emplace_back((posIt - 1)/2);             

            if ((posIt & 1) == 0)
                result.auxiliaryHashes.emplace_back(this->_nodes[posIt - 1]);
            else
            {
                /* We are the left node. Two cases: */
                if (nextIt == sortedPositions.end() || *nextIt != posIt + 1)
                    result.auxiliaryHashes.emplace_back(this->_nodes[posIt + 1]);
                else
                    ++nextIt;
            }
            it = nextIt;
        }

        std::swap(sortedPositions, newsortedPositions);
    }

    return result;
}

template<typename FieldT>
bool MT_CommitmentScheme<FieldT>::verify(const hash::digest &root, const hash::sparseMultipleDecommitment<FieldT> &decommitment, const std::vector<size_t> &positions)
{
    auto randomnessIt = decommitment.randomnessHashes.begin();
    auto auxiliaryIt = decommitment.auxiliaryHashes.begin();

    //assert(this->_hashInnerNodeF(result.contentHashes.back(), result.randomnessHashes.back(), this->_digestLength) == this->_nodes[pos + this->_leavesNumber - 1]);

    /* transform to sorted set of indices */
    std::vector<size_t> sortedPositions(positions); 
    std::sort(sortedPositions.begin(), sortedPositions.end());

    /* remove possible duplicates */
    sortedPositions.erase(std::unique(sortedPositions.begin(), sortedPositions.end()), sortedPositions.end()); 

    if (sortedPositions.size()!=decommitment.contents.size()) {
        std::cout << "positions queried and content inconsistent, Fiat-Shamir inconsistent between party and verifier." << std::endl;
        exit(1);
    }

    // std::cerr << "digest1:" << (int)root[0] << std::endl;
    // std::cerr << "digest2:" << (int)decommitment.randomnessHashes[0][0] << std::endl;
    // std::cerr << "digest3:" << (int)decommitment.auxiliaryHashes[0][0] << std::endl;
    // std::cerr << "digest4:" << (int)decommitment.contentHashes[0][0] << std::endl;
    // std::cerr << "digest5:" << (int)sortedPositions[0] << std::endl;

    /* Check contentHash against leaves */
    for (size_t idx = 0 ; idx < sortedPositions.size(); idx++) {
        const hash::digest contentDigest = this->_hashLeavesContentF(&decommitment.contents[idx][0], _leavesEntropy, this->_digestLength);
        if (contentDigest != decommitment.contentHashes[idx]) {
            DBG("Content/Hash Integrity Failure, Early False Response.");
            return false;
        }
    }

    /* one time position adjustment */
    for (auto &pos : sortedPositions) pos += this->_leavesNumber - 1;

    /* preparing and populating the hash table containing each position and digest */
    std::vector<hash::PosAndDigest> sortedHashTable;
    sortedHashTable.reserve(sortedPositions.size());

    std::vector<hash::digest> combinedDigests;
    for (auto &contentDigest : decommitment.contentHashes)
    {
        const hash::digest randomDigest = *randomnessIt++;
        combinedDigests.emplace_back(this->_hashInnerNodeF(contentDigest, randomDigest, this->_digestLength));
    }

    std::transform(
                sortedPositions.begin(), 
                sortedPositions.end(), 
                combinedDigests.begin(),
                std::back_inserter(sortedHashTable),
                [](const size_t pos, const hash::digest &hash) {return std::make_pair(pos, hash);});

    /* now proceed to calculating the root based on decommitment provided */
    while (true)
    {
        auto it = sortedHashTable.begin();

        /* have we arrived at the root yet? */
        if (it->first == 0 && it == --sortedHashTable.end()) break;

        std::vector<hash::PosAndDigest> newSortedHashTable;

        while (it != sortedHashTable.end())
        {
            const size_t posIt = it->first;
            const hash::digest hashIt = it->second;

            auto nextIt = ++it;

            hash::digest leftHash;
            hash::digest rightHash;

            if ((posIt & 1) == 0)
            {
                leftHash = *auxiliaryIt++;
                rightHash = hashIt;
            }
            else
            {
                leftHash = hashIt;

                if (nextIt == sortedHashTable.end() || nextIt->first != posIt + 1) rightHash = *auxiliaryIt++;
                else
                {
                    rightHash = nextIt->second;
                    ++nextIt;
                }
            }

            const size_t parentPos = (posIt - 1)/2;
            const hash::digest parentHash = this->_hashInnerNodeF(leftHash,rightHash,this->_digestLength);
            newSortedHashTable.emplace_back(std::make_pair(parentPos, parentHash));

            it = nextIt;
        }

        std::swap(sortedHashTable, newSortedHashTable);
    }

    return (sortedHashTable.begin()->second == root);
}

}