#pragma once

/* Third Party */
#include <sodium/crypto_generichash_blake2b.h>
#include <sodium/randombytes.h>

/* Third Party */
#include <sodium/crypto_generichash_blake2b.h>
#include <sodium/randombytes.h>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

/* Ligero */
#include "Common.hpp"
#include "Math.hpp"

constexpr unsigned short iSeedSize = 32;

/** ============================================================= */
/**  Hash definitions */
namespace hash {
typedef std::string digest;
typedef std::pair<size_t, hash::digest> PosAndDigest;

/** Structure contains the decommitment of a set of positions
 * in the Merkle Tree
 */
template <typename FieldT>
struct sparseMultipleDecommitment {
  std::vector<digest> randomnessHashes;
  std::vector<digest> auxiliaryHashes;
  std::vector<digest> contentHashes;
  std::vector<std::vector<FieldT>> contents;
  std::vector<size_t> positions;

  friend class boost::serialization::access;

 private:
  /** Implements serialization and deserialization
   * @param ar Archive Used
   * @param version Version of the implementation for the class serialized
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &randomnessHashes;
    ar &auxiliaryHashes;
    ar &contentHashes;
    ar &contents;
  }
};

/** Forward declarations Hashing Functions
 */
template <typename FieldT>
using HashLeafContentF =
    std::function<digest(const FieldT *, const size_t, const size_t)>;

typedef std::function<digest(const digest &, const digest &, size_t)>
    HashInnerNodeF;
typedef std::function<digest(const uint8_t *, const size_t, const size_t)>
    HashZKRandomnessF;

template <typename FieldT>
using HashIntegerSeedGeneratorF =
    std::function<std::vector<unsigned char>(const hash::digest &)>;
using HashStringSeedGeneratorF = std::function<std::vector<unsigned char>(
    const char *data, const size_t length)>;

/** Wrappers for Libsodium Hashing Functions
 */

/** Computes the Blake2b hash of an existing digest (hash).
 * @param root A hashed value, for instance the root of a merkle tree
 * @return A vector of unsigned chars (length set by constant iSeedSize)
 * containing the hash
 */
template <typename FieldT>
std::vector<unsigned char> blake2bIntegerSeedGenerator(
    const hash::digest &root) {
  std::vector<unsigned char> result(iSeedSize);
  size_t scaledResult;
  size_t index = 0;

  if (crypto_generichash_blake2b((unsigned char *)&result[0], sizeof(result),
                                 (unsigned char *)&root[0], root.size(),
                                 (unsigned char *)&index, sizeof(index)) != 0)
    throw std::runtime_error(
        "Libsodium Encryption Error Randomness Generator "
        "(crypto_generichash_blake2b)");

  return result;
}

/** Computes the Blake2b hash of generic content.
 * @param data a pointer to an array of char containing the data to be hashed.
 * @param length the length of the array
 * @return A vector of unsigned chars (length set by constant iSeedSize)
 * containing the hash
 */
std::vector<unsigned char> blake2bStringHash(const char *data,
                                             const size_t length) {
  std::vector<unsigned char> result(iSeedSize);

  if (crypto_generichash_blake2b(
          (unsigned char *)&result[0], iSeedSize,
          (result.empty() ? NULL : (unsigned char *)data), length, NULL,
          0) != 0)
    throw std::runtime_error(
        "Libsodium Encryption Error Field Element "
        "(crypto_generichash_blake2b)");

  return result;
}

/** Computes a Blake2b hash of desired size from values in a specified field.
 * @param data a pointer to an array of char containing the data to be hashed
 * @param length the length of the array
 * @param digestLengthBytes desired length of the output hash, in bytes
 * @return A hash of length digestLengthBytes
 */
template <typename FieldT>
digest blake2bFieldElementHash(const FieldT *data, const size_t length,
                               const size_t digestLengthBytes) {
  hash::digest result(digestLengthBytes, '-');

  if (crypto_generichash_blake2b(
          (unsigned char *)&result[0], digestLengthBytes,
          (result.empty() ? NULL : (unsigned char *)data),
          sizeof(FieldT) * length, NULL, 0) != 0)
    throw std::runtime_error(
        "Libsodium Encryption Error Field Element "
        "(crypto_generichash_blake2b)");

  return result;
}

/** Computes a Blake2b hash of desired size from content stored as an array of
 * unsigned char (uint8_t).
 * @param data a pointer to an array of char containing the data to be hashed
 * @param length the length of the array
 * @param digestLengthBytes desired length of the output hash, in bytes
 * @return A hash of length digestLengthBytes
 */
digest blake2bZKElementHash(const uint8_t *data, const size_t length,
                            const size_t digestLengthBytes) {
  hash::digest result(digestLengthBytes, '-');

  if (crypto_generichash_blake2b((unsigned char *)&result[0], digestLengthBytes,
                                 (unsigned char *)data,
                                 length * sizeof(uint8_t), NULL, 0) != 0) {
    std::cout << length * sizeof(uint8_t) << std::endl;
    throw std::runtime_error(
        "Libsodium Encryption Error ZK Element (crypto_generichash_blake2b)");
  }
  return result;
}

/** Generates a single output hash of desired size from two concatenated input
 * hashes
 * @param first first input hash, provided as a hash::digest
 * @param second second input hash, provided as a hash::digest
 * @param digestLengthBytes desired length of the output hash, in bytes
 * @return A hash of length digestLengthBytes
 */
digest blake2bTwoToOneHash(const hash::digest &first,
                           const hash::digest &second,
                           const size_t digestLengthBytes) {
  const hash::digest firstPlusSecond = first + second;
  hash::digest result(digestLengthBytes, '-');

  if (crypto_generichash_blake2b((unsigned char *)&result[0], digestLengthBytes,
                                 (unsigned char *)&firstPlusSecond[0],
                                 firstPlusSecond.size(), NULL, 0))
    throw std::runtime_error(
        "Libsodium Encryption Error Inner Node (crypto_generichash_blake2b)");

  return result;
}

/** Hashing of the combined Public Data and Sigma Public Data
 * @param hash the hash function to be used
 * @param pdata the structure containing all generic Public Data
 * @param spdata the structure containing Public Data tied to the Sigma Protocol
 * (Sigma Public Data)
 * @return the hash is returned as a string
 */
std::string hashPublicData(const hash::HashStringSeedGeneratorF &hash,
                           PublicData &pdata, SigmaProtocolPublicData &spdata) {
  std::string msgPdata = serialize<PublicData>(pdata);
  std::string msgSpdata = serialize<SigmaProtocolPublicData>(spdata);
  std::string msg = msgPdata + msgSpdata;

  std::vector<unsigned char> outputVector = hash(msg.c_str(), msg.size());
  std::string out(outputVector.begin(), outputVector.end());

  return out;
}

}  // namespace hash
