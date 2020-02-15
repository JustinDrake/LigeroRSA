#pragma once

/* Third Party */
#include <sodium/randombytes.h>
#include <sodium/crypto_generichash_blake2b.h>

/* Third Party */
#include <sodium/randombytes.h>
#include <sodium/crypto_generichash_blake2b.h>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

/* Ligero */
#include "Common.hpp"
#include "Math.hpp"

constexpr unsigned short iSeedSize = 32;

//==============================================================================
// Hash definitions
namespace hash {
    typedef std::string digest;
    typedef std::pair<size_t, hash::digest> PosAndDigest;

    // template<class Archive, typename dataType>
    // void archive(Archive& ar, std::vector<dataType> v) {
    //     for (auto h : v) {ar & h;}
    // }

    template<typename FieldT>
    struct sparseMultipleDecommitment {
        std::vector<digest>                     randomnessHashes;
        std::vector<digest>                     auxiliaryHashes;
        std::vector<digest>                     contentHashes;
        std::vector<std::vector<FieldT>>        contents;
        std::vector<size_t>                     positions;
        // int a;
        friend class boost::serialization::access;

        private:
            // Implements serialization and deserialization
            template <class Archive>
            void serialize(Archive& ar, const unsigned int version)
            {
                // archive<Archive, digest>(ar,randomnessHashes);
                // archive<Archive, digest>(ar,auxiliaryHashes);
                // archive<Archive, digest>(ar,contentHashes);
                // archive<Archive, std::vector<std::vector<FieldT>>>(ar,contents);
                ar & randomnessHashes;
                ar & auxiliaryHashes;
                ar & contentHashes;
                ar & contents;
                // ar & a;
            }
    };

    /* Hashing Functions */
    template<typename FieldT>
    using HashLeafContentF = std::function<digest(const FieldT *, const size_t, const size_t)>;

    typedef std::function<digest(const digest&, const digest&, size_t)> HashInnerNodeF;
    typedef std::function<digest(const uint8_t *,const size_t, const size_t)> HashZKRandomnessF;

    template<typename FieldT>
    using HashIntegerSeedGeneratorF = std::function<std::vector<unsigned char>(const hash::digest &)>;
    using HashStringSeedGeneratorF = std::function<std::vector<unsigned char>(const char *data, const size_t length)>;

/* 
    Wrappers for Hash Functions 
    */

    // For now our random generator is limited to 64 bits number
    template<typename FieldT>
    std::vector<unsigned char> blake2bIntegerSeedGenerator(const hash::digest &root)
    {
        std::vector<unsigned char> result(iSeedSize);
        size_t scaledResult;
        size_t index = 0;

        if (crypto_generichash_blake2b(
            (unsigned char*)&result[0],
            sizeof(result),
            (unsigned char*)&root[0],
            root.size(),
            (unsigned char*)&index, 
            sizeof(index)) !=0) throw std::runtime_error("Libsodium Encryption Error Randomness Generator (crypto_generichash_blake2b)");

        return result;
    }

    std::vector<unsigned char> blake2bStringHash(const char *data, const size_t length)
    {
        std::vector<unsigned char> result(iSeedSize);

        if (crypto_generichash_blake2b(
            (unsigned char*)&result[0],
            iSeedSize,
            (result.empty() ? NULL : (unsigned char*)data),
            length,
            NULL, 0) != 0) throw std::runtime_error("Libsodium Encryption Error Field Element (crypto_generichash_blake2b)");

        return result;
    }

    template<typename FieldT>
    digest blake2bFieldElementHash(const FieldT *data, const size_t length, const size_t digestLengthBytes)
    {
        hash::digest result(digestLengthBytes, '-');

        if (crypto_generichash_blake2b(
            (unsigned char*)&result[0],
            digestLengthBytes,
            (result.empty() ? NULL : (unsigned char*)data),
            sizeof(FieldT) * length,
            NULL, 0) != 0) throw std::runtime_error("Libsodium Encryption Error Field Element (crypto_generichash_blake2b)");

        return result;
    }

    digest blake2bZKElementHash(const uint8_t *data, const size_t length, const size_t digestLengthBytes)
    {
        hash::digest result(digestLengthBytes, '-');

        if (crypto_generichash_blake2b(
            (unsigned char*)&result[0],
            digestLengthBytes,
            (unsigned char *)data,
            length * sizeof(uint8_t),
            NULL, 0) != 0) {
                std::cout << length * sizeof(uint8_t) << std::endl;
                throw std::runtime_error("Libsodium Encryption Error ZK Element (crypto_generichash_blake2b)");
            }
        return result;
    }

    digest blake2bTwoToOneHash(const hash::digest &first, const hash::digest &second, const size_t digestLengthBytes)
    {
        const hash::digest firstPlusSecond = first + second;
        hash::digest result(digestLengthBytes, '-');

        if (crypto_generichash_blake2b(
            (unsigned char*)&result[0],
            digestLengthBytes,
            (unsigned char*)&firstPlusSecond[0],
            firstPlusSecond.size(),
            NULL, 0)) throw std::runtime_error("Libsodium Encryption Error Inner Node (crypto_generichash_blake2b)");

        return result;
    }
} // namespace hash
