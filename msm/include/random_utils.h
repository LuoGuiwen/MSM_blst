/**
 * @file    random_utils.h
 * @brief   Cryptographically secure random number generation using SHA256.
 *
 * Uses OpenSSL's SHA256 hash function and RAND_bytes to generate
 * uniformly distributed random scalars less than the BLS12-381 group order r.
 *
 * @author  Guiwen Luo
 * @date    Oct 2022
 */

#ifndef MSM_RANDOM_UTILS_H
#define MSM_RANDOM_UTILS_H

#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>

#include <openssl/rand.h>
#include <openssl/sha.h>

#include "external/uint256/uint256.hpp"

/* ---- Group order of BLS12-381 ---- */
const uint256_t r_GROUP_ORDER = {
    uint64_t(0xffffffff00000001), uint64_t(0x53bda402fffe5bfe),
    uint64_t(0x3339d80809a1d805), uint64_t(0x73eda753299d7d48)
};

/**
 * Generate a 255-bit random scalar less than r using SHA256-based RNG.
 * Suitable for experimental / benchmarking purposes.
 */
inline uint256_t random_scalar_less_than_r_SHA256()
{
    uint256_t ret = r_GROUP_ORDER;

    const size_t BUFFER_LENGTH = 32;  // 256 bits
    unsigned char buffer[BUFFER_LENGTH];
    uint64_t randomNumbers[4];

    while (ret >= r_GROUP_ORDER) {
        unsigned char hash[SHA256_DIGEST_LENGTH];

        // Initialize the buffer with random data
        RAND_bytes(buffer, sizeof(buffer));

        // Compute the SHA256 hash of the buffer
        SHA256(buffer, sizeof(buffer), hash);

        // Store the hash in the array of uint64_t
        memcpy(randomNumbers, hash, sizeof(uint64_t) * 4);

        ret.data[3] = randomNumbers[3] >> 1;
        ret.data[2] = randomNumbers[2];
        ret.data[1] = randomNumbers[1];
        ret.data[0] = randomNumbers[0];
    }

    return ret;
}

/** Generate a reproducible benchmark scalar from a monotonically increasing counter. */
inline uint256_t deterministic_scalar_less_than_r_SHA256(uint64_t counter)
{
    unsigned char input[16];
    unsigned char hash[SHA256_DIGEST_LENGTH];

    for (uint64_t attempt = 0;; ++attempt) {
        for (int i = 0; i < 8; ++i) {
            input[i] = static_cast<unsigned char>(counter >> (8 * i));
            input[8 + i] = static_cast<unsigned char>(attempt >> (8 * i));
        }
        SHA256(input, sizeof(input), hash);

        uint256_t ret = 0;
        for (int limb = 0; limb < 4; ++limb) {
            uint64_t value = 0;
            for (int byte_index = 0; byte_index < 8; ++byte_index)
                value |= uint64_t(hash[8 * limb + byte_index]) << (8 * byte_index);
            ret.data[limb] = value;
        }
        ret.data[3] >>= 1;
        if (ret < r_GROUP_ORDER)
            return ret;
    }
}

/**
 * Generate 4×64-bit random numbers via SHA256(rand_bytes).
 */
inline void generateRandomNumber(uint64_t* numbers)
{
    const size_t BUFFER_LENGTH = 32;  // 256 bits
    unsigned char buffer[BUFFER_LENGTH];
    unsigned char hash[SHA256_DIGEST_LENGTH];

    // Initialize the buffer with random data
    RAND_bytes(buffer, sizeof(buffer));

    // Compute the SHA256 hash of the buffer
    SHA256(buffer, sizeof(buffer), hash);

    // Store the hash in the array of uint64_t
    memcpy(numbers, hash, sizeof(uint64_t) * 4);
}

/**
 * Print a 256-bit random number in hex.
 */
inline void printRandomNumber(const uint64_t* numbers)
{
    std::cout << "Random Number: ";
    for (int i = 0; i < 4; ++i) {
        std::cout << std::hex << std::setw(16) << std::setfill('0') << numbers[i] << " ";
    }
    std::cout << std::dec << std::endl;
}

/**
 * Generate a random scalar less than r using Mersenne Twister (non-cryptographic).
 * Only for experimental comparison — prefer random_scalar_less_than_r_SHA256.
 */
inline uint256_t random_scalar_less_than_r()
{
    uint256_t ret = r_GROUP_ORDER;

    std::random_device rd;
    std::mt19937 gen(rd());

    // Generate pseudo-random uint32_t uniformly distributed in [0, 2^32-1]
    std::uniform_int_distribution<> dis(0, 0xffffffff);

    while (ret >= r_GROUP_ORDER) { // ensure scalar < r
        ret.data[3] = (uint64_t(dis(gen)) >> 1) + (uint64_t(dis(gen)) << 31);  // 255 bit
        ret.data[2] = uint64_t(dis(gen)) + (uint64_t(dis(gen)) << 32);
        ret.data[1] = uint64_t(dis(gen)) + (uint64_t(dis(gen)) << 32);
        ret.data[0] = uint64_t(dis(gen)) + (uint64_t(dis(gen)) << 32);
    }

    return ret;
}

#endif /* MSM_RANDOM_UTILS_H */
