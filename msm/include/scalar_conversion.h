/**
 * @file    scalar_conversion.h
 * @brief   Scalar conversion between integer, q-ary, MB, and q-half representations.
 *
 * **Important**: This header must be included AFTER config/config_file.h and AFTER
 * the global `digit_decomposition* DIGIT_CONVERSION_HASH_TABLE` is declared.
 * The functions use constexpr values (h_LEN_SCALAR, EXPONENT_OF_q, etc.) and the
 * global DIGIT_CONVERSION_HASH_TABLE directly.
 *
 * @author  Guiwen Luo
 * @date    Oct 2022
 */

#ifndef MSM_SCALAR_CONVERSION_H
#define MSM_SCALAR_CONVERSION_H

#include <array>
#include <cstdint>

#include "external/uint256/uint256.hpp"

// digit_decomposition is defined in bindings/blst.h

/** MB expression: array of (m, b) pairs, length = h_LEN_SCALAR (from config). */
typedef std::array<std::array<int, 2>, h_LEN_SCALAR> scalar_MB_expr;

/* ---- Scalar conversion functions ---- */

/**
 * Convert a uint256_t scalar to standard q-ary representation.
 * Each digit is in [0, q-1], obtained by masking the lower EXPONENT_OF_q bits.
 */
inline void trans_uint256_t_to_standard_q_ary_expr(
    std::array<int, h_LEN_SCALAR>& ret_std_expr, const uint256_t& a)
{
    uint256_t tmp = a;
    uint32_t mask = (1 << EXPONENT_OF_q) - 1;
    for (int i = 0; i < h_LEN_SCALAR; ++i) {
        ret_std_expr[i] = tmp.data[0] & mask;
        tmp = tmp >> EXPONENT_OF_q;
    }
}

/**
 * Convert a uint256_t scalar to M×B radix-q expression (CHES method).
 * Uses DIGIT_CONVERSION_HASH_TABLE for the conversion.
 */
inline void trans_uint256_t_to_MB_radixq_expr(
    std::array<std::array<int, 2>, h_LEN_SCALAR>& ret_MB_expr, const uint256_t& a)
{
    // Convert to standard q-ary expression first
    std::array<int, h_LEN_SCALAR> tmp_std_expr;
    uint256_t tmp = a;
    uint32_t mask = (1 << EXPONENT_OF_q) - 1;

    for (int i = 0; i < h_LEN_SCALAR; ++i) {
        tmp_std_expr[i] = tmp.data[0] & mask;
        tmp = tmp >> EXPONENT_OF_q;
    }

    digit_decomposition tmp_tri;
    for (int i = 0; i < h_LEN_SCALAR; ++i) {
        tmp_tri = DIGIT_CONVERSION_HASH_TABLE[tmp_std_expr[i]];
        if (tmp_tri.alpha == 0) {
            ret_MB_expr[i][0] = tmp_tri.m;
            ret_MB_expr[i][1] = tmp_tri.b;
        } else {
            ret_MB_expr[i][0] = -tmp_tri.m;
            ret_MB_expr[i][1] = tmp_tri.b;
            tmp_std_expr[i + 1] += 1;
        }
    }
}

/**
 * Convert a uint256_t scalar to standard q-ary representation for BGMW95.
 */
inline void trans_uint256_t_to_standard_q_ary_expr_BGMW95(
    std::array<int, h_BGMW95>& ret_std_expr, const uint256_t& a)
{
    uint256_t tmp = a;
    uint32_t mask = (1 << EXPONENT_OF_q_BGMW95) - 1;
    for (int i = 0; i < h_BGMW95; ++i) {
        ret_std_expr[i] = tmp.data[0] & mask;
        tmp = tmp >> EXPONENT_OF_q_BGMW95;
    }
}

/**
 * Convert a uint256_t scalar to signed q-half representation for BGMW95.
 * Digits are in [-q/2, q/2].
 */
inline void trans_uint256_t_to_qhalf_expr(
    std::array<int, h_BGMW95>& ret_qhalf_expr, const uint256_t& a)
{
    uint256_t tmp = a;
    int qhalf = int(q_RADIX_PIPPENGER_VARIANT >> 1);
    uint32_t mask = uint32_t(q_RADIX_PIPPENGER_VARIANT - 1);
    for (int i = 0; i < h_BGMW95; ++i) {
        ret_qhalf_expr[i] = tmp.data[0] & mask;
        tmp = tmp >> EXPONENT_OF_q_BGMW95;
    }
    for (int i = 0; i < h_BGMW95 - 1; ++i) {
        if (ret_qhalf_expr[i] > qhalf) {
            ret_qhalf_expr[i] -= q_RADIX_PIPPENGER_VARIANT;
            ret_qhalf_expr[i + 1] += 1;
            // System parameters guarantee ret_qhalf_expr[h-1] <= q/2.
        }
    }
}

#endif /* MSM_SCALAR_CONVERSION_H */
