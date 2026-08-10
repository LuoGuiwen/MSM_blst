/**
 * @file    msm_common.h
 * @brief   Master include for the MSM library — types, I/O, utilities, and
 *          sub-module headers.
 *
 * Include order matters:
 *   1. blst bindings
 *   2. config/config_file.h
 *   3. Declare global `digit_decomposition* DIGIT_CONVERSION_HASH_TABLE`
 *   4. This header (msm_common.h)
 *
 * This header replaces the old monolithic `auxiliaryfunc.h`.
 *
 * @author  Guiwen Luo
 * @date    Oct 2022
 */

#ifndef MSM_COMMON_H
#define MSM_COMMON_H

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>

#include "bindings/blst.h"
#include "external/uint256/uint256.hpp"
#include "msm/include/bucket_index.h"
#include "msm/include/bucket_set.h"
#include "msm/include/random_utils.h"
#include "msm/include/scalar_conversion.h"

/* ===================================================================
 *  Ostream operators for blst types
 * =================================================================== */

inline std::ostream& operator<<(std::ostream& os, const blst_fp& b)
{
    os << std::hex << "0x" << std::setfill('0')
       << std::setw(16) << b.l[5] << " " << std::setw(16) << b.l[4]
       << " " << std::setw(16) << b.l[3] << " "
       << std::setw(16) << b.l[2] << " " << std::setw(16) << b.l[1]
       << " " << std::setw(16) << b.l[0];
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const blst_fp2& b)
{
    os << b.fp[0] << " + " << std::endl << b.fp[1] << "*i  ";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const blst_fr& b)
{
    os << std::hex << b.l[3] << " " << b.l[2] << " " << b.l[1] << " " << b.l[0] << " ";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const uint8_t bt)
{
    return os << unsigned(bt);
}

inline std::ostream& operator<<(std::ostream& os, const blst_scalar& scalar)
{
    os << std::hex << "0x" << " " << std::setfill('0')
       << std::setw(2) << scalar.b[31] << " " << std::setw(2) << scalar.b[30] << " "
       << std::setw(2) << scalar.b[29] << " " << std::setw(2) << scalar.b[28] << " "
       << std::setw(2) << scalar.b[27] << " " << std::setw(2) << scalar.b[26] << " "
       << std::setw(2) << scalar.b[25] << " " << std::setw(2) << scalar.b[24] << ", "
       << std::setw(2) << scalar.b[23] << ' ' << std::setw(2) << scalar.b[22] << ' '
       << std::setw(2) << scalar.b[21] << " " << std::setw(2) << scalar.b[20] << " "
       << std::setw(2) << scalar.b[19] << " " << std::setw(2) << scalar.b[18] << " "
       << std::setw(2) << scalar.b[17] << " " << std::setw(2) << scalar.b[16] << ", "
       << std::setw(2) << scalar.b[15] << ' ' << std::setw(2) << scalar.b[14] << ' '
       << std::setw(2) << scalar.b[13] << " " << std::setw(2) << scalar.b[12] << " "
       << std::setw(2) << scalar.b[11] << " " << std::setw(2) << scalar.b[10] << " "
       << std::setw(2) << scalar.b[9]  << " " << std::setw(2) << scalar.b[8]  << ", "
       << std::setw(2) << scalar.b[7]  << ' ' << std::setw(2) << scalar.b[6]  << ' '
       << std::setw(2) << scalar.b[5]  << " " << std::setw(2) << scalar.b[4]  << " "
       << std::setw(2) << scalar.b[3]  << " " << std::setw(2) << scalar.b[2]  << " "
       << std::setw(2) << scalar.b[1]  << " " << std::setw(2) << scalar.b[0]  << " ";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const uint64_t* b)
{
    os << std::hex << b[3] << " " << b[2] << " " << b[1] << " " << b[0] << " ";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const digit_decomposition b)
{
    os << std::dec << "[" << b.m << " " << b.b << " " << b.alpha << "]";
    return os;
}

/* ===================================================================
 *  Utility functions
 * =================================================================== */

/** Prefetch hint (no-op on most platforms but prevents compiler reordering). */
inline void vec_prefetch(const void* ptr, size_t len)
{
    (void)ptr;
    (void)len;
}

/** Convert a uint32_t to a 4-byte little-endian representation. */
inline void byte_str_from_uint32(uint8_t ret[4], const uint32_t a)
{
    uint32_t w = a;
    ret[0] = (uint8_t)w;
    ret[1] = (uint8_t)(w >> 8);
    ret[2] = (uint8_t)(w >> 16);
    ret[3] = (uint8_t)(w >> 24);
}

/** Zero-initialize a memory region (volatile to prevent optimization). */
inline void vec_zero(void* ret, size_t num)
{
    volatile limb_t* rp = (volatile limb_t*)ret;
    size_t i;

    num /= sizeof(limb_t);

    for (i = 0; i < num; i++)
        rp[i] = 0;
}

/* ---- Functions from blst library (used by the MSM implementation) ---- */

/** Determine the optimal Pippenger window size for a given number of points. */
inline size_t pippenger_window_size(size_t npoints)
{
    size_t wbits;

    for (wbits = 0; npoints >>= 1; wbits++)
        ;

    return wbits > 12   ? wbits - 3
           : (wbits > 4 ? wbits - 2
                        : (wbits ? 2 : 1));
}

/**
 * Booth encoding: window value encoding using the fact that -P is trivially
 * calculated, halving the size of the pre-computed table.
 * Attributed to A. D. Booth.
 */
inline limb_t booth_encode(limb_t wval, size_t sz)
{
    limb_t mask = 0 - (wval >> sz); /* "sign" bit -> mask */

    wval = (wval + 1) >> 1;
    wval = (wval & ~mask) | ((0 - wval) & mask);

    /* &0x1f, but <=0x10, is index in table, rest is extended "sign" bit */
    return wval;
}

/** Extract a window value from a byte array (works up to 25 bits). */
inline limb_t get_wval_limb(const uint8_t* d, size_t off, size_t bits)
{
    size_t i, top = (off + bits - 1) / 8;
    limb_t ret, mask = (limb_t)0 - 1;

    d += off / 8;
    top -= off / 8 - 1;

    for (ret = 0, i = 0; i < 4;) {
        ret |= (*d & mask) << (8 * i);
        mask = (limb_t)0 - ((++i - top) >> (8 * sizeof(top) - 1));
        d += 1 & mask;
    }

    return ret >> (off % 8);
}

#endif /* MSM_COMMON_H */
