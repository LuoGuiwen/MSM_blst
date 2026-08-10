/**
 * @file    bucket_set.h
 * @brief   Bucket set construction for the CHES MSM algorithm.
 *
 * Implements the bucket set construction described in the CHES paper.
 * The construction uses the property of 2- and 3-adic valuations (omega2, omega3)
 * to determine which integers belong to the bucket set.
 *
 * Two overloads of construct_bucket_set are provided:
 *  - int version (used by MSM benchmarks): fills a pre-sized array
 *  - long version (standalone utility): also returns the computed size
 *
 * @author  Guiwen Luo
 * @date    Oct 2022
 */

#ifndef MSM_BUCKET_SET_H
#define MSM_BUCKET_SET_H

#include <iostream>
#include <set>

/**
 * 2-adic valuation: returns the exponent of the highest power of 2 dividing n.
 */
inline int omega2(int n)
{
    int rem = n % 2;
    int exponent = 0;
    while (rem == 0) {
        exponent++;
        n >>= 1;
        rem = n % 2;
    }
    return exponent;
}

/**
 * 3-adic valuation: returns the exponent of the highest power of 3 dividing n.
 */
inline int omega3(int n)
{
    int rem = n % 3;
    int exponent = 0;
    while (rem == 0) {
        exponent++;
        n = n / 3;
        rem = n % 3;
    }
    return exponent;
}

/**
 * Construct the bucket set for the CHES MSM algorithm.
 *
 * Fills the pre-allocated bucket_set[] array (must be at least q/4 elements)
 * with the constructed values in ascending order.
 *
 * @param bucket_set  Output array, pre-allocated to sufficient size.
 * @param q           Radix (power of 2).
 * @param ah          Leading term parameter.
 */
inline void construct_bucket_set(int bucket_set[], int q, int ah)
{
    std::set<int> B = {0, 1};

    for (int i = 2; i <= q / 2; ++i) {
        if (((omega2(i) + omega3(i)) % 2) == 0) {
            B.insert(i);
        }
    }

    for (int i = q / 4; i < q / 2; ++i) {
        if ((B.find(i) != B.end()) && (B.find(q - 2 * i) != B.end())) {
            B.erase(q - 2 * i);
        }
    }
    for (int i = q / 6; i < q / 4; ++i) {
        if ((B.find(i) != B.end()) && (B.find(q - 3 * i) != B.end())) {
            B.erase(q - 3 * i);
        }
    }

    for (int i = 1; i <= ah + 1; ++i) {
        if (((omega2(i) + omega3(i)) % 2) == 0) {
            B.insert(i);
        }
    }

    int index = 0;
    for (auto b : B) {
        bucket_set[index] = b;
        ++index;
    }
}

/**
 * Construct the bucket set and return its actual size (long version).
 * Provided for standalone analysis / debugging.
 */
inline void construct_bucket_set(long bucket_set[], long& bsize, const long q, const long ah)
{
    std::set<long> B = {0, 1};

    for (long i = 2; i <= q / 2; ++i) {
        if (((omega2(i) + omega3(i)) % 2) == 0) {
            B.insert(i);
        }
    }

    for (long i = q / 4; i < q / 2; ++i) {
        if ((B.find(i) != B.end()) && (B.find(q - 2 * i) != B.end())) {
            B.erase(q - 2 * i);
        }
    }
    for (long i = q / 6; i < q / 4; ++i) {
        if ((B.find(i) != B.end()) && (B.find(q - 3 * i) != B.end())) {
            B.erase(q - 3 * i);
        }
    }

    for (long i = 1; i <= ah + 1; ++i) {
        if (((omega2(i) + omega3(i)) % 2) == 0) {
            B.insert(i);
        }
    }

    long index = 0;
    for (auto b : B) {
        bucket_set[index] = b;
        ++index;
    }

    bsize = index;
}

/**
 * Verify that a constructed bucket set satisfies the required covering properties.
 *
 * @return true if the bucket set is valid, false otherwise.
 */
inline bool check_bucket_set_validity(const long bucket_set[], const long bsize,
                                      const long q, const long ah)
{
    std::set<long> C = {0};
    long i;

    for (long j = 0; j < bsize; ++j) {
        i = bucket_set[j];
        if (i <= (ah + 1)) C.insert(i);
        if (2 * i <= (ah + 1)) C.insert(2 * i);
        if (3 * i <= (ah + 1)) C.insert(3 * i);
    }

    if ((C.size() != size_t(ah + 2)) || (*C.rbegin() != (ah + 1))) {
        std::cout << "leading_c check false!" << std::endl;
        return false;
    }

    C.clear();
    C.insert(0);

    for (long j = 0; j < bsize; ++j) {
        i = bucket_set[j];
        C.insert(i);
        C.insert(q - i);

        if (2 * i <= q) { C.insert(2 * i); C.insert(q - 2 * i); }
        if (3 * i <= q) { C.insert(3 * i); C.insert(q - 3 * i); }
    }

    if ((C.size() != size_t(q + 1)) || (*C.rbegin() != q)) {
        std::cout << C.size() << std::endl;
        std::cout << *C.rbegin() << std::endl;
        std::cout << "q check false!" << std::endl;
        return false;
    }

    return true;
}

/**
 * Compute the maximum gap between consecutive elements in the bucket set.
 */
inline long max_gap_in_bucket_set(const long bucket_set[], const long bsize)
{
    long tmp, d_max = 0;
    for (long j = 0; j < bsize - 1; ++j) {
        tmp = bucket_set[j + 1] - bucket_set[j];
        if (tmp > d_max) d_max = tmp;
    }
    return d_max;
}

#endif /* MSM_BUCKET_SET_H */
