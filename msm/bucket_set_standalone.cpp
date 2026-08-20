/**

compile with: g++ -std=c++17 -o bucket_set_test -g -O2 msm/bucket_set_standalone.cpp
run with: ./bucket_set_test <exponent_of_q>
example: ./bucket_set_test 13 constructs the bucket set for q = 2^13

use the following code in command line to unleash the stack restriction:

ulimit -s unlimited

Standalone utility to construct and validate a bucket set for given q and ah.

-- Guiwen. Oct, 2022.

*/

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "../external/uint256/uint256.hpp"
#include "include/bucket_set.h"

namespace {

const uint256_t BLS12_381_SCALAR_FIELD_ORDER = {
    uint64_t(0xffffffff00000001), uint64_t(0x53bda402fffe5bfe),
    uint64_t(0x3339d80809a1d805), uint64_t(0x73eda753299d7d48)
};

} // namespace

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <exponent_of_q>" << std::endl;
        std::cerr << "Example: " << argv[0] << " 13  (q = 2^13)" << std::endl;
        return EXIT_FAILURE;
    }

    char* end = nullptr;
    const long requested_exponent = std::strtol(argv[1], &end, 10);
    if (argv[1][0] == '\0' || *end != '\0' || requested_exponent <= 0 ||
        requested_exponent >= std::numeric_limits<long>::digits) {
        std::cerr << "exponent_of_q must be an integer in [1, "
                  << (std::numeric_limits<long>::digits - 1) << "]"
                  << std::endl;
        return EXIT_FAILURE;
    }

    const long exponent_of_q = requested_exponent;
    const long q_RADIX = 1L << exponent_of_q;

    // h = ceil(log_q(r)). Since q = 2^exponent_of_q and r is not a
    // power of two, this is floor(log2(r) / exponent_of_q) + 1.
    const long h =
        static_cast<long>(BLS12_381_SCALAR_FIELD_ORDER.get_msb()) /
            exponent_of_q + 1;

    // r_{h-1} = floor(r / q^(h-1)).
    const long a_LEADING_TERM = static_cast<long>(
        BLS12_381_SCALAR_FIELD_ORDER >>
        uint256_t(exponent_of_q * (h - 1)));

    long* bucket_set;
    const long bucket_set_capacity =
        std::max(q_RADIX / 2, a_LEADING_TERM + 1) + 1;
    bucket_set = new long[bucket_set_capacity];

    long bsize;

    construct_bucket_set(bucket_set, bsize, q_RADIX, a_LEADING_TERM);

    std::cout << std::dec << bsize << std::endl;
    std::cout << max_gap_in_bucket_set(bucket_set, bsize) << std::endl;
    std::cout << float(bsize) / q_RADIX << std::endl;

    std::cout << check_bucket_set_validity(bucket_set, bsize, q_RADIX, a_LEADING_TERM) << std::endl;

    delete[] bucket_set;
    return 0;
}
