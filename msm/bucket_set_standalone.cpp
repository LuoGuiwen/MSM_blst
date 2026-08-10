/**

compile with: g++ -std=c++17 -o bucket_set_test -g -O2 msm/bucket_set_standalone.cpp

use the following code in command line to unleash the stack restriction:

ulimit -s unlimited

Standalone utility to construct and validate a bucket set for given q and ah.

-- Guiwen. Oct, 2022.

*/

#include <iostream>
#include "msm/include/bucket_set.h"

int main()
{
    long q_RADIX = (long)1 << 31;
    long a_LEADING_TERM = 115;

    long* bucket_set;
    bucket_set = new long[q_RADIX / 4];

    long bsize;

    construct_bucket_set(bucket_set, bsize, q_RADIX, a_LEADING_TERM);

    std::cout << std::dec << bsize << std::endl;
    std::cout << max_gap_in_bucket_set(bucket_set, bsize) << std::endl;
    std::cout << float(bsize) / q_RADIX << std::endl;

    std::cout << check_bucket_set_validity(bucket_set, bsize, q_RADIX, a_LEADING_TERM) << std::endl;

    delete[] bucket_set;
    return 0;
}
