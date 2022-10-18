/*

compile with: g++ -std=c++17 -o main_test -g -O2 main_bucket_set_construction.cpp 

use the following code in command line to unleash the stack restriction:

ulimit -s unlimited

*/

#include <algorithm>
#include <set>
#include <array>
#include <iostream>

long omega2( long n){
    long rem = n % 2;
    long exponent = 0;
    while( rem == 0){
        exponent ++;
        n >>= 1;
        rem = n % 2;
    }
    return exponent;
}

long omega3( long n){
    long rem = n % 3;
    long exponent = 0;
    while( rem == 0){
        exponent ++;
        n = n/3;
        rem = n % 3;
    }
    return exponent;
}

// Correctness is checked in sagemath
void construct_bucket_set(long bucket_set[], long& bsize, const long q, const long ah){

    std::set<long> B = {0,1};
    
    for(long i = 2; i <= q/2; ++i){
        if (((omega2(i) + omega3(i))%2) == 0){
            B.insert(i);
        }
    }
    
    for(long i = q/4; i < q/2; ++i){
        if ((B.find(i) != B.end()) && (B.find(q - 2*i) != B.end()) ) // if i is in B and q-3*i is in B
        {
            B.erase(q - 2*i);
        }
    }
    for(long i = q/6; i < q/4; ++i){
        if ((B.find(i) != B.end()) && (B.find(q - 3*i) != B.end()) ) // if i is in B and q-3*i is in B
        {
            B.erase(q - 3*i);
        }
    }
    
    for(long i = 1; i <= ah+1; ++i){
        if (((omega2(i) + omega3(i))%2) == 0){
            B.insert(i);
        }
    }
    
    long index = 0;
    for(auto b: B){bucket_set[index] = b; ++index;};
    
    bsize = index;
}

bool check_bucket_set_validity(const long bucket_set[], const long bsize, const long q, const long ah){

   std::set<long> C = {0};

    long i;

    for(long j = 0; j < bsize; ++j){
        i = bucket_set[j];
        if ( i <= (ah + 1) ) C.insert(i);
        if ( 2*i <= (ah + 1) ) C.insert(2*i);
        if ( 3*i <= (ah + 1) ) C.insert(3*i);
    }

    if ((C.size() != size_t(ah + 2)) || (*C.rbegin() != (ah+1))){
        std::cout << "leading_c check false!" << std::endl;
        return 0;
    }

    C.clear();
    C.insert(0);

    for(long j = 0; j < bsize; ++j){
        i = bucket_set[j];
        C.insert(i);
        C.insert(q - i);
       
        if ( 2*i <= q ){ C.insert(2*i); C.insert(q - 2*i);}
        if ( 3*i <= q ){ C.insert(3*i); C.insert(q - 3*i);}
    }

    if ((C.size() != size_t(q + 1)) || (*C.rbegin() != q)){
        std::cout << C.size() << std::endl;
        std::cout << *C.rbegin() << std::endl;
        std::cout << "q check false!" << std::endl;
        return 0;
    }
    

    return 1;
}

long max_gap_in_bucket_set(const long bucket_set[], const long bsize){
    long tmp, d_max = 0;
    for(long j = 0; j < bsize - 1; ++j){
        tmp = bucket_set[j+1] -bucket_set[j];
        if(tmp > d_max) d_max = tmp;
    }
    return d_max;
}

// #include "ches_config_files/config_file_n_exp_16.h" 

int main(){
 
    long q_RADIX = (long) 1 << 31; 
    long a_LEADING_TERM = 115;

    long * bucket_set;
    bucket_set = new long [q_RADIX/4];

    long bsize;

    construct_bucket_set(bucket_set, bsize, q_RADIX, a_LEADING_TERM);

    std::cout << std::dec << bsize << std::endl;
    std::cout << max_gap_in_bucket_set(bucket_set, bsize) << std::endl;
    std::cout << float(bsize)/q_RADIX << std::endl;

    std::cout << check_bucket_set_validity(bucket_set, bsize, q_RADIX, a_LEADING_TERM) << std::endl;


    delete[] bucket_set;
    return 0;
}


