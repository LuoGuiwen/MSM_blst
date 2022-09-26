// Compile: g++ -std=c++17 -o main_test -g -O2 main_test_for_ubuntu.cpp libblst.a

#include "bindings/blst.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>

// using namespace std;

/* Define ostreams */
std::ostream& operator<<(std::ostream& os, const blst_fp& b)
{
    os << std::hex << b.l[0] << " "<< b.l[1]<<  " "<< b.l[2]<< " "<< b.l[3] <<" "<<  b.l[4] <<" "<<  b.l[5] <<" ";
    return os;
}


/* Define global variables */

const size_t N_POINTS = (size_t) (1 << 20);  // I cannot test 2**19 and 2**20, it will cause segmentation fault. 
const blst_p1 G1_BLST_DEFAULT_GENERATOR = *blst_p1_generator(); // Default generator in blst_p1/

blst_fr FR_ONE;
blst_fp FP_ONE;
blst_p1 G1_GENERATOR; 
blst_p1_affine G1_GENERATOR_AFFINE; 


void init(){

    // initialize the identity, i.e., one, in fr.
    uint64_t blst_fr_one_vec[] = {uint64_t(1),uint64_t(0),uint64_t(0),uint64_t(0)}; 
    blst_fr_from_uint64(&FR_ONE, blst_fr_one_vec); 

    // initialize the identity, i.e., one, in fp.
    uint64_t fp_one_vec[] = {uint64_t(0x1), uint64_t(0x0), uint64_t(0x0), uint64_t(0x0), uint64_t(0x0), uint64_t(0x0)};
    blst_fp_from_uint64(&FP_ONE, fp_one_vec);   

    /* 
    initialize the generator in blst_p1 by Guiwen. G1_GENERATOR = 11 * 10177 * 859267 * 52437899* (Point whose x-coordinate is 4).
    
    G1_GENERATOR =  
    {0x7f127a0a6f06434698a6b6598fc6d8bd7e8482362c69b416d8640c18c1caec0ab474874acad9e91be475966f7413a26, 
    0x5be03c7afc54a0b30376055f27a4ff60e8ca9060651b98fa6caa6937bed9116b52ad54fbc4e22cd69b8519cb9bfd662}
    */   
    uint64_t G1x_vec[] = {uint64_t(0xbe475966f7413a26), uint64_t(0xab474874acad9e91), uint64_t(0x6d8640c18c1caec0), uint64_t(0xd7e8482362c69b41), uint64_t(0x698a6b6598fc6d8b), uint64_t(0x07f127a0a6f06434)} ;
    uint64_t G1y_vec[] = {uint64_t(0x69b8519cb9bfd662), uint64_t(0xb52ad54fbc4e22cd), uint64_t(0xa6caa6937bed9116), uint64_t(0x0e8ca9060651b98f), uint64_t(0x30376055f27a4ff6), uint64_t(0x05be03c7afc54a0b)} ;
    blst_fp G1x, G1y;
    blst_fp_from_uint64(&G1x, G1x_vec);
    blst_fp_from_uint64(&G1y, G1y_vec);   
    G1_GENERATOR = {G1x, G1y, FP_ONE};  // integer a -> a << 384 mod p in montgomery form
    G1_GENERATOR_AFFINE = {G1x, G1y};
    std::cout << "Check G1_generator is in G1: " << blst_p1_in_g1(&G1_GENERATOR) <<std::endl;
}


/*  */

blst_scalar random_blst_scalar(){

    // credit to url: https://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library // this doesn't work on my computer
    // std::random_device rd;
    // std::mt19937 mt(rd());
    // std::uniform_real_distribution<int> dist(1, (int) 0xffff); 
    uint64_t scalar_vec[4];
    scalar_vec[0] = (uint64_t) rand();
    scalar_vec[1] = (uint64_t) rand();
    scalar_vec[2] = (uint64_t) rand();
    scalar_vec[3] = (uint64_t) 0xfffffff;
    blst_scalar scalar;
    blst_scalar_from_uint64(&scalar, scalar_vec);
    return scalar;
}

void test(){

    size_t nbits = 256;
    
    uint8_t* scalars[N_POINTS];

    for(size_t i = 0; i < N_POINTS; ++i){
            scalars[i] = random_blst_scalar().b;
    }

    blst_p1_affine** points;
    points = new blst_p1_affine*[N_POINTS];

    blst_p1 tmp_P, tmp2_P = G1_GENERATOR;
    blst_p1_affine tmp_P_affine;

    for(size_t i = 0; i < N_POINTS; ++i){
            blst_p1_double(&tmp_P, &tmp2_P);
            tmp2_P = tmp_P;
            blst_p1_to_affine(&tmp_P_affine, &tmp_P);
            points[i] = & tmp_P_affine;
        }

    blst_p1 ret_P; // Mont coordinates
    blst_p1_affine ret_P_affine;

    /*blst pippenger*/

    std::cout<< "blst pippenger test: " << std::endl;

    size_t TEST_NUM = 1;
    
    limb_t* scratch;
    scratch = new limb_t[blst_p1s_mult_pippenger_scratch_sizeof(N_POINTS)/sizeof(limb_t)];
    auto st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< TEST_NUM; ++i)
    blst_p1s_mult_pippenger(&ret_P, points, N_POINTS, scalars, nbits, scratch);
    auto ed = std::chrono::steady_clock::now();   
    delete[] scratch;
    std::chrono::microseconds diff = std::chrono::duration_cast<std::chrono::microseconds>(ed -st); 
    std::cout << "blst pippenger Wall clock time elapse is: " << std::dec << diff.count()/TEST_NUM << " us "<< std::endl;
    blst_p1_to_affine(&ret_P_affine, &ret_P);
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;

    std::cout << std::dec << N_POINTS <<std::endl;
    std::cout << "TEST END" <<std::endl;


}

int main(){

    init();
    test();

    return 0;
}

