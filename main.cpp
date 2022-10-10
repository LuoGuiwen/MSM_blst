/*

compile with: g++ -std=c++17 -o main_test -g -O2 main.cpp libblst.a

use the following code in command line to unleash the stack restriction:

ulimit -s unlimited


All functions must be invoked after init().
*/


#include "bindings/blst.h"
// #include "src/consts.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <iomanip>

#include <set>
#include <array>

// using namespace std; // don't use this line, since there is a uint_8 and byte definition that are not compatible with std.


/***----***/
/* Define global variables and their initialization*/
/***----***/


#include "config_radix13.h" //define configuration in a seperate file.

/*
One only needs to change this config_radixnn.h file to do other computation.
*/



#include "auxiliaryfunc.h"

std::set<int> MULTI_SET = {1, 2, 3};

/* This once leads to a big BUG. */
std::array<int, q_RADIX/2 +1> BUCKET_VALUE_TO_ITS_INDEX; //Consider the maximum element in the BUCKET_SET, which is very very important

std::array<std::array<int,3>, q_RADIX+1>  DIGIT_CONVERSION_HASH_TABLE;

std::array< blst_p1_affine, N_POINTS> *FIX_POINTS_LIST; 
std::array< std::array<blst_p1_affine, h_LEN_SCALAR >, N_POINTS> *PRECOMPUTATION_POINTS_LIST_nh; // Define the pointer then use new to allocate memory in heap, since this array is too big.
std::array< std::array<blst_p1_affine, h_LEN_SCALAR >, N_POINTS> *PRECOMPUTATION_POINTS_LIST_nh_2;
std::array< std::array<blst_p1_affine, h_LEN_SCALAR >, N_POINTS> *PRECOMPUTATION_POINTS_LIST_nh_3;

uint256_t* SCALARS_ARRAY;

// std::array<blst_p1_affine, N_POINTS*h_LEN_SCALAR> *POINT_INTERMEDIATE;

/* initialization later on using init() */
blst_fr FR_ONE;
blst_fp FP_ONE, FP_MONT_ONE, FP_Z;

blst_p1 G1_GENERATOR, G1_INFINITY; 
blst_p1_affine G1_GENERATOR_AFFINE, G1_AFFINE_INFINITY; 

typedef std::array<std::array< int, 2>, h_LEN_SCALAR> scalar_MB_expr;
void print( const scalar_MB_expr &expr){

    std::cout << std::dec << "{";
    for(auto a : expr){
        std::cout <<"[ "<<a[0]<<", "<<a[1]<<" ]";
    }
    std::cout <<"}"<< std::endl;
}

void print( const blst_p1 &P){

    std::cout << std::hex << "[";
    std::cout << P.x <<", "<<std::endl;
    std::cout << P.y <<", "<<std::endl;
    std::cout << P.z <<std::endl;
    std::cout <<"]"<< std::endl;
}

int omega2( int n){
    int rem = n % 2;
    int exponent = 0;
    while( rem == 0){
        exponent ++;
        n >>= 1;
        rem = n % 2;
    }
    return exponent;
}

int omega3( int n){
    int rem = n % 3;
    int exponent = 0;
    while( rem == 0){
        exponent ++;
        n = n/3;
        rem = n % 3;
    }
    return exponent;
}


// Correctness is checked in sagemath
std::set<int> construct_bucket_set( int q, int ah){

    std::set<int> B = {0, 1};
    
    for(int i = 2; i <= q/2; ++i){
        if (((omega2(i) + omega3(i))%2) == 0){
            B.insert(i);
        }
    }
    
    for(int i = 1; i < (q -ah)/2; ++i){
        if ((B.find(i) != B.end()) && (B.find(q - 2*i) != B.end()) ) // if i is in B and q-3*i is in B
        {
            B.erase(q - 2*i);
        }
    }
    for(int i = 1; i < q/4; ++i){
        if ((B.find(i) != B.end()) && (B.find(q - 3*i) != B.end()) ) // if i is in B and q-3*i is in B
        {
            B.erase(q - 3*i);
        }
    }

    return B;
}

void init(){

    // initialize the identity, i.e., one, in fr.
    uint64_t blst_fr_one_vec[] = {uint64_t(1),uint64_t(0),uint64_t(0),uint64_t(0)}; 
    blst_fr_from_uint64(&FR_ONE, blst_fr_one_vec); 

    // initialize the identity, i.e., one, in fp.
    uint64_t fp_one_vec[] = {uint64_t(0x1), uint64_t(0x0), uint64_t(0x0), uint64_t(0x0), uint64_t(0x0), uint64_t(0x0)};
    blst_fp_from_uint64(&FP_ONE, fp_one_vec);   


    uint64_t fp_mont_one_vec[] ={uint64_t(0x760900000002fffd), uint64_t(0xebf4000bc40c0002), uint64_t(0x5f48985753c758ba), uint64_t(0x77ce585370525745), uint64_t(0x5c071a97a256ec6d), uint64_t(0x15f65ec3fa80e493)};
    blst_fp_from_uint64(&FP_MONT_ONE, fp_mont_one_vec);   

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
    // G1_GENERATOR = {G1x, G1y, FP_ONE};  // integer a -> a << 384 mod p in montgomery form
    G1_GENERATOR = {G1x, G1y, FP_ONE}; 
    G1_GENERATOR_AFFINE = {G1x, G1y};
    std::cout << "Check G1_generator is in G1: " << blst_p1_in_g1(&G1_GENERATOR) <<std::endl;

    FP_Z = {0,0,0,0,0,0};
     // We set the value of x equal to modulus to represent inifinty
    G1_INFINITY = { FP_Z, FP_ONE, FP_Z} ;
    blst_p1_to_affine(&G1_AFFINE_INFINITY, &G1_INFINITY);

    std::cout << "Check G1_INFINITY is in G1: " << blst_p1_in_g1(&G1_INFINITY) <<std::endl;
    std::cout << "Check G1_INFINITY is 0 in G1: " << blst_p1_is_inf(&G1_INFINITY) <<std::endl;

    std::cout << "Check G1_AFFINE_INFINITY is in G1: " << blst_p1_affine_in_g1(&G1_AFFINE_INFINITY) <<std::endl;
    std::cout << "Check G1_AFFINE_INFINITY is 0 in G1: " << blst_p1_affine_is_inf(&G1_AFFINE_INFINITY) <<std::endl;


    //Initialize BUCKET_SET

    // BUCKET_SET = construct_bucket_set(q_RADIX, ah_LEADING);
    std::cout<< "BUCKET_SET directly read from file. The size of BUCKET_SET is: " << B_SIZE << std::endl;

    for(size_t i = 0; i < B_SIZE; ++i){
        BUCKET_VALUE_TO_ITS_INDEX[BUCKET_SET[i]] = i;
    }
    
    // Initialize  DIGIT_CONVERSION_HASH_TABLE;
    
    for(int m: MULTI_SET){
        for (int b: BUCKET_SET){
            if (m*b <= q_RADIX) DIGIT_CONVERSION_HASH_TABLE[q_RADIX - m*b] = { m, b, 1};
        }
    }
    for(int m: MULTI_SET){
        for (int b: BUCKET_SET){
            if (m*b <= q_RADIX) DIGIT_CONVERSION_HASH_TABLE[m*b] = { m, b, 0};
        }
    }
    std::cout<< "Example element in DIGIT_CONVERSION_HASH_TABLE: " << DIGIT_CONVERSION_HASH_TABLE[6577][0] <<" " << DIGIT_CONVERSION_HASH_TABLE[6577][1] <<" " <<DIGIT_CONVERSION_HASH_TABLE[6577][2] <<" " <<std::endl;

    // Initialize *FIX_POINTS_LIST

    FIX_POINTS_LIST = new std::array< blst_p1_affine, N_POINTS>;
    blst_p1 tmp_P, tmp2_P = G1_GENERATOR;
    blst_p1_affine tmp_P_affine;

    for(size_t i = 0; i < N_POINTS; ++i){
            blst_p1_double(&tmp_P, &tmp2_P);
            tmp2_P = tmp_P;  // carefully use the assign operator, it may not well-defined.
            blst_p1_to_affine(&tmp_P_affine, &tmp_P);
            // std::memcpy(&points[i], &tmp_P_affine, sizeof(blst_p1_affine));
            (*FIX_POINTS_LIST)[i] = tmp_P_affine;
        }
    std::cout<< "FIX_POINTS_LIST Generated" <<std::endl;

    // POINT_INTERMEDIATE = new std::array< blst_p1_affine, N_POINTS*h_LEN_SCALAR>;
    // std::cout<< "POINT_INTERMEDIATE allocated" <<std::endl;

    // Initialize SCALARS_ARRAY
    SCALARS_ARRAY = new uint256_t[N_POINTS];
    for(size_t i = 0; i < N_POINTS; ++i)\
        SCALARS_ARRAY[i] = random_scalar_less_than_r();
    std::cout<< "SCALARS_ARRAY Generated" <<std::endl;
}


blst_p1_affine single_scalar_multiplication(uint256_t scalar, blst_p1_affine Q){


    blst_p1_affine aret; 
    blst_p1 ret = G1_INFINITY; // ret = G1_INFINITY; 
    // blst_p1 xyzQ = {Q.x, Q.y, FP_ONE}; // Don't initialize with  xyzQ = {Q.x, Q.y, FP_MONT_ONE}, I don't know why but it will lead to wrong result.

    blst_p1 xyzQ;
    blst_p1_from_affine(&xyzQ, &Q);

    while (scalar > 0){
        if ( scalar.data[0] & 1 ){
            blst_p1_add_or_double(&ret, &ret, &xyzQ);
        }
        
        blst_p1_add_or_double(&xyzQ, &xyzQ, &xyzQ);    // tested. No need to use temp variable.
        scalar = scalar >> 1; 
    }

    blst_p1_to_affine(&aret, &ret);
    return aret;
} 


blst_p1 single_scalar_multiplication(uint256_t scalar, blst_p1 Q){

    blst_p1 ret = G1_INFINITY; // ret = G1_INFINITY; 

    while (scalar > 0){
        if ( scalar.data[0] & 1 ){
            blst_p1_add_or_double(&ret, &ret, &Q);
        }
        
        blst_p1_add_or_double(&Q, &Q, &Q);    // tested. No need to use temp variable.
        scalar = scalar >> 1; 
    }
    
    return ret;
} 


blst_p1_affine trivial_mult_scalar_multiplication_2(){
    blst_p1 ret1 = G1_INFINITY;
    blst_p1 ret2 = G1_INFINITY;
    blst_p1 tmp;
    blst_scalar scalar;

    for(size_t i = 0; i < N_POINTS; ++i){
        blst_p1_from_affine(&tmp, &(*FIX_POINTS_LIST)[i]);
        blst_scalar_from_uint32( &scalar, SCALARS_ARRAY[i].data);
        blst_p1_mult( &ret1, &tmp, scalar.b, (size_t) 255);
        blst_p1_add_or_double(&ret2, &ret2, &ret1);
        }    
    blst_p1_affine ret_a;
    blst_p1_to_affine(&ret_a, &ret2);
    return ret_a;           
}

void trivial_mult_scalar_multiplication(blst_p1 *ret, const blst_p1_affine *const points[], size_t npoints, const byte *const scalars[] , size_t nbits){
    blst_p1 ret1 = G1_INFINITY;
    blst_p1 ret2 = G1_INFINITY;
    blst_p1 tmp_P;

    for(size_t i = 0; i < npoints; ++i){

        blst_p1_from_affine(&tmp_P, points[i]);
        blst_p1_mult( &ret1, &tmp_P, scalars[i], nbits);
        blst_p1_add_or_double(&ret2, &ret2, &ret1);
        }          

    ret = &ret2;                 
}


void init_PRECOMPUTATION_POINTS_LIST_nh(){
    /*
    ### Initialize the precomputation ###
    */
    PRECOMPUTATION_POINTS_LIST_nh   = new std::array< std::array< blst_p1_affine, h_LEN_SCALAR >, N_POINTS>;
    PRECOMPUTATION_POINTS_LIST_nh_2 = new std::array< std::array< blst_p1_affine, h_LEN_SCALAR >, N_POINTS>;
    PRECOMPUTATION_POINTS_LIST_nh_3 = new std::array< std::array< blst_p1_affine, h_LEN_SCALAR >, N_POINTS>;

    auto st = std::chrono::steady_clock::now();
    blst_p1_affine Pt;

    for(int i = 0; i< N_POINTS; ++i){
        blst_p1_affine qjQi = (*FIX_POINTS_LIST)[i];
        for(uint j = 0; j< h_LEN_SCALAR; ++j){
            (*PRECOMPUTATION_POINTS_LIST_nh)[i][j] = qjQi;
            (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j] = single_scalar_multiplication(2, qjQi);
            (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j] = single_scalar_multiplication(3, qjQi);
            // auto tmp_Pa = (*PRECOMPUTATION_POINTS_LIST_nh)[i][j];
            // blst_fp_sub(&tmp_Pa.y, &FP_Z, &tmp_Pa.y);
            // (*PRECOMPUTATION_POINTS_LIST_nh_n1)[i][j] = tmp_Pa;

            // tmp_Pa = (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j];
            // blst_fp_sub(&tmp_Pa.y, &FP_Z, &tmp_Pa.y);
            // (*PRECOMPUTATION_POINTS_LIST_nh_n2)[i][j] = tmp_Pa;

            // tmp_Pa = (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j];
            // blst_fp_sub(&tmp_Pa.y, &FP_Z, &tmp_Pa.y);
            // (*PRECOMPUTATION_POINTS_LIST_nh_n3)[i][j] = tmp_Pa;

            qjQi = single_scalar_multiplication(q_RADIX, qjQi);    
        }
    }
    auto ed = std::chrono::steady_clock::now();   

    std::chrono::microseconds diff = std::chrono::duration_cast<std::chrono::microseconds>(ed -st);
    std::cout<< "PRECOMPUTATION_POINTS_LIST_nh SUCCESSFULLY CONSTRUCTED" << std::endl;
    std::cout << "PRECOMPUTATION Wall clock time elapse is: " << diff.count() << " us "<< std::endl;
}

/*  */
void trans_uint256_t_to_standard_q_ary_expr( std::array<int, h_LEN_SCALAR> &ret_std_expr, const uint256_t &a){
    uint256_t tmp = a;
    uint32_t mask = (1 << EXPONENT_OF_q) - 1;
    for (int i=0; i< h_LEN_SCALAR; ++i){
        ret_std_expr[i] = tmp.data[0] & mask;// we only need the bit-wise xor with the last 32-bit of tmp.
        tmp = tmp >> EXPONENT_OF_q;
    }
}

void trans_uint256_t_to_MB_radixq_expr(std::array<std::array< int, 2>, h_LEN_SCALAR>  &ret_MB_expr, const uint256_t &a){
    
    // convert to the standard q ary expression first
    std::array<int, h_LEN_SCALAR> tmp_std_expr;
    uint256_t tmp = a;
    uint32_t mask = (1<<EXPONENT_OF_q) -1;

    for (int i=0; i< h_LEN_SCALAR; ++i){
        tmp_std_expr[i] = tmp.data[0] & mask;// we only need the bit-wise xor with the last 32-bit of tmp.
        tmp = tmp >> EXPONENT_OF_q;
    }

    std::array<int, 3> tmp_tri;
    for (int i = 0; i< h_LEN_SCALAR -1; ++i){
        tmp_tri = DIGIT_CONVERSION_HASH_TABLE[tmp_std_expr[i]];
        if(tmp_tri[2]==0){
            ret_MB_expr[i][0] = tmp_tri[0];
            ret_MB_expr[i][1] = tmp_tri[1];
        }

        else{
            ret_MB_expr[i][0] = - tmp_tri[0];
            ret_MB_expr[i][1] = tmp_tri[1];
            tmp_std_expr[i+1] += 1; 
        }

        tmp_tri = DIGIT_CONVERSION_HASH_TABLE[tmp_std_expr[h_LEN_SCALAR-1]];
        ret_MB_expr[h_LEN_SCALAR - 1][0] = tmp_tri[0];
        ret_MB_expr[h_LEN_SCALAR - 1][1] = tmp_tri[1];
    }

}

/*  divmod_2008_08654 */
void divmod_2008_08654( uint256_t &q, uint32_t &r, const uint256_t &scalar, const int &m, const int &b, const uint32_t &c, uint32_t &mask){
    uint256_t v1 = scalar + c; // It is better to choose c = 1 if possible.
    q = v1 >> b;
    for (int i = 0; i < m; ++i){
        q = ((q + v1) >> b);
    }
    uint256_t tmp = scalar + q;
    r = tmp.data[0] & mask;
}

void trans_uint256_t_to_MB_radixq_expr_2008_08654(std::array<std::array< int, 2>, h_LEN_SCALAR> &ret_MB_expr, const uint256_t &a, uint32_t &mask){
    
    // convert to the standard q ary expression first
    std::array<int, h_LEN_SCALAR> tmp_std_expr;
    uint256_t tmp = a;

    uint256_t quotient;
    uint32_t reminder;
    for (int i=0; i< h_LEN_SCALAR; ++i){
        divmod_2008_08654( quotient, reminder, tmp , 12 - i, EXPONENT_OF_q, 5, mask);
        tmp = quotient;
        tmp_std_expr[i] = reminder;// we only need the bit-wise xor with the last 32-bit of tmp.
    }

    std::array<int, 3> tmp_tri;
    for (int i = 0; i< h_LEN_SCALAR -1; ++i){
        tmp_tri = DIGIT_CONVERSION_HASH_TABLE[tmp_std_expr[i]];
        if(tmp_tri[2]==0){
            ret_MB_expr[i][0] = tmp_tri[0];
            ret_MB_expr[i][1] = tmp_tri[1];
        }

        else{
            ret_MB_expr[i][0] = - tmp_tri[0];
            ret_MB_expr[i][1] = tmp_tri[1];
            tmp_std_expr[i+1] += 1; 
        }

        tmp_tri = DIGIT_CONVERSION_HASH_TABLE[tmp_std_expr[h_LEN_SCALAR-1]];
        ret_MB_expr[h_LEN_SCALAR - 1][0] = tmp_tri[0];
        ret_MB_expr[h_LEN_SCALAR - 1][1] = tmp_tri[1];
    }

}

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



// The accumulation algorithms depicted in our CHES submission

blst_p1 accumulation_consecutive( const std::vector<int> & consecutiveBucketSetList, const std::vector<blst_p1> & intermediateSumList){
    blst_p1 tmp = G1_INFINITY;
    blst_p1 tmp1 = G1_INFINITY;

    for(auto i = consecutiveBucketSetList.size() -1; i > 0; --i ){
        blst_p1_add_or_double(&tmp, &tmp, &intermediateSumList[i]);
        blst_p1_add_or_double(&tmp1, &tmp1, &tmp);
    }
    return tmp1;
}


//
blst_p1 accumulation_d(int bucketSetList[], const std::vector<blst_p1>& intermediateSumList, const size_t bucket_set_size, const int d_max){

    blst_p1 tmp = G1_INFINITY;

    //# We don't use the tmp_d[0]
    std::vector<blst_p1> tmp_d (d_max + 1,  G1_INFINITY); 
    // for(uint i = 0; i<= max_d; ++i) tmp_d.push_back(tmp); // those two initializations are equivalent

    for( auto i = bucket_set_size -1; i>0; --i){
        blst_p1_add_or_double(&tmp, &tmp, &(intermediateSumList[i]));
        int differ = bucketSetList[i] -  bucketSetList[i-1];
        blst_p1_add_or_double(&(tmp_d[differ]), &(tmp_d[differ]), &tmp);
    }


    tmp = G1_INFINITY;
    blst_p1 tmp1 = G1_INFINITY;

    for(int i = d_max; i > 0; --i ){
        blst_p1_add_or_double(&tmp, &tmp, &(tmp_d[i]));
        blst_p1_add_or_double(&tmp1, &tmp1, &tmp);
    }

    return tmp1;    
}





 /* Calculate sum of bucket_set_ascend[i]*buckets[i] for i=0 to i= bucket_set_size - 1. 0 is in bucket_set */\
/*Correctness verified*/
void blst_p1_integrate_buckets(blst_p1 *out, blst_p1xyzz buckets[], int bucket_set_ascend[], size_t bucket_set_size, int d_max){
    
    blst_p1xyzz tmp, tmp1, tmp_d[d_max+1];
    vec_zero(&tmp, sizeof(tmp));
    vec_zero(tmp_d, sizeof(tmp_d[0])*(d_max+1));
    for(size_t i = bucket_set_size -1; i> 0; --i){
        blst_p1xyzz_dadd(&tmp, &tmp, &buckets[i]);
        int differ =  bucket_set_ascend[i] -  bucket_set_ascend[i-1];
        blst_p1xyzz_dadd(&tmp_d[differ], &tmp_d[differ], &tmp);   
    }
    
    vec_zero(&tmp, sizeof(tmp));
    vec_zero(&tmp1, sizeof(tmp1));

    for(int i = d_MAX_DIFF; i > 0; --i ){
        blst_p1xyzz_dadd(&tmp, &tmp, &(tmp_d[i]));
        blst_p1xyzz_dadd(&tmp1, &tmp1, &tmp);
    }

    blst_p1xyzz_to_Jacobian(out, &tmp1); 
}

/*Correctness has been verified*/
blst_p1_affine pippenger_variant_submission_CHES_backup(){
    
    std::cout <<"pippenger_variant_submission_CHES() with accumulation_d" << std::endl;

    blst_p1xyzz* buckets;
    buckets = new blst_p1xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 
  

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;
   
    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){

            int m = ret_MB_expr[j][0];
            int b = ret_MB_expr[j][1];
            int booth_idx = BUCKET_VALUE_TO_ITS_INDEX[b];
            unsigned char booth_sign;
            blst_p1_affine tmp_Pa;

            if (m> 0) {
                tmp_Pa = (m == 1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m==2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
                booth_sign = 0;
            }
            else{
                tmp_Pa = (m == -1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m==-2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
                booth_sign = 1;
            }
            blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
        
        }
    }
    
    blst_p1 ret;
    blst_p1_affine res_affine;

    blst_p1_integrate_buckets(&ret, buckets, BUCKET_SET, B_SIZE, d_MAX_DIFF);
    blst_p1_to_affine( &res_affine, &ret);

    delete[] buckets;

    return res_affine;
}

/*Correctness has been verified*/
blst_p1_affine pippenger_variant_submission_CHES(){
    
    std::cout <<"pippenger_variant_submission_CHES() with accumulation_d" << std::endl;

    blst_p1xyzz* buckets;
    buckets = new blst_p1xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 
  

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;
   
    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);

        int j = 0;

        int booth_idx = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j][1]];
        int booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j+1][1]];
        unsigned char booth_sign;
        blst_p1_affine tmp_Pa;

        int m = ret_MB_expr[j][0];
        if (m > 0) {
            tmp_Pa = (m == 1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][0]: \
            ((m==2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][0]:\
            (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][0]);
            booth_sign = 0;
        }
        else{
            tmp_Pa = (m == -1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][0]: \
            ((m==-2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][0]:\
            (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][0]);
            booth_sign = 1;
        }
        blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
        
        ++j;

        for( ; j< h_LEN_SCALAR-1; ++j){
            m = ret_MB_expr[j][0];
            booth_idx = booth_idx_nxt;
            booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j+1][1]];
            blst_p1_prefetch(buckets, booth_idx_nxt);

            unsigned char booth_sign;
            blst_p1_affine tmp_Pa;

            if (m> 0) {
                tmp_Pa = (m == 1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m==2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
                booth_sign = 0;
            }
            else{
                tmp_Pa = (m == -1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m==-2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
                booth_sign = 1;
            }
            blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
        }

        m = ret_MB_expr[j][0];
        booth_idx = booth_idx_nxt;
        booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j][1]];

        if (m> 0) {
            tmp_Pa = (m == 1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
            ((m==2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
            (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
            booth_sign = 0;
        }
        else{
            tmp_Pa = (m == -1)? (*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
            ((m==-2)? (*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
            (*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
            booth_sign = 1;
        }
        blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    }
    
    blst_p1 ret;
    blst_p1_affine res_affine;

    blst_p1_integrate_buckets(&ret, buckets, BUCKET_SET, B_SIZE, d_MAX_DIFF);
    blst_p1_to_affine( &res_affine, &ret);

    delete[] buckets;

    return res_affine;
}






void blst_p1_bucket_CHES(blst_p1xyzz buckets[], int booth_idx, const blst_p1_affine *p, unsigned char booth_sign){
    blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], p, booth_sign);
}

void blst_p1_tile_pippenger_CHES_q_over_5(blst_p1 *ret, \
                                    const blst_p1_affine *const points[], \
                                    size_t npoints, \
                                    const int scalars[], const unsigned char booth_signs[], \
                                    blst_p1xyzz buckets[], int bucket_set_ascend[], size_t bucket_set_size){
    
    // Initialization    
    vec_zero(buckets, sizeof(buckets[0])*bucket_set_size); \
    vec_zero(ret, sizeof(*ret)); \

    int i, scalar, booth_idx, booth_idx_nxt;
    unsigned char booth_sign;
    
    const blst_p1_affine *point = *points++;
    
    booth_idx = BUCKET_VALUE_TO_ITS_INDEX[*scalars++];
    booth_sign = *booth_signs++;

    booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[*scalars++];
    --npoints;
    if(booth_idx) blst_p1_bucket_CHES(buckets, booth_idx, point, booth_sign);
    for(i = 1; i < npoints; ++i){
        booth_idx = booth_idx_nxt;
        booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[*scalars++];
        blst_p1_prefetch(buckets, booth_idx_nxt);

        point = *points++; 
        booth_sign = *booth_signs++;
        if(booth_idx) blst_p1_bucket_CHES(buckets, booth_idx, point,  booth_sign);
    }
    point = *points;
    booth_sign = *booth_signs;
    blst_p1_bucket_CHES(buckets, booth_idx_nxt, point, booth_sign);// Carefully, it must be booth_idx_nxt
    blst_p1_integrate_buckets(ret, buckets, bucket_set_ascend, bucket_set_size, d_MAX_DIFF);

}

blst_p1_affine pippenger_variant_submission_CHES_blst_tile(){
    
    std::cout <<"pippenger_variant_submission_CHES_blst_tile() with custom-tailored tile accumulation" << std::endl;


    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;

    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints];
    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];

    blst_p1_affine** points_ptr;
    points_ptr = new blst_p1_affine* [npoints]; // points_ptr is an array of pointers that point to blst_p1_affine points.

    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){
            size_t idx = i*h_LEN_SCALAR + j;
            int m = ret_MB_expr[j][0];
            scalars[idx]  = ret_MB_expr[j][1];

            if (m> 0) {
                points_ptr[idx] = (m == 1)? &(*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m==2)? &(*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                &(*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
                booth_signs[idx] = 0; 
            }
            else{
                points_ptr[idx] = (m == -1)? &(*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m== -2)? &(*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                &(*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
                booth_signs[idx] = 1; 
            }  
        }
    }
    blst_p1 ret; // Mont coordinates

    blst_p1xyzz* buckets;
    buckets = new blst_p1xyzz [B_SIZE];
    blst_p1_tile_pippenger_CHES_q_over_5(&ret, points_ptr, npoints, scalars, booth_signs,\
                                         buckets, BUCKET_SET, B_SIZE);
    delete[] buckets;
    delete[] points_ptr;
    delete[] booth_signs;    
    delete[] scalars;    
    blst_p1_affine res_affine;
    blst_p1_to_affine( &res_affine, &ret);

    return res_affine;
}


// Use  blst_p1s_tile_pippenger to implement the CHES 
blst_p1_affine pippenger_variant_submission_CHES_blst_tile_copy(){
    
    std::cout <<"pippenger_variant_submission_CHES_blst_tile() with blst_p1s_tile_pippenger" << std::endl;

    // blst_p1_affine tmp_Pa = G1_AFFINE_INFINITY;
    // blst_p1 tmp_P = G1_INFINITY;

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;

    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    uint8_t* INTERMEDIATE_SCALARS;
    INTERMEDIATE_SCALARS = new  uint8_t[npoints*4];

    uint8_t** scalars_ptr;
    scalars_ptr = new uint8_t* [npoints];

    blst_p1_affine** points_ptr;
    points_ptr = new blst_p1_affine* [npoints]; // points_ptr is an array of pointers that point to blst_p1_affine points.
    
    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){
            size_t idx = i*h_LEN_SCALAR + j;
            int m = ret_MB_expr[j][0];

            uint32_t b_value = uint32_t(ret_MB_expr[j][1]);
            // blst_scalar_from_uint32(&INTERMEDIATE_SCALARS[idx], b_value);
            uint8_t tmp_uint8_t_array[4];
            byte_str_from_uint32(tmp_uint8_t_array, b_value);
            for(int i = 0; i<4; ++i){
                INTERMEDIATE_SCALARS[idx*4 +i] = tmp_uint8_t_array[i];
            }
            scalars_ptr[idx] = &INTERMEDIATE_SCALARS[idx*4];
            // scalars_ptr[idx] = INTERMEDIATE_SCALARS[1].b;
            if (m> 0) {
                points_ptr[idx] = (m == 1)? &(*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m==2)? &(*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                &(*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
            }
            else{
                // points_ptr[idx] = (m == -1)? &(*PRECOMPUTATION_POINTS_LIST_nh_n1)[i][j]: \
                // ((m== -2)? &(*PRECOMPUTATION_POINTS_LIST_nh_n2)[i][j]:\
                // &(*PRECOMPUTATION_POINTS_LIST_nh_n3)[i][j]);
               
                points_ptr[idx] = (m == -1)? &(*PRECOMPUTATION_POINTS_LIST_nh)[i][j]: \
                ((m== -2)? &(*PRECOMPUTATION_POINTS_LIST_nh_2)[i][j]:\
                &(*PRECOMPUTATION_POINTS_LIST_nh_3)[i][j]);
            }  
        }
    }
    blst_p1 ret_P; // Mont coordinates
    size_t window = EXPONENT_OF_q;

    limb_t* SCRATCH;
    SCRATCH = new limb_t[blst_p1s_mult_pippenger_scratch_sizeof_CHES_q_over_5(size_t (1<< EXPONENT_OF_q))/sizeof(limb_t)];
    
    // blst_p1s_tile_pippenger_CHES(&ret_P, points_ptr,\
    //                         npoints, scalars_ptr,\
    //                         window,  SCRATCH,\
    //                         0, window);
    
    delete[] SCRATCH;
    delete[] points_ptr;
    delete[] scalars_ptr;
    delete[] INTERMEDIATE_SCALARS;

    blst_p1_affine res_affine;
    blst_p1_to_affine( &res_affine, &ret_P);

    return res_affine;
}



void test(){

    std::cout << "Test starts.\n";
    blst_scalar scalars[N_POINTS];
    uint8_t* scalars_ptr[N_POINTS];

    for(size_t i = 0; i < N_POINTS; ++i){
            blst_scalar_from_uint32( &scalars[i], SCALARS_ARRAY[i].data);
            scalars_ptr[i] = (scalars[i].b);
    }

    blst_p1_affine** points_ptr;
    points_ptr = new blst_p1_affine*[N_POINTS]; // points_ptr is an array of pointers that point to blst_p1_affine points.


    for(size_t i = 0; i < N_POINTS; ++i){
            points_ptr[i] = &((*FIX_POINTS_LIST)[i]);
        }

    blst_p1 ret_P; // Mont coordinates
    blst_p1_affine ret_P_affine;

    /*blst pippenger*/
    size_t nbits = 255;
    
    std::cout<< "blst pippenger test: " << std::endl;

    size_t TEST_NUM = 5;

    limb_t* SCRATCH;

    auto st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< TEST_NUM; ++i)
    {
        SCRATCH = new limb_t[blst_p1s_mult_pippenger_scratch_sizeof(N_POINTS)/sizeof(limb_t)];
        blst_p1s_mult_pippenger(&ret_P, points_ptr, N_POINTS, scalars_ptr, nbits, SCRATCH);   
        delete[] SCRATCH;
    }
    auto ed = std::chrono::steady_clock::now();   
 

    std::chrono::microseconds diff1 = std::chrono::duration_cast<std::chrono::microseconds>(ed -st); 
    std::cout << "blst pippenger Wall clock time elapse is: " << std::dec << diff1.count()/TEST_NUM << " us "<< std::endl;
    blst_p1_to_affine(&ret_P_affine, &ret_P);
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;

    /*nh + q/5 method with tile*/
    st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< TEST_NUM; ++i)
    {
        ret_P_affine = pippenger_variant_submission_CHES_blst_tile();
    }
    ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff2 = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    std::cout << "CHES 'nh+ q/5' with blst tile Wall clock time elapse is: " << std::dec << diff2.count()/TEST_NUM << " us "<< std::endl;
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;

    std::cout << "Improvement: " << float(diff1.count() - diff2.count())/float(diff1.count()) <<std::endl;

    /*nh + q/5 method*/
    st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< TEST_NUM; ++i)
    {
        ret_P_affine = pippenger_variant_submission_CHES();
    }
    ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff3 = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    std::cout << "CHES 'nh+ q/5' method Wall clock time elapse is: " << std::dec << diff3.count()/TEST_NUM << " us "<< std::endl;
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;
    
    std::cout << "Improvement: " << float(diff1.count() - diff3.count())/float(diff1.count()) <<std::endl;





    // /*trivial method 2 */
    // st = std::chrono::steady_clock::now();
    // ret_P_affine  =  trivial_mult_scalar_multiplication_2();

    // ed = std::chrono::steady_clock::now();   
    // diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    // std::cout << "trivial method 2 Wall clock time elapse is: " << std::dec << diff.count() << " us "<< std::endl;
    // std::cout << ret_P_affine.x <<std::endl;
    // std::cout << ret_P_affine.y <<std::endl;
 

    //  /*trivial method*/

    // st = std::chrono::steady_clock::now();
    // trivial_mult_scalar_multiplication(&ret_P, points_ptr, N_POINTS, scalars_ptr, 32);
    // ed = std::chrono::steady_clock::now();   
    // diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    // std::cout << "trivial method Wall clock time elapse is: " << std::dec << diff.count() << " us "<< std::endl;
    // blst_p1_to_affine(&ret_P_affine, &ret_P);
    // std::cout << ret_P_affine.x <<std::endl;
    // std::cout << ret_P_affine.y <<std::endl;

    std::cout << "TEST END" <<std::endl;

}

void test_scalar_conversion_bench(){

    scalar_MB_expr ret_MB_expr;
    uint32_t mask = (1 << EXPONENT_OF_q) - 1;
    auto tmp = random_scalar_less_than_r();

    auto st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< N_POINTS; ++i)
    {
    trans_uint256_t_to_MB_radixq_expr(ret_MB_expr,SCALARS_ARRAY[i]);   
    }
    auto ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff = std::chrono::duration_cast<std::chrono::microseconds>(ed -st); 
    std::cout << "scalars conversion clock time elapse is: " << std::dec << diff.count() << " us "<< std::endl;
    print(ret_MB_expr);

    st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< N_POINTS; ++i)
    {
    trans_uint256_t_to_MB_radixq_expr_2008_08654(ret_MB_expr, tmp, mask);   
    }
    ed = std::chrono::steady_clock::now();   
    diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st); 
    std::cout << "scalars conversion clock time elapse is: " << std::dec << diff.count() << " us "<< std::endl;
    
}

void test_three_scalar_multiplication_methods(){// Correctness has been tested.
    
    uint256_t a_int = random_scalar_less_than_r();
    blst_p1 xyzQ, xyzret2,  xyzret3;
    blst_p1_affine ret1, ret2, ret3, Q = *blst_p1_affine_generator();
    
    blst_p1_from_affine(&xyzQ, &Q);
    std::cout << "Base point in projective coordinates is: "<< std::endl;
    std::cout << xyzQ.x << std::endl;
    std::cout << xyzQ.y << std::endl;
    std::cout << xyzQ.z << std::endl;  

    ret1 =  single_scalar_multiplication(a_int, Q);

    blst_scalar s;
    blst_scalar_from_uint32( &s, a_int.data);
    blst_p1_mult( &xyzret2, &xyzQ, s.b, (size_t) 255);
    blst_p1_to_affine(&ret2, &xyzret2);

    xyzret3  = single_scalar_multiplication( a_int, xyzQ);
    blst_p1_to_affine(&ret3, &xyzret3);
    
    std::cout << "ret1 is: "<< std::endl;
    std::cout << ret1.x << std::endl;
    std::cout << ret1.y << std::endl;

    std::cout << "ret2 is: "<< std::endl;
    std::cout << ret2.x << std::endl;
    std::cout << ret2.y << std::endl;

    std::cout << "ret3 is: "<< std::endl;
    std::cout << ret3.x << std::endl;
    std::cout << ret3.y << std::endl;

    uint32_t intlist[8] = { 0, 0, 0, 0, 0, 0, 0, 0};
    intlist[0] = 7891;
    blst_scalar_from_uint32( &s, intlist);
    std::cout << "small integer conversion test is: "<< std::endl;
    std::cout << s << std::endl;
}


void test_xyzz(){
    std::cout << " test_xyzz() " << std::endl;
    blst_p1_affine P = *blst_p1_affine_generator();
    blst_p1 Q = *blst_p1_generator();
    blst_p1 ret = Q;
    blst_p1_affine reta;
    blst_p1 ret2 = Q;
    blst_p1_affine ret2a;

    unsigned char booth_sign = 0; // if booth_sign == 1, it is a subtraction.
    blst_p1xyzz Pxyzz;
    vec_zero(&Pxyzz, sizeof(Pxyzz));
    blst_p1xyzz_to_Jacobian(&ret, &Pxyzz);
    std::cout <<"xyzz initialization is inf test: "<<  blst_p1_is_inf(&ret) << std::endl;

    blst_p1xyzz_dadd_affine(&Pxyzz, &Pxyzz, &P, booth_sign);
    blst_p1_to_affine(&reta, &ret);

    blst_p1_add_or_double(&ret2, &Q, &Q);
    blst_p1_to_affine(&ret2a, &ret2);

    std::cout << Pxyzz.x << std::endl;
    std::cout << Pxyzz.y << std::endl;
    std::cout << Pxyzz.zzz << std::endl;
    std::cout << Pxyzz.zz << std::endl;

    std::cout << P.x << std::endl;
    std::cout << P.y << std::endl;
}

int main(){


    init();

    init_PRECOMPUTATION_POINTS_LIST_nh();
    test();

    // test_scalar_conversion_bench();
    // test_three_scalar_multiplication_methods();// Correct!


    // blst_p1_affine Q = *blst_p1_affine_generator();
    // blst_p1_affine Qn = Q;

    // blst_fp_sub(&Qn.y, &FP_Z, &Qn.y);
    // blst_p1 ret, xyzQ, xyzQn;

    // blst_p1_from_affine(&xyzQ, &Q);

    // blst_p1_from_affine(&xyzQn, &Qn);   
    // blst_p1_add_or_double(&ret, &xyzQn, &xyzQ);


    // std::cout << ret.x << std::endl;
    // std::cout << ret.y << std::endl;
    // std::cout << ret.z << std::endl;

    // for(int e = 10; e<= 30; ++e)
    // std::cout << std::dec << e<< "  "<< pippenger_window_size( (size_t) (1<< e) -8) << std::endl;

    // byte ret1[4];
    // uint a = 0x12343a79;
    // byte_str_from_uint32( ret1, a);

    // std::cout << ret1[0] << " "<< ret1[1] << " "<< ret1[2] << " "<< ret1[3] << " "<<std::endl;
    // std::cout << "booth_encode test" << std::endl;
    // for(int i = 0; i< 20; ++i){
    //    std::cout << i << " " << booth_encode(i, 4) << std::endl;
    // }

    // std::cout << "get_wval_limb test" << std::endl;
    // for(int i = 0; i< (1<<13); ++i){
    //     uint8_t d[4];
    //     byte_str_from_uint32( d, (uint32_t) i);
    //     limb_t wval = (get_wval_limb(d, 0 , 13) << 1);
    //    std::cout << i << " " << wval << std::endl;

    // //    booth_encode(wval, 13)
    // }

    // int nn = 4;
    // while(nn) std::cout <<  SCALARS_ARRAY[nn--] << std::endl; 

    test_xyzz();


    return 0;
}

