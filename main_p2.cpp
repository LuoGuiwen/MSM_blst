/***----***
 
compile with: g++ -std=c++17 -o main_test -g -O2 main_p2.cpp libblst.a

If segmentation fault occurs, possibly it can be tentatively circumvented by using the following code in command line to unleash the stack restriction:

ulimit -s unlimited

All functions must be invoked after init().

-- Guiwen. Oct, 2022.

***----***/


#include "bindings/blst.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <iomanip>

#include <set>
#include <array>
// don't use namespace std, there is a uint_8 and byte definition that are not compatible with std.

/***----***
A) Define global variables and their initializations
***----***/

#include "config_file_n_exp_10.h" //define configuration in a seperate file.

std::set<int> MULTI_SET = {1, 2, 3};
int* BUCKET_SET;
int* BUCKET_VALUE_TO_ITS_INDEX; 

typedef struct {int m; int b; int alpha;} digit_decomposition;
digit_decomposition* DIGIT_CONVERSION_HASH_TABLE;

#include "auxiliaryfunc.h"

uint256_t* SCALARS_ARRAY;
blst_p2_affine *FIX_POINTS_LIST;
blst_p2_affine* PRECOMPUTATION_POINTS_LIST_3nh;
blst_p2_affine* PRECOMPUTATION_POINTS_LIST_BGMW95;

void init(){

    // Initialize FIX_POINTS_LIST
    FIX_POINTS_LIST = new blst_p2_affine [N_POINTS];
    blst_p2 tmp_P, tmp2_P = *blst_p2_generator();
    blst_p2_affine tmp_P_affine;

    for(size_t i = 0; i < N_POINTS; ++i){
            blst_p2_double(&tmp_P, &tmp2_P);
            tmp2_P = tmp_P;  
            blst_p2_to_affine(&tmp_P_affine, &tmp_P);
            (FIX_POINTS_LIST)[i] = tmp_P_affine;
        }
    std::cout<< "FIX_POINTS_LIST Generated" <<std::endl;

    // Initialize SCALARS_ARRAY
    SCALARS_ARRAY = new uint256_t[N_POINTS];
    for(size_t i = 0; i < N_POINTS; ++i)\
        SCALARS_ARRAY[i] = random_scalar_less_than_r();
    std::cout<< "SCALARS_ARRAY Generated" <<std::endl;
}

void free_init(){
    delete[] SCALARS_ARRAY;
    delete[] FIX_POINTS_LIST;
}

blst_p2_affine single_scalar_multiplication(uint256_t scalar, blst_p2_affine Q){

    blst_p2_affine aret; 
    blst_p2 ret = {0,1,0}; // ret = INFINITY; 

    blst_p2 xyzQ;
    blst_p2_from_affine(&xyzQ, &Q);

    while (scalar > 0){
        if ( scalar.data[0] & 1 ){
            blst_p2_add_or_double(&ret, &ret, &xyzQ);
        }
        
        blst_p2_add_or_double(&xyzQ, &xyzQ, &xyzQ);    // tested. No need to use temp variable.
        scalar = scalar >> 1; 
    }

    blst_p2_to_affine(&aret, &ret);
    return aret;
} 


void init_pippenger_BGMW95(){
    PRECOMPUTATION_POINTS_LIST_BGMW95 = new blst_p2_affine [h_BGMW95*N_POINTS];
    auto st = std::chrono::steady_clock::now();

    for(int i = 0; i< N_POINTS; ++i){
        blst_p2_affine qjQi = FIX_POINTS_LIST[i];
        for(int j = 0; j< h_BGMW95; ++j){
            auto idx = i*h_BGMW95 +j;
            PRECOMPUTATION_POINTS_LIST_BGMW95[idx] = qjQi;
            qjQi = single_scalar_multiplication(q_RADIX_PIPPENGER_VARIANT, qjQi);    
        }
    }
    auto ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff = std::chrono::duration_cast<std::chrono::microseconds>(ed -st);
    std::cout<< "PRECOMPUTATION_POINTS_LIST_BGMW95 SUCCESSFULLY CONSTRUCTED" << std::endl;
    std::cout << "PRECOMPUTATION Wall clock time elapse is: " << diff.count() << " us "<< std::endl;
}   

void free_init_pippenger_BGMW95(){
    delete[] PRECOMPUTATION_POINTS_LIST_BGMW95;
}

void init_pippenger_CHES_q_over_5(){

    //Initialize BUCKET_SET and BUCKET_VALUE_TO_ITS_INDEX 
    BUCKET_SET = new int[B_SIZE];
    construct_bucket_set(BUCKET_SET, q_RADIX, a_LEADING_TERM);
    std::cout<< "BUCKET_SET constructed. The size of BUCKET_SET is: " << B_SIZE << std::endl;

    BUCKET_VALUE_TO_ITS_INDEX = new int[q_RADIX/2 +1];
    for(size_t i = 0; i < B_SIZE; ++i){
        BUCKET_VALUE_TO_ITS_INDEX[BUCKET_SET[i]] = i;
    }
    
    // Initialize  DIGIT_CONVERSION_HASH_TABLE;
    DIGIT_CONVERSION_HASH_TABLE = new digit_decomposition[q_RADIX+1];
    for(int m: MULTI_SET){
        for (int i = 0; i < B_SIZE; ++i){
            int b = BUCKET_SET[i];
            if (m*b <= q_RADIX) DIGIT_CONVERSION_HASH_TABLE[q_RADIX - m*b] = {m,b,1};
        }
    }
    for(int m: MULTI_SET){
        for (int i = 0; i < B_SIZE; ++i){
            int b = BUCKET_SET[i];
            if (m*b <= q_RADIX) DIGIT_CONVERSION_HASH_TABLE[m*b] = {m,b,0};
        }
    }
    std::cout<< "DIGIT_CONVERSION_HASH_TABLE constructed." <<std::endl;
    
// ### Initialize the precomputation ###
   PRECOMPUTATION_POINTS_LIST_3nh = new blst_p2_affine [3*N_POINTS*h_LEN_SCALAR];

    auto st = std::chrono::steady_clock::now();
    blst_p2_affine Pt;

    for(int i = 0; i< N_POINTS; ++i){
        blst_p2_affine qjQi = FIX_POINTS_LIST[i];
        for(int j = 0; j< h_LEN_SCALAR; ++j){
            for(int m = 1; m <=3; ++m){
            size_t idx_i_j_m =  3*(i*h_LEN_SCALAR +j) + m-1;  
            if(m==1) PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m]  = qjQi;
            else{PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m] = single_scalar_multiplication(m, qjQi);}
            }
            qjQi = single_scalar_multiplication(q_RADIX, qjQi);    
        }
    }
    auto ed = std::chrono::steady_clock::now();   

    std::chrono::microseconds diff = std::chrono::duration_cast<std::chrono::microseconds>(ed -st);
    std::cout<< "PRECOMPUTATION_POINTS_LIST_3nh SUCCESSFULLY CONSTRUCTED" << std::endl;
    std::cout << "PRECOMPUTATION Wall clock time elapse is: " << diff.count() << " us "<< std::endl;
}

void free_init_pippenger_CHES_q_over_5(){

    delete[] PRECOMPUTATION_POINTS_LIST_3nh;
    delete[] DIGIT_CONVERSION_HASH_TABLE;
    delete[] BUCKET_VALUE_TO_ITS_INDEX;
    delete[] BUCKET_SET;
}


/***----***
B) Define pippenger's bucket method and its variants
***----***/

/*Correctness has been verified*/
blst_p2_affine pippenger_variant_CHES_primitive_implementation(){
    
    blst_p2xyzz* buckets;
    buckets = new blst_p2xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 
  

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;
   
    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);

        int j = 0;

        int booth_idx = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j][1]];
        int booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j+1][1]];
        unsigned char booth_sign;
        blst_p2_affine tmp_Pa;
        size_t idx_i_j_m;

        int m = ret_MB_expr[j][0];
        if (m> 0) {               
            idx_i_j_m =  3*(i*h_LEN_SCALAR+j) + m-1; 
            tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
            booth_sign = 0;
        }
        else{
            idx_i_j_m =  3*(i*h_LEN_SCALAR+j)  - m-1; 
            tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
            booth_sign = 1;
        }
        blst_p2xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
        
        ++j;
   
 
        for( ; j< h_LEN_SCALAR-1; ++j){
            m = ret_MB_expr[j][0];
            booth_idx = booth_idx_nxt;
            booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j+1][1]];
            blst_p2_prefetch_CHES(buckets, booth_idx_nxt);

            if (m> 0) {               
                idx_i_j_m =  3*(i*h_LEN_SCALAR+j) + m-1; 
                tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_sign = 0;
            }
            else{
                idx_i_j_m =  3*(i*h_LEN_SCALAR+j)  - m-1; 
                tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_sign = 1;
            }
            blst_p2xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
        }

        m = ret_MB_expr[j][0];
        booth_idx = booth_idx_nxt;
        booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j][1]];

            if (m> 0) {               
                idx_i_j_m =  3*(i*h_LEN_SCALAR+j) + m-1; 
                tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_sign = 0;
            }
            else{
                idx_i_j_m =  3*(i*h_LEN_SCALAR+j)  - m-1; 
                tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_sign = 1;
            }
        blst_p2xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    }
    
    blst_p2 ret;
    blst_p2_affine res_affine;

    blst_p2_integrate_buckets_accumulation_d_CHES(&ret, buckets, BUCKET_SET, B_SIZE, d_MAX_DIFF);
    blst_p2_to_affine( &res_affine, &ret);

    delete[] buckets;

    return res_affine;
}

/* Correctess verified*/
blst_p2_affine pippenger_variant_q_over_5_CHES_prefetch_1step_ahead(){
    
    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints+2];
    int* scalar_p = scalars;

    std::array< int, h_LEN_SCALAR> ret_std_expr;
    size_t s_idx = 0;
    for(int i = 0; i< N_POINTS; ++i){
        trans_uint256_t_to_standard_q_ary_expr(ret_std_expr, SCALARS_ARRAY[i]);
        for(int j = 0; j< h_LEN_SCALAR; ++j){
            // std::cout << *scalars++ << std::endl;
            scalars[s_idx++] = ret_std_expr[j];
        }
    }

    int scalar_now, point_idx, point_idx_nxt, scalar_nxt, booth_idx, booth_idx_nxt, mul, mul_nxt, b, b_nxt;
    unsigned char booth_sign, booth_sign_nxt;
    blst_p2_affine tmp_Pa;
    blst_p2xyzz* buckets;
    buckets = new blst_p2xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 
  
    digit_decomposition tmp_tri, tmp_tri_nxt;
    size_t size_tri = sizeof(tmp_tri);
    size_t size_point = sizeof(blst_p2_affine);

    size_t i = 0;

    tmp_tri = DIGIT_CONVERSION_HASH_TABLE[scalars[i]];
    booth_idx = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri.b];
    booth_sign = tmp_tri.alpha;
    if(tmp_tri.alpha) scalars[i+1] +=1; // if tmp_tri[2] == 1
    // i = 0
    point_idx = 3*i + tmp_tri.m - 1;
    tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx];

    i = 1;
    tmp_tri_nxt = DIGIT_CONVERSION_HASH_TABLE[scalars[i]];
    booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt.b];  
    booth_sign_nxt = tmp_tri_nxt.alpha;  
    if(tmp_tri_nxt.alpha) scalars[i+1] +=1;
    // i = 1
    point_idx_nxt = 3*i + tmp_tri_nxt.m - 1;
    vec_prefetch(&PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt], size_point);
    blst_p2xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);

    for( i = 2; i < npoints; ++i){

        tmp_tri = tmp_tri_nxt;
        booth_idx = booth_idx_nxt;
        booth_sign = booth_sign_nxt;
        point_idx = point_idx_nxt;
       
        tmp_tri_nxt = DIGIT_CONVERSION_HASH_TABLE[scalars[i]];
        //i == 2
        booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt.b];  
        booth_sign_nxt = tmp_tri_nxt.alpha;  
        if(tmp_tri_nxt.alpha) scalars[i+1] +=1;

        point_idx_nxt = 3*i + tmp_tri_nxt.m - 1;

        vec_prefetch(&PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt], size_point);
        vec_prefetch(&buckets[booth_idx_nxt], size_point);

        tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx];
        if(booth_idx) blst_p2xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    }
    tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt];
    if(booth_idx_nxt) blst_p2xyzz_dadd_affine(&buckets[booth_idx_nxt], &buckets[booth_idx_nxt], &tmp_Pa, booth_sign_nxt);

    blst_p2 ret; // Mont coordinates
    blst_p2_affine res_affine;
    blst_p2_integrate_buckets_accumulation_d_CHES(&ret, buckets, BUCKET_SET, B_SIZE, d_MAX_DIFF);
    blst_p2_to_affine( &res_affine, &ret);
    
    delete[] buckets; 
    delete[] scalars;    

    blst_p2_to_affine( &res_affine, &ret);
    return res_affine;
}

/* Correctess verified*/
blst_p2_affine pippenger_variant_q_over_5_CHES_prefetch_2step_ahead(){
    
    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints+2];

    int* scalar_p = scalars;

    std::array< int, h_LEN_SCALAR> ret_std_expr;
    size_t s_idx = 0;
    for(int i = 0; i< N_POINTS; ++i){
        trans_uint256_t_to_standard_q_ary_expr(ret_std_expr, SCALARS_ARRAY[i]);
        for(int j = 0; j< h_LEN_SCALAR; ++j){
            scalars[s_idx++] = ret_std_expr[j];
        }
    }

    int  point_idx, point_idx_nxt, booth_idx, booth_idx_nxt;
    unsigned char booth_sign, booth_sign_nxt;
    
    blst_p2_affine tmp_Pa;
    blst_p2xyzz* buckets;
    buckets = new blst_p2xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 
  
    digit_decomposition tmp_tri, tmp_tri_nxt, tmp_tri_nxt2;
    size_t size_tri = sizeof(tmp_tri);
    size_t size_point = sizeof(blst_p2_affine);

    size_t i = 0;

    tmp_tri = DIGIT_CONVERSION_HASH_TABLE[scalars[i]];
    booth_idx = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri.b];
    booth_sign = tmp_tri.alpha;
    if(tmp_tri.alpha) scalars[i+1] +=1; // if tmp_tri[2] == 1
    // i = 0
    point_idx = 3*i + tmp_tri.m - 1;
    tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx];

    i = 1;
    tmp_tri_nxt = DIGIT_CONVERSION_HASH_TABLE[scalars[i]];
    booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt.b];  
    booth_sign_nxt = tmp_tri_nxt.alpha;  
    if(tmp_tri_nxt.alpha) scalars[i+1] +=1;
    // i = 1
    point_idx_nxt = 3*i + tmp_tri_nxt.m - 1;

    vec_prefetch(&PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt], size_point);
    blst_p2xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    
    i = 2;
    tmp_tri_nxt2 = DIGIT_CONVERSION_HASH_TABLE[scalars[i]];

    while(i < npoints){

        tmp_tri = tmp_tri_nxt;
        booth_idx = booth_idx_nxt;
        booth_sign = booth_sign_nxt;
        point_idx = point_idx_nxt;

        // i == 2
        tmp_tri_nxt = tmp_tri_nxt2;
        booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt.b];  
        booth_sign_nxt = tmp_tri_nxt.alpha;  
        if(tmp_tri_nxt.alpha) scalars[i+1] +=1;
        point_idx_nxt = 3*i + tmp_tri_nxt.m - 1;

        ++i;
        //i == 3
        tmp_tri_nxt2 = DIGIT_CONVERSION_HASH_TABLE[scalars[i]];

        vec_prefetch(&DIGIT_CONVERSION_HASH_TABLE[scalars[i+1]], size_tri);
        vec_prefetch(&BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt2.b], 4);
        vec_prefetch(&PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt], size_point);
        vec_prefetch(&buckets[booth_idx_nxt], size_point);

        tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx];
        if(booth_idx) blst_p2xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    }
    tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt];
    if(booth_idx_nxt) blst_p2xyzz_dadd_affine(&buckets[booth_idx_nxt], &buckets[booth_idx_nxt], &tmp_Pa, booth_sign_nxt);

    blst_p2 ret; // Mont coordinates
    blst_p2_affine res_affine;
    blst_p2_integrate_buckets_accumulation_d_CHES(&ret, buckets, BUCKET_SET, B_SIZE, d_MAX_DIFF);
    blst_p2_to_affine( &res_affine, &ret);
    
    delete[] buckets; 
    delete[] scalars;    

    blst_p2_to_affine( &res_affine, &ret);
    return res_affine;
}

blst_p2_affine pippenger_variant_q_over_5_CHES_3nh(){
    
    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;

    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints+1]; // add 1 slot redundancy for prefetch, 20221010 very important
    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];

    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [npoints]; 

    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){
            size_t idx = i*h_LEN_SCALAR + j;
            int m = ret_MB_expr[j][0];
            scalars[idx]  = ret_MB_expr[j][1];

            if (m> 0) {
                size_t idx_i_j_m =  3*idx + m-1; 
                points_ptr[idx] = &PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_signs[idx] = 0; 
            }
            else{
                size_t idx_i_j_m =  3*idx  -m - 1; 
                points_ptr[idx] = &PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_signs[idx] = 1; 
            }  
        }
    }
    blst_p2 ret; // Mont coordinates

    blst_p2xyzz* buckets;

    buckets = new blst_p2xyzz [B_SIZE];

    blst_p2_tile_pippenger_d_CHES(&ret, points_ptr, npoints, scalars, booth_signs,\
                                         buckets, BUCKET_SET, BUCKET_VALUE_TO_ITS_INDEX , B_SIZE, d_MAX_DIFF);
    delete[] buckets;
    delete[] points_ptr;
    delete[] booth_signs;    
    delete[] scalars;    
    blst_p2_affine res_affine;
    blst_p2_to_affine( &res_affine, &ret);

    return res_affine;
}

blst_p2_affine pippenger_variant_q_over_5_CHES_noindexhash(){
    
    // try to use direct hash to obtain the booth_idxx, buckets take two much spaces.

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;

    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints];
    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];

    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [npoints]; 

    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){
            size_t idx = i*h_LEN_SCALAR + j;
            int m = ret_MB_expr[j][0];
            scalars[idx]  = ret_MB_expr[j][1];

            if (m> 0) {
                size_t idx_i_j_m =  3*idx + m-1; 
                points_ptr[idx] = &PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_signs[idx] = 0; 
            }
            else{
                size_t idx_i_j_m =  3*idx  -m - 1; 
                points_ptr[idx] = &PRECOMPUTATION_POINTS_LIST_3nh[idx_i_j_m];
                booth_signs[idx] = 1; 
            }  
        }
    }
    blst_p2 ret; // Mont coordinates

    blst_p2xyzz* buckets;
    buckets = new blst_p2xyzz [q_RADIX/2 +1];

    blst_p2_tile_pippenger_d_CHES_noindexhash(&ret, \
                                    points_ptr, \
                                    npoints, \
                                    scalars, booth_signs, \
                                    buckets, BUCKET_SET,\
                                    B_SIZE, d_MAX_DIFF);

    delete[] buckets;
    delete[] points_ptr;
    delete[] booth_signs;    
    delete[] scalars;    
    blst_p2_affine res_affine;
    blst_p2_to_affine( &res_affine, &ret);

    return res_affine;
}

blst_p2_affine pippenger_variant_BGMW95(){
    
    std::array< int, h_BGMW95> ret_qhalf_expr;

    uint64_t npoints = N_POINTS*h_BGMW95;

    int* scalars;
    scalars = new int [npoints];
    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];

    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [npoints]; 

    for(int i = 0; i< N_POINTS; ++i){
        trans_uint256_t_to_qhalf_expr(ret_qhalf_expr, SCALARS_ARRAY[i]);

        for(int j = 0; j< h_BGMW95; ++j){
            size_t idx = i*h_BGMW95 + j;
            scalars[idx]  = ret_qhalf_expr[j];
            points_ptr[idx] =  &PRECOMPUTATION_POINTS_LIST_BGMW95[idx];
            if ( scalars[idx] >= 0) {
                booth_signs[idx] = 0; 
            }
            else{
                scalars[idx] = -scalars[idx];
                booth_signs[idx] = 1; 
            }  
        }
    }
    blst_p2 ret; // Mont coordinates

    blst_p2xyzz* buckets;
    int qhalf = int(q_RADIX_PIPPENGER_VARIANT>>1);
    buckets = new blst_p2xyzz [qhalf + 1];

    blst_p2_tile_pippenger_BGMW95(&ret, \
                                    points_ptr, \
                                    npoints, \
                                    scalars, booth_signs,\
                                    buckets,\
                                    EXPONENT_OF_q_BGMW95);
    delete[] buckets;
    delete[] points_ptr;
    delete[] booth_signs;    
    delete[] scalars;  

    blst_p2_affine res_affine;
    blst_p2_to_affine( &res_affine, &ret);
    return res_affine;
}

blst_p2_affine pippenger_blst_built_in(){

    blst_scalar* scalars;
    scalars = new blst_scalar [N_POINTS];

    uint8_t** scalars_ptr;
    scalars_ptr = new uint8_t* [N_POINTS];

    for(size_t i = 0; i < N_POINTS; ++i){
            blst_scalar_from_uint32( &scalars[i], SCALARS_ARRAY[i].data);
            scalars_ptr[i] = (scalars[i].b);
    }

    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [N_POINTS]; 

    for(size_t i = 0; i < N_POINTS; ++i){
            points_ptr[i] = &(FIX_POINTS_LIST)[i];
        }

    limb_t* SCRATCH;
    SCRATCH = new limb_t[blst_p2s_mult_pippenger_scratch_sizeof(N_POINTS)/sizeof(limb_t)];

    blst_p2 ret; // Mont coordinates
    size_t nbits = 255;

    blst_p2s_mult_pippenger(&ret, points_ptr, N_POINTS, scalars_ptr, nbits, SCRATCH);   


    delete[] scalars;
    delete[] scalars_ptr;
    delete[] points_ptr;
    delete[] SCRATCH;

    blst_p2_affine res_affine;
    blst_p2_to_affine( &res_affine, &ret);
    return res_affine;
}

void test_pippengers(){
    std::cout << "PIPPENGERS TEST OVER G2 for NPOINTS:  2**" << N_EXP << std::endl;

    size_t TEST_NUM = 10;
    if(N_EXP <= 16) TEST_NUM = 20;

    blst_p2_affine ret_P_affine;

    /*nh + q/5 method */
    auto st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< TEST_NUM; ++i)
    {
        ret_P_affine = pippenger_variant_q_over_5_CHES_3nh();
    }
    auto ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff1 = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    std::cout << "1. CHES 'nh+ q/5'.. Wall clock time elapse is: " << std::dec << diff1.count()/TEST_NUM << " us "<< std::endl;
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;
    std::cout << std::endl;

    /* nh + q/5 method 2 */
    st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< TEST_NUM; ++i)
    {
        ret_P_affine = pippenger_variant_q_over_5_CHES_prefetch_2step_ahead();
    }
    ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff2 = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    std::cout << "2. CHES 'nh+ q/5' prfetch 2step ahead. Wall clock time elapse is: " << std::dec << diff2.count()/TEST_NUM << " us "<< std::endl;
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;
    std::cout << std::endl;

    /*nh + q/2 method*/
    st = std::chrono::steady_clock::now();
    for(size_t i = 0; i< TEST_NUM; ++i)
    {
        ret_P_affine = pippenger_variant_BGMW95();
    }
    ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff3 = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    std::cout << "3. pippenger_variant_BGMW95. Wall clock time elapse is: " << std::dec << diff3.count()/TEST_NUM << " us "<< std::endl;
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;
    std::cout << std::endl;   

    /*blst pippenger h(n+q/2) method*/
    st = std::chrono::steady_clock::now();    
    for(size_t i = 0; i< TEST_NUM; ++i) 
    {
        ret_P_affine = pippenger_blst_built_in();
    }
    ed = std::chrono::steady_clock::now(); 

    std::chrono::microseconds diff4 = std::chrono::duration_cast<std::chrono::microseconds>(ed -st); 
    std::cout << "4. pippenger_blst_built_in. Wall clock time elapse is: " << std::dec << diff4.count()/TEST_NUM << " us "<< std::endl;
    std::cout << ret_P_affine.x <<std::endl;
    std::cout << ret_P_affine.y <<std::endl;
    std::cout << std::endl;

    std::cout << "Improvement, blst BGMW95 vs pipp: " << float(diff4.count() - diff3.count())/float(diff4.count()) <<std::endl;
    std::cout << "Improvement, blst CHES_q_over_5 vs pipp: " << float(diff4.count() - diff1.count())/float(diff4.count()) <<std::endl;
    std::cout << "Improvement, blst CHES_q_over_5 vs BGMW95: " << float(diff3.count() - diff1.count())/float(diff3.count()) <<std::endl;

    std::cout << "TEST END\n" <<std::endl;

}

void test_scalar_conversion_bench(){

    size_t TEST_NUM = 10;

    if(N_EXP <= 16) TEST_NUM = 50;

    scalar_MB_expr ret_MB_expr;

    auto st = std::chrono::steady_clock::now();
    for(int j=0; j< TEST_NUM; ++j){
        for(size_t i = 0; i< N_POINTS; ++i)
            {
            trans_uint256_t_to_MB_radixq_expr(ret_MB_expr,SCALARS_ARRAY[i]);   
            }
    }
    auto ed = std::chrono::steady_clock::now();   
    std::chrono::microseconds diff = std::chrono::duration_cast<std::chrono::microseconds>(ed -st); 
    std::cout << "CHES q_over_5 scalars conversion clock time elapse is: " << std::dec << diff.count()/TEST_NUM << " us "<< std::endl;
    // print(ret_MB_expr);

    std::array< int, h_BGMW95> ret_qhalf_expr;

    st = std::chrono::steady_clock::now();
    for(int j=0; j< TEST_NUM; ++j){
        for(size_t i = 0; i< N_POINTS; ++i)
            {
            trans_uint256_t_to_qhalf_expr(ret_qhalf_expr, SCALARS_ARRAY[i]);   
            }
    }
    ed = std::chrono::steady_clock::now();   
    diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st); 
    std::cout << "BGMW95 scalars conversion clock time elapse is: " << std::dec << diff.count()/TEST_NUM << " us\n"<< std::endl;  
}



int main(){
 
    init();
    init_pippenger_CHES_q_over_5();
    init_pippenger_BGMW95();
    // Test should be down between init() and free_init();

    test_pippengers();
    test_scalar_conversion_bench();

    std::cout <<  "First several scalars used in our computation" << std::endl; 
    int nn = 4;
    while(nn) std::cout <<  SCALARS_ARRAY[nn--] << std::endl; 

    free_init_pippenger_BGMW95();
    free_init_pippenger_CHES_q_over_5();
    free_init();

    return 0;
}

