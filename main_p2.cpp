/***----***
 
Compile with: g++ -std=c++17 -o main_test -g -O2 main_p2.cpp libblst.a

If segmentation fault occurs, possibly it can be tentatively circumvented by using the following code in command line to unleash the stack restriction:

ulimit -s unlimited

All functions must be invoked after init_xx().

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
Define configuration
***----***/
#include "ches_config_files/config_file.h" //define configuration in a seperate file.
bool TEST_PIPPENGER_Q_OVER_5_CHES = 1;
bool TEST_PIPPENGER_BGMW95 = 1; // Due to the momory size limitation, some time we need to switch off one of the test.


/***----***
A) Define global variables and their initializations
***----***/
digit_decomposition* DIGIT_CONVERSION_HASH_TABLE;
#include "auxiliaryfunc.h"

std::set<int> MULTI_SET = {1, 2, 3};
int* BUCKET_SET;
int* BUCKET_VALUE_TO_ITS_INDEX; 

blst_p2_affine* FIX_POINTS_LIST;
blst_p2_affine* PRECOMPUTATION_POINTS_LIST_3nh;
blst_p2_affine* PRECOMPUTATION_POINTS_LIST_BGMW95;

void init_fix_point_list(){

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
}

void free_init_fix_point_list(){
    delete[] FIX_POINTS_LIST;
}

blst_p2_affine single_scalar_multiplication(uint256_t scalar, const blst_p2_affine Q){

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

        // We do not need to do computation in this scenario
    if(TEST_PIPPENGER_Q_OVER_5_CHES && (q_RADIX == q_RADIX_PIPPENGER_VARIANT)){
        for(int i = 0; i< N_POINTS; ++i){
            for(int j = 0; j< h_BGMW95; ++j){
                auto idx = i*h_BGMW95 +j;
                PRECOMPUTATION_POINTS_LIST_BGMW95[idx] = PRECOMPUTATION_POINTS_LIST_3nh[3*idx];
            }
        }
    }
    else{
        for(int i = 0; i< N_POINTS; ++i){
            blst_p2_affine qjQi = FIX_POINTS_LIST[i];
            for(int j = 0; j< h_BGMW95; ++j){
                auto idx = i*h_BGMW95 +j;
                PRECOMPUTATION_POINTS_LIST_BGMW95[idx] = qjQi;
                qjQi = single_scalar_multiplication(q_RADIX_PIPPENGER_VARIANT, qjQi);    
            }
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
            // ((m-1)*h_LEN_SCALAR + j)*N_POINTS + i;
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

blst_p2_affine pippenger_variant_q_over_5_CHES(uint256_t scalars_array[]){
    
    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;

    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints+1]; // add 1 slot redundancy for prefetch, 20221010 very important
    scalars[npoints] = 0;  //

    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];

    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [npoints];

    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, scalars_array[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){
            size_t idx = i*h_LEN_SCALAR + j;
            int m = ret_MB_expr[j][0];
            scalars[idx]  = ret_MB_expr[j][1];

            if (m> 0) {
                size_t idx_i_j_m =  3*idx + m-1; 
                points_ptr[idx] = PRECOMPUTATION_POINTS_LIST_3nh + idx_i_j_m;
                booth_signs[idx] = 0; 
            }
            else{
                size_t idx_i_j_m =  3*idx  -m - 1; 
                points_ptr[idx] = PRECOMPUTATION_POINTS_LIST_3nh + idx_i_j_m;
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


blst_p2_affine pippenger_variant_q_over_5_CHES_integral_scalar_conversion(uint256_t scalars_array[]){
    
    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints+2]; // add 2 slots redundancy for prefetch
    scalars[npoints] = 0; 
    scalars[npoints+1] = 0; 
    // Convert a scalar to its standard q-ary form
    std::array< int, h_LEN_SCALAR> ret_std_expr;
    int* scalars_p = scalars;
     for(int i = 0; i< N_POINTS; ++i){
        trans_uint256_t_to_standard_q_ary_expr(ret_std_expr, scalars_array[i]);
        for(int j = 0; j< h_LEN_SCALAR; ++j){
            *scalars_p++ = ret_std_expr[j];
        }
    }
    
    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];

    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [npoints];

    // convert scalar from standard q-ary form to MxB representation, and obtain the scalars, booth_signs, and points_ptr.
    blst_p2_construct_nh_scalars_nh_points(scalars, booth_signs, points_ptr, npoints,\
                            PRECOMPUTATION_POINTS_LIST_3nh, DIGIT_CONVERSION_HASH_TABLE);

    blst_p2xyzz* buckets;
    buckets = new blst_p2xyzz [B_SIZE];
    blst_p2 ret;

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

blst_p2_affine pippenger_variant_BGMW95(uint256_t scalars_array[]){
    
    std::array< int, h_BGMW95> ret_qhalf_expr;

    uint64_t npoints = N_POINTS*h_BGMW95;
    
    int* scalars;
    scalars = new int [npoints];
    
    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];
    
    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [npoints]; 

    // This is only for BLS12-381 curve
    if (N_EXP == 13 || N_EXP == 14 || N_EXP == 16 || N_EXP == 17){

        uint64_t  tt = uint64_t(1) << 62;

        for(int i = 0; i< N_POINTS; ++i){

            uint256_t aa = scalars_array[i];
           
            bool condition =  (aa.data[3] > tt); // a > 0.5*q*q**(h-1)
            if (condition == true) {
                aa = r_GROUP_ORDER - aa;
            }

            trans_uint256_t_to_qhalf_expr(ret_qhalf_expr, aa);
            
            if (condition == true) {
                for(int j = 0; j< h_BGMW95; ++j){
                    size_t idx = i*h_BGMW95 + j;
                    scalars[idx]  = ret_qhalf_expr[j];
                    points_ptr[idx] =  PRECOMPUTATION_POINTS_LIST_BGMW95 + idx;

                    // std::cout << "YES" << std::endl;
                    if ( scalars[idx] > 0) {
                        booth_signs[idx] = 1; 
                    }
                    else{
                        scalars[idx] = - scalars[idx];
                        booth_signs[idx] = 0; 
                    } 
                } 
            }
            else{
                for(int j = 0; j< h_BGMW95; ++j){
                    size_t idx = i*h_BGMW95 + j;
                    scalars[idx]  = ret_qhalf_expr[j];
                    points_ptr[idx] =  PRECOMPUTATION_POINTS_LIST_BGMW95 + idx;

                    if ( scalars[idx] > 0) {
                        booth_signs[idx] = 0; 
                    }
                    else{
                        scalars[idx] = - scalars[idx];
                        booth_signs[idx] = 1; 
                    } 
                }
            }    
        }
    }
    else{
        for(int i = 0; i< N_POINTS; ++i){
            trans_uint256_t_to_qhalf_expr(ret_qhalf_expr, scalars_array[i]);

            for(int j = 0; j< h_BGMW95; ++j){
                size_t idx = i*h_BGMW95 + j;
                scalars[idx]  = ret_qhalf_expr[j];
                points_ptr[idx] =  PRECOMPUTATION_POINTS_LIST_BGMW95 + idx;
                if ( scalars[idx] > 0) {
                    booth_signs[idx] = 0; 
                }
                else{
                    scalars[idx] = -scalars[idx];
                    booth_signs[idx] = 1; 
                }  
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

blst_p2_affine pippenger_blst_built_in(uint256_t scalars_array[]){

    blst_scalar* scalars;
    scalars = new blst_scalar [N_POINTS];

    uint8_t** scalars_ptr;
    scalars_ptr = new uint8_t* [N_POINTS];

    for(size_t i = 0; i < N_POINTS; ++i){
            blst_scalar_from_uint64( &scalars[i], scalars_array[i].data);
            scalars_ptr[i] = (scalars[i].b);
    }

    blst_p2_affine** points_ptr;
    points_ptr = new blst_p2_affine* [N_POINTS]; // points_ptr is an array of pointers that point to blst_p2_affine points.

    for(size_t i = 0; i < N_POINTS; ++i){
            points_ptr[i] = FIX_POINTS_LIST + i;
        }

    limb_t* SCRATCH;
    SCRATCH = new limb_t[blst_p2s_mult_pippenger_scratch_sizeof(N_POINTS)/sizeof(limb_t)];

    blst_p2 ret; // Mont coordinates
    size_t nbits = 255;

    blst_p2s_mult_pippenger(&ret, points_ptr, N_POINTS, scalars_ptr, nbits, SCRATCH);   

    delete[] SCRATCH;
    delete[] points_ptr;
    delete[] scalars_ptr;
    delete[] scalars;

    blst_p2_affine res_affine;
    blst_p2_to_affine( &res_affine, &ret);
    return res_affine;
}

void test_pippengers(){
    std::cout << "\nPIPPENGERS TEST OVER G2 for NPOINTS:  2**" << N_EXP << std::endl;

    size_t TEST_NUM = 5;
    size_t LOOP_NUM;

    if(N_EXP <= 8) LOOP_NUM = 40;
    else if(N_EXP <= 12) LOOP_NUM = 10;
    else if(N_EXP<= 16) LOOP_NUM = 5;
    else LOOP_NUM = 1;
    
    std::chrono::microseconds acc_t1, acc_t2, acc_t3, acc_t4, acc_t5, acc_conver_q_over_5, acc_conver_bgmw, diff, min_t12; // time accumulation
    acc_t1 = std::chrono::microseconds::zero();
    acc_t2 = std::chrono::microseconds::zero();
    acc_t3 = std::chrono::microseconds::zero();
    acc_t4 = std::chrono::microseconds::zero();
    acc_t5 = std::chrono::microseconds::zero();

    acc_conver_q_over_5 = std::chrono::microseconds::zero();
    acc_conver_bgmw = std::chrono::microseconds::zero();

    auto st = std::chrono::steady_clock::now();
    auto ed = std::chrono::steady_clock::now();   
    
    blst_p2_affine ret_P_affine_1, ret_P_affine_2, ret_P_affine_3, ret_P_affine_4, ret_P_affine_5;


    // Initialize SCALARS_ARRAY

    scalar_MB_expr ret_MB_expr;
    std::array< int, h_BGMW95> ret_qhalf_expr;

    for( int idx = 1; idx <= TEST_NUM; ++idx){

        uint256_t* SCALARS_ARRAY;
        SCALARS_ARRAY = new uint256_t[N_POINTS];
        std::cout << "This is No." << idx << " SCALARS_ARRAY." << std::endl;
        for(size_t i = 0; i < N_POINTS; ++i)\
            SCALARS_ARRAY[i] = random_scalar_less_than_r_SHA256();

        if(TEST_PIPPENGER_Q_OVER_5_CHES){
            /*nh + q/5 method */
            st = std::chrono::steady_clock::now();
            for(size_t i = 0; i< LOOP_NUM; ++i)
            {
                ret_P_affine_1 = pippenger_variant_q_over_5_CHES(SCALARS_ARRAY);
            }
            ed = std::chrono::steady_clock::now();   
            diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
            acc_t1 += diff;

            /* nh + q/5 method 2 */

            st = std::chrono::steady_clock::now();
            for(size_t i = 0; i< LOOP_NUM; ++i)
            {
                ret_P_affine_2 = pippenger_variant_q_over_5_CHES_integral_scalar_conversion(SCALARS_ARRAY);
            }
            ed = std::chrono::steady_clock::now();   
            diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
            acc_t2 += diff;
        }

        if(TEST_PIPPENGER_BGMW95){
            /*nh + q/2 method BGMW95*/
            st = std::chrono::steady_clock::now();
            for(size_t i = 0; i< LOOP_NUM; ++i)
            {
                ret_P_affine_3 = pippenger_variant_BGMW95(SCALARS_ARRAY);
            }
            ed = std::chrono::steady_clock::now();   
            diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
            acc_t3 += diff;
        }

        /*blst pippenger h(n+q/2) method*/
        st = std::chrono::steady_clock::now();    
        for(size_t i = 0; i< LOOP_NUM; ++i) 
        {
            ret_P_affine_4 = pippenger_blst_built_in(SCALARS_ARRAY);
        }
        ed = std::chrono::steady_clock::now(); 
        diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
        acc_t4 += diff;

        /* scalar conversion benchmark*/
        if(TEST_PIPPENGER_BGMW95){ 
            st = std::chrono::steady_clock::now();
            for(int j=0; j< LOOP_NUM; ++j){
                for(size_t i = 0; i< N_POINTS; ++i)
                    {
                    trans_uint256_t_to_qhalf_expr(ret_qhalf_expr, SCALARS_ARRAY[i]);     
                    }
            }
            ed = std::chrono::steady_clock::now();   
            diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
            acc_conver_bgmw += diff;
        }

        if(TEST_PIPPENGER_Q_OVER_5_CHES){ 
            st = std::chrono::steady_clock::now();
            for(int j=0; j< LOOP_NUM; ++j){
                for(size_t i = 0; i< N_POINTS; ++i)
                    {
                    trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, SCALARS_ARRAY[i]);   
                    }
            }
            ed = std::chrono::steady_clock::now();   
            diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
            acc_conver_q_over_5 += diff;
        }

        std::cout <<  "First scalar: " << SCALARS_ARRAY[0] << std::endl; 
        delete[] SCALARS_ARRAY;
    }

    if(TEST_PIPPENGER_Q_OVER_5_CHES){
        std::cout << "\n1. CHES 'nh+ q/5'. Wall clock time elapse is: " << std::dec << acc_t1.count()/(TEST_NUM*LOOP_NUM) << " us "<< std::endl;
        std::cout << ret_P_affine_1.x <<std::endl;
        std::cout << ret_P_affine_1.y <<std::endl;
        std::cout << std::endl;

        std::cout << "2. CHES 'nh+ q/5' integral scalar conversion. Wall clock time elapse is: " << std::dec << acc_t2.count()/(TEST_NUM*LOOP_NUM) << " us "<< std::endl;
        std::cout << ret_P_affine_2.x <<std::endl;
        std::cout << ret_P_affine_2.y <<std::endl;
        std::cout << std::endl;
    }

    if(TEST_PIPPENGER_BGMW95){
        std::cout << "3. pippenger_variant_BGMW95. Wall clock time elapse is: " << std::dec << acc_t3.count()/(TEST_NUM*LOOP_NUM) << " us "<< std::endl;
        std::cout << ret_P_affine_3.x <<std::endl;
        std::cout << ret_P_affine_3.y <<std::endl;
        std::cout << std::endl;  
    }

    std::cout << "4. pippenger_blst_built_in. Wall clock time elapse is: " << std::dec << acc_t4.count()/(TEST_NUM*LOOP_NUM) << " us "<< std::endl;
    std::cout << ret_P_affine_4.x <<std::endl;
    std::cout << ret_P_affine_4.y <<std::endl;
    std::cout << std::endl; 

    min_t12 = (acc_t1> acc_t2)? acc_t2 : acc_t1; 
    // min_t12 = (min_t12 > acc_t5)? acc_t5 : min_t12; // min_t12 = min(a1,a2,a5)

    if(TEST_PIPPENGER_BGMW95){
    std::cout << "Improvement, blst BGMW95 vs pipp: " <<\
    std::fixed << std::setprecision(3) << 100*float(acc_t4.count() - acc_t3.count())/float(acc_t4.count()) << '%' <<std::endl;
    }

    if(TEST_PIPPENGER_Q_OVER_5_CHES){
    std::cout << "Improvement, blst CHES_q_over_5 vs pipp: " <<\
    std::fixed << std::setprecision(3) << 100*float(acc_t4.count() - min_t12.count())/float(acc_t4.count()) << '%' <<std::endl;
    }

    if(TEST_PIPPENGER_Q_OVER_5_CHES && TEST_PIPPENGER_BGMW95){
    std::cout << "Improvement, blst CHES_q_over_5 vs BGMW95: " <<\
    std::fixed << std::setprecision(3) << 100*float(acc_t3.count() - min_t12.count())/float(acc_t3.count()) << '%' <<std::endl;
    }

    if(TEST_PIPPENGER_BGMW95){ 
    std::cout << "\nBGMW95 scalars conversion clock time elapse is: " <<\
    std::dec << acc_conver_bgmw.count()/(TEST_NUM*LOOP_NUM) << " us"<< std::endl;  
    std::cout << "It takes up what percentage of BGMW95 time: " <<\
    std::fixed << std::setprecision(3) << 100*float(acc_conver_bgmw.count())/float(acc_t3.count()) << '%' <<std::endl;
    }

    if(TEST_PIPPENGER_Q_OVER_5_CHES){
    std::cout << "CHES q_over_5 scalars conversion clock time elapse is: " <<\
    std::dec << acc_conver_q_over_5.count()/(TEST_NUM*LOOP_NUM) << " us"<< std::endl;
    std::cout << "It takes up what percentage of q_over_5 time: " << \
    std::fixed << std::setprecision(3) << 100*float(acc_conver_q_over_5.count())/float(min_t12.count()) << '%' <<std::endl;
    }
    std::cout << "\nTEST END" <<std::endl;    
}

int main(){
 
    init_fix_point_list();
    if(TEST_PIPPENGER_Q_OVER_5_CHES) init_pippenger_CHES_q_over_5();
    if(TEST_PIPPENGER_BGMW95) init_pippenger_BGMW95();
    
    // Test should be down between init_xx() and free_init_xx();
    test_pippengers();

    if(TEST_PIPPENGER_BGMW95) free_init_pippenger_BGMW95();
    if(TEST_PIPPENGER_Q_OVER_5_CHES) free_init_pippenger_CHES_q_over_5();
    free_init_fix_point_list();

    return 0;
}

