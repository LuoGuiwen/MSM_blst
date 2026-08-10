//



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


/*Correctness has been verified*/
blst_p1_affine pippenger_variant_CHES_primitive_implementation(uint256_t scalars_array[]){
    
    blst_p1xyzz* buckets;
    buckets = new blst_p1xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 
  

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;
   
    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, scalars_array[i]);

        int j = 0;

        int booth_idx = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j][1]];
        int booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j+1][1]];
        unsigned char booth_sign;
        blst_p1_affine tmp_Pa;
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
        blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
        
        ++j;
   
 
        for( ; j< h_LEN_SCALAR-1; ++j){
            m = ret_MB_expr[j][0];
            booth_idx = booth_idx_nxt;
            booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[ret_MB_expr[j+1][1]];
            blst_p1_prefetch_CHES(buckets, booth_idx_nxt);

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
            blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
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
        blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    }
    
    blst_p1 ret;
    blst_p1_affine res_affine;

    blst_p1_integrate_buckets_accumulation_d_CHES(&ret, buckets, BUCKET_SET, B_SIZE, d_MAX_DIFF);
    blst_p1_to_affine( &res_affine, &ret);

    delete[] buckets;

    return res_affine;
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
        blst_p1_prefetch_CHES(buckets, booth_idx_nxt);

        point = *points++; 
        booth_sign = *booth_signs++;
        if(booth_idx) blst_p1_bucket_CHES(buckets, booth_idx, point,  booth_sign);
    }
    point = *points;
    booth_sign = *booth_signs;
    blst_p1_bucket_CHES(buckets, booth_idx_nxt, point, booth_sign);// Carefully, it must be booth_idx_nxt
    blst_p1_integrate_buckets_accumulation_d_CHES(ret, buckets, bucket_set_ascend, bucket_set_size, d_MAX_DIFF);

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


blst_p1_affine pippenger_variant_BGMW95_blst_tile_function(uint256_t scalars_array[]){
    
    std::array< int, h_BGMW95> ret_std_expr;

    uint64_t npoints = N_POINTS*h_BGMW95;

    uint8_t* scalars;
    scalars = new uint8_t [npoints*4];// it is required that radix is less than 2^30
    uint8_t** scalars_ptr;
    scalars_ptr = new uint8_t* [npoints];

    blst_p1_affine** points_ptr;
    points_ptr = new blst_p1_affine* [npoints]; 

    for(int i = 0; i< N_POINTS; ++i){
        trans_uint256_t_to_standard_q_ary_expr_BGMW95(ret_std_expr, scalars_array[i]);

        for(int j = 0; j< h_BGMW95; ++j){
            size_t idx = i*h_BGMW95 + j;
            byte_str_from_uint32(&scalars[idx*4],ret_std_expr[j]);
            scalars_ptr[idx] = &scalars[idx*4];
            points_ptr[idx] =  &PRECOMPUTATION_POINTS_LIST_BGMW95[idx];
        }
    }

    blst_p1 ret; // Mont coordinates

    size_t window_width = EXPONENT_OF_q_BGMW95;
    limb_t* scratch;
    scratch = new limb_t [(q_RADIX_PIPPENGER_VARIANT)*sizeof(blst_p1xyzz)/sizeof(limb_t)];

    blst_p1s_tile_pippenger(&ret, points_ptr, npoints, scalars_ptr, window_width, scratch, 0, window_width+1);

    delete[] scratch;
    delete[] points_ptr;  
    delete[] scalars_ptr;  
    delete[] scalars;  

    blst_p1_affine res_affine;
    blst_p1_to_affine( &res_affine, &ret);
    return res_affine;
}



blst_p1_affine pippenger_variant_q_over_5_CHES_noindexhash(uint256_t scalars_array[]){
    
    // try to use direct hash to obtain the booth_idxx. 

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;

    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints];
    unsigned char* booth_signs; // it acts as a bool type
    booth_signs = new unsigned char [npoints];

    blst_p1_affine** points_ptr;
    points_ptr = new blst_p1_affine* [npoints]; 

    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, scalars_array[i]);

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
    blst_p1 ret; // Mont coordinates

    blst_p1xyzz* buckets;
    buckets = new blst_p1xyzz [q_RADIX/2 +1];

    blst_p1_tile_pippenger_d_CHES_noindexhash(&ret, \
                                    points_ptr, \
                                    npoints, \
                                    scalars, booth_signs, \
                                    buckets, BUCKET_SET,\
                                    B_SIZE, d_MAX_DIFF);

    delete[] buckets;
    delete[] points_ptr;
    delete[] booth_signs;    
    delete[] scalars;    
    blst_p1_affine res_affine;
    blst_p1_to_affine( &res_affine, &ret);

    return res_affine;
}



/* Correctess verified*/
blst_p1_affine pippenger_variant_q_over_5_CHES_prefetch_1step_ahead(uint256_t scalars_array[]){
    
    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints+2];
    int* scalar_p = scalars;

    std::array< int, h_LEN_SCALAR> ret_std_expr;
    size_t s_idx = 0;
    for(int i = 0; i< N_POINTS; ++i){
        trans_uint256_t_to_standard_q_ary_expr(ret_std_expr, scalars_array[i]);
        for(int j = 0; j< h_LEN_SCALAR; ++j){
            // std::cout << *scalars++ << std::endl;
            scalars[s_idx++] = ret_std_expr[j];
        }
    }

    int scalar_now, point_idx, point_idx_nxt, scalar_nxt, booth_idx, booth_idx_nxt, mul, mul_nxt, b, b_nxt;
    unsigned char booth_sign, booth_sign_nxt;
    blst_p1_affine tmp_Pa;
    blst_p1xyzz* buckets;
    buckets = new blst_p1xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 
  
    digit_decomposition tmp_tri, tmp_tri_nxt;
    size_t size_tri = sizeof(tmp_tri);
    size_t size_point = sizeof(blst_p1_affine);

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
    blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);

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
        if(booth_idx) blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    }
    tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt];
    if(booth_idx_nxt) blst_p1xyzz_dadd_affine(&buckets[booth_idx_nxt], &buckets[booth_idx_nxt], &tmp_Pa, booth_sign_nxt);

    blst_p1 ret; // Mont coordinates
    blst_p1_affine res_affine;
    blst_p1_integrate_buckets_accumulation_d_CHES(&ret, buckets, BUCKET_SET, B_SIZE, d_MAX_DIFF);
    blst_p1_to_affine( &res_affine, &ret);
    
    delete[] buckets; 
    delete[] scalars;    

    return res_affine;
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



void pippenger_tile_CHES_prefetch_2step_ahead_input_std_scalar(blst_p1 *ret, \
                                    const blst_p1_affine PRECOMPUTATION_POINTS_LIST_3nh[], \
                                    size_t npoints, \
                                    int scalars[], digit_decomposition DIGIT_CONVERSION_HASH_TABLE[], blst_p1xyzz buckets[],\
                                    int bucket_set_ascend[], int BUCKET_VALUE_TO_ITS_INDEX[],\
                                    size_t bucket_set_size, int d_max){
    
    int  point_idx, point_idx_nxt, booth_idx, booth_idx_nxt;
    unsigned char booth_sign, booth_sign_nxt;
    
    blst_p1_affine tmp_Pa;
  
    digit_decomposition tmp_tri, tmp_tri_nxt, tmp_tri_nxt2;
    size_t size_tri = sizeof(tmp_tri);
    size_t size_point = sizeof(blst_p1_affine);

    // point to the beginning of scalars array.
    int* scalars_p = scalars;
    int i = 0;

    tmp_tri = DIGIT_CONVERSION_HASH_TABLE[*scalars_p++];
    booth_idx = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri.b];
    booth_sign = tmp_tri.alpha;
    if(tmp_tri.alpha) ++(*scalars_p); // if tmp_tri[2] == 1
    // i = 0
    point_idx = 3*i + tmp_tri.m - 1;
    tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx];

    i = 1;
    tmp_tri_nxt = DIGIT_CONVERSION_HASH_TABLE[*scalars_p++];
    booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt.b];  
    booth_sign_nxt = tmp_tri_nxt.alpha;  
    if(tmp_tri_nxt.alpha) ++(*scalars_p);
    // i = 1
    point_idx_nxt = 3*i + tmp_tri_nxt.m - 1;

    vec_prefetch(&PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt], size_point);
    blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    
    i = 2;
    tmp_tri_nxt2 = DIGIT_CONVERSION_HASH_TABLE[*scalars_p++];

    while(i < npoints){

        tmp_tri = tmp_tri_nxt;
        booth_idx = booth_idx_nxt;
        booth_sign = booth_sign_nxt;
        point_idx = point_idx_nxt;

        // i == 2
        tmp_tri_nxt = tmp_tri_nxt2;
        booth_idx_nxt = BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt.b];  
        booth_sign_nxt = tmp_tri_nxt.alpha;  
        if(tmp_tri_nxt.alpha)  ++(*scalars_p);
        point_idx_nxt = 3*i + tmp_tri_nxt.m - 1;

        ++i;
        //i == 3
        tmp_tri_nxt2 = DIGIT_CONVERSION_HASH_TABLE[*scalars_p++];

        vec_prefetch(&BUCKET_VALUE_TO_ITS_INDEX[tmp_tri_nxt2.b], 4);
        vec_prefetch(&PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt], size_point);
        vec_prefetch(&buckets[booth_idx_nxt], size_point);
        vec_prefetch(&DIGIT_CONVERSION_HASH_TABLE[*scalars_p], size_tri);

        tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx];
        if(booth_idx) blst_p1xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], &tmp_Pa, booth_sign);
    }
    tmp_Pa = PRECOMPUTATION_POINTS_LIST_3nh[point_idx_nxt];
    if(booth_idx_nxt) blst_p1xyzz_dadd_affine(&buckets[booth_idx_nxt], &buckets[booth_idx_nxt], &tmp_Pa, booth_sign_nxt);

    blst_p1_integrate_buckets_accumulation_d_CHES(ret, buckets, bucket_set_ascend, bucket_set_size, d_max);
}



void construct_nh_scalars_nh_points(int nh_scalars[], unsigned char booth_signs[], blst_p1_affine* nh_points_ptr[], const size_t npoints, blst_p1_affine PRECOMPUTATION_POINTS_LIST_3nh[], const digit_decomposition DIGIT_CONVERSION_HASH_TABLE[]){

    int* scalars_p = nh_scalars;
    unsigned char* booth_signs_p = booth_signs;
    blst_p1_affine** points_p = nh_points_ptr;
    
    digit_decomposition tmp_tri;

    size_t point_idx;
    int i = 0;
    for(  ; i< npoints -1; ++i){
        tmp_tri = DIGIT_CONVERSION_HASH_TABLE[*scalars_p];
        *scalars_p++ = tmp_tri.b;
        *booth_signs_p++ = tmp_tri.alpha;
        if(tmp_tri.alpha) ++(*scalars_p); // if tmp_tri[2] == 1
        point_idx = 3*i + tmp_tri.m - 1; // m is always greater than 0.
        *points_p++ = &PRECOMPUTATION_POINTS_LIST_3nh[point_idx];
    }

    tmp_tri = DIGIT_CONVERSION_HASH_TABLE[*scalars_p];
    *scalars_p = tmp_tri.b;
    *booth_signs_p = tmp_tri.alpha;
    point_idx = 3*i + tmp_tri.m - 1; 
    *points_p = &PRECOMPUTATION_POINTS_LIST_3nh[point_idx];
    }


/* Correctess verified*/
blst_p1_affine pippenger_variant_q_over_5_CHES_prefetch_2step_ahead(uint256_t scalars_array[]){
    
    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    int* scalars;
    scalars = new int [npoints+2]; // add 2 slot redundancy for prefetch
    scalars[npoints] = 0; // initialization
    scalars[npoints+1] = 0;

    std::array< int, h_LEN_SCALAR> ret_std_expr;

    int idx;
    for(int i = 0; i< N_POINTS; ++i){
        trans_uint256_t_to_standard_q_ary_expr(ret_std_expr, scalars_array[i]);
        for(int j = 0; j< h_LEN_SCALAR; ++j){
            idx = i*h_LEN_SCALAR +j;
            scalars[idx] = ret_std_expr[j];
        }
    }

    blst_p1xyzz* buckets;
    buckets = new blst_p1xyzz [B_SIZE];
    vec_zero(buckets, sizeof(buckets[0])*B_SIZE); 

    blst_p1 ret;

    // here the scalar's standard q-ary representation int scalars[], and the precomputation array PRECOMPUTATION_POINTS_LIST_3nh is
    // directly input to the function, the MB conversion and booth_sign are dealt within the function. See multi_scalar.c for the code.
    blst_p1_tile_pippenger_CHES_prefetch_2step_ahead_input_std_scalar(&ret, \
                                    PRECOMPUTATION_POINTS_LIST_3nh, \
                                    npoints, \
                                    scalars,  DIGIT_CONVERSION_HASH_TABLE,\
                                    buckets,\
                                    BUCKET_SET, BUCKET_VALUE_TO_ITS_INDEX,\
                                    B_SIZE, d_MAX_DIFF);
  
    delete[] buckets; 
    delete[] scalars;    

    blst_p1_affine res_affine;
    blst_p1_to_affine( &res_affine, &ret);
    return res_affine;
}

/* initialization later on using init() */
blst_fr FR_ONE;
blst_fp FP_ONE, FP_MONT_ONE, FP_Z;

blst_p1 G1_GENERATOR, G1_INFINITY; 
blst_p1_affine G1_GENERATOR_AFFINE, G1_AFFINE_INFINITY; 


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
