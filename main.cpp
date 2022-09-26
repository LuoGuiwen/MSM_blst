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

/* 
Define 256 bit scalar with shift operations Guiwen, 
WARNING: shift should be no more than 32 !!! 
*/

class uint256_t {
    public:

    uint32_t data[8]; // the scalar is equal to data[0] + data[1]<< 32 + data[2] << 64 + ...
    
    constexpr uint256_t(const uint32_t a = 0)
        : data{ a, 0, 0, 0, 0, 0, 0, 0 }
    {}

    constexpr uint256_t(const uint32_t a0, const uint32_t a1, const uint32_t a2, const uint32_t a3, const uint32_t a4, const uint32_t a5, const uint32_t a6, const uint32_t a7)
        :data{a0, a1, a2, a3, a4, a5, a6, a7}
    {}

    constexpr uint256_t(const uint256_t& other)
        : data{ other.data[0], other.data[1], other.data[2], other.data[3], other.data[4], other.data[5], other.data[6], other.data[7] }
    {}

    constexpr uint256_t& operator=(const uint256_t& other) = default;
    explicit constexpr operator bool() const { return static_cast<bool>(data[0]); };

    // define right shift operation
    constexpr uint256_t operator>>(const uint32_t & total_shift) const // total_shift <= 31, because our radix is no more than 2**31
    {
        if (total_shift == 0) {
            return *this;
        }

        uint32_t limb_shift = total_shift;

        uint32_t shifted_limbs[8] = { 0 };

        uint32_t remainder_shift = 32UL - limb_shift;

        shifted_limbs[7] = data[7] >> limb_shift;

        uint32_t remainder = (data[7]) << remainder_shift; // elegant !!!!

        shifted_limbs[6] = (data[6] >> limb_shift) + remainder;

        remainder = (data[6]) << remainder_shift;

        shifted_limbs[5] = (data[5] >> limb_shift) + remainder;

        remainder = (data[5]) << remainder_shift;

        shifted_limbs[4] = (data[4] >> limb_shift) + remainder;

        remainder = (data[4]) << remainder_shift;

        shifted_limbs[3] = (data[3] >> limb_shift) + remainder;
        
        remainder = (data[3]) << remainder_shift;

        shifted_limbs[2] = (data[2] >> limb_shift) + remainder;
        
        remainder = (data[2]) << remainder_shift;

        shifted_limbs[1] = (data[1] >> limb_shift) + remainder;
        
        remainder = (data[1]) << remainder_shift;

        shifted_limbs[0] = (data[0] >> limb_shift) + remainder;
        
        uint256_t result(0);

        for (int i = 0; i < 8; ++i) {
            result.data[i] = shifted_limbs[i];
        }

        return result;
    };

    // define left shift operation
    constexpr uint256_t operator<<(const uint32_t & total_shift) const // total_shift <= 31, because our radix is no more than 2**31
    {
        if (total_shift == 0) {
            return *this;
        }

        uint32_t limb_shift = total_shift;

        uint32_t shifted_limbs[8] = { 0 };

        uint32_t remainder_shift = 32UL - limb_shift;

        shifted_limbs[0] = data[0] << limb_shift;

        uint32_t remainder = (data[0]) >> remainder_shift; // elegant !!!!

        shifted_limbs[1] = (data[1] << limb_shift) + remainder;

        remainder = (data[1]) >> remainder_shift;

        shifted_limbs[2] = (data[2] << limb_shift) + remainder;

        remainder = (data[2]) >> remainder_shift;

        shifted_limbs[3] = (data[3] << limb_shift) + remainder;

        remainder = (data[3]) >> remainder_shift;

        shifted_limbs[4] = (data[4] << limb_shift) + remainder;

        remainder = (data[4]) >> remainder_shift;
       
        shifted_limbs[5] = (data[5] << limb_shift) + remainder;

        remainder = (data[5]) >> remainder_shift;

        shifted_limbs[6] = (data[6] << limb_shift) + remainder;

        remainder = (data[6]) >> remainder_shift;

        shifted_limbs[7] = (data[7] << limb_shift) + remainder;
        
        uint256_t result(0);

        for (int i = 0; i < 8; ++i) {
            result.data[i] = shifted_limbs[i];
        }

        return result;
    };

    // define addition with a uint256_t.

    // compute a + b + carry, returning the carry
    constexpr std::pair<uint32_t, uint32_t> addc(const uint32_t a,
                                                        const uint32_t b,
                                                        const uint32_t carry_in) const
    {
        const uint32_t sum = a + b;
        const uint32_t carry_temp = sum < a;
        const uint32_t r = sum + carry_in;
        const uint32_t carry_out = carry_temp + (r < carry_in);
        return { r, carry_out };
    };

    constexpr uint32_t addc_discard_hi(const uint32_t a, const uint32_t b, const uint32_t carry_in) const
    {
        return a + b + carry_in;
    };

    constexpr uint256_t operator+(const uint256_t& other) const
    {
        const auto [r0, t0] = addc(data[0], other.data[0], 0);
        const auto [r1, t1] = addc(data[1], other.data[1], t0);
        const auto [r2, t2] = addc(data[2], other.data[2], t1);
        const auto [r3, t3] = addc(data[3], other.data[3], t2);
        const auto [r4, t4] = addc(data[4], other.data[4], t3);
        const auto [r5, t5] = addc(data[5], other.data[5], t4);
        const auto [r6, t6] = addc(data[6], other.data[6], t5);

        const auto r7 = addc_discard_hi(data[7], other.data[7], t6);
        return { r0, r1, r2, r3, r4, r5, r6, r7};
    };    

    constexpr bool operator>(const uint256_t& other) const
    {
    bool t0 = data[7] > other.data[7];
    bool t1 = data[7] == other.data[7] && data[6] > other.data[6];
    bool t2 = data[7] == other.data[7] && data[6] == other.data[6] && data[5] > other.data[5];
    bool t3 = data[7] == other.data[7] && data[6] == other.data[6] && data[5] == other.data[5] && data[4] > other.data[4];
    bool t4 = data[7] == other.data[7] && data[6] == other.data[6] && data[5] == other.data[5] && data[4] == other.data[4] && data[3] > other.data[3];
    bool t5 = data[7] == other.data[7] && data[6] == other.data[6] && data[5] == other.data[5] && data[4] == other.data[4] && data[3] == other.data[3] &&  data[2] > other.data[2];
    bool t6 = data[7] == other.data[7] && data[6] == other.data[6] && data[5] == other.data[5] && data[4] == other.data[4] && data[3] == other.data[3] &&  data[2] == other.data[2]  &&  data[1] > other.data[1];
    bool t7 = data[7] == other.data[7] && data[6] == other.data[6] && data[5] == other.data[5] && data[4] == other.data[4] && data[3] == other.data[3] &&  data[2] == other.data[2]  &&  data[1] == other.data[1]&&  data[0] > other.data[0];
    return t0 || t1 || t2 || t3 || t4 || t5 || t6 || t7;
    }
};

inline std::ostream& operator<<(std::ostream& os, uint256_t const& a)
{
    std::ios_base::fmtflags f(os.flags());
    os  << std::hex << "0x" << std::setfill('0') << std::setw(8) << a.data[7] << std::setw(8) << a.data[6]
        << std::setw(8) << a.data[5] << std::setw(8) << a.data[4] << std::setw(8) << a.data[3] << std::setw(8) << a.data[2]
        << std::setw(8) << a.data[1] << std::setw(8) << a.data[0];
    os.flags(f);
    return os;
}

/* Define ostreams for other related class*/
std::ostream& operator<<(std::ostream& os, const blst_fp& b)
{
    os << std::hex << "0x" << std::setfill('0') \
    << std::setw(16) << b.l[5] << " "<< std::setw(16) << b.l[4]<<  " "<< std::setw(16) << b.l[3] << " "\
    << std::setw(16) << b.l[2] <<" "<< std::setw(16) << b.l[1] <<" "<< std::setw(16) << b.l[0];
    return os;
}

std::ostream& operator<<(std::ostream& os, const blst_fr& b)
{
    os << std::hex << b.l[3] << " "<< b.l[2]<<  " "<< b.l[1]<< " "<< b.l[0] <<" ";
    return os;
}

std::ostream& operator<< (std::ostream& os, const uint8_t bt) {
    return os << unsigned(bt);
}

std::ostream& operator<<(std::ostream& os, const blst_scalar& scalar)
{
    os << std::hex << "0x" << " " << std::setfill('0') \
    << std::setw(2) << scalar.b[31] << " "<< std::setw(2) << scalar.b[30] << " " << std::setw(2) << scalar.b[29] << " "<< std::setw(2) <<scalar.b[28] <<" "\
    << std::setw(2)  <<scalar.b[27] << " "<< std::setw(2) <<scalar.b[26] <<" "<< std::setw(2) <<scalar.b[25] << " " << std::setw(2) <<scalar.b[24] << ", " \
    << std::setw(2) << scalar.b[23] << ' '<<std::setw(2) <<  scalar.b[22] << ' ' << std::setw(2) << scalar.b[21] << " "<< std::setw(2) << scalar.b[20] <<" "\
    << std::setw(2) << scalar.b[19] << " "<< std::setw(2) << scalar.b[18] <<" "<< std::setw(2) << scalar.b[17] << " " << std::setw(2) << scalar.b[16] << ", " \
    << std::setw(2) << scalar.b[15] << ' '<< std::setw(2) << scalar.b[14] << ' ' << std::setw(2) << scalar.b[13] << " "<< std::setw(2) << scalar.b[12] <<" "\
    << std::setw(2) << scalar.b[11] << " "<< std::setw(2) << scalar.b[10] <<" "<< std::setw(2) << scalar.b[9] << " " << std::setw(2) << scalar.b[8] << ", " \
    << std::setw(2) << scalar.b[7] << ' '<< std::setw(2) << scalar.b[6] << ' ' << std::setw(2) << scalar.b[5] << " "<< std::setw(2) << scalar.b[4] <<" "\
    << std::setw(2) << scalar.b[3] << " "<< std::setw(2) << scalar.b[2] <<" "<< std::setw(2) << scalar.b[1] << " " << std::setw(2) << scalar.b[0] << " ";
    return os;
}

std::ostream& operator<<(std::ostream& os, const uint64_t* b)
{
    os << std::hex << b[3] << " "<< b[2]<<  " "<< b[1]<< " "<< b[0] <<" ";
    return os;
}

/***----***/
/* Define global variables and their initialization*/
/***----***/

const size_t N_POINTS = (size_t) (1 << 18);  
const blst_p1 G1_BLST_DEFAULT_GENERATOR = *blst_p1_generator(); // Default generator in blst_p1
const int EXPONENT_OF_q = 14;
const int q_RADIX = (int) (1 << EXPONENT_OF_q);
const int h_LEN_SCALAR = 19;
const int ah_LEADING = 7;
const int d_MAX_DIFF = 60;

std::set<int> MULTI_SET = {1, 2, 3};
std::set<int> BUCKET_SET;
std::vector<int> BUCKET_SET_ASCEND_VEC;

/* This once leads to a big BUG. */
std::array<int, q_RADIX/2 +1> BUCKET_VALUE_TO_ITS_INDEX; //Consider the maximum element in the BUCKET_SET, which is very very important

std::array<std::array<int,3>, q_RADIX+1>  DIGIT_CONVERSION_HASH_TABLE;

std::array< blst_p1_affine, N_POINTS> *FIX_POINTS_LIST; 
std::array< std::array< std::array< blst_p1_affine, 6>, h_LEN_SCALAR >, N_POINTS> *PRECOMPUTATION_POINTS_LIST_6nh; // Define the pointer then use new to allocate memory in heap, since this array is too big.

std::array< uint256_t, N_POINTS> *SCALARS_LIST;

// std::array<blst_p1_affine, N_POINTS*h_LEN_SCALAR> *POINT_INTERMEDIATE;

/* initialization later on using init() */
blst_fr FR_ONE;
blst_fp FP_ONE, FP_MONT_ONE, FP_Z;

blst_p1 G1_GENERATOR, G1_INFINITY; 
blst_p1_affine G1_GENERATOR_AFFINE, G1_AFFINE_INFINITY; 



/*
Functions
*/


// uint256_t random_scalar_less_than_r(){
//     uint256_t ret;
//     ret.data[7] = uint32_t(rand() % (0x73eda753)); // make sure the scalar is less than r.
//     ret.data[6] = uint32_t(rand());
//     ret.data[5] = uint32_t(rand());
//     ret.data[4] = uint32_t(rand());
//     ret.data[3] = uint32_t(rand());
//     ret.data[2] = uint32_t(rand());
//     ret.data[1] = uint32_t(rand());
//     ret.data[0] = uint32_t(rand());

//     return ret;
// }


uint256_t random_scalar_less_than_r(){

    uint256_t ret;
    auto seed = time(NULL);
    srand((unsigned) seed);
    ret.data[7] = uint32_t(rand() % (0x73eda753)); // make sure the scalar is less than r.

    srand((unsigned) seed+1);
    ret.data[6] = uint32_t(rand());

    srand((unsigned) seed+2);
    ret.data[5] = uint32_t(rand());

    srand((unsigned) seed+3);
    ret.data[4] = uint32_t(rand());

    srand((unsigned) seed+4);
    ret.data[3] = uint32_t(rand());

    srand((unsigned) seed+5);
    ret.data[2] = uint32_t(rand());

    srand((unsigned) seed+6);
    ret.data[1] = uint32_t(rand());

    srand((unsigned) seed+7);
    ret.data[0] = uint32_t(rand());

    return ret;
}

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

    BUCKET_SET = construct_bucket_set(q_RADIX, ah_LEADING);
    std::copy(BUCKET_SET.begin(), BUCKET_SET.end(), std::back_inserter(BUCKET_SET_ASCEND_VEC));

    std::cout<< "BUCKET_SET initialized. The size of BUCKET_SET is: " << BUCKET_SET.size() << std::endl;

    for(size_t i = 0; i < BUCKET_SET.size(); ++i){
        BUCKET_VALUE_TO_ITS_INDEX[BUCKET_SET_ASCEND_VEC[i]] = i;
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

    // Initialize SCALARS_LIST
    SCALARS_LIST = new std::array< uint256_t, N_POINTS>;
        for(size_t i = 0; i < N_POINTS; ++i){
            (*SCALARS_LIST)[i] = random_scalar_less_than_r();
            //  std::cout<< (*SCALARS_LIST)[i] <<std::endl;// it always generates the same psudo-random number.
        }
      std::cout<< "SCALARS_LIST Generated" <<std::endl;


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
        // ret1 = single_scalar_multiplication((*SCALARS_LIST)[i], tmp);
        blst_scalar_from_uint32( &scalar, (*SCALARS_LIST)[i].data);
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


void init_PRECOMPUTATION_POINTS_LIST_6nh(){
    /*
    ### Initialize the precomputation ###
    */
    PRECOMPUTATION_POINTS_LIST_6nh = new std::array< std::array< std::array< blst_p1_affine, 6>, h_LEN_SCALAR >, N_POINTS>;

    auto st = std::chrono::steady_clock::now();
    blst_p1_affine Pt;

    for(int i = 0; i< N_POINTS; ++i){
        blst_p1_affine qjQi = (*FIX_POINTS_LIST)[i];
        for(uint j = 0; j< h_LEN_SCALAR; ++j){
            for(uint m= 0;  m < 3; ++m){
                (*PRECOMPUTATION_POINTS_LIST_6nh)[i][j][m] = single_scalar_multiplication( m+1, qjQi);
                auto tmp_Pa = (*PRECOMPUTATION_POINTS_LIST_6nh)[i][j][m];
                blst_fp_sub(&tmp_Pa.y, &FP_Z, &tmp_Pa.y);
                (*PRECOMPUTATION_POINTS_LIST_6nh)[i][j][5 - m] = tmp_Pa;
            }
            qjQi = single_scalar_multiplication(q_RADIX, qjQi);    
        }
    }
    auto ed = std::chrono::steady_clock::now();   

    std::chrono::microseconds diff = std::chrono::duration_cast<std::chrono::microseconds>(ed -st);
    std::cout<< "PRECOMPUTATION_POINTS_LIST_6nh SUCCESSFULLY CONSTRUCTED" << std::endl;
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
blst_p1 accumulation_d(const std::vector<int>& bucketSetList, const std::vector<blst_p1>& intermediateSumList, const int d_max){
    if (bucketSetList.size()!= intermediateSumList.size()){
        std::cout << "ERROR: lengths do not match" <<std::endl;
        return intermediateSumList[0];
    }
    blst_p1 tmp = G1_INFINITY;

    //# We don't use the tmp_d[0]
    std::vector<blst_p1> tmp_d (d_max + 1,  G1_INFINITY); 
    // for(uint i = 0; i<= max_d; ++i) tmp_d.push_back(tmp); // those two initializations are equivalent

    for( auto i = bucketSetList.size()-1; i>0; --i){
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

blst_p1_affine pippenger_variant_submission_CHES(){
    
    std::cout <<"pippenger_variant_submission_CHES() with accumulation_d" << std::endl;

    blst_p1_affine tmp_Pa = G1_AFFINE_INFINITY;
    blst_p1 tmp_P = G1_INFINITY;

    std::vector<blst_p1> intermediate_sum_list (BUCKET_SET.size(), G1_INFINITY); 

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;
   
    for(int i = 0; i< N_POINTS; ++i){

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, (*SCALARS_LIST)[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){

            int m = ret_MB_expr[j][0];
            int b = ret_MB_expr[j][1];
            int index_b = BUCKET_VALUE_TO_ITS_INDEX[b];

            if (m>= 0) {
                tmp_Pa = (*PRECOMPUTATION_POINTS_LIST_6nh)[i][j][ m-1 ];
            }
            else{
                tmp_Pa = (*PRECOMPUTATION_POINTS_LIST_6nh)[i][j][ -m-1 ];
                blst_fp_sub(&tmp_Pa.y, &FP_Z, &tmp_Pa.y);
            }
           
            blst_p1_from_affine(&tmp_P, &tmp_Pa);
            blst_p1_add_or_double(&(intermediate_sum_list[index_b]), &(intermediate_sum_list[index_b]), &tmp_P);
        }
    }
    auto res = accumulation_d(BUCKET_SET_ASCEND_VEC, intermediate_sum_list, d_MAX_DIFF);
    blst_p1_affine res_affine;
    blst_p1_to_affine( &res_affine, &res);

    return res_affine;
}




// Use  blst_p1s_tile_pippenger to implement the CHES 
blst_p1_affine pippenger_variant_submission_CHES_blst_tile(){
    
    std::cout <<"pippenger_variant_submission_CHES_blst_tile() with blst_p1s_tile_pippenger" << std::endl;

    // blst_p1_affine tmp_Pa = G1_AFFINE_INFINITY;
    // blst_p1 tmp_P = G1_INFINITY;

    std::array<std::array< int, 2>, h_LEN_SCALAR> ret_MB_expr;

    uint64_t npoints = N_POINTS*h_LEN_SCALAR;

    blst_scalar* INTERMEDIATE_SCALARS;
    INTERMEDIATE_SCALARS = new  blst_scalar[npoints];

    uint8_t** scalars_ptr;
    scalars_ptr = new uint8_t* [npoints];

    blst_p1_affine** points_ptr;
    points_ptr = new blst_p1_affine* [npoints]; // points_ptr is an array of pointers that point to blst_p1_affine points.
    
    for(int i = 0; i< N_POINTS; ++i){

        // std::cout << (*SCALARS_LIST)[i]<< std::endl;

        trans_uint256_t_to_MB_radixq_expr(ret_MB_expr, (*SCALARS_LIST)[i]);

        for(int j = 0; j< h_LEN_SCALAR; ++j){
            size_t idx = i*h_LEN_SCALAR + j;
            int m = ret_MB_expr[j][0];

            uint32_t b_value[8] = { 0, 0, 0, 0, 0, 0, 0, 0};
            b_value[0] = uint32_t(ret_MB_expr[j][1]);
            blst_scalar_from_uint32(&INTERMEDIATE_SCALARS[idx], b_value);
            scalars_ptr[idx] = INTERMEDIATE_SCALARS[idx].b;
            // scalars_ptr[idx] = INTERMEDIATE_SCALARS[1].b;
          
            points_ptr[idx] = (m > 0)? &(*PRECOMPUTATION_POINTS_LIST_6nh)[i][j][ m-1 ]: \
            &(*PRECOMPUTATION_POINTS_LIST_6nh)[i][j][ 6 + m ];
        }
    }
    blst_p1 ret_P; // Mont coordinates
    size_t window = EXPONENT_OF_q;

    limb_t* SCRATCH;
    SCRATCH = new limb_t[blst_p1s_mult_pippenger_scratch_sizeof(npoints*8)/sizeof(limb_t)];
    blst_p1s_tile_pippenger(&ret_P, points_ptr,\
                            npoints, scalars_ptr,\
                            window,  SCRATCH,\
                            0, window);
    
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
            blst_scalar_from_uint32( &scalars[i], (*SCALARS_LIST)[i].data);
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

    // std::cout << float(diff1.count() - 2*diff2.count())/float(diff1.count() - diff2.count()) <<std::endl;
    // /*nh + q/5 method*/
    // st = std::chrono::steady_clock::now();
    // for(size_t i = 0; i< TEST_NUM; ++i)
    // {
    //     ret_P_affine = pippenger_variant_submission_CHES();
    // }
    // ed = std::chrono::steady_clock::now();   
    // diff = std::chrono::duration_cast<std::chrono::microseconds>(ed - st);
    // std::cout << "CHES 'nh+ q/5' method Wall clock time elapse is: " << std::dec << diff.count()/TEST_NUM << " us "<< std::endl;
    // std::cout << ret_P_affine.x <<std::endl;
    // std::cout << ret_P_affine.y <<std::endl;






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
    trans_uint256_t_to_MB_radixq_expr(ret_MB_expr,(*SCALARS_LIST)[i]);   
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

size_t pippenger_window_size(size_t npoints)
{
    size_t wbits;

    for (wbits=0; npoints>>=1; wbits++) ;

    return wbits>12 ? wbits-3 : (wbits>4 ? wbits-2 : (wbits ? 2 : 1));
}

int main(){


    init();
    init_PRECOMPUTATION_POINTS_LIST_6nh();
    
    test();

    test_scalar_conversion_bench();
    // test_three_scalar_multiplication_methods();// Correct!


    blst_p1_affine Q = *blst_p1_affine_generator();
    blst_p1_affine Qn = Q;

    blst_fp_sub(&Qn.y, &FP_Z, &Qn.y);
    blst_p1 ret, xyzQ, xyzQn;

    blst_p1_from_affine(&xyzQ, &Q);

    blst_p1_from_affine(&xyzQn, &Qn);   
    blst_p1_add_or_double(&ret, &xyzQn, &xyzQ);


    std::cout << ret.x << std::endl;
    std::cout << ret.y << std::endl;
    std::cout << ret.z << std::endl;

    for(int e = 10; e<= 30; ++e)
    std::cout << std::dec << e<< "  "<< pippenger_window_size( (size_t) (1<< e) -8) << std::endl;

    return 0;
}

