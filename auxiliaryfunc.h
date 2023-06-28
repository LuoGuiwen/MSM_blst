#include "src_from_aztec/numeric/uint256/uint256.hpp"
#include <openssl/rand.h>
#include <openssl/sha.h>

const uint256_t r_GROUP_ORDER = {
    uint64_t(0xffffffff00000001), uint64_t(0x53bda402fffe5bfe), uint64_t(0x3339d80809a1d805), uint64_t(0x73eda753299d7d48)
};

// a simple yet useful prefetch function
void vec_prefetch(const void *ptr, size_t len)
{   (void)ptr; (void)len;   }


/* Define ostreams for other related class*/
std::ostream& operator<<(std::ostream& os, const blst_fp& b)
{
    os << std::hex << "0x" << std::setfill('0') \
    << std::setw(16) << b.l[5] << " "<< std::setw(16) << b.l[4]<<  " "<< std::setw(16) << b.l[3] << " "\
    << std::setw(16) << b.l[2] <<" "<< std::setw(16) << b.l[1] <<" "<< std::setw(16) << b.l[0];
    return os;
}

std::ostream& operator<<(std::ostream& os, const blst_fp2& b)
{
    os << b.fp[0] <<" + " << std::endl\
    << b.fp[1] << "*i  ";
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

std::ostream& operator<<(std::ostream& os, const digit_decomposition b)
{
    os << std::dec << "[" << b.m << " "<< b.b<<  " "<< b.alpha<< "]";
    return os;
}


/*
Auxiliary Functions
*/

typedef std::array<std::array< int, 2>, h_LEN_SCALAR> scalar_MB_expr;

void print( const scalar_MB_expr &expr){

    std::cout << std::dec << "{";
    for(auto a : expr){
        std::cout <<"[ "<<a[0]<<", "<<a[1]<<" ]";
    }
    std::cout <<"}"<< std::endl;
}

/* scalar conversion */
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

    digit_decomposition tmp_tri;
    for (int i = 0; i< h_LEN_SCALAR; ++i){
        digit_decomposition tmp_tri = DIGIT_CONVERSION_HASH_TABLE[tmp_std_expr[i]];
        if(tmp_tri.alpha ==0){
            ret_MB_expr[i][0] = tmp_tri.m;
            ret_MB_expr[i][1] = tmp_tri.b;
        }

        else{
            ret_MB_expr[i][0] = - tmp_tri.m;
            ret_MB_expr[i][1] = tmp_tri.b;
            tmp_std_expr[i+1] += 1; 
        }
    }
}


void trans_uint256_t_to_standard_q_ary_expr_BGMW95( std::array<int, h_BGMW95> &ret_std_expr, const uint256_t &a){
    uint256_t tmp = a;
    uint32_t mask = (1 << EXPONENT_OF_q_BGMW95) - 1;
    for (int i=0; i< h_BGMW95; ++i){
        ret_std_expr[i] = tmp.data[0] & mask;// we only need the bit-wise xor with the last 32-bit of tmp.
        tmp = tmp >> EXPONENT_OF_q_BGMW95;
    }
}

void trans_uint256_t_to_qhalf_expr( std::array<int, h_BGMW95> &ret_qhalf_expr, const uint256_t &a){
    uint256_t tmp = a;
    int qhalf = int (q_RADIX_PIPPENGER_VARIANT>>1);
    uint32_t mask = uint32_t (q_RADIX_PIPPENGER_VARIANT - 1);
    for (int i=0; i< h_BGMW95; ++i){
        ret_qhalf_expr[i] = tmp.data[0] & mask;// we only need the bit-wise xor with the last 32-bit of tmp.
        tmp = tmp >> EXPONENT_OF_q_BGMW95;
    }
    for (int i=0; i< h_BGMW95 - 1; ++i){
            if(ret_qhalf_expr[i] > qhalf){
            ret_qhalf_expr[i] -= q_RADIX_PIPPENGER_VARIANT;
            ret_qhalf_expr[i+1] +=1;
            // system parameter makes sure ret_qhalf_expr[h-1]<= q/2.
        }
    }
}



/*Random Number Generator*/


// The following random function is only for experiment purpose.
// It is not cryptographically secure.

uint256_t random_scalar_less_than_r(){

    uint256_t ret = r_GROUP_ORDER;

    std::random_device rd;
    // Initialize Mersenne Twister pseudo-random number generator
    std::mt19937 gen(rd());

    // Generate pseudo-random uint32_t
    // uniformly distributed in range
    std::uniform_int_distribution<> dis(0, 0xffffffff); // '0xffffffff' = 2**32 -1

    while (ret >= r_GROUP_ORDER){ // make sure the scalar is less than r_GROUP_ORDER.
        ret.data[3] = (uint64_t(dis(gen))>>1) + (uint64_t(dis(gen))<<31);  // random number is 255 bit
        ret.data[2] = uint64_t(dis(gen)) + (uint64_t(dis(gen))<<32);
        ret.data[1] = uint64_t(dis(gen)) + (uint64_t(dis(gen))<<32);
        ret.data[0] = uint64_t(dis(gen)) + (uint64_t(dis(gen))<<32);
    }
    
    return ret;
}


uint256_t random_scalar_less_than_r_SHA256(){

    uint256_t ret = r_GROUP_ORDER;
    
    const size_t BUFFER_LENGTH = 32;  // 256 bits
    unsigned char buffer[BUFFER_LENGTH];

    uint64_t randomNumbers[4];

    while (ret >= r_GROUP_ORDER){ 

        unsigned char hash[SHA256_DIGEST_LENGTH];

        // Initialize the buffer with random data
        RAND_bytes(buffer, sizeof(buffer));

        // Compute the SHA256 hash of the buffer
        SHA256(buffer, sizeof(buffer), hash);

        // Store the hash in the array of uint64_t
        memcpy(randomNumbers, hash, sizeof(uint64_t) * 4);

        ret.data[3] = randomNumbers[3] >> 1;
        ret.data[2] = randomNumbers[2];
        ret.data[1] = randomNumbers[1];
        ret.data[0] = randomNumbers[0];
    }

    return ret;
}


void generateRandomNumber(uint64_t* numbers) {
    const size_t BUFFER_LENGTH = 32;  // 256 bits
    unsigned char buffer[BUFFER_LENGTH];
    unsigned char hash[SHA256_DIGEST_LENGTH];

    // Initialize the buffer with random data
    RAND_bytes(buffer, sizeof(buffer));

    // Compute the SHA256 hash of the buffer
    SHA256(buffer, sizeof(buffer), hash);

    // Store the hash in the array of uint64_t
    memcpy(numbers, hash, sizeof(uint64_t) * 4);
}

void printRandomNumber(const uint64_t* numbers) {
    std::cout << "Random Number: ";
    for (int i = 0; i < 4; ++i) {
        std::cout << std::hex << std::setw(16) << std::setfill('0') << numbers[i] << " ";
    }
    std::cout << std::dec << std::endl;
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
void construct_bucket_set(int bucket_set[], int q, int ah){

    std::set<int> B = {0, 1};
    
    for(int i = 2; i <= q/2; ++i){
        if (((omega2(i) + omega3(i))%2) == 0){
            B.insert(i);
        }
    }
    
    for(int i = q/4; i < q/2; ++i){
        if ((B.find(i) != B.end()) && (B.find(q - 2*i) != B.end()) ) // if i is in B and q-3*i is in B
        {
            B.erase(q - 2*i);
        }
    }
    for(int i = q/6; i < q/4; ++i){
        if ((B.find(i) != B.end()) && (B.find(q - 3*i) != B.end()) ) // if i is in B and q-3*i is in B
        {
            B.erase(q - 3*i);
        }
    }
    
    for(int i = 1; i <= ah+1; ++i){
        if (((omega2(i) + omega3(i))%2) == 0){
            B.insert(i);
        }
    }
    
    int index = 0;
    for(auto b: B){bucket_set[index] = b; ++index;};
}


void byte_str_from_uint32(uint8_t ret[4], const uint32_t a)
{
    uint32_t w = a;
    ret[0] = (byte)w;
    ret[1] = (byte)(w >> 8);
    ret[2] = (byte)(w >> 16);
    ret[3] = (byte)(w >> 24);
}



void vec_zero(void *ret, size_t num)
{
    volatile limb_t *rp = (volatile limb_t *)ret;
    size_t i;

    num /= sizeof(limb_t);

    for (i = 0; i < num; i++)
        rp[i] = 0;
}


/*functions originated from blst library*/

size_t pippenger_window_size(size_t npoints)
{
    size_t wbits;

    for (wbits=0; npoints>>=1; wbits++) ;

    return wbits>12 ? wbits-3 : (wbits>4 ? wbits-2 : (wbits ? 2 : 1));
}

/*
 * Window value encoding that utilizes the fact that -P is trivially
 * calculated, which allows to halve the size of pre-computed table,
 * is attributed to A. D. Booth, hence the name of the subroutines...
 */
limb_t booth_encode(limb_t wval, size_t sz)
{
    limb_t mask = 0 - (wval >> sz);     /* "sign" bit -> mask */

    wval = (wval + 1) >> 1;
    wval = (wval & ~mask) | ((0-wval) & mask);

    /* &0x1f, but <=0x10, is index in table, rest is extended "sign" bit */
    return wval;
}

/* Works up to 25 bits. */
limb_t get_wval_limb(const byte *d, size_t off, size_t bits)
{
    size_t i, top = (off + bits - 1)/8;
    limb_t ret, mask = (limb_t)0 - 1;

    d   += off/8;
    top -= off/8-1;

    /* this is not about constant-time-ness, but branch optimization */
    for (ret=0, i=0; i<4;) {
        ret |= (*d & mask) << (8*i);
        mask = (limb_t)0 - ((++i - top) >> (8*sizeof(top)-1));
        d += 1 & mask;
    }

    return ret >> (off%8);
}

