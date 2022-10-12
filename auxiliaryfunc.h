// a simple yet useful prefetch function
void vec_prefetch(const void *ptr, size_t len)
{   (void)ptr; (void)len;   }


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


uint256_t random_scalar_less_than_r(){
    
    uint256_t ret;

    std::random_device rd;
    // Initialize Mersenne Twister pseudo-random number generator
    std::mt19937 gen(rd());

    // Generate pseudo-random numbers
    // uniformly distributed in range
    std::uniform_int_distribution<> dis(0, 2147483647); // 2147483647 = 2**31 -1

    ret.data[7] = uint32_t(dis(gen)) % (0x73eda753); // make sure the scalar is less than r.
    ret.data[6] = uint32_t(dis(gen));
    ret.data[5] = uint32_t(dis(gen));
    ret.data[4] = uint32_t(dis(gen));
    ret.data[3] = uint32_t(dis(gen));
    ret.data[2] = uint32_t(dis(gen));
    ret.data[1] = uint32_t(dis(gen));
    ret.data[0] = uint32_t(dis(gen));

    return ret;
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

