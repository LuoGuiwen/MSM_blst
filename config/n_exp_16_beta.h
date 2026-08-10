/*
n = 2**16
*/

const int N_EXP = 16; 
const size_t N_POINTS = (size_t) 1<< N_EXP; 
constexpr int EXPONENT_OF_q = 18;
constexpr int q_RADIX = (int) (1 << EXPONENT_OF_q);
constexpr int h_LEN_SCALAR = 15;
constexpr int a_LEADING_TERM = 7;
constexpr int d_MAX_DIFF = 6;
constexpr int B_SIZE = 54618;

//parameters for q/2 (BGMW95) 
const int EXPONENT_OF_q_BGMW95 = 17;
const int q_RADIX_PIPPENGER_VARIANT = (int) (1 << EXPONENT_OF_q_BGMW95);
const int h_BGMW95 = 15;