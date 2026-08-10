/*
n = 2**10
*/

const int N_EXP = 10; 
const size_t N_POINTS = (size_t) 1<< N_EXP; 
constexpr int EXPONENT_OF_q = 13;
constexpr int q_RADIX = (int) (1 << EXPONENT_OF_q);
constexpr int h_LEN_SCALAR = 20;
constexpr int a_LEADING_TERM = 231;
constexpr int d_MAX_DIFF = 6;
constexpr int B_SIZE = 1725;

//parameters for q/2 (BGMW95) 
const int EXPONENT_OF_q_BGMW95 = 12;
const int q_RADIX_PIPPENGER_VARIANT = (int) (1 << EXPONENT_OF_q_BGMW95);
const int h_BGMW95 = 22;
 