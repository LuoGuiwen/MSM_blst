/*
n = 2**17
*/

// parameters for q/5
const int N_EXP = 17; 
const size_t N_POINTS = (size_t) 1<< N_EXP; 
constexpr int EXPONENT_OF_q = 19;
constexpr int q_RADIX = (int) (1 << EXPONENT_OF_q);
constexpr int h_LEN_SCALAR = 14;
constexpr int a_LEADING_TERM = 231;
constexpr int d_MAX_DIFF = 6;
constexpr int B_SIZE = 109244;

//parameters for q/2 (BGMW95) 
const int EXPONENT_OF_q_BGMW95 = 17;
const int q_RADIX_PIPPENGER_VARIANT = (int) (1 << EXPONENT_OF_q_BGMW95);
const int h_BGMW95 = 15;