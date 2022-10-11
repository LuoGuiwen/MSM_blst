/*
n = 2**12
*/

const int N_EXP = 12; 
const size_t N_POINTS = (size_t) 1<< N_EXP; 
constexpr int EXPONENT_OF_q = 14;
constexpr int q_RADIX = (int) (1 << EXPONENT_OF_q);
constexpr int h_LEN_SCALAR = 19;
constexpr int a_LEADING_TERM = 7;
constexpr int d_MAX_DIFF = 6;
constexpr int B_SIZE = 3417;



//parameters for q/2 (BGMW95) 
const int EXPONENT_OF_q_BGMW95 = 13;
const int q_RADIX_PIPPENGER_VARIANT = (int) (1 << EXPONENT_OF_q_BGMW95);
const int h_BGMW95 = 20;
 