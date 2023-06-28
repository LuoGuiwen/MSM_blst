/*
n = 2**24
*/

const int N_EXP = 24; 
const size_t N_POINTS = (size_t) 1<< N_EXP; 
constexpr int EXPONENT_OF_q = 26;
constexpr int q_RADIX = (int) (1 << EXPONENT_OF_q);
constexpr int h_LEN_SCALAR = 10;
constexpr int a_LEADING_TERM = 1899369;
constexpr int d_MAX_DIFF = 6;
constexpr int B_SIZE = 14139299;

//parameters for q/2 (BGMW95) 
const int EXPONENT_OF_q_BGMW95 = 24;
const int q_RADIX_PIPPENGER_VARIANT = (int) (1 << EXPONENT_OF_q_BGMW95);
const int h_BGMW95 = 11;