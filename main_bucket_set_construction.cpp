/*

compile with: g++ -std=c++17 -o main_test -g -O2 main_p1.cpp libblst.a

use the following code in command line to unleash the stack restriction:

ulimit -s unlimited

*/


#include <algorithm>
#include <set>
#include <array>



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
void construct_bucket_set(int bucket_set[], int& bsize, const int q, const int ah){

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
    
    bsize = index;
}

int check_bucket_set_validity(int bucket_set[], const int bsize, const int q, const int ah){
    
}
