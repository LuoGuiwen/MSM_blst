
# MSM_blst

Credit: This is developed on top of a modified blst library. Check the original library at [blst library](https://github.com/supranational/blst).

## Compilation
In the terminal under 
<code>
MSM_blst 
</code>
folder,
one first runs 
<code>
./build.sh
</code>
to build the modified blst library, then runs

<pre><code>
g++ -std=c++17 -o main_test -g -O2 main_p1.cpp libblst.a
</code></pre> 
or 
<pre><code>
g++ -std=c++17 -o main_test -g -O2 main_p2.cpp libblst.a
</code></pre> 

to complie the corresponding benchmark over
<code>
G_1 
</code>
or
<code>
G_2
</code>
respectively.

Type in
<pre><code>
./main_test
</code></pre> 
to run the benchmark.

## Configuration
In 
<code>
main_p1.cpp
</code>
or 
<code>
main_p2.cpp
</code>
there is a 
<pre><code>
/***----***
Configuration
***----***/
</code></pre> 
snippet. One can adjust the integer 
<code>
xx (8<= xx <= 24)
</code>
in the line 
<pre><code> 
#include "ches_config_files/config_file_n_exp_11.h" 
</code></pre> 
to run the code for different number of points
<code> 
n = 2^{xx}
</code>.
One can decide whether to run the test for a specific algorithm by assigning bool values to
<code>
TEST_PIPPENGER_Q_OVER_5_CHES
</code>
and
<code>
TEST_PIPPENGER_BGMW95
</code>.
