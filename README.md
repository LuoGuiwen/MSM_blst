
# MSM_blst

Credit: This is developed on top of a modified blst library. Check the original library at [blst library](https://github.com/supranational/blst).


## Overview
The arithmetic of computing multiple scalar multiplications in an elliptic
curve group then adding them together is called multi-scalar multiplication (MSM).
MSM over fixed points dominates the time consumption in the pairing-based trusted
setup zero-knowledge succinct non-interactive argument of knowledge (zkSNARK),
thus for practical applications we would appreciate fast algorithms to compute it.
This paper proposes a bucket set construction that can be utilized in the context
of Pippenger’s bucket method to speed up MSM over fixed points with the help of
precomputation. If instantiating the proposed construction over BLS12-381 curve,
when computing n-scalar multiplications for n = 2e
(10 ≤ e ≤ 21), theoretical analysis
indicates that the proposed construction saves more than 21% computational cost
compared to Pippenger’s bucket method, and that it saves 2.6% to 9.6% computational
cost compared to the most popular variant of Pippenger’s bucket method. Finally, our
experimental result demonstrates the feasibility of accelerating the computation of
MSM over fixed points using large precomputation tables as well as the effectiveness
of our new construction.


## Compilation
The code is tested on intel Mac OS, on M1 Mac OS, and on intel unbuntu. They can be compiled and run using the same commands as explained below. Note MSM_blst is not compatible with the original blst library since some of the source code in blst has been modified. 

In the terminal under 
<code>
MSM_blst 
</code>
folder,
one first runs 
<code>
run.sh
</code>
to build the modified blst library and complie and run the corresponding benchmark over
<code>
G_1 
</code>
or
<code>
G_2
</code>
respectively.


### bash script parameters:
`benchmark`: to complie the corresponding benchmark over `G1` or `G2` respectively.

`config`: to select a list of config files to be included in the algorithm. 
- Defualt for `G1` is only the `10` config file.
- Defualt for `G2` is only the `16` config file

Example:

compile and run 3 times, each for config 8,9,10 with all their benchmark to be 1
<pre><code>
./run.sh benchmark=1 config=8,9,10 
</code></pre>
compile and run 2 times, each for config 15,16 with all their benchmark to be 2
<pre><code>
./run.sh benchmark=2 config=15,16 
</code></pre>
default config=10
<pre><code>
./run.sh benchmark=1 
</code></pre>
default config=16
<pre><code>
./run.sh benchmark=2      
</code></pre>


## Expected Output
<!-- difference between the two .cpp files -->
```main_p1.cpp``` and ```main_p1.cpp``` each contains 4 methods to do the computation. Which are
1. Our Construction (CHES 'nh+ q/5')
2. Our Construction with integral scalar conversion
3. Pippenger Variant (pippenger_variant_BGMW95)
4. Pippenger (pippenger_blst_built_in)

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
#include "ches_config_files/config_file_n_exp_xx.h" 
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
