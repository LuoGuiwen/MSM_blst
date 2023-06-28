
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
The code is tested on M1 Mac OS. 
<!---
It relies on SHA256 hash function in openssl library to generate the random scalars. Before running the following code, one should first install openssl ( For example, by running <code>
brew install openssl
</code>).
-->
In the terminal under 
<code>
MSM_blst 
</code>
folder,
one directly runs 
<code>
./run.sh group=1
</code>
or
<code>
./run.sh group=2
</code>
to build the modified blst library, complie and run the corresponding test in
<code>
\mathbb{G}_1
</code>
or
<code>
\mathbb{G}_2
</code>
respectively. The number of points $n$ is $2^{10}$ by defualt. Check the next part for setting configuration.

Note MSM_blst is not compatible with the original blst library since some of the source code in blst has been modified. 

## Bash script parameters:
`group`: Complie the corresponding group over in $\mathbb{G}_1$ or $\mathbb{G}_2$ respectively.

`config`: Select a list of config files to be included in the algorithm. 


Examples:

The following command
```
./run.sh group=1 config=10,11,12 
```
execute the test for $n$-scalar multiplications in $\mathbb{G}_1$, and it will run through $n = 2^{10}, 2^{11}, 2^{12}$ sequentially.

The following command
```
./run.sh group=2 config=15,17
```
execute the test for $n$-scalar multiplications in $\mathbb{G}_2$, and it will run through $n = 2^{15}, 2^{17}$ sequentially.

We have `config=10` by defualt.


## Expected output 
<!-- difference between the two .cpp files -->
```main_p1.cpp``` and ```main_p2.cpp```  have 4 methods to compute $n$-scalar multiplications. They are
1. Our Construction (CHES 'nh+ q/5'),
2. Our Construction with integral scalar conversion,
3. Pippenger Variant (pippenger_variant_BGMW95),
4. Pippenger (pippenger_blst_built_in).

In 
<code>
main_p1.cpp
</code>
and
<code>
main_p2.cpp
</code>
there is a 
```
/***----***
Configuration
***----***/
```
snippet. One can decide whether to run the test for a specific algorithm by assigning bool values to
`TEST_PIPPENGER_Q_OVER_5_CHES`
and
`TEST_PIPPENGER_BGMW95`.
