# MSM_blst — Multi-Scalar Multiplication over Fixed Points on BLS12-381

Fast multi-scalar multiplication (MSM) over fixed points on the BLS12-381 elliptic
curve, using a novel **bucket set construction** that accelerates Pippenger's bucket
method via precomputation.

> Built on top of a modified [blst](https://github.com/supranational/blst) library.

---

## Directory Structure

```
MSM_blst/
├── README.md
├── LICENSE
├── Makefile                  # Build system
├── build.sh                  # Build libblst.a (modified blst)
├── run.sh                    # One-step build & run script
│
├── src/                      # Modified blst library source
│   ├── multi_scalar.c        #   Core MSM routines (Pippenger + CHES variants)
│   ├── server.c, exp.c, ...
│   └── asm/                  #   Assembly optimizations
│
├── bindings/                 # blst language bindings (C, Rust, Go, Python, …)
│   └── blst.h
│
├── build/                    # Build helper scripts & platform support
│
├── msm/                      # ★ MSM application layer (this project)
│   ├── include/
│   │   ├── msm_common.h      #   Master header: types, ostream ops, utilities
│   │   ├── bucket_set.h      #   Bucket set construction & validation
│   │   ├── scalar_conversion.h # Scalar ↔ q-ary / MB / q-half conversion
│   │   └── random_utils.h    #   SHA256-based random scalar generation
│   ├── bench_g1.cpp          #   Benchmark harness for G₁
│   ├── bench_g2.cpp          #   Benchmark harness for G₂
│   └── bucket_set_standalone.cpp  # Standalone bucket set construction utility
│
├── config/                   # Pre-computed configuration files (n = 2⁸ … 2²¹)
│   ├── config_file.h         #   Currently active config (copied by Makefile)
│   ├── n_exp_10.h
│   ├── n_exp_11.h
│   └── ...
│
├── external/                 # Third-party dependencies
│   └── uint256/              #   uint256 type (from Aztec Protocol)
│
└── archive/                  # Discarded / experimental code kept for reference
    └── discarded_function.h
```

---

## Dependencies

| Dependency  | Purpose                        | Install                                  |
|-------------|--------------------------------|------------------------------------------|
| **OpenSSL** | SHA256 hash for random scalars | `brew install openssl` (macOS)           |
| **g++**     | C++17 compiler                 | Xcode CLI Tools / `brew install gcc`     |
| **make**    | Build system                   | Pre-installed on macOS                   |

The modified **blst** library is bundled in `src/` and built locally — no external
installation needed.

---

## Quick Start

### 1. Install dependencies and configure paths

```bash
# macOS — install OpenSSL if not already present
brew install openssl

# Find your OpenSSL library path (version number may differ)
ls /opt/homebrew/Cellar/openssl@3/
# Example output: 3.6.3
```

Then edit `Makefile` to match your system. Use `brew --prefix openssl` to check:

```bash
brew --prefix openssl
# → /opt/homebrew/opt/openssl@3
```

Update the two paths in `Makefile`:
```makefile
openssl_include_directory = /opt/homebrew/include/
openssl_library_directory = /opt/homebrew/Cellar/openssl@3/<your-version>/lib
```

### 2. One-command: build & run

`run.sh` handles everything — building the blst library, copying the config,
compiling, and executing:

```bash
# G₁ benchmark, n = 2¹⁰ (default)
./run.sh group=1

# G₂ benchmark, n = 2¹⁰ (default)
./run.sh group=2

# Multiple n values in one go
./run.sh group=1 config=10,11,12     # n = 2¹⁰, 2¹¹, 2¹² sequentially
./run.sh group=2 config=15,17        # n = 2¹⁵, 2¹⁷ sequentially
```

**`config` parameter:** `8`–`21` (maps to `config/n_exp_*.h`), plus `_beta` variants
like `16_beta`, `17_beta`, `20_beta`. Default is `10`.

### 3. Manual step-by-step (useful for debugging)

```bash
# Step 1 — build the modified blst library
make libblst.a

# Step 2 — select config and compile
make group=1 config=10
# Produces: ./main_test_p1

# Step 3 — prevent stack overflow for large n, then run
ulimit -s unlimited
./main_test_p1

# Clean up build artifacts
make clean
```

> **Important for large n (≥ 16):** Precomputation tables can consume several GB
> of memory. Always run `ulimit -s unlimited` first, and ensure your machine has
> sufficient RAM (> 8 GB for n = 2¹⁶, > 32 GB for n = 2²¹).

> **Note:** `MSM_blst` is **not** compatible with the upstream blst library —
> some source files in `src/` have been modified to support the CHES bucket method.

---

## Configuration

Configuration files in `config/` define parameters for different numbers of points
(`n = 2^N_EXP`). Each file sets:

| Parameter              | Description                                   |
|------------------------|-----------------------------------------------|
| `N_EXP`                | Exponent: n = 2^N_EXP                         |
| `EXPONENT_OF_q`        | log₂(radix) for CHES method                   |
| `h_LEN_SCALAR`         | Number of radix-q digits                      |
| `B_SIZE`               | Pre-computed bucket set size                  |
| `d_MAX_DIFF`           | Maximum gap in bucket set                     |
| `EXPONENT_OF_q_BGMW95` | log₂(radix) for BGMW95 variant                |
| `q_RADIX_PIPPENGER_VARIANT` | Radix for Pippenger variant              |
| `h_BGMW95`             | Number of digits for BGMW95                   |

Available configs: `n_exp_8` through `n_exp_21` (n = 2⁸ … 2²¹), plus `_beta` variants.

---

## Algorithms Tested

Each benchmark compares **4 methods** for computing n-scalar multiplications:

| # | Method                                   | Description                              |
|---|------------------------------------------|------------------------------------------|
| 1 | **CHES** (`nh + q/5`)                    | Our construction (paper)                 |
| 2 | **CHES** integral scalar conversion      | Variant with different conversion path   |
| 3 | **BGMW95** (`nh + q/2`)                  | Pippenger variant from BGMW95            |
| 4 | **blst built-in** Pippenger              | Original blst implementation             |

Enable/disable methods by editing the flags at the top of `msm/bench_g1.cpp`
or `msm/bench_g2.cpp`:

```cpp
bool TEST_PIPPENGER_Q_OVER_5_CHES = 1;   // Enable CHES methods (1 & 2)
bool TEST_PIPPENGER_BGMW95 = 1;          // Enable BGMW95 method (3)
```

---

## Expected Output

The benchmark prints, for each method:
- Average wall-clock time (µs)
- Coordinates of the computed point (x, y)
- Percentage improvement over baselines (blst built-in and BGMW95)
- Scalar conversion overhead

Example:
```
1. CHES 'nh+ q/5'. Wall clock time elapse is: 12345 us
2. CHES 'nh+ q/5' integral scalar conversion. Wall clock time elapse is: 12000 us
3. pippenger_variant_BGMW95. Wall clock time elapse is: 15000 us
4. pippenger_blst_built_in. Wall clock time elapse is: 18000 us
Improvement, blst CHES_q_over_5 vs pipp: 33.333%
...
```

---

## Standalone Utilities

### Bucket Set Construction

```bash
g++ -std=c++17 -o bucket_set_test -g -O2 msm/bucket_set_standalone.cpp
./bucket_set_test
```

Prints the bucket set size, maximum gap, density ratio, and validity check.

---

## Troubleshooting

**Segmentation fault?** Increase the stack size before running:
```bash
ulimit -s unlimited
./main_test_p1
```
This is often needed for n ≥ 2¹² due to large stack-allocated arrays.

**OpenSSL not found at link time?** Check that the path in `Makefile` matches your
installed version:
```bash
# Show installed OpenSSL version and path
brew --prefix openssl
ls /opt/homebrew/Cellar/openssl@3/

# Then update Makefile:
#   openssl_library_directory = /opt/homebrew/Cellar/openssl@3/<version>/lib
```

**Build fails with `'stdlib.h' file not found`?** Your Xcode SDK path may have
changed. Verify with:
```bash
xcrun --show-sdk-path
```
The `Makefile` auto-detects this — if it still fails, try:
```bash
SDKROOT=$(xcrun --show-sdk-path) make group=1 config=10
```

**Clean rebuild:**
```bash
make clean
rm -f config/config_file.h   # also reset the active config
```

---

## Credit

This project is developed on top of a modified [blst](https://github.com/supranational/blst)
library. The MSM algorithm and bucket set construction are from the CHES paper by
Guiwen Luo.

---

## License

See [LICENSE](./LICENSE) for details.
