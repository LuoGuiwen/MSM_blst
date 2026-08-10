# Bucket value-to-index benchmark

Date: 2026-08-10  
Platform: Apple M1 Pro (6 performance + 2 efficiency cores), 16 GB RAM  
Compiler: Apple Clang 17.0.0, `-O2`  

## Implementations

- `direct` (default): `int[q_RADIX / 2 + 1]`
- `rank`: one 16-byte `bucket_rank_word` per 64 possible bucket values

The rank table stores a 64-bit presence word and the number of set bits before
that word. A lookup performs one table read, one mask, one population count,
and one addition. The original direct implementation remains unchanged and is
still the default.

## Reproduction

Default direct implementation:

```sh
make group=1 config=12 bucket_index=direct
./main_test_p1
```

Rank implementation with reproducible inputs and explicit repetition counts:

```sh
make group=1 config=12 bucket_index=rank \
  bench_test_num=3 bench_loop_num=5 bench_deterministic=1
./main_test_p1
```

The same options work with `group=2`. `bench_test_num`, `bench_loop_num`, and
`bench_deterministic` affect only the benchmark harness; their defaults retain
the original behavior.

## G1 results

Times are average microseconds per complete CHES MSM. Positive deltas mean rank
is slower. Each run also compared both CHES outputs with BGMW95 and built-in
blst; all coordinates matched.

| Config | q | Input | Direct table | Rank table | CHES direct | CHES rank | Delta | Integral direct | Integral rank | Delta |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| N_EXP=12 | 2^14 | fixed | 32 KiB | 2 KiB | 27,787 | 28,402 | +2.21% | 27,616 | 28,259 | +2.33% |
| N_EXP=14 | 2^16 | random | 128 KiB | 8 KiB | 98,269 | 101,405 | +3.19% | 98,728 | 98,063 | -0.67% |
| N_EXP=16 beta | 2^18 | random | 512 KiB | 32 KiB | 362,975 | 376,622 | +3.76% | 356,929 | 365,204 | +2.32% |
| N_EXP=16 | 2^19 | random | 1,024 KiB | 64 KiB | 374,814 | 380,284 | +1.46% | 371,361 | 377,066 | +1.54% |
| N_EXP=17 | 2^20 | random | 2,048 KiB | 128 KiB | 732,342 | 736,920 | +0.63% | 705,097 | 717,726 | +1.79% |
| N_EXP=20 | 2^22 | fixed, 1 sample | 8,192 KiB | 512 KiB | 5,881,164 | 6,222,089 | +5.80% | 5,175,948 | 5,054,298 | -2.35% |

Only the q selected by the repository configuration for each N_EXP is relevant
to an end-to-end MSM decision. An earlier artificial N_EXP=16/q=2^22 cache
experiment has therefore been excluded. The N_EXP=20 row is a single paired
sample: its two CHES paths moved in opposite directions, while unrelated
BGMW95/blst timings also moved by about 2%, so it is not sufficient evidence
of a speedup or slowdown.

## Memory

The rank representation is approximately 16 times smaller for every q. At the
N_EXP=20 optimum q=2^22, the mapping table itself shrinks from 8,388,612 bytes
to 524,304 bytes. Whole-process memory is dominated by fixed-point data, and
precomputation is outside the blockchain business-time objective.

## Decision

- Keep `direct` as the default.
- Treat `rank` only as a known low-memory candidate, not as a demonstrated MSM
  time optimization.
- Do not select q to favor the lookup structure. Benchmark only the q already
  chosen as optimal for the corresponding N_EXP.
