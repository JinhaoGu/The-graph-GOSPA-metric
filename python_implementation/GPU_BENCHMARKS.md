# GPU backend: time and memory by graph size

All numbers in this file come from the GPU backend (`gospa_gpu.py`). Each row
is a single call to `LP_graph_GOSPA` or `LP_graph_GOSPA_directed` with the
default settings: PDLP, and `presolve="auto"`, which turns presolve off at
these sizes.

## Setup

| | |
|---|---|
| GPU | NVIDIA A100-SXM4-80GB (Google Colab), driver 580.82.07 |
| Software | cuOpt 26.10.00, SciPy 1.16.3, NumPy 2.1.3, Python 3.13 |
| Parameters | c = 3, p = 1, epsilon = 1, 2-D node attributes |
| Graphs | X is random with edge density 0.5. Y is a noisy copy of X: 5% of the nodes are dropped, positions get Gaussian noise (sigma 0.3) and 5% of the edges are flipped. |

Columns:

- **build**: assembling the LP with SciPy on the CPU (`graphGOSPA._build_lp`).
- **solve**: the call to `gospa_gpu.solve_lp`. This includes building the cuOpt model.
- **RAM** and **VRAM**: peak process RSS and peak GPU memory during the call.

## Undirected graphs

| n | variables | nonzeros | build | solve | total | RAM | VRAM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 28,395 | 0.9 M | 0.04 s | 0.54 s | 0.6 s | 0.9 GB | 0.5 GB |
| 200 | 113,189 | 7.4 M | 0.29 s | 1.14 s | 1.4 s | 1.1 GB | 0.9 GB |
| 300 | 258,888 | 25.5 M | 1.13 s | 2.00 s | 3.1 s | 1.5 GB | 2.1 GB |
| 400 | 466,389 | 61.5 M | 2.70 s | 3.78 s | 6.5 s | 2.7 GB | 4.3 GB |
| 500 | 711,975 | 116.2 M | 5.09 s | 4.96 s | 10.1 s | 4.7 GB | 7.7 GB |
| 600 | 1,012,763 | 196.0 M | 8.53 s | 11.13 s | 19.7 s | 7.5 GB | 12.6 GB |
| 700 | 1,399,967 | 320.2 M | 14.18 s | 18.08 s | 32.3 s | 11.8 GB | 20.1 GB |
| 800 | 1,808,754 | 469.7 M | 20.40 s | 20.18 s | 40.6 s | 17.1 GB | 29.3 GB |
| 900 | 2,331,864 | 686.1 M | 29.50 s | 32.99 s | 62.5 s | 24.7 GB | 42.4 GB |
| 1000 | 2,842,948 | 925.6 M | 40.21 s | 45.65 s | 85.9 s | 33.1 GB | 57.0 GB |
| 1100 | 3,420,937 | 1220.4 M | 53.04 s | 65.78 s | 118.9 s | 43.5 GB | 74.9 GB |

## Directed graphs

A directed graph has two edge-mismatch blocks (one per direction). The LP
therefore has about twice the nonzeros of an undirected graph of the same size.

| n | variables | nonzeros | build | solve | total | RAM | VRAM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 47,195 | 1.9 M | 0.09 s | 0.57 s | 0.7 s | 0.9 GB | 0.6 GB |
| 200 | 188,389 | 14.7 M | 0.66 s | 1.34 s | 2.0 s | 1.2 GB | 1.4 GB |
| 300 | 431,088 | 50.8 M | 2.44 s | 2.85 s | 5.3 s | 2.4 GB | 3.7 GB |
| 400 | 776,789 | 122.9 M | 5.78 s | 6.33 s | 12.1 s | 4.9 GB | 8.1 GB |
| 500 | 1,185,975 | 232.5 M | 10.99 s | 10.56 s | 21.6 s | 8.8 GB | 14.8 GB |
| 600 | 1,687,163 | 392.5 M | 18.54 s | 18.34 s | 36.9 s | 14.4 GB | 24.6 GB |
| 700 | 2,332,367 | 638.5 M | 30.31 s | 39.40 s | 69.7 s | 23.1 GB | 39.6 GB |
| 800 | 3,013,554 | 939.8 M | 44.79 s | 51.79 s | 96.6 s | 33.7 GB | 57.9 GB |

## Complete graphs

Undirected, edge density 1. Y is X with jittered positions and the same edges.

| n | total | RAM | VRAM |
|---:|---:|---:|---:|
| 100 | 1.0 s | 0.9 GB | 0.6 GB |
| 200 | 4.0 s | 1.3 GB | 1.5 GB |
| 300 | 8.3 s | 2.5 GB | 3.9 GB |
| 500 | 50.7 s | 9.4 GB | 15.8 GB |
| 700 | 156.4 s | 24.7 GB | 42.2 GB |

Before the `kron` fix, the n = 100 case did not finish in 600 s. This was
true for the `main` branch and for the L1 formulation without the fix. When
the adjacency density is 0.5 or more, `scipy.sparse.kron` returned BSR
blocks that stored the zeros too. With the L1 formulation, the LP matrix then had about 50x more stored
entries than real nonzeros (101 M instead of 2.0 M at n = 100), and a
directed n = 500 graph crashed SciPy with a segmentation fault. Every `kron`
call now asks for CSR output.

## Where the memory and time go

The number of nonzeros of the constraint matrix determines both time and
memory. It grows as n^3, and linearly with the edge density rho:

```text
nonzeros ~= 1.9 * rho * n^3     undirected
nonzeros ~= 3.7 * rho * n^3     directed
```

For n >= 300 each nonzero needs about **62 bytes of VRAM** and **36 bytes of
RAM**. RAM peaks while SciPy assembles the LP, before cuOpt is called. VRAM
peaks during the solve.

Fitted VRAM rules for n >= 300 (GB):

```text
undirected, density 0.5:  VRAM ~= 0.7 + 5.6e-8  * n^3
directed,   density 0.5:  VRAM ~= 0.8 + 1.12e-7 * n^3
complete (undirected):    VRAM ~= 0.6 + 1.21e-7 * n^3
```

The total time also grows roughly as n^3: about 90 ns * n^3 for undirected
graphs at density 0.5. The CPU build takes about as long as the GPU solve.
During the build the GPU is idle, so moving the assembly to the GPU is the
next obvious speed-up.

## Largest n that fits in GPU memory

These values come from the VRAM rules above. The measured range is up to
n = 1100 (undirected), 800 (directed) and 700 (complete). Values beyond those
are extrapolations.

| GPU memory | undirected, density 0.5 | directed, density 0.5 | complete |
|---:|---:|---:|---:|
| 16 GB | 650 | 510 | 500 |
| 24 GB | 745 | 590 | 580 |
| 40 GB | 890 | 705 | 685 |
| 48 GB | 945 | 750 | 730 |
| 80 GB | 1120 | 890 | 870 |

Host RAM must also hold the LP while it is built. That is about 0.6x the VRAM
figure.

## Solver choice

`method="pdlp"` (the default) runs only cuOpt's GPU PDLP solver. cuOpt's
`"concurrent"` mode also starts dual simplex and barrier on the CPU, and it
returns only after they stop. On these LPs PDLP always finishes first. Dual
simplex alone did not finish in 900 s at n = 300, and barrier asked for
66.7 GB of GPU memory at n = 300 and failed. At n = 300 the solve took
6.4 s in concurrent mode and 1.9 s with PDLP only. From n = 700 upwards the
two modes take the same time.

## Comparison with the previous LP formulation

The previous formulation stored the edge-mismatch matrix twice
(`h >= G w`, `h >= -G w`) and added a dense `e >= sum(h)` row. The
previous code also used cuOpt's concurrent mode. It was measured on the same
GPU with the same graphs:

| n | VRAM before → now | RAM before → now | total before → now |
|---:|---:|---:|---:|
| 500 | 14.6 → 7.7 GB | 10.7 → 4.7 GB | 18.2 → 10.1 s |
| 700 | 39.4 → 20.1 GB | 28.3 → 11.8 GB | 50.6 → 32.3 s |

The largest undirected graph that fits in 80 GB went from about 845 to about
1120 nodes.

The original measurements of the previous formulation are kept below for
reference. They used 2-D node attributes, adjacency density 0.5, c=3, p=1,
epsilon=1, `presolve=0`, and cuOpt's default concurrent mode. The `model` column is
the conversion to cuOpt's `DataModel`.

| n | variables | nonzeros | LP build | model | solve | RAM peak | VRAM peak |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 20,202 | 2.0 M | 0.08 s | 0.02 s | 0.74 s | 1.3 GB | 1.4 GB |
| 200 | 80,402 | 16.0 M | 0.58 s | 0.20 s | 4.42 s | 3.3 GB | 2.3 GB |
| 300 | 180,602 | 54.0 M | 2.08 s | 0.69 s | 2.36 s | 3.2 GB | 4.6 GB |
| 400 | 320,802 | 128.6 M | 4.92 s | 1.63 s | 5.24 s | 6.7 GB | 9.2 GB |
| 500 | 501,002 | 250.2 M | 9.57 s | 3.18 s | 8.43 s | 12.4 GB | 16.6 GB |
| 600 | 721,202 | 430.7 M | 16.37 s | 5.42 s | 15.85 s | 20.9 GB | 27.6 GB |
| 700 | 981,402 | 685.4 M | 25.86 s | 8.64 s | 23.20 s | 32.8 GB | 43.0 GB |
| 800 | 1,281,602 | 1025.1 M | 38.88 s | 12.88 s | 34.67 s | 48.8 GB | 63.6 GB |


## CPU vs GPU

From the original measurements: on the same graphs the CPU implementation
(`highs-ipm`) takes about 4.6 s at n = 50, while the GPU path takes 0.13 s.
Above n = 100 the GPU path is the only practical option. The CPU solver needs
minutes at n = 300 and does not scale to n = 500.
