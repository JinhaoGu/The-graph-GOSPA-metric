# GPU backend: time and memory by graph size

All numbers in this file come from the GPU backend (`gospa_gpu.py`). Each row
is a single call to `LP_graph_GOSPA` or `LP_graph_GOSPA_directed` with the
default settings: PDLP; `presolve="auto"`, which turns presolve off at these
sizes; and `build="auto"`, which assembles the LP with CuPy on the GPU for
LPs with at least 10 M stored entries and with SciPy on the CPU below that.

## Setup

| | |
|---|---|
| GPU | NVIDIA A100-SXM4-80GB (Google Colab), driver 580.82.07 |
| Software | cuOpt 26.10.00, CuPy 14.0.1, SciPy 1.16.3, NumPy 2.1.3, Python 3.13 |
| Parameters | c = 3, p = 1, epsilon = 1, 2-D node attributes |
| Graphs | X is random with edge density 0.5. Y is a noisy copy of X: 5% of the nodes are dropped, positions get Gaussian noise (sigma 0.3) and 5% of the edges are flipped. |

Columns:

- **build**: assembling the LP. The label says where: GPU (`gospa_gpu._build_lp_gpu`)
  or CPU (`graphGOSPA._build_lp`).
- **solve**: building the cuOpt model and solving it.
- **RAM** and **VRAM**: peak process RSS and peak GPU memory during the call.

## Undirected graphs

| n | variables | nonzeros | build | solve | total | RAM | VRAM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 28,395 | 0.9 M | 0.04 s (CPU) | 0.57 s | 0.6 s | 0.9 GB | 0.5 GB |
| 200 | 113,189 | 7.4 M | 0.29 s (CPU) | 1.16 s | 1.4 s | 1.0 GB | 0.9 GB |
| 300 | 258,888 | 25.5 M | 0.82 s (GPU) | 1.76 s | 2.6 s | 1.4 GB | 2.2 GB |
| 400 | 466,389 | 61.5 M | 1.41 s (GPU) | 3.52 s | 4.9 s | 2.3 GB | 4.6 GB |
| 500 | 711,975 | 116.2 M | 2.29 s (GPU) | 4.72 s | 7.0 s | 3.6 GB | 8.3 GB |
| 600 | 1,012,763 | 196.0 M | 3.63 s (GPU) | 10.98 s | 14.6 s | 5.5 GB | 13.3 GB |
| 700 | 1,399,967 | 320.2 M | 5.64 s (GPU) | 17.88 s | 23.6 s | 8.4 GB | 20.2 GB |
| 800 | 1,808,754 | 469.7 M | 8.07 s (GPU) | 20.04 s | 28.1 s | 11.9 GB | 29.3 GB |
| 900 | 2,331,864 | 686.1 M | 11.46 s (GPU) | 32.87 s | 44.4 s | 17.0 GB | 42.4 GB |
| 1000 | 2,842,948 | 925.6 M | 15.42 s (GPU) | 44.82 s | 60.3 s | 22.7 GB | 57.0 GB |
| 1100 | 3,420,937 | 1220.4 M | 19.99 s (GPU) | 66.46 s | 86.5 s | 29.6 GB | 74.9 GB |

## Directed graphs

A directed graph has two edge-mismatch blocks (one per direction). The LP
therefore has about twice the nonzeros of an undirected graph of the same size.

| n | variables | nonzeros | build | solve | total | RAM | VRAM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 47,195 | 1.9 M | 0.09 s (CPU) | 0.62 s | 0.7 s | 0.9 GB | 0.6 GB |
| 200 | 188,389 | 14.7 M | 0.64 s (GPU) | 1.09 s | 1.7 s | 1.2 GB | 1.4 GB |
| 300 | 431,088 | 50.8 M | 1.24 s (GPU) | 2.62 s | 3.9 s | 2.0 GB | 3.7 GB |
| 400 | 776,789 | 122.9 M | 2.38 s (GPU) | 6.11 s | 8.5 s | 3.7 GB | 8.1 GB |
| 500 | 1,185,975 | 232.5 M | 4.19 s (GPU) | 10.33 s | 14.5 s | 6.3 GB | 14.8 GB |
| 600 | 1,687,163 | 392.5 M | 6.92 s (GPU) | 18.15 s | 25.1 s | 10.1 GB | 24.6 GB |
| 700 | 2,332,367 | 638.5 M | 10.85 s (GPU) | 39.26 s | 50.2 s | 15.9 GB | 39.6 GB |
| 800 | 3,013,554 | 939.8 M | 15.67 s (GPU) | 51.38 s | 67.1 s | 23.0 GB | 57.9 GB |

## Complete graphs

Undirected, edge density 1. Y is X with jittered positions and the same edges.

| n | total | RAM | VRAM |
|---:|---:|---:|---:|
| 100 | 1.0 s | 0.9 GB | 0.6 GB |
| 200 | 3.7 s | 1.2 GB | 1.6 GB |
| 300 | 7.0 s | 2.1 GB | 4.1 GB |
| 500 | 43.8 s | 6.7 GB | 15.8 GB |
| 700 | 136.6 s | 17.0 GB | 42.2 GB |

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

For n >= 300 each nonzero needs about **62 bytes of VRAM** and **24 bytes of
RAM**. RAM peaks when cuOpt's model is created, because for a moment both our
copy of the matrix and cuOpt's copy exist (12 bytes per nonzero each). VRAM
peaks during the solve.

Fitted VRAM rules for n >= 300 (GB):

```text
undirected, density 0.5:  VRAM ~= 0.7 + 5.6e-8  * n^3
directed,   density 0.5:  VRAM ~= 0.8 + 1.12e-7 * n^3
complete (undirected):    VRAM ~= 0.6 + 1.21e-7 * n^3
```

The total time also grows roughly as n^3: about 65 ns * n^3 for undirected
graphs at density 0.5. With the GPU build, assembling the LP takes about a
quarter of the total; the rest is the PDLP solve.

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

Host RAM must also hold the LP for a moment. That is about 0.4x the VRAM
figure.

## LP assembly on the GPU

`build="gpu"` generates the constraint matrix on the GPU with CuPy, in chunks
of mismatch columns. It converts each chunk to CSR on the GPU and copies it
into preallocated host arrays. The result is the same matrix as the SciPy
build, entry for entry. It was checked on 560 random cases: directed and
undirected; int, bool, weighted, complete and sparse adjacency matrices; with
and without self-loops; with the default chunk size and with a tiny one
that forces many chunks. Unlike the SciPy build, there
are no intermediate copies on the host. The host copy is also released as
soon as cuOpt has its own copy. Undirected graphs, density 0.5:

| n | build CPU → GPU | total | RAM | VRAM |
|---:|---:|---:|---:|---:|
| 500 | 5.09 → 2.29 s | 10.1 → 7.0 s | 4.7 → 3.6 GB | 7.7 → 8.3 GB |
| 700 | 14.18 → 5.64 s | 32.3 → 23.6 s | 11.8 → 8.4 GB | 20.1 → 20.2 GB |
| 1000 | 40.21 → 15.42 s | 85.9 → 60.3 s | 33.1 → 22.7 GB | 57.0 → 57.0 GB |
| 1100 | 53.04 → 19.99 s | 118.9 → 86.6 s | 43.5 → 29.6 GB | 74.9 → 74.9 GB |

On mid-sized graphs the GPU build adds up to about 0.7 GB of VRAM. On large
graphs the solve sets the VRAM peak, so the peak does not change. For small
LPs, the first CuPy calls cost more than they save: at n = 100 the GPU build
takes 0.41 s and the CPU build 0.04 s. For this reason `build="auto"` keeps
those LPs on the CPU.

CuPy 14.0.1 note: `cupyx.scipy.sparse.csr_matrix.eliminate_zeros()` also
drops negative entries. The GPU build filters zeros on the COO arrays instead.

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
| 500 | 14.6 → 8.3 GB | 10.7 → 3.6 GB | 18.2 → 7.0 s |
| 700 | 39.4 → 20.2 GB | 28.3 → 8.4 GB | 50.6 → 23.6 s |

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
