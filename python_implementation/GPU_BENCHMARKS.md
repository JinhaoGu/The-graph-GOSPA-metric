# GPU backend: time and memory by graph size

Environment: NVIDIA A100-SXM4-80GB, CUDA 13.3, cuOpt 26.10.0, float64,
`presolve=0`, 2-D node attributes, adjacency density 0.5, c=3, p=1, epsilon=1.
Times are for one call of `LP_graph_GOSPA`.

| n | variables | nonzeros | LP build | model | solve | RAM peak | VRAM peak |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 20,202 | 2.0 M | 0.08 s | 0.02 s | 0.74 s | 1.3 GB | 1.4 GB |
| 200 | 80,402 | 16.0 M | 0.58 s | 0.20 s | 4.42 s | 3.3 GB | 2.3 GB |
| 300 | 180,602 | 54.0 M | 2.08 s | 0.69 s | 2.36 s | 3.2 GB | 4.6 GB |
| 400 | 320,802 | 128.6 M | 4.92 s | 1.63 s | 5.24 s | 6.7 GB | 9.2 GB |
| 500 | 501,002 | 250.2 M | 9.57 s | 3.18 s | 8.43 s | 12.4 GB | 16.6 GB |

`LP build` is the Python/scipy model assembly, `model` is the conversion to
cuOpt's `DataModel`, and `solve` is `Solve()`.

## Where the memory goes

- **RAM peak** happens while the LP is assembled in Python, before cuOpt is
  called. `hstack`/`vstack` keep several copies of the constraint matrix
  alive at once.
- **VRAM peak** happens during the solve. It grows with the number of
  nonzeros, which is O(n^3) at fixed density.

## Rules of thumb

For density 0.5 and n >= 200:

```text
RAM  (GB) ~= 0.5 + 9.5e-8 * n^3
VRAM (GB) ~= 1.3 + 1.22e-7 * n^3
```

Equivalently, at n >= 300: about **50 bytes of RAM** and **70 bytes of VRAM**
per nonzero of the constraint matrix. Both scale roughly linearly with the
adjacency density.

## Largest n that fits in GPU memory

Solving the VRAM rule for n:

| GPU memory | approximate max n |
|---:|---:|
| 12 GB | 450 |
| 24 GB | 565 |
| 40 GB | 670 |
| 80 GB | 845 |

These are extrapolations from the measured points; at n = 1000 the LP would
need roughly 120 GB of VRAM, so it does not fit on a single 80 GB device.

## CPU vs GPU

On the same graphs the CPU implementation (`highs-ipm`) takes about 4.6 s at
n = 50, while the GPU path takes 0.13 s. Above n = 100 the GPU path is the
only practical option: the CPU side needs minutes at n = 300 and does not
scale to n = 500.
