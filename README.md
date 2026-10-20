# The graph GOSPA metric

This repository contains Matlab and Python code with the linear programming implementation of the graph generalised optimal subpattern assignment (GOSPA) metric proposed in [1]. The graph GOSPA metric is a mathematically principled metric for graphs. Given two graphs, it penalises node attribute error for properly assigned nodes, the number of missed nodes, false nodes, and edge mismatches. The repository contains the implementations for undirected graphs and for directed graphs.

The graph GOSPA metric is an extension of the GOSPA metric for sets of objects proposed in [2], which was also extended for sets of trajectories in [3].

[1] J. Gu, A. F. García-Fernández, Robert E. Firth, L. Svensson, “Graph GOSPA metric: a metric to measure the discrepancy between graphs of different sizes” in IEEE Transactions on Signal Processing, 2024 (https://arxiv.org/abs/2311.07596)

[2] A. S. Rahmathullah, Á. F. García-Fernández and L. Svensson, "Generalized optimal sub-pattern assignment metric," 2017 20th International Conference on Information Fusion (Fusion), Xi'an, China, 2017, pp. 1-8, doi: 10.23919/ICIF.2017.8009645.

[3] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A Metric on the Space of Finite Sets of Trajectories for Evaluation of Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020, doi: 10.1109/TSP.2020.3005309.”

## Usage:
Below are usage examples of the graph GOSPA metric.
### Python:
```python
import numpy as np
from graphGOSPA import LP_graph_GOSPA,LP_graph_GOSPA_directed

# define graph X and Y 
X_attr= np.array([[0,0],[10,10],[10,20]])
X_adj=np.array([[0,1,1],[1,0,1],[1,1,0]])

Y_attr= np.array([[0,0],[10,10],[10,20]])
Y_adj=np.array([[0,1,1],[1,0,1],[1,1,0]])

# choose parameters
c=3 # penalty for missing or false nodes 
p=1 # p-norm
epsilon=1 # penalty for edge mismatch

# for undirected graphs
dxy,loc_cost,false_cost,miss_cost,edge_cost=LP_graph_GOSPA(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon)


# for directed graphs
dxy,loc_cost,false_cost,miss_cost,edge_cost=LP_graph_GOSPA_directed(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon)

```

`python_implementation/graphGOSPA_original.py` keeps the original CPU implementation, with the earlier LP formulation, for reference and comparison. It has the same functions and returns the same values. The full repository state before that change is tagged `original-cpu-version`.

### MATLAB:

```MATLAB
% define graphs X and Y
X = struct('xState', [], 'adj', []);
Y = struct('xState', [], 'adj', []);
nx = 3;
X.xState = [0,0;10,10;10,20];
X.adj = [0,1,1;1,0,1;1,1,0];
ny = 3;
Y.xState = [0,0;10,10;10,20]; 
Y.adj = [0,1,1;1,0,1;1,1,0];

% choose the parameters
c=3; %penalty for missing or false nodes 
p=1; %p-norm
epsilon=1; %penalty for edge mismatch



% graph GOSPA for undirected graphs 
[dxy,loc_cost,fa_cost,miss_cost,edge_cost]=LPgraphGOSPA(X,Y,c,p,epsilon);

% graph GOSPA for directed graphs
[dxy,loc_cost,fa_cost,miss_cost,edge_cost]=LPgraphGOSPA_directed(X,Y,c,p,epsilon);

```


## GPU backend (optional)

`python_implementation/gospa_gpu.py` solves the same LP with the
[NVIDIA cuOpt](https://github.com/NVIDIA/cuopt) GPU solver. The functions have
the same names, arguments and return values as in `graphGOSPA.py`. The LP is
built by the same code, so only the solver changes.

Use the GPU backend when the graphs have more than about 50 nodes. At that size
the CPU solver already takes seconds, and from a few hundred nodes it needs
minutes or does not finish (see [GPU_BENCHMARKS.md](python_implementation/GPU_BENCHMARKS.md)).

### Requirements and installation

- Linux, and an NVIDIA GPU with a recent driver. Google Colab GPU runtimes work.
- numpy and scipy, plus cuOpt (this also installs CuPy, which the GPU LP assembly uses):

```bash
pip install --extra-index-url https://pypi.nvidia.com cuopt-cu13
```

Check the installation:

```python
import gospa_gpu
print(gospa_gpu.gpu_status())   # {'available': True, ...}
```

### Usage

Put `gospa_gpu.py` next to `graphGOSPA.py` and call it like the CPU version:

```python
import gospa_gpu as M

# undirected graphs
dxy,loc_cost,false_cost,miss_cost,edge_cost=\
    M.LP_graph_GOSPA(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon)

# directed graphs
dxy,loc_cost,false_cost,miss_cost,edge_cost=\
    M.LP_graph_GOSPA_directed(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon)
```

Optional keyword arguments:

| argument | default | meaning |
|---|---|---|
| `method` | `"pdlp"` | cuOpt algorithm: `"pdlp"`, `"concurrent"`, `"barrier"` or `"dual_simplex"`. Keep `"pdlp"`: the other solvers run on the CPU or run out of memory on these LPs, and they never finish first. |
| `presolve` | `"auto"` | `"auto"` enables presolve up to 10,000 LP variables (about 57 nodes per graph for undirected graphs, 44 for directed ones) and disables it above, where it is slow. Pass `True` or `False` to force it. |
| `tolerance` | `1e-8` | Relative optimality tolerance. Results match the CPU solver to about 1e-8 relative. |
| `time_limit` | `1800` | Time limit in seconds. |
| `build` | `"auto"` | Where the LP is assembled: `"gpu"` (CuPy), `"cpu"` (SciPy) or `"auto"`. `"auto"` uses the GPU from 10 M stored entries on (about 215 nodes per graph at density 0.5) and the CPU below that. Both give the same matrix. The GPU build is about 2.5x faster and needs less host RAM. |

If cuOpt stops before reaching an optimal solution (time limit or numerical
trouble), a `RuntimeError` is raised. A non-optimal solution is never
returned, because it would overestimate the metric.

### Memory and graph size

GPU memory grows with the number of nonzeros of the constraint matrix. For
dense graphs that is O(n^3), where n is the number of nodes.

Largest graph that fits, from measurements on an A100 80 GB:

| GPU memory | undirected, density 0.5 | directed, density 0.5 | complete graphs |
|---:|---:|---:|---:|
| 16 GB | ~650 nodes | ~510 | ~500 |
| 24 GB | ~745 | ~590 | ~580 |
| 40 GB | ~890 | ~705 | ~685 |
| 80 GB | ~1120 | ~890 | ~870 |

For example, an undirected graph pair with n = 1000 at density 0.5 takes 60 s
and 57 GB of GPU memory. The machine also needs host RAM of about 0.4x the GPU
memory, because cuOpt keeps a copy of the LP on the host.

Full measurements, time scaling and a comparison with the previous LP
formulation are in [GPU_BENCHMARKS.md](python_implementation/GPU_BENCHMARKS.md).
