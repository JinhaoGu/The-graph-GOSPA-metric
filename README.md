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
import numpy
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


