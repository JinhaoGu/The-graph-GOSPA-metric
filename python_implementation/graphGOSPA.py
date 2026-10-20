#Author: Jinhao Gu
#This code is a python implementation of the graph GOSPA metric proposed in the paper
# "Graph GOSPA metric: a metric to measure the discrepancy between graphs of different sizes"
# by Jinhao Gu, Á. F. García-Fernández, Robert E. Firth, Lennart Svensson
#
#Overview of the linear program (LP)
#-----------------------------------
#Graph X has n_x nodes and graph Y has n_y nodes. The node assignment is an
#(n_x+1)x(n_y+1) matrix W:
#  W(i,j)    for i<n_x, j<n_y : weight of assigning node i of X to node j of Y
#  W(i,n_y)  (last column)    : node i of X is left unassigned (missed node)
#  W(n_x,j)  (last row)       : node j of Y is left unassigned (false node)
#  W(n_x,n_y)                 : dummy-to-dummy entry. It has no cost and is
#                               not constrained, so it is always free.
#Every real row and every real column of W sums to one. W is relaxed to
#continuous values in [0,1], which is the LP relaxation used in the paper above.
#
#Objective = node cost + edge cost:
#  node cost : sum_ij D(i,j) W(i,j), where D is the cost matrix built by locCostComp
#  edge cost : eps^p/2 * ||A_X W_r - W_r A_Y||_1 for undirected graphs, with
#              W_r = W(:n_x,:n_y). Directed graphs add a second block in the
#              reverse direction and weight each block by eps^p/4.
#
#The LP vector of unknowns is
#  x = [ vec(W) , s_plus , s_minus ]
#Here vec() stacks the columns of W (column-major, MATLAB order), so W(i,j)
#is at index j*(n_x+1)+i. The slacks turn the L1 norm into linear
#constraints (see _build_lp).
#
#The kron identity used below to vectorise matrix products is
#  vec(A W B) = (B^T kron A) vec(W)
#Every sps.kron call passes format='csr'. Without it, scipy returns BSR when
#the second factor is dense enough (2*nnz >= size, e.g. adjacency density
#>= 0.5). BSR stores whole n x (n+1) blocks, zeros included, which made the
#LP matrix about 30x larger for complete graphs and crashed scipy at n ~ 500.
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog


def computeLocCostPerTime(x,y,c,p):
    #Scalar version of the cost between one node attribute x and one node
    #attribute y. Not used by the LP code, which calls the vectorised
    #locCostComp; kept for backward compatibility.
    if np.all(~np.isnan(x)) & np.all(~np.isnan(y)):
        #neither x nor y has nan
        return np.linalg.norm(x-y)**p # Use Euclidean distance as the cost function

    elif np.any(np.isnan(x) & ~np.isnan(y)) | np.any(~np.isnan(x) & np.isnan(y)):
        #exactly one of x or y has nan
        return c**p/2

    else:
        #both x and y have nan
        return 0


def locCostComp(X_attr,Y_attr,c,p):
    '''
    Build the (n_x+1)x(n_y+1) matrix D of node costs, which multiplies W in the objective.

    D(i,j), i<n_x, j<n_y : ||x_i - y_j||^p when neither attribute has a NaN;
                           c^p/2 when exactly one of them has a NaN in some
                           coordinate; 0 when both are NaN in the same coordinates.
    D(i,n_y), D(n_x,j)   : c^p/2, the cost of a missed or a false node.
    D(n_x,n_y)           : 0, the dummy-to-dummy entry.
    '''
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]
    tmpCost=c**p/2
    #Every entry starts as the missed/false node cost; the real block is overwritten below
    locCostMat=np.full((n_x+1,n_y+1),tmpCost)
    locCostMat[-1,-1]=0

    if n_x>0 and n_y>0:
        #Broadcast to n_x x n_y x D so that all pairs are handled in one go
        x_nan=np.isnan(X_attr)[:,None,:]
        y_nan=np.isnan(Y_attr)[None,:,:]
        #no_hole: neither node has a NaN anywhere
        no_hole=~(x_nan.any(axis=2) | y_nan.any(axis=2))
        #one_hole: some coordinate is NaN in exactly one of the two nodes
        one_hole=(x_nan & ~y_nan).any(axis=2) | (~x_nan & y_nan).any(axis=2)
        #Pairwise Euclidean distances (NaN where a NaN is present; those are replaced by np.where)
        dist=np.linalg.norm(np.asarray(X_attr,dtype=float)[:,None,:]-np.asarray(Y_attr,dtype=float)[None,:,:],axis=2)
        locCostMat[:n_x,:n_y]=np.where(no_hole,dist**p,np.where(one_hole,tmpCost,0.0))

    return locCostMat


def _is_numeric_dtype(dtype):
    #Accept integer, floating point and boolean arrays
    return (np.issubdtype(dtype,np.integer) or
            np.issubdtype(dtype,np.floating) or
            np.issubdtype(dtype,np.bool_))


def _as_attribute_array(attr,name):
    #Convert node attributes to a 2-D numeric array (nodes x attribute dimension).
    #An empty input such as [] becomes a 0x0 array, which means an empty graph.
    attr=np.asarray(attr)
    if attr.size==0 and attr.ndim<2:
        attr=attr.reshape(0,0)
    if attr.ndim!=2:
        raise ValueError(f'{name} must be a 2-D array of node attributes')
    if not _is_numeric_dtype(attr.dtype):
        raise ValueError(f'{name} must contain numeric node attributes')
    return attr


def _as_adjacency_array(adj,size,name):
    #Convert an adjacency matrix to a numeric array and check that it is size x size,
    #where size is the number of nodes given by the attribute array
    adj=np.asarray(adj)
    if adj.size==0 and adj.ndim<2 and size==0:
        adj=adj.reshape(0,0)
    if adj.shape!=(size,size):
        raise ValueError(f'{name} must have shape ({size}, {size}), got {adj.shape}')
    if not _is_numeric_dtype(adj.dtype):
        raise ValueError(f'{name} must contain numeric entries')
    return adj


def _check_inputs(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon):
    #Validate and normalise all inputs; returns the arrays in canonical form.
    #c>0, epsilon>0 and p>=1 are needed for the result to be a metric.
    X_attr=_as_attribute_array(X_attr,'X_attr')
    Y_attr=_as_attribute_array(Y_attr,'Y_attr')
    X_adj=_as_adjacency_array(X_adj,X_attr.shape[0],'X_adj')
    Y_adj=_as_adjacency_array(Y_adj,Y_attr.shape[0],'Y_adj')
    if not c>0:
        raise ValueError('c must be positive')
    if not epsilon>0:
        raise ValueError('epsilon must be positive')
    if np.ndim(p)!=0:
        raise ValueError('p must be a scalar')
    if not (np.isfinite(p) and p>=1):
        raise ValueError('p must be finite and greater than or equal to 1')
    return X_attr,Y_attr,X_adj,Y_adj


def _empty_result(n_x,n_y,c,p):
    #Closed form when one graph is empty: every node of X is missed and every
    #node of Y is false, each costing c^p/2. There are no assigned pairs and no
    #edge mismatch is counted, so no LP is needed.
    miss_cost_p=n_x*c**p/2
    false_cost_p=n_y*c**p/2
    dxy=(miss_cost_p+false_cost_p)**(1/p)
    return dxy,0.0,false_cost_p,miss_cost_p,0.0


def _assignment_constraints(n_x,n_y,nParam):
    '''
    Equality constraints that make W a valid (relaxed) assignment:
      for every real column j<n_y : sum_{i=0..n_x} W(i,j) = 1
      for every real row    i<n_x : sum_{j=0..n_y} W(i,j) = 1
    The dummy row and dummy column are not constrained themselves.
    The matrix has nParam columns: zero coefficients are added for the
    variables after vec(W), i.e. the slacks.
    '''
    WLen=(n_x+1)*(n_y+1)
    #Constraint 1: sum_i W(i,j)=1 for every real column j
    #vec(W) stores column j as one contiguous block of n_x+1 entries, so the
    #rows are kron(I_{n_y}, ones(1,n_x+1)). The dummy column (the last n_x+1
    #entries of vec(W)) gets zeros.
    Aeq1=sps.hstack((
        sps.kron(sps.eye(n_y,format='csr'),np.ones((1,n_x+1)),format='csr'),
        sps.csr_matrix((n_y,n_x+1)),
        sps.csr_matrix((n_y,nParam-WLen))),format='csr')
    #Constraint 2: sum_j W(i,j)=1 for every real row i
    #keep_x=[I_{n_x}, 0] selects the real rows from one column of W. Repeating
    #it once per column (kron with ones(1,n_y+1)) sums row i across all
    #columns, including the dummy column.
    keep_x=sps.hstack((sps.eye(n_x,format='csr'),
                       sps.csr_matrix((n_x,1))),format='csr')
    Aeq2=sps.hstack((
        sps.kron(np.ones((1,n_y+1)),keep_x,format='csr'),
        sps.csr_matrix((n_x,nParam-WLen))),format='csr')
    beq=np.ones(n_x+n_y)
    return sps.vstack((Aeq1,Aeq2),format='csr'),beq


def _edge_difference_matrix(n_x,n_y,X_adj,Y_adj):
    '''
    Sparse matrix G with G @ vec(W)[:n_y*(n_x+1)] = vec(A_X W_r - W_r A_Y),
    where W_r = W(:n_x,:n_y) is the block of real-to-real assignments.
    It has n_x*n_y rows, one per entry of the edge mismatch matrix.
    Only the first n_y columns of W are used, so the matrix has n_y*(n_x+1)
    columns. The caller pads the dummy column of W with zeros.

    Size: about n_y*nnz(A_X) + n_x*nnz(A_Y) nonzeros. This is O(n^3) for dense
    graphs and dominates the memory needed by the LP.
    '''
    #A_X W - W A_Y for the first n_x rows and n_y columns of W
    #X_pad=[A_X, 0] multiplies a full column of W (n_x+1 entries) and ignores
    #its dummy entry, so X_pad @ W(:,:n_y) = A_X @ W_r
    X_pad=sps.hstack((sps.csr_matrix(X_adj),
                      sps.csr_matrix((n_x,1))),format='csr')
    #keep_x=[I, 0] drops the dummy row: keep_x @ W(:,:n_y) = W_r
    keep_x=sps.hstack((sps.eye(n_x,format='csr'),
                       sps.csr_matrix((n_x,1))),format='csr')
    #vec(X_pad W I)     = (I kron X_pad) vec(W)       -> vec(A_X W_r)
    MX=sps.kron(sps.eye(n_y,format='csr'),X_pad,format='csr')
    #vec(keep_x W A_Y)  = (A_Y^T kron keep_x) vec(W)  -> vec(W_r A_Y)
    MY=sps.kron(sps.csr_matrix(Y_adj.T),keep_x,format='csr')
    return MX-MY


def _reverse_edge_difference_matrix(n_x,n_y,X_adj,Y_adj):
    '''
    Edge mismatch in the reverse direction, used for directed graphs:
    G2 @ vec(W)[:n_y*(n_x+1)] = vec(W_r A_Y^T - A_X^T W_r).
    That matrix is the transpose of A_Y W_r^T - W_r^T A_X, so it has the same
    entries and the same L1 norm. For symmetric (undirected) adjacency
    matrices it equals -(A_X W_r - W_r A_Y) and adds no new information.
    '''
    #A_Y W^T - W^T A_X for the first n_x rows and n_y columns of W
    keep_x=sps.hstack((sps.eye(n_x,format='csr'),
                       sps.csr_matrix((n_x,1))),format='csr')
    X_pad_t=sps.hstack((sps.csr_matrix(X_adj.T),
                        sps.csr_matrix((n_x,1))),format='csr')
    #vec(keep_x W A_Y^T) = (A_Y kron keep_x) vec(W)    -> vec(W_r A_Y^T)
    M1=sps.kron(sps.csr_matrix(Y_adj),keep_x,format='csr')
    #vec(X_pad_t W)      = (I kron X_pad_t) vec(W)     -> vec(A_X^T W_r)
    M2=sps.kron(sps.eye(n_y,format='csr'),X_pad_t,format='csr')
    return M1-M2


def _build_lp(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon,directed):
    '''
    Assemble the graph GOSPA LP as min f'x s.t. Aeq x = beq, x >= 0.

    Variables are [W (vectorised, column-major), s_plus, s_minus]. Every edge
    mismatch block G w is written as G w - s_plus + s_minus = 0 and penalised
    through the slacks, so G is stored once and there is no dense row.

    Why the slacks give the L1 norm: any vector g can be written as
    g = s_plus - s_minus with s_plus, s_minus >= 0. Both slacks have the
    positive cost edge_weight, so at the optimum at most one of
    s_plus[k], s_minus[k] is nonzero for each k. Then
    s_plus + s_minus = |g| entrywise, and edge_weight*sum(s_plus+s_minus)
    equals edge_weight*||G w||_1.

    The previous formulation used h >= G w, h >= -G w and e >= sum(h), which
    stores G twice and adds a dense row. This one has about half the nonzeros
    and the same optimal value.

    Returns f, Aeq, beq, the node cost matrix DAB, WLen = number of entries
    of W, and the per-slack weight edge_weight. The caller needs the last
    three to split the optimum into its cost components.
    '''
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]
    DAB=locCostComp(X_attr,Y_attr,c,p)

    nxny=n_x*n_y
    WLen=(n_x+1)*(n_y+1)
    #One mismatch block for undirected graphs, two (forward and reverse) for directed ones.
    #The weights follow the definition of the metric in [1]: eps^p/2 for the
    #undirected block, eps^p/4 for each directed block.
    blocks=[_edge_difference_matrix(n_x,n_y,X_adj,Y_adj)]
    edge_weight=epsilon**p/2
    if directed:
        blocks.append(_reverse_edge_difference_matrix(n_x,n_y,X_adj,Y_adj))
        edge_weight=epsilon**p/4
    #Number of mismatch entries; there is one s_plus and one s_minus per entry
    nSlack=len(blocks)*nxny
    nParam=WLen+2*nSlack

    ############# Objective function ################
    #Node costs on vec(W) (same column-major order as the variables), then
    #edge_weight on every slack
    f=np.empty(nParam)
    f[:WLen]=DAB.reshape(WLen,order='F')
    f[WLen:]=edge_weight

    ###########Equality constraints############
    #Assignment constraints on W
    Aeq_W,beq_W=_assignment_constraints(n_x,n_y,nParam)

    #Edge mismatch: G w - s_plus + s_minus = 0
    #Column layout: [ G | zeros for the dummy column of W | -I (s_plus) | +I (s_minus) ]
    I=sps.eye(nSlack,format='csr')
    G=blocks[0] if len(blocks)==1 else sps.vstack(blocks,format='csr')
    #blocks would keep the mismatch matrices alive, so drop it before the copies below
    del blocks
    Aeq_E=sps.hstack((G,sps.csr_matrix((nSlack,n_x+1)),-I,I),format='csr')
    #Aeq_E holds a copy of G; free G and I before the final vstack to lower the peak memory
    del G,I

    Aeq=sps.vstack((Aeq_W,Aeq_E),format='csr')
    beq=np.concatenate((beq_W,np.zeros(nSlack)))
    return f,Aeq,beq,DAB,WLen,edge_weight


def _solve_lp(f,Aeq,beq):
    #CPU solver: HiGHS interior point through scipy. Variables are >= 0 by default.
    #Returns (x, objective) so that GPU solvers can be swapped in with the same interface.
    res=linprog(f,A_eq=Aeq,b_eq=beq,method='highs-ipm')
    if not res.success:
        raise RuntimeError(f'LP solver failed: {res.message}')
    return res.x,res.fun


def _metric_components(dxy,loc_cost,false_cost,miss_cost,edge_cost,p):
    #Interior point solvers can return values slightly below zero (e.g. -1e-12)
    #for costs that are really zero. Clamp them so that the p-th root and the
    #returned components are well defined.
    def clamp(value):
        return 0.0 if value<=0 else value
    return (clamp(dxy)**(1/p),clamp(loc_cost),clamp(false_cost),
            clamp(miss_cost),clamp(edge_cost))


def _components_from_solution(x,obj,DAB,n_x,n_y,WLen,edge_weight,p):
    #Split the optimal objective (which is dxy^p) into its four parts
    #Undo the column-major vectorisation to recover W
    Wx=x[:WLen].reshape((n_x+1,n_y+1),order='F')
    #Assigned pairs: real rows and real columns of W
    loc_cost=np.sum(DAB[:n_x,:n_y]*Wx[:n_x,:n_y])
    #False nodes: nodes of Y assigned to the dummy row
    false_cost=np.sum(DAB[n_x,:n_y]*Wx[n_x,:n_y])
    #Missed nodes: nodes of X assigned to the dummy column
    miss_cost=np.sum(DAB[:n_x,n_y]*Wx[:n_x,n_y])
    #Edge mismatch: all slacks share the same weight
    edge_cost=edge_weight*np.sum(x[WLen:])
    return _metric_components(obj,loc_cost,false_cost,miss_cost,edge_cost,p)


def _graph_GOSPA(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon,directed,solve):
    '''
    Shared driver for the CPU and GPU backends.
    solve(f, Aeq, beq) -> (x, objective) solves min f'x s.t. Aeq x = beq, x >= 0.
    The CPU backend passes _solve_lp; gospa_gpu passes a cuOpt solver.
    '''
    X_attr,Y_attr,X_adj,Y_adj=_check_inputs(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon)
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]

    #An empty graph has a closed-form result and would give an LP with no assignment variables
    if n_x==0 or n_y==0:
        return _empty_result(n_x,n_y,c,p)

    f,Aeq,beq,DAB,WLen,edge_weight=_build_lp(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon,directed)
    x,obj=solve(f,Aeq,beq)
    return _components_from_solution(x,obj,DAB,n_x,n_y,WLen,edge_weight,p)


def LP_graph_GOSPA(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon):
    '''
    This function calculates the GOSPA metric between two undirected graphs using linear programming.
    Input:
    X_attr: NxD array of node attributes for graph X
    Y_attr: MxD array of node attributes for graph Y
    X_adj: NxN symmetric adjacency matrix for graph X
    Y_adj: MxM symmetric adjacency matrix for graph Y
    c: penalty for missing or false nodes
    p: p-norm
    epsilon: penalty for edge mismatch

    Returns:
    graph GOSPA cost, localisation cost, false node cost, miss node cost, edge mismatch cost.
    The four component costs are returned to the p-th power, so dxy^p equals their sum.
    '''
    return _graph_GOSPA(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon,False,_solve_lp)


def LP_graph_GOSPA_directed(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon):
    '''
    This function calculates the GOSPA metric between two directed graphs using linear programming.
    Input:

    X_attr: NxD array of node attributes for graph X
    Y_attr: MxD array of node attributes for graph Y
    X_adj: NxN adjacency matrix for graph X
    Y_adj: MxM adjacency matrix for graph Y
    c: penalty for missing or false nodes
    p: p-norm
    epsilon: penalty for edge mismatch

    Returns:
    graph GOSPA cost, localisation cost, false node cost, miss node cost, edge mismatch cost.
    The four component costs are returned to the p-th power, so dxy^p equals their sum.
    '''
    return _graph_GOSPA(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon,True,_solve_lp)
