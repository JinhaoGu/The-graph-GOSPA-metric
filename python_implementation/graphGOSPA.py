#Author: Jinhao Gu
#This code is a python implementation of the graph GOSPA metric proposed in the paper 
# "Graph GOSPA metric: a metric to measure the discrepancy between graphs of different sizes"
# by Jinhao Gu, Á. F. García-Fernández, Robert E. Firth, Lennart Svensson
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog


def computeLocCostPerTime(x,y,c,p):
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
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]
    tmpCost=c**p/2
    locCostMat=np.full((n_x+1,n_y+1),tmpCost)
    locCostMat[-1,-1]=0

    if n_x>0 and n_y>0:
        x_nan=np.isnan(X_attr)[:,None,:]
        y_nan=np.isnan(Y_attr)[None,:,:]
        no_hole=~(x_nan.any(axis=2) | y_nan.any(axis=2))
        one_hole=(x_nan & ~y_nan).any(axis=2) | (~x_nan & y_nan).any(axis=2)
        dist=np.linalg.norm(X_attr[:,None,:]-Y_attr[None,:,:],axis=2)
        locCostMat[:n_x,:n_y]=np.where(no_hole,dist**p,np.where(one_hole,tmpCost,0.0))

    return locCostMat


def _is_numeric_dtype(dtype):
    return (np.issubdtype(dtype,np.integer) or
            np.issubdtype(dtype,np.floating) or
            np.issubdtype(dtype,np.bool_))


def _as_attribute_array(attr,name):
    attr=np.asarray(attr)
    if attr.size==0 and attr.ndim<2:
        attr=attr.reshape(0,0)
    if attr.ndim!=2:
        raise ValueError(f'{name} must be a 2-D array of node attributes')
    if not _is_numeric_dtype(attr.dtype):
        raise ValueError(f'{name} must contain numeric node attributes')
    return attr


def _as_adjacency_array(adj,size,name):
    adj=np.asarray(adj)
    if adj.shape!=(size,size):
        raise ValueError(f'{name} must have shape ({size}, {size}), got {adj.shape}')
    if not _is_numeric_dtype(adj.dtype):
        raise ValueError(f'{name} must contain numeric entries')
    return adj


def _check_inputs(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon):
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
    miss_cost_p=n_x*c**p/2
    false_cost_p=n_y*c**p/2
    dxy=(miss_cost_p+false_cost_p)**(1/p)
    return dxy,0.0,false_cost_p,miss_cost_p,0.0


def _assignment_constraints(n_x,n_y,nParam):
    WLen=(n_x+1)*(n_y+1)
    #Constraint 1: sum_i W(i,j)=1 for every real column j
    Aeq1=sps.hstack((
        sps.kron(sps.eye(n_y,format='csr'),np.ones((1,n_x+1))),
        sps.csr_matrix((n_y,n_x+1)),
        sps.csr_matrix((n_y,nParam-WLen))),format='csr')
    #Constraint 2: sum_j W(i,j)=1 for every real row i
    keep_x=sps.hstack((sps.eye(n_x,format='csr'),
                       sps.csr_matrix((n_x,1))),format='csr')
    Aeq2=sps.hstack((
        sps.kron(np.ones((1,n_y+1)),keep_x),
        sps.csr_matrix((n_x,nParam-WLen))),format='csr')
    beq=np.ones(n_x+n_y)
    return sps.vstack((Aeq1,Aeq2),format='csr'),beq


def _edge_difference_matrix(n_x,n_y,X_adj,Y_adj):
    #A_X W - W A_Y for the first n_x rows and n_y columns of W
    X_pad=sps.hstack((sps.csr_matrix(X_adj),
                      sps.csr_matrix((n_x,1))),format='csr')
    keep_x=sps.hstack((sps.eye(n_x,format='csr'),
                       sps.csr_matrix((n_x,1))),format='csr')
    MX=sps.kron(sps.eye(n_y,format='csr'),X_pad)
    MY=sps.kron(sps.csr_matrix(Y_adj.T),keep_x)
    return MX-MY


def _reverse_edge_difference_matrix(n_x,n_y,X_adj,Y_adj):
    #A_Y W^T - W^T A_X for the first n_x rows and n_y columns of W
    keep_x=sps.hstack((sps.eye(n_x,format='csr'),
                       sps.csr_matrix((n_x,1))),format='csr')
    X_pad_t=sps.hstack((sps.csr_matrix(X_adj.T),
                        sps.csr_matrix((n_x,1))),format='csr')
    M1=sps.kron(sps.csr_matrix(Y_adj),keep_x)
    M2=sps.kron(sps.eye(n_y,format='csr'),X_pad_t)
    return M1-M2


def _solve_lp(f,A,b,Aeq,beq):
    res=linprog(f,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=beq,method='highs-ipm')
    if not res.success:
        raise RuntimeError(f'LP solver failed: {res.message}')
    return res


def _metric_components(dxy,loc_cost,false_cost,miss_cost,edge_cost,p):
    def clamp(value):
        return 0.0 if value<=0 else value
    return (clamp(dxy)**(1/p),clamp(loc_cost),clamp(false_cost),
            clamp(miss_cost),clamp(edge_cost))


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
    X_attr,Y_attr,X_adj,Y_adj=_check_inputs(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon)
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]

    if n_x==0 or n_y==0:
        return _empty_result(n_x,n_y,c,p)

    DAB=locCostComp(X_attr,Y_attr,c,p)

    nxny=n_x*n_y
    WLen=(n_x+1)*(n_y+1)
    nParam=WLen+1+nxny

    ############# Objective function ################
    f=np.zeros(nParam)
    f[:WLen]=DAB.reshape(WLen,order='F')
    f[WLen]=epsilon**p/2

    ###########Equality constraints############
    Aeq,beq=_assignment_constraints(n_x,n_y,nParam)

    ###########Inequality constraints############
    #Constraint 1: e1 >= sum(h1)
    A1=sps.hstack((
        sps.csr_matrix((1,WLen)),
        sps.csr_matrix(np.array([[-1.0]])),
        np.ones((1,nxny))),format='csr')

    #Constraints 2 and 3: h1 >= +/- (A_X W - W A_Y)
    G=sps.hstack((_edge_difference_matrix(n_x,n_y,X_adj,Y_adj),
                  sps.csr_matrix((nxny,n_x+1))),format='csr')
    neg_I=-sps.eye(nxny,format='csr')
    A2=sps.hstack((G,sps.csr_matrix((nxny,1)),neg_I),format='csr')
    A3=sps.hstack((-G,sps.csr_matrix((nxny,1)),neg_I),format='csr')

    A=sps.vstack((A1,A2,A3),format='csr')
    b=np.zeros(A.shape[0])

    #Solve the LP
    res=_solve_lp(f,A,b,Aeq,beq)

    Wx=res.x[:WLen].reshape((n_x+1,n_y+1),order='F')
    loc_cost=np.sum(DAB[:n_x,:n_y]*Wx[:n_x,:n_y])
    false_cost=np.sum(DAB[n_x,:n_y]*Wx[n_x,:n_y])
    miss_cost=np.sum(DAB[:n_x,n_y]*Wx[:n_x,n_y])
    edge_cost=epsilon**p/2*res.x[WLen]

    return _metric_components(res.fun,loc_cost,false_cost,miss_cost,edge_cost,p)


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
    X_attr,Y_attr,X_adj,Y_adj=_check_inputs(X_attr,Y_attr,X_adj,Y_adj,c,p,epsilon)
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]

    if n_x==0 or n_y==0:
        return _empty_result(n_x,n_y,c,p)

    DAB=locCostComp(X_attr,Y_attr,c,p)

    nxny=n_x*n_y
    WLen=(n_x+1)*(n_y+1)
    nParam=WLen+2+2*nxny

    ############# Objective function ################
    f=np.zeros(nParam)
    f[:WLen]=DAB.reshape(WLen,order='F')
    f[WLen]=epsilon**p/4
    f[WLen+1]=epsilon**p/4

    ###########Equality constraints############
    Aeq,beq=_assignment_constraints(n_x,n_y,nParam)

    ###########Inequality constraints############
    #Constraint 1: e1 >= sum(h1) and e2 >= sum(h2)
    minus_one=sps.csr_matrix(np.array([[-1.0]]))
    A1_1=sps.hstack((
        sps.csr_matrix((1,WLen)),minus_one,sps.csr_matrix((1,1)),
        np.ones((1,nxny)),sps.csr_matrix((1,nxny))),format='csr')
    A1_2=sps.hstack((
        sps.csr_matrix((1,WLen)),sps.csr_matrix((1,1)),minus_one,
        sps.csr_matrix((1,nxny)),np.ones((1,nxny))),format='csr')
    A1=sps.vstack((A1_1,A1_2),format='csr')

    #Constraints 2 and 3: h1 >= +/- (A_X W - W A_Y)
    G1=sps.hstack((_edge_difference_matrix(n_x,n_y,X_adj,Y_adj),
                   sps.csr_matrix((nxny,n_x+1))),format='csr')
    #Constraints 4 and 5: h2 >= +/- (A_Y W^T - W^T A_X)
    G2=sps.hstack((_reverse_edge_difference_matrix(n_x,n_y,X_adj,Y_adj),
                   sps.csr_matrix((nxny,n_x+1))),format='csr')

    zeros_e=sps.csr_matrix((nxny,2))
    neg_I_h1=sps.hstack((-sps.eye(nxny,format='csr'),
                         sps.csr_matrix((nxny,nxny))),format='csr')
    neg_I_h2=sps.hstack((sps.csr_matrix((nxny,nxny)),
                         -sps.eye(nxny,format='csr')),format='csr')
    A2=sps.hstack((G1,zeros_e,neg_I_h1),format='csr')
    A3=sps.hstack((-G1,zeros_e,neg_I_h1),format='csr')
    A4=sps.hstack((G2,zeros_e,neg_I_h2),format='csr')
    A5=sps.hstack((-G2,zeros_e,neg_I_h2),format='csr')

    A=sps.vstack((A1,A2,A3,A4,A5),format='csr')
    b=np.zeros(A.shape[0])

    #Solve the LP
    res=_solve_lp(f,A,b,Aeq,beq)

    Wx=res.x[:WLen].reshape((n_x+1,n_y+1),order='F')
    loc_cost=np.sum(DAB[:n_x,:n_y]*Wx[:n_x,:n_y])
    false_cost=np.sum(DAB[n_x,:n_y]*Wx[n_x,:n_y])
    miss_cost=np.sum(DAB[:n_x,n_y]*Wx[:n_x,n_y])
    edge_cost=epsilon**p/4*(res.x[WLen]+res.x[WLen+1])

    return _metric_components(res.fun,loc_cost,false_cost,miss_cost,edge_cost,p)
