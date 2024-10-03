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
        return np.linalg.norm(x-y)**p

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
    locCostMat=np.zeros((n_x+1,n_y+1))

    for i in range(n_x+1):
        if i<=n_x-1:# x not dummy
            for j in range(n_y+1):
                if j<=n_y-1: # y not dummy
                    locCostMat[i,j]=computeLocCostPerTime(X_attr[i,:],Y_attr[j,:],c,p)
                else:
                    locCostMat[i,j]=tmpCost
        
        else:
            for j in range(n_y): # x is dummy
                locCostMat[i,j]=tmpCost

    return locCostMat


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
    graph GOSPA cost, localisation cost, false node cost, miss node cost, edge mismatch cost
    '''
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]
    DAB=locCostComp(X_attr,Y_attr,c,p)

    nxny=n_x*n_y
    nxny2=(n_x+1)*(n_y+1)
    WLen=nxny2
    h1Len=nxny

    eLen=1
    nParam=WLen+eLen+h1Len
    WPos=np.arange(WLen)
    e1Pos=np.arange(WLen,WLen+1)
    h1Pos=np.arange(WLen+eLen,WLen+eLen+h1Len)
    # print(nParam)
    
    ############# Objective function ################
    f=np.zeros([nParam,1])
    f[WPos]=np.reshape(DAB,(WLen,1),order='F')
    f[e1Pos]=epsilon**p/2

    ###########Equality constraint############
    #Constraint 1
    index_x=np.tile(np.arange(n_y),(n_x+1,1))
    index_x=index_x.flatten(order='F')
    # index_y=np.zeros([n_y*(n_x+1),1])
    index_y=np.arange(n_y*(n_x+1))
    index_y=index_y.flatten(order='F')
    

    Aeq1=sps.coo_matrix((np.ones(len(index_x)),(index_x,index_y)),shape=(n_y,nParam))
    beq1=np.ones([n_y,1])


    #Constraint 2
    index_x=np.tile(np.arange(n_x),(n_y+1,1)).T
    index_x=index_x.flatten(order='F')

    index_y=np.tile(np.arange(1,(n_x+1)/n_x*len(index_x),step=n_x+1),(n_x,1))-1
    index_y2=np.tile(np.arange(n_x),(np.size(index_y,1),1)).T
    index_y=index_y+index_y2
    index_y=index_y.flatten(order='F')
    Aeq2=sps.coo_matrix((np.ones(len(index_x)),(index_x,index_y)),shape=(n_x, nParam))
    beq2=np.ones([n_x,1])
    Aeq=sps.vstack((Aeq1,Aeq2))
    beq=np.vstack((beq1,beq2))

    #inequality constraint
    #Constraint 1
    index_minus_x=0
    index_minus_y=WLen+index_minus_x
    value_minus=-1
    index_one_x=np.hstack((index_minus_x,np.zeros(nxny)))
    index_one_y=WLen+eLen+np.arange(nxny)
    index_one_y=np.hstack((index_minus_y,index_one_y))
    value_one=np.hstack((value_minus,np.ones(nxny)))
    A1=sps.coo_matrix((value_one,(index_one_x,index_one_y)),shape=(1,nParam))

    #Constraint 2
    index_1_x=np.tile(np.arange(nxny),(n_x,1))

    index_1_x=index_1_x.flatten(order='F')
    index_1_y=np.tile(np.arange(1,(n_x+1)/n_x*nxny,step=n_x+1),(n_x,1))-1 # todo
    index_1_y2=np.tile(np.arange(n_x),(np.size(index_1_y,1),1)).T
    index_1_y=index_1_y+index_1_y2
    index_1_y=np.tile(index_1_y,(1,n_x))
    index_1_y=np.reshape(index_1_y,(nxny*n_x),order='F')
    mask=sps.coo_matrix((np.ones(len(index_1_x)),(index_1_x,index_1_y)),shape=(nxny,(n_x+1)*n_y))

    A_adj1=np.tile(np.hstack((X_adj,np.zeros([n_x,1]))),(1,n_y))
    ind=np.tile(np.arange(n_x),(n_y,1))
    ind=ind.flatten(order='F')

    new_adj=np.zeros([nxny,(n_x+1)*n_y])

    for i in range(nxny):
        new_adj[i,:]=A_adj1[ind[i],:]


    new_adj=sps.coo_matrix(new_adj)

    A_adj_1=sps.hstack((new_adj.multiply(mask),np.zeros([nxny,n_x+1])))


    index_2_x=np.tile(np.arange(nxny),(n_y,1))
    index_2_x=index_2_x.flatten(order='F')
    index_2_y=np.tile(np.arange(1,(n_x+1)/n_x*nxny,step=n_x+1),(n_x,1))-1
    index_2_y2=np.tile(np.arange(n_x),(np.size(index_2_y,1),1)).T
    index_2_y=index_2_y+index_2_y2
    index_2_y=np.tile(index_2_y,(1,n_y)).T
    index_2_y=index_2_y.flatten(order='F')
    mask=sps.coo_matrix((np.ones(len(index_2_x)),(index_2_x,index_2_y)),shape=(nxny,(n_x+1)*n_y))

    A_adj2=np.tile(Y_adj.T,(n_x,1))
    ind=np.tile(np.arange(n_y),(n_x+1,1))
    ind=ind.flatten(order='F')
    new_adj=np.zeros([nxny,(n_x+1)*n_y])

    for i in range((n_x+1)*n_y):
        new_adj[:,i]=A_adj2[:,ind[i]]

    new_adj=sps.coo_matrix(new_adj)
    A_adj_2=sps.hstack((new_adj.multiply(mask),np.zeros([nxny,n_x+1])))

    A2_adj=sps.hstack(
                        (A_adj_1-A_adj_2,np.zeros([nxny,eLen]),
                        -1*(np.zeros([nxny,h1Len])+
                        np.eye(h1Len))
                        )
                    )

    A3_adj=sps.hstack(
                        (A_adj_2-A_adj_1,np.zeros([nxny,eLen]),
                        -1*(np.zeros([nxny,h1Len])+np.eye(h1Len))
                        )
                    )


    A=sps.vstack((A1,A2_adj,A3_adj))
    lenb,_=A.shape
    b=np.zeros([lenb,1])

    # lb=np.zeros([nParam,1])
    # ub=np.full([nParam,1],None)
    #Solve the LP
    res=linprog(f,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=beq,method='highs-ipm')

    W=res.x
    Wx=np.reshape(W[0:nxny2],(n_x+1,n_y+1),order='F')
    loc_cost=np.sum(np.multiply(DAB[0:n_x,0:n_y],Wx[0:n_x,0:n_y]))
    false_cost=np.sum(np.multiply(DAB[n_x,0:n_y],Wx[n_x,0:n_y]))
    miss_cost=np.sum(np.multiply(DAB[0:n_x,n_y],Wx[0:n_x,n_y]))
    edge_cost=epsilon**p/2 * W[e1Pos].item()
    
    
    dxy=res.fun+1e-9
    loc_cost=loc_cost+1e-9
    false_cost=loc_cost+1e-9
    miss_cost=miss_cost+1e-9
    edge_cost=edge_cost+1e-9

    return dxy**(1/p),loc_cost**(1/p),false_cost**(1/p),miss_cost**(1/p),edge_cost**(1/p)



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
    graph GOSPA cost, localisation cost, false node cost, miss node cost, edge mismatch cost
    '''
    
    n_x=X_attr.shape[0]
    n_y=Y_attr.shape[0]
    DAB=locCostComp(X_attr,Y_attr,c,p)

    nxny=n_x*n_y
    nxny2=(n_x+1)*(n_y+1)
    WLen=nxny2
    h1Len=nxny
    h2Len=nxny
    eLen=2
    nParam=WLen+eLen+h1Len+h2Len
    WPos=np.arange(WLen)
    e1Pos=np.arange(WLen,WLen+1)
    e2Pos=np.arange(WLen+1,WLen+2)
    h1Pos=np.arange(WLen+eLen,WLen+eLen+h1Len) 
    h2Pos=np.arange(WLen+eLen+h1Len,WLen+eLen+h1Len+h2Len)
    
    ############# Objective function ################
    f=np.zeros([nParam,1])
    f[WPos]=np.reshape(DAB,(WLen,1),order='F')
    f[e1Pos]=epsilon**p/4
    f[e2Pos]=epsilon**p/4

    ###########Equality constraint############
    #Constraint 1
    index_x=np.tile(np.arange(n_y),(n_x+1,1))
    index_x=index_x.flatten(order='F')
    # index_y=np.zeros([n_y*(n_x+1),1])
    index_y=np.arange(n_y*(n_x+1))
    index_y=index_y.flatten(order='F')
    # index_y[index_rep]=index_rep

    Aeq1=sps.coo_matrix((np.ones(len(index_x)),(index_x,index_y)),shape=(n_y,nParam))
    beq1=np.ones([n_y,1])


    #Constraint 2
    index_x=np.tile(np.arange(n_x),(n_y+1,1)).T
    index_x=index_x.flatten(order='F')

    index_y=np.tile(np.arange(1,(n_x+1)/n_x*len(index_x),step=n_x+1),(n_x,1))-1
    index_y2=np.tile(np.arange(n_x),(np.size(index_y,1),1)).T
    index_y=index_y+index_y2
    index_y=index_y.flatten(order='F')
    Aeq2=sps.coo_matrix((np.ones(len(index_x)),(index_x,index_y)),shape=(n_x, nParam))
    beq2=np.ones([n_x,1])
    Aeq=sps.vstack((Aeq1,Aeq2))
    beq=np.vstack((beq1,beq2))

    #inequality constraint
    #Constraint 1
    index_minus_x=0
    index_minus_y=WLen+index_minus_x
    value_minus=-1
    index_one_x=np.hstack((index_minus_x,np.zeros(nxny+nxny)))
    index_one_y=WLen+eLen+np.arange(nxny+nxny)
    index_one_y=np.hstack((index_minus_y,index_one_y))
    value_one=np.hstack((value_minus,np.ones(nxny),np.zeros(nxny)))
    Ae1=sps.coo_matrix((value_one,(index_one_x,index_one_y)),shape=(1,nParam))


    index_minus_x=0
    index_minus_y=WLen+index_minus_x+1
    value_minus=-1
    index_one_x=np.hstack((index_minus_x,np.zeros(nxny)))
    index_one_y=WLen+eLen+nxny+np.arange(nxny)
    index_one_y=np.hstack((index_minus_y,index_one_y))
    value_one=np.hstack((value_minus,np.ones(nxny)))
    Ae2=sps.coo_matrix((value_one,(index_one_x,index_one_y)),shape=(1,nParam))

    A1=sps.vstack((Ae1,Ae2))
    
    #Constraint 2
    index_1_x=np.tile(np.arange(nxny),(n_x,1))

    index_1_x=index_1_x.flatten(order='F')
    index_1_y=np.tile(np.arange(1,(n_x+1)/n_x*nxny,step=n_x+1),(n_x,1))-1 # todo
    index_1_y2=np.tile(np.arange(n_x),(np.size(index_1_y,1),1)).T
    index_1_y=index_1_y+index_1_y2
    index_1_y=np.tile(index_1_y,(1,n_x))
    index_1_y=np.reshape(index_1_y,(nxny*n_x),order='F')
    mask=sps.coo_matrix((np.ones(len(index_1_x)),(index_1_x,index_1_y)),shape=(nxny,(n_x+1)*n_y))

    A_adj1=np.tile(np.hstack((X_adj,np.zeros([n_x,1]))),(1,n_y))
    ind=np.tile(np.arange(n_x),(n_y,1))
    ind=ind.flatten(order='F')

    new_adj=np.zeros([nxny,(n_x+1)*n_y])

    for i in range(nxny):
        new_adj[i,:]=A_adj1[ind[i],:]


    new_adj=sps.coo_matrix(new_adj)

    A_adj_1=sps.hstack((new_adj.multiply(mask),np.zeros([nxny,n_x+1])))


    index_2_x=np.tile(np.arange(nxny),(n_y,1))
    index_2_x=index_2_x.flatten(order='F')
    index_2_y=np.tile(np.arange(1,(n_x+1)/n_x*nxny,step=n_x+1),(n_x,1))-1
    index_2_y2=np.tile(np.arange(n_x),(np.size(index_2_y,1),1)).T
    index_2_y=index_2_y+index_2_y2
    index_2_y=np.tile(index_2_y,(1,n_y)).T
    index_2_y=index_2_y.flatten(order='F')
    mask=sps.coo_matrix((np.ones(len(index_2_x)),(index_2_x,index_2_y)),shape=(nxny,(n_x+1)*n_y))

    A_adj2=np.tile(Y_adj.T,(n_x,1))
    ind=np.tile(np.arange(n_y),(n_x+1,1))
    ind=ind.flatten(order='F')
    new_adj=np.zeros([nxny,(n_x+1)*n_y])

    for i in range((n_x+1)*n_y):
        new_adj[:,i]=A_adj2[:,ind[i]]

    new_adj=sps.coo_matrix(new_adj)
    A_adj_2=sps.hstack((new_adj.multiply(mask),np.zeros([nxny,n_x+1])))

    A2_adj=sps.hstack(
                        (A_adj_1-A_adj_2,np.zeros([nxny,eLen]),
                        -1*(np.zeros([nxny,h1Len+h2Len])+
                        np.hstack((np.eye(h1Len),np.zeros([h2Len,h2Len]))))
                        )
                    )

    A3_adj=sps.hstack(
                        (A_adj_2-A_adj_1,np.zeros([nxny,eLen]),
                        -1*(np.zeros([nxny,h1Len+h2Len])+
                        np.hstack((np.eye(h1Len),np.zeros([h2Len,h2Len]))))
                        )
                    )


    #Contraint 3

    index_1_x=np.tile(np.arange(nxny),(n_y,1))
    index_1_x=index_1_x.flatten(order='F')

    index_1_y=np.tile(np.arange(1,(n_x+1)/n_x*nxny,step=n_x+1),(n_x,1))-1
    index_1_y2=np.tile(np.arange(n_x),(np.size(index_1_y,1),1)).T
    index_1_y=index_1_y+index_1_y2
    index_1_y=np.tile(index_1_y.T,(1,n_y))
    index_1_y=index_1_y.flatten(order='F')
    mask=sps.coo_matrix((np.ones(len(index_1_x)),(index_1_x,index_1_y)),shape=(nxny,(n_x+1)*n_y))

    A_adj3=np.repeat(np.repeat(Y_adj,n_x,axis=1),n_x,axis=0)
    ind=np.tile(np.arange(1,(n_x+1)/n_x*nxny,step=n_x+1),(n_x,1))-1
    ind2=np.tile(np.arange(n_x),(np.size(ind,1),1)).T
    ind=ind+ind2
    ind=ind.flatten(order='F')
    ind=ind.astype(int)

    new_adj=np.zeros([nxny,(n_x+1)*n_y])

    for i in range(nxny):
        new_adj[:,ind[i]]=A_adj3[:,i]

    new_adj=sps.coo_matrix(new_adj)
    A_adj_3=sps.hstack((new_adj.multiply(mask),np.zeros([nxny,n_x+1])))



    index_2_x=np.tile(np.arange(nxny),(n_x,1))
    index_2_x=index_2_x.flatten(order='F')
    index_2_y=np.tile(np.arange((n_x+1)*n_y,step=n_x+1),(n_x,1))
    index_2_y2=np.tile(np.arange(n_x),(n_y,1)).T
    index_2_y=index_2_y+index_2_y2
    index_2_y=np.tile(index_2_y.T,(1,n_x)).T
    index_2_y=index_2_y.flatten(order='F')

    mask=sps.coo_matrix((np.ones(len(index_2_x)),(index_2_x,index_2_y)),shape=(nxny,(n_x+1)*n_y))
    A_adj4=np.tile(np.hstack((X_adj.T,np.zeros([n_x,1]))),(n_y,n_y))

    new_adj=sps.coo_matrix(A_adj4)

    A_adj_4=sps.hstack((new_adj.multiply(mask),np.zeros([nxny,n_x+1])))


    A4_adj=sps.hstack(
                        (A_adj_3-A_adj_4,np.zeros([nxny,eLen]),
                        -1*(np.zeros([nxny,h1Len+h2Len])+np.hstack((np.zeros([h1Len,h1Len]),np.eye(h2Len))))
                        )
                    )

    A5_adj=sps.hstack(
                        (A_adj_4-A_adj_3,np.zeros([nxny,eLen]),
                        -1*(np.zeros([nxny,h1Len+h2Len])+np.hstack((np.zeros([h1Len,h1Len]),np.eye(h2Len))))
                        )
                    )



    A=sps.vstack((A1,A2_adj,A3_adj,A4_adj,A5_adj))
    lenb,_=A.shape
    b=np.zeros([lenb,1])
    
    # lb=np.zeros([nParam,1])
    # ub=np.full([nParam,1],None)
    #Solve the LP
    res=linprog(f,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=beq,method='highs-ipm')

    W=res.x
    Wx=np.reshape(W[0:nxny2],(n_x+1,n_y+1),order='F')
    loc_cost=np.sum(np.multiply(DAB[0:n_x,0:n_y],Wx[0:n_x,0:n_y]))
    false_cost=np.sum(np.multiply(DAB[n_x,0:n_y],Wx[n_x,0:n_y]))
    miss_cost=np.sum(np.multiply(DAB[0:n_x,n_y],Wx[0:n_x,n_y]))
    edge_cost=epsilon**p/4 * (W[e1Pos].item()+W[e2Pos].item())
    
    dxy=res.fun+1e-9
    loc_cost=loc_cost+1e-9
    false_cost=loc_cost+1e-9
    miss_cost=miss_cost+1e-9
    edge_cost=edge_cost+1e-9
    
    return dxy**(1/p),loc_cost**(1/p),false_cost**(1/p),miss_cost**(1/p),edge_cost**(1/p)
