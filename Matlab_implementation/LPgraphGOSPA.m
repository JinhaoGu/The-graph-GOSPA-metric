function [dxy,loc_cost,fa_cost,miss_cost,edge_cost]=LPgraphMetric(X, Y, c, p,epsilon)

nx = size(X.xState, 1);
ny = size(Y.xState, 1);

% disp(isdirected);
if(nx==0 && ny==0)
    dxy=0;
    loc_cost=0;
    miss_cost=0;
    fa_cost=0;
    edge_cost=0;
    return;
end

% if (isdirected=="undirected" || isdirected=="directed")
%     epsilon=epsilon/2;
% end


%%%%%%%%%% localisation cost computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAB = locCostComp_v2(X, Y, c, p);

[Wx,dxy,loc_cost, miss_cost, fa_cost, edge_cost]=LP_graph_metric(X,Y,DAB,nx,ny,c,p,epsilon);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% utils %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Wx,dxy,loc_cost, miss_cost, fa_cost, edge_cost]=LP_graph_metric(X,Y,DAB,nx,ny,c,p,epsilon)
% x = [W_1,1 W_2,1 .. W_nx+1,1, .. W_1,ny+1 W_2,ny+1 ...W_nx+1,ny+1,
% e,
% h_1,1 .. h_nx,ny ... h_1,1 ... h_nx,ny]'

%%% Length of the variable components in x
nxny=nx*ny;
nxny2=(nx+1)*(ny+1);
WLen=nxny2;
hLen=nxny;
eLen=1;
nParam = WLen + eLen + hLen;%total number of variables

%%% Position of the variable components in x
WPos = (1:WLen);
ePos=WLen+1;
hPos = WLen + eLen + (1:hLen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  objective function f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = zeros(nParam, 1);
f(WPos) = reshape(DAB, [WLen,1]); % for vec(W(1)) to vec(W(T)), loc cost
f(ePos) = epsilon^p/2; % edge cost
%%% TO DO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%  equality constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constraint 1

index_x=repmat(1:ny,nx+1,1);
index_x=index_x(:);
index_y=zeros(ny*(nx+1),1);
index_rep=1:ny*(nx+1);
index_y(index_rep)=index_rep;
Aeq1=sparse(index_x,index_y,1,ny,nParam);
beq1 = ones(ny, 1);


% %%%% Constraint 2 %%%%

index_x=repmat(reshape(1:nx,nx,1),ny+1,1);
index_x=index_x(:);
index_y=repmat(1:nx+1:(nx+1)/nx*length(index_x),nx,1);
index_y2=repmat((0:nx-1)',1,size(index_y,2));
index_y=index_y2+index_y;
index_y=index_y(:);
Aeq2=sparse(index_x,index_y,1,nx,nParam);
beq2 = ones(nx, 1);
Aeq = [Aeq1; Aeq2]; beq = [beq1; beq2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  upper and lower bound constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 <= W < inf, 0 <= e < inf, 0 <= h < inf
lb = zeros(nParam, 1);
ub = inf(nParam, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%  inequality constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Constraint 1 %%%%

index_minus_x=1;
index_minus_y=WLen+index_minus_x;
value_minus=-1;
index_one_x=ones(nxny,1);
index_one_y=WLen + eLen+(1:nxny) ;
value_one=ones(1,length(index_one_y));
A1=sparse([index_minus_x';index_one_x],[index_minus_y';index_one_y'],[value_minus';value_one']);


%%%% Constraint 2 %%%%
index_1_x=repmat(1:nxny,nx,1);
index_1_x=index_1_x(:);

index_1_y=repmat(1:nx+1:(nx+1)/nx*nxny,nx,1);
index_1_y2=repmat((0:nx-1)',1,size(index_1_y,2));
index_1_y=index_1_y2+index_1_y;
index_1_y=repmat(index_1_y,1,nx);
index_1_y=index_1_y(:);
mask=sparse(index_1_x,index_1_y,1,nxny,(nx+1)*ny);


A_adj1=repmat([X.adj,zeros(nx,1)],1,ny);
ind=repmat(1:nx,ny,1);
ind=ind(:);
new_adj=zeros(nxny,(nx+1)*ny);
% new_adj=[];
for i =1:nx*ny
    new_adj(i,:)=A_adj1(ind(i),:);
%     new_adj=[new_adj;A_adj1(ind(i),:)];
end
new_adj=sparse(new_adj);
A_adj_1=[mask.*new_adj,zeros(nxny,nx+1)];


index_2_x=repmat(1:nxny,ny,1);
index_2_x=index_2_x(:);
index_2_y=repmat(1:nx+1:(nx+1)/nx*nxny,nx,1);
index_2_y2=repmat((0:nx-1)',1,size(index_1_y,2));
% index_2_y3=repmat((0:ny-1),1,size(index_1_y,2));
index_2_y=index_2_y2+index_2_y;
index_2_y=repmat(index_2_y',ny,1);
index_2_y=index_2_y(:);

mask=sparse(index_2_x,index_2_y,1,nxny,ny*(nx+1));
A_adj2=repmat(Y.adj',nx,1); %Need a transpose for adjacency matrix 
ind=repmat(1:ny,nx+1,1);
ind=ind(:);
new_adj=zeros(nxny,(nx+1)*ny);
% new_adj=[];
for i =1:(nx+1)*ny
    new_adj(:,i)=A_adj2(:,ind(i));
%     new_adj=[new_adj,A_adj2(:,ind(i))];
end

A_adj_2=[mask.*new_adj,zeros(nxny,nx+1)];


A2_adj=[A_adj_1-A_adj_2,zeros(nxny,eLen),-1*(zeros(nxny,hLen)+eye(hLen))];

A3_adj=[A_adj_2-A_adj_1,zeros(nxny,eLen),-1*(zeros(nxny,hLen)+eye(hLen))];

A=[A1;A2_adj;A3_adj];

b=sparse(eLen+nxny+nxny,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  optimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linProgOptions = optimoptions('linprog','Display','off','Algorithm','interior-point'); %If this line returns an error, 
% it may be required to install Matlab optimization toolbox
[x, dxy] = linprog(f, A, b, Aeq, beq, lb, ub, linProgOptions);%% deleted []
dxy=dxy^(1/p);

Wx = reshape(x(1:WLen),[nx+1,ny+1]);
% disp(Wx);
loc_cost=sum(DAB(1:nx,1:ny).*Wx(1:nx,1:ny),"all")^(1/p);
fa_cost=sum(DAB(nx+1,1:ny).*Wx(nx+1,1:ny))^(1/p);
miss_cost=sum(DAB(1:nx,ny+1).*Wx(1:nx,1+ny))^(1/p);
% edge_cost=epsilon^p*sum(abs(X.adj*Wx(1:nx,1:ny)-Wx(1:nx,1:ny)*Y.adj),'all');
edge_cost=double(epsilon^p/2*x(WLen+eLen)).^double(1/p);
%%% TO DO
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function locCostMat = locCostComp_v2(X, Y, c, p)
% function locCostMat = locCostComp(stMat, X, Y, c, p)
% computing the localisation cost for every (i,j)
nx = size(X.xState, 1);
ny = size(Y.xState, 1);
tmpCost = c^p / 2; % cost for being unassigned
locCostMat  = zeros(nx+1, ny+1);
for xind = 1:nx+1
    if (xind<=nx)% xind not dummy
        for yind = 1:ny+1
            if (yind<=ny)% yind not dummy
            locCostMat(xind,yind)=computeLocCostPerTime( ...
                                X.xState(xind,:), Y.xState(yind,:), c, p);
            else
            locCostMat(xind,yind) = tmpCost;
            end
        end
    else % xind is dummy
        for yind=1:ny 
            locCostMat(xind,yind) = tmpCost;
        end
    end
end
end


function d = computeLocCostPerTime(x, y, c, p)
if all(~isnan(x)) && all(~isnan(y))
    % neither x nor y has hole
%     d = min(norm(x-y, p)^p,c^p);
      d = norm(double(x-y))^p; 
elseif any(isnan(x) & ~isnan(y)) || any(~isnan(x) & isnan(y))
    % exactly one of x and y has hole
    d = c^p/2;
else
    d = 0;
end
end
