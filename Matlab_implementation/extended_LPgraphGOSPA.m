function [dxy,loc_cost,fa_cost,miss_cost,edge_cost]=extended_LPgraphGOSPA(X, Y, c, p,epsilon)

nx = size(X.xState, 1);
ny = size(Y.xState, 1);

if(nx==0 && ny==0)
    dxy=0;
    loc_cost=0;
    miss_cost=0;
    fa_cost=0;
    edge_cost=0;
    return;
end


%%%%%%%%%% localisation cost computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAB = locCostComp_v2(X, Y, c, p);

[dxy,loc_cost, miss_cost, fa_cost, edge_cost]=LP_graph_metric(X,Y,DAB,nx,ny,c,p,epsilon);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% utils %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Wx,dxy,loc_cost, miss_cost, fa_cost, edge_cost]=LP_graph_metric(X,Y,DAB,nx,ny,c,p,epsilon)
% x = [W_1,1 W_2,1 .. W_nx+1,1, .. W_1,ny+1 W_2,ny+1 ...W_nx+1,ny+1,
% e1,e2,
% h1_1,1 .. h1_nx,1 ... h1_nx,ny]',
% h2_1,1 .. h2_nx,1 ... h2_nx,ny]'

%%% Length of the variable components in x
nxny=nx*ny;
nxny2=(nx+1)*(ny+1);
WLen=nxny2;
h1Len=nxny;
h2Len=nxny;
eLen=2;

nParam = WLen + eLen+ h1Len+h2Len;%total number of variables

%%% Position of the variable components in x
WPos = (1:WLen);
e1Pos=WLen+1;
e2Pos=e1Pos+1;
h1Pos = WLen + eLen+  (1:h1Len);
h2Pos =WLen + eLen+h1Len+(1:h2Len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%  objective function f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = zeros(nParam, 1);
f(WPos) = reshape(DAB, [WLen,1]); % for vec(W(1)) to vec(W(T)), loc cost
f(e1Pos) = epsilon^p/4; % edge cost
f(e2Pos) = epsilon^p/4;
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
index_one_x=ones(nxny+nxny,1);
index_one_y=WLen + eLen+(1:nxny+nxny) ;
value_one=[ones(1,nxny),zeros(1,nxny)];
A_ie1=sparse([index_minus_x';index_one_x],[index_minus_y';index_one_y'],[value_minus';value_one']);

index_minus_x=1;
index_minus_y=WLen+index_minus_x+1;
value_minus=-1;
index_one_x=ones(nxny,1);
index_one_y=WLen + eLen+nxny+(1:nxny) ;
value_one=ones(1,length(index_one_y));
A_ie2=sparse([index_minus_x';index_one_x],[index_minus_y';index_one_y'],[value_minus';value_one']);

A1=[A_ie1;A_ie2];

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

for i =1:nx*ny
    new_adj(i,:)=A_adj1(ind(i),:);

end
new_adj=sparse(new_adj);
A_adj_1=[mask.*new_adj,zeros(nxny,nx+1)];


index_2_x=repmat(1:nxny,ny,1);
index_2_x=index_2_x(:);
index_2_y=repmat(1:nx+1:(nx+1)/nx*nxny,nx,1);
index_2_y2=repmat((0:nx-1)',1,size(index_1_y,2));

index_2_y=index_2_y2+index_2_y;
index_2_y=repmat(index_2_y',ny,1);
index_2_y=index_2_y(:);

mask=sparse(index_2_x,index_2_y,1,nxny,ny*(nx+1));
A_adj2=repmat(Y.adj',nx,1); %Need a transpose for adjacency matrix 
ind=repmat(1:ny,nx+1,1);
ind=ind(:);
new_adj=zeros(nxny,(nx+1)*ny);

for i =1:(nx+1)*ny
    new_adj(:,i)=A_adj2(:,ind(i));

end

A_adj_2=[mask.*new_adj,zeros(nxny,nx+1)];


A2_adj=[A_adj_1-A_adj_2,zeros(nxny,eLen),-1*(zeros(nxny,h1Len+h2Len) ...
    +[eye(h1Len),zeros(h2Len,h2Len)])];

A3_adj=[A_adj_2-A_adj_1,zeros(nxny,eLen),-1*(zeros(nxny,h1Len+h2Len) ...
    +[eye(h1Len),zeros(h2Len,h2Len)])];


%%%% Constraint 3 %%%%
index_1_x=repmat(1:nxny,ny,1);
index_1_x=index_1_x(:);

index_1_y=repmat(1:nx+1:(nx+1)/nx*nxny,nx,1);
index_1_y2=repmat((0:nx-1)',1,size(index_1_y,2));
index_1_y=index_1_y2+index_1_y;
index_1_y=repmat(index_1_y',1,ny);
index_1_y=index_1_y(:);
mask=sparse(index_1_x,index_1_y,1,nxny,ny*(nx+1));

A_adj3=repelem(Y.adj,nx,nx);
ind=repmat(1:nx+1:(nx+1)/nx*nxny,nx,1);
ind2=repmat((0:nx-1)',1,size(ind,2));
ind=ind2+ind;
ind=ind(:);
new_adj=zeros(nxny,ny*(nx+1));

for i =1:nx*ny
    new_adj(:,ind(i))=A_adj3(:,i);
end
new_adj=sparse(new_adj);
A_adj_3=[mask.*new_adj,zeros(nxny,nx+1)];



index_2_x=repmat(1:nxny,nx,1);
index_2_x=index_2_x(:);

index_2_y=repmat(1:nx+1:(nx+1)*ny,nx,1);
index_2_y2=repmat((0:nx-1)',1,ny);
index_2_y=index_2_y2+index_2_y;
index_2_y=repmat(index_2_y',1,nx)';
index_2_y=index_2_y(:);


mask=sparse(index_2_x,index_2_y,1,nxny,ny*(nx+1));
A_adj4=repmat([X.adj',zeros(nx,1)],ny,ny); %Need a transpose for adjacency matrix 
new_adj=A_adj4;
A_adj_4=[mask.*new_adj,zeros(nxny,nx+1)];

A4_adj=[A_adj_3-A_adj_4,zeros(nxny,eLen),-1*(zeros(nxny,h1Len+h2Len) ...
    +[zeros(h1Len,h1Len),eye(h2Len)])];

A5_adj=[A_adj_4-A_adj_3,zeros(nxny,eLen),-1*(zeros(nxny,h1Len+h2Len) ...
    +[zeros(h1Len,h1Len),eye(h2Len)])];

A=[A1;A2_adj;A3_adj;A4_adj;A5_adj];
[lenb,~]= size(A);
b=sparse(lenb,1);

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
edge_cost=epsilon^p/4*(x(e1Pos)+x(e2Pos))^(1/p);
%%% TO DO

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function locCostMat = locCostComp_v2(X, Y, c, p)

% computing the localisation cost for every (i,j)
nx = size(X.xState, 1);
ny = size(Y.xState, 1);
tmpCost = c.^p / 2; % cost for being unassigned
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
