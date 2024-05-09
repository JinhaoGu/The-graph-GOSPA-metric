clc;
clear;
%% parameters pre-setting

example_id=11;%0-3, example 1; 4-7,example2; 8-11,example3.

c=3; 
p=1; 
epsilon=1;


%%
% structures for X and Y. X will define the grounth truth set of
% graphs. Y will define the estimated sets of graphs.
% Meaning of the struct fields is explained in LPGraphMetric.m

X = struct('xState', [], 'adj', []);
Y = struct('xState', [], 'adj', []);

%Parameters to define sets of graphs
 delVal = 1;
 DelVal = 0.2;
 switch example_id
    case 0
        % Example 1a
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 3;
        Y.xState = [0,0;10,10;10,20]; %0,0;1,1.5;2,5
        Y.adj = [0,1,1;1,0,1;1,1,0];
    case 1
        % Example 1b
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; %0,0;1,1.5;2,5
        Y.adj = [0,1,0;1,0,1;0,1,0];
    case 2
        % Example 1c
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 2;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal]; %0,0;1,1.5;2,5
        Y.adj = [0,1;1,0];
            

    case 3
        % Example 1d
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20]; 
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,15]; 
        Y.adj = [0,1,0;1,0,1;0,1,0];

    case 4
        % Example 2a
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20]; 
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; 
        Y.adj = [0,0.3,0.4;0.3,0,0.5;0.4,0.5,0];


    case 5

        % Example 2b
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; %0,0;1,1.5;2,5
        Y.adj = [0,0.3,0;0.3,0,0.5;0,0.5,0];

   case 6

        % Example 2c
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal]; %0,0;1,1.5;2,5
        Y.adj = [0,0.3;0.3,0,];
        
      case 7

        % Example 2d
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,15]; %0,0;1,1.5;2,5
        Y.adj = [0,0.3,0;0.3,0,0.6;0,0.6,0];

      case 8
        % Example 3a
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,0;0,0,1;1,0,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; 
        Y.adj = [0,1,1;0,0,1;0,0,0];

        case 9
        % Example 3b
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,0;0,0,1;1,0,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; 
        Y.adj = [0,1,0;0,0,1;0,0,0];


        case 10
        % Example 3c
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,0;0,0,1;1,0,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal]; 
        Y.adj = [0,1;0,0];
        

        case 11
        % Example 3d
        %stDim Target state dimension
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        stDim = 1; nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,0;0,0,1;1,0,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,15]; 
        Y.adj = [0,1,0;0,0,1;0,0,0];
 end
%%
if (example_id<=7)
    [dxy,loc_cost,fa_cost,miss_cost,edge_cost]=LPgraphGOSPA(X,Y,c,p,epsilon);

else
    [dxy,loc_cost,fa_cost,miss_cost,edge_cost]=extended_LPgraphGOSPA(X,Y,c,p,epsilon);

end
 %%
disp('Metric Value')
disp(dxy)
disp('Localisation Cost')
disp(loc_cost)
disp('False Detection Cost')
disp(fa_cost)
disp('Miss Detection Cost')
disp(miss_cost)
disp('Edge Mismatch Cost')
disp(edge_cost)



 