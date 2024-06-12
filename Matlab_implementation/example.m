clc;
clear;
%% parameters pre-setting

example_id=1;% 0-3, example 1 undirected graphs;
             % 4-7,example2 undirected graphs; 
             % 8-11,example3 directed graphs.

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
 DelVal = 0.7;
 switch example_id
    case 0
        % Example 1a
        
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 3;
        Y.xState = [0,0;10,10;10,20]; %0,0;1,1.5;2,5
        Y.adj = [0,1,1;1,0,1;1,1,0];
    case 1
        % Example 1b

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; %0,0;1,1.5;2,5
        Y.adj = [0,1,0;1,0,1;0,1,0];
    case 2
        % Example 1c

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 2;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal]; %0,0;1,1.5;2,5
        Y.adj = [0,1;1,0];
            

    case 3
        % Example 1d
   
        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20]; 
        X.adj = [0,1,1;1,0,1;1,1,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,15]; 
        Y.adj = [0,1,0;1,0,1;0,1,0];

    case 4
        % Example 2a

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20]; 
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; 
        Y.adj = [0,0.3,0.4;0.3,0,0.5;0.4,0.5,0];


    case 5

        % Example 2b

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; %0,0;1,1.5;2,5
        Y.adj = [0,0.3,0;0.3,0,0.5;0,0.5,0];

   case 6

        % Example 2c

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal]; %0,0;1,1.5;2,5
        Y.adj = [0,0.3;0.3,0,];
        
      case 7

        % Example 2d

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,0.3,0.7;0.3,0,0.5;0.7,0.5,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,15]; %0,0;1,1.5;2,5
        Y.adj = [0,0.3,0;0.3,0,0.6;0,0.6,0];

      case 8
        % Example 3a

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,0;0,0,1;1,0,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; 
        Y.adj = [0,1,1;0,0,1;0,0,0];

        case 9
        % Example 3b

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,0;0,0,1;1,0,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal;10,20]; 
        Y.adj = [0,1,0;0,0,1;0,0,0];


        case 10
        % Example 3c

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
          nx = 3;
        X.xState = [0,0;10,10;10,20];
        X.adj = [0,1,0;0,0,1;1,0,0];
        ny = 3;
        Y.xState = [0,0-DelVal*delVal;10,10-DelVal*delVal]; 
        Y.adj = [0,1;0,0];
        

        case 11
        % Example 3d

        %nx: number of nodes in graph X
        %ny: number of nodes in graph Y
        %adj_x: adjacency matrix of x
        %adj_y: adjacency matrix of y
        nx = 3;
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
    [dxy,loc_cost,fa_cost,miss_cost,edge_cost]=LPgraphGOSPA_directed(X,Y,c,p,epsilon);

end
 %%
disp('Graph GOSPA Metric Value')
disp(dxy)
disp(‘Node attribute cost’)
disp(loc_cost)
disp('False node cost')
disp(fa_cost)
disp('Missed node cost')
disp(miss_cost)
disp('Edge Mismatch Cost')
disp(edge_cost)
%%
if (example_id<=7)
    Gx=graph(X.adj);
    Gy=graph(Y.adj);
    figure(1);
    clf
    hold on
    plot(Gx,'XData',X.xState(:,1),'YData',X.xState(:,2),'Marker','x','MarkerSize',12,'NodeColor','b',LineWidth=2)
    plot(Gy,'XData',Y.xState(:,1),'YData',Y.xState(:,2),'Marker','o','MarkerSize',12,'NodeColor','r',LineWidth=2)
    
    hold off
    grid on
    xlabel('Position X')
    ylabel('Position Y')
    title('Graphs (X in blue, Y in red)')
else
    Gx=digraph(X.adj);
    Gy=digraph(Y.adj);
    figure(1);
    clf
    hold on
    plot(Gx,'XData',X.xState(:,1),'YData',X.xState(:,2),'Marker','x','MarkerSize',12,'NodeColor','b',LineWidth=2)
    plot(Gy,'XData',Y.xState(:,1),'YData',Y.xState(:,2),'Marker','o','MarkerSize',12,'NodeColor','r',LineWidth=2)
    
    hold off
    grid on
    xlabel('Position X')
    ylabel('Position Y')
    title('Graphs (X in blue, Y in red)')
end




 