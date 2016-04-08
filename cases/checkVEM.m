clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/vem/mat/VEM2D/stable/')
addpath('../../vem/vem/mat/')
run('../../project-mechanics-fractures/mystartup.m')

% G = pebiGrid(1/10, [1,1]);
G = unitSquare(10,10);
G = sortEdges(G);
G = computeVEM2DGeometry(G);

% f = @(X) sin(X(:,1));
% gD = @(X) sin(X(:,1));
f = @(X) -2*ones(size(X,1),1);
gD = @(X) X(:,1).^2;

boundaryEdges = find(any(G.faces.neighbors == 0,2));

bc = VEM_addBC(G, [], boundaryEdges, 'pressure', f);

sol1 = VEM2D_v3(G,f,1,bc);
sol1 = cellAverages(G,sol1);
U1 = sol1.nodeValues;
figure;
plotCellData(G,sol1.cellAverages);

sol2 = VEM2D_v3(G,f,2,bc);
U2 = [sol2.nodeValues; sol2.edgeValues; sol2.cellMoments];
figure;
plotCellData(G,sol2.cellMoments);

u = [gD(G.nodes.coords); gD(G.faces.centroids); polygonInt_v2(G, 1:G.cells.num, gD, 7)];

err1 = norm(U1-u(1:G.nodes.num),2)

err2 = norm(U2 - u,2)