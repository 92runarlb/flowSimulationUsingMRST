clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/vem/mat/VEM2D/stable/')
addpath('../../vem/vem/mat/')
addpath('../../vem/vem/mat/VEM2D/')
run('../../project-mechanics-fractures/mystartup.m')

% G = pebiGrid(1/10, [1,1]);
G = unitSquare(10,10);
G = sortEdges(G);
G = computeVEM2DGeometry(G);


f = @(X) zeros(size(X,1),1);
C = -[.2,.2];
gD = @(X) -log(1./(sqrt(sum(bsxfun(@minus, X, C).^2,2)))) + 1000;

% f = @(X) -2*ones(size(X,1),1);
% gD = @(X) X(:,1).^2;

% f = @(X) zeros(size(X,1),1);
% gD = @(X) X(:,1);

% f = @(X) -6*X(:,2);
% gD = @(X) X(:,1) + X(:,2).^3;



boundaryEdges = find(any(G.faces.neighbors == 0,2));

bc_VEM = VEM_addBC(G, [], boundaryEdges, 'pressure', gD);
% bc_MRST = addBC([], boundaryEdges, 'pressure',gD(G.faces.centroids(boundaryEdges,:)));

% gravity reset off 
% 
% fluid = initSingleFluid('mu' , 1, 'rho', 1);
% rock.poro = ones(G.cells.num,1);
% rock.perm = ones([G.cells.num,1]);
% 
% %% Initialize state
% sInit = initState(G, [], 0, [0.0,1]);
% S     = computeMimeticIP(G, rock, 'Verbose', true);
% trans = computeTrans(G,rock);
% %% Solve Laplace
% 
% sTPFA = incompTPFA(sInit, G, trans, fluid, 'bc', bc_MRST);
% sMIM  = solveIncompFlow(sInit, G, S, fluid, 'bc', bc_MRST);
sVEM1 = VEM2D_v3(G,0,1,bc_VEM, 'findCellAverages', true);
sVEM2 = VEM2D_v3(G,0,2,bc_VEM);

% subplot(2,2,1)
% plotCellData(G,sTPFA.pressure);
% title('TPFA');
% colorbar;
% axis equal
% subplot(2,2,2)
% plotCellData(G,sMIM.pressure);
% title('Mimetic');
% colorbar;
% axis equal
% subplot(2,2,3)
% plotCellData(G,sVEM1.cellMoments);
% title('VEM 1st order');
% colorbar;
% axis equal
% subplot(2,2,4)
% plotCellData(G,sVEM2.cellMoments);
% title('VEM 2nd order');
% colorbar;
% axis equal