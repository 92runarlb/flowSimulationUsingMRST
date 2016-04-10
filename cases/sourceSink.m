clc; clear all; close all;

mrstModule add mimetic
% addpath('../');
% addpath('../../');
% addpath('~/Documents/master/pebiGridding/voronoi2D/')
% run('../../startup.m');


run('~/NTNU/5/master/project-mechanics-fractures/mystartup.m')
addpath('~/NTNU/5/master/flowSimulationUsingMRST/pebi/')
addpath('~/NTNU/5/master/pebiGridding/voronoi2D/')
addpath('~/NTNU/5/master/vem/vem/mat/')
addpath('~/NTNU/5/master/vem/vem/mat/VEM2D/')
addpath('~/NTNU/5/master/vem/vem/mat/VEM2D/stable/')


n = 5;

% x     = linspace(0.2, 0.8, 10);
% y     = 0.8 - 0.5*x - 0.05* sin(6*pi*x);
% fault = {[x' , y']};
% G = compositePebiGrid(1/20, [1, 1], ...
%                        'faultLines', fault, 'faultGridFactor', 1/sqrt(2));

% x = linspace(.2,.8,10);
% y = 1-x;
% fault = {[x' , y']};
% G = pebiGrid(1/10, [1, 1], ...
%                   'faultLines', fault);

xMax = 1; yMax = 1;
G = cartGrid([n,100*n],[xMax, yMax]);
G = twister(G);

G = computeGeometry(G);

% faultFaces = 1:G.faces.num;
% faultFaces = faultFaces(G.faces.centroids(:,1) == xMax/2);
% faultFaces = 1:G.faces.num;
% faultFaces = faultFaces(G.faces.tag);
G = mrstGridWithFullMappings(G);
G = sortEdges(G);

% G = VEM2D_makeInternalBoundary(G, faultFaces);
f = zeros(G.cells.num,1);
G = computeVEM2DGeometry(G,f,1,1);
 
% sourceCoords = [.2,.2];
% source = sum(bsxfun(@minus, G.cells.centroids, sourceCoords).^2,2);
% source = find(source == min(source));
% source = source(1);
% 
% sinkCoords = [.8,.8];
% sink = sum(bsxfun(@minus, G.cells.centroids, sinkCoords).^2,2);
% sink = find(sink == min(sink));
% sink = sink(1);
% 
% Q = 10;
% src = addSource([], source, Q);
% src = addSource(src, sink, -Q);

tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
isExternal = abs(G.faces.centroids(boundaryEdges,1)) < tol | ...
             abs(G.faces.centroids(boundaryEdges,1) - xMax) < tol | ...
             abs(G.faces.centroids(boundaryEdges,2)) < tol |...
             abs(G.faces.centroids(boundaryEdges,2) - yMax) < tol;
isInternal = ~isExternal;

gD = @(X) X(:,1).^2 -X(:,2).^2 + yMax^2;

bc_VEM = VEM_addBC(G, [], boundaryEdges(isExternal), 'pressure', gD);
bc_VEM = VEM_addBC(G, bc_VEM, boundaryEdges(isInternal), 'flux', gD);
bc_MRST = addBC([], boundaryEdges(isExternal), 'pressure', gD(G.faces.centroids(boundaryEdges(isExternal),:)));
bc_MRST = addBC(bc_MRST, boundaryEdges(isInternal), 'flux', gD(G.faces.centroids(boundaryEdges(isExternal),:)));

            
%% Set fluid and rock properties
gravity reset off 

fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);


%% add Sources
pv = sum(poreVolume(G, rock));

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);
%% Solve Laplace

sTPFA = incompTPFA(sInit, G, trans, fluid, 'bc', bc_MRST);
sMIM  = solveIncompFlow(sInit, G, S, fluid, 'bc', bc_MRST);
sVEM1 = VEM2D_v3(G,0,1,bc_VEM, 'findCellAverages', true);
sVEM2 = VEM2D_v3(G,0,2,bc_VEM);

G.cells.centroids(max(sTPFA.pressure) == sTPFA.pressure,:)
G.cells.centroids(max(sMIM.pressure) == sMIM.pressure,:)
G.nodes.coords(max(sVEM1.nodeValues) == sVEM1.nodeValues,:)
G.nodes.coords(max(sVEM2.nodeValues) == sVEM2.nodeValues,:)

subplot(2,2,1)
plotCellData(G,sTPFA.pressure,'edgeColor', 'none');
title('TPFA');
colorbar;
axis equal
subplot(2,2,2)
plotCellData(G,sMIM.pressure,'edgeColor', 'none');
title('Mimetic');
colorbar;
axis equal
subplot(2,2,3)
plotCellData(G,sVEM1.cellMoments,'edgeColor', 'none');
title('VEM 1st order');
colorbar;
axis equal
subplot(2,2,4)
plotCellData(G,sVEM2.cellMoments,'edgeColor', 'none');
title('VEM 2nd order');
colorbar;
axis equal