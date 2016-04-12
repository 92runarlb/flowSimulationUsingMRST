clc;

clear; close all; %addpath ../../voronoi2D/

mrstModule add mimetic
addpath('../../vem/mat/VEM2D/')
addpath('../../pebiGridding/voronoi2D/')

% this script performs a simple 2 phase(oil and water) flow simulation in
% Using a composite grid. The case consist of two wells (one injector and
% one producer) and a linear no-flow fault seperating them.

%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid

gT = 2;


%% Grid parameters
% shared parameters
xmax = 20;                              % Set grid dimentions
ymax = 20;

faultLine = {[16, 5; 3,10.1]};
wellLine = {[5,5], [15.0,15.0]};                % Set source center
switch gT
  case 1
    nx = 10;
    G = cartGrid([nx,nx],[xmax,ymax]);
    G = computeGeometry(G);
    % Find wells and fault, OBS!!! THIS IS VERY BUGGY!!
    w1 = wellLine{1};
    w2 = wellLine{2};
    w = [w1;w2];
    D = pdist2(G.cells.centroids, w);
    [~, I] = min(D, [], 1);
    G.cells.tag = false(G.cells.num,1);
    G.cells.tag(I') = true(size(I'));
    % Find faults
    n = nx*100;
    fault = faultLine{1};
    fault(:,2) = fault(:,2);    
    dx = fault(2,1) - fault(1,1);
    dy = fault(2,2) - fault(1,2);
    
    vx = fault(2,:) - fault(1,:);
    spacing = linspace(0,1,n)';
    liney = fault(1,2) +ceil((dy*spacing- mod(dy*spacing, 0.5/nx))*nx)/nx;
    linex = fault(1,1) + dx*spacing- mod(dx*spacing, 1/nx);
    line = [linex,liney];
    [line, ~, IC] = uniquetol(line,1/nx*1e-6, 'ByRows', true);
    line = line(IC,:);
    line = unique(line,'rows','stable');
    line = 0.5*(line(1:end-1,:)+line(2:end,:));
    
    D = pdist2(G.faces.centroids, line);
    [~, I] = min(D, [], 1);
    G.faces.tag = false(G.faces.num,1);
    G.faces.tag(I') = true(size(I'));
    
  case 2 
    gridSize = xmax*1/10;                   % Size of gridcells

    mlqtMax = 2;                            % Set number of reminement levels
    wellGridSize = 0.75/2^mlqtMax;          % Set relative well grid size
    mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                            % Size around wells to be refined 
    G = compositePebiGrid(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                         'wellGridFactor', wellGridSize, ...
                         'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes,...
                         'faultLines', faultLine);
    G = computeGeometry(G);
  case 3
    gridSize = xmax*1/7.5;   
    wellGridSize = 0.75/2^2;
    epsilon = gridSize*2;
    G = pebiGrid(gridSize, [xmax, ymax], 'wellLines', wellLine,   ...
                'wellGridFactor', wellGridSize, 'wellRefinement',true,...
                'epsilon',epsilon,'faultLines', faultLine);
    G = computeGeometry(G);
  otherwise
    error('unknown grid case')
end                  

G = sortEdges(G);
G = mrstGridWithFullMappings(G);
G = VEM2D_makeInternalBoundary(G, find(G.faces.tag));
plotGrid(G);
figure;
G = computeVEM2DGeometry(G);

%%  Set BC
tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
isExternal = abs(G.faces.centroids(boundaryEdges,1)) < tol | ...
             abs(G.faces.centroids(boundaryEdges,1) - xmax) < tol | ...
             abs(G.faces.centroids(boundaryEdges,2)) < tol |...
             abs(G.faces.centroids(boundaryEdges,2) - ymax) < tol;
isInternal = ~isExternal;

bc_VEM = VEM2D_addBC([], G, boundaryEdges(isExternal), 'pressure', 0);
bc_VEM = VEM2D_addBC(bc_VEM, G, boundaryEdges(isInternal), 'flux', 0);

bc_MRST = addBC([], boundaryEdges(isExternal), 'pressure', 0);
bc_MRST = addBC(bc_MRST, boundaryEdges(isInternal), 'flux', 0);

            
%% Set fluid and rock properties
gravity reset off 

fluid = initSingleFluid('mu' , 1, 'rho', 1);

rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);


%% add Sources
srcCells = find(G.cells.tag);
pv = sum(poreVolume(G, rock));
src = addSource([],srcCells(1),5);
src = addSource(src, srcCells(2), -5);

% eps = .019;
% srcCentroid = G.cells.centroids(srcCells(1));
% fsrc = @(X) src.rate(1)/(2*pi*eps)*exp(-sum((bsxfun(@minus, X, srcCentroid).^2),2)/(2*eps));
% sinkCentroid = G.cells.centroids(srcCells(2));
% fsink = @(X) src.rate(2)/(2*pi*eps)*exp(-sum((bsxfun(@minus, X, sinkCentroid).^2),2)/(2*eps));
% f = @(X) fsrc(X) + fsink(X);

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);
%% Solve Laplace

sTPFA = incompTPFA(sInit, G, trans, fluid, 'src', src, 'bc', bc_MRST);
sMIM  = solveIncompFlow(sInit, G, S, fluid,'src', src, 'bc', bc_MRST);
sVEM1 = VEM2D(G,0,1,bc_VEM, 'src', src, 'cellAverages', true);
sVEM2 = VEM2D(G,0,2,bc_VEM, 'src', src);

subplot(2,2,1)
plotCellData(G,sTPFA.pressure);
title('TPFA');
colorbar;
axis equal
subplot(2,2,2)
plotCellData(G,sMIM.pressure);
title('Mimetic');
colorbar;
axis equal
subplot(2,2,3)
plotCellData(G,sVEM1.cellMoments);
title('VEM 1st order');
colorbar;
axis equal
subplot(2,2,4)
plotCellData(G,sVEM2.cellMoments);
title('VEM 2nd order');
colorbar;
axis equal