clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/mat/VEM2D/')

% this script performs a simple 2 phase(oil and water) flow simulation in
% Using a composite grid. The case consist of two wells (one injector and
% one producer) and a linear no-flow fault seperating them.

%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid



load('../data/grids/pressureDropGridCartFine.mat');

%%  Set BC
tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);

bc_VEM = VEM2D_addBC([], G, boundaryEdges, 'pressure', 0);
bc_MRST = addBC([], boundaryEdges, 'pressure', 0);

            
%% Set fluid and rock properties
gravity reset off 

fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);


%% add Sources
Q = 100;
srcCells = find( ...
                 abs(G.cells.centroids(:,1)-G.cells.centroids(G.cells.tag,1)) < 1/26 ...
             &   abs(G.cells.centroids(:,2)-G.cells.centroids(G.cells.tag,2)) < 1/26);
Qscaled = Q/numel(srcCells);
src = addSource([],srcCells,Qscaled*ones(numel(srcCells),1));

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);
%% Solve Laplace for coarse Cartesian

sTPFAcart = incompTPFA(sInit, G, trans, fluid, 'src', src, 'bc', bc_MRST);
sMIMcart  = solveIncompFlow(sInit, G, S, fluid,'src', src, 'bc', bc_MRST);
sVEM1cart = VEM2D(G,0,1,bc_VEM, 'src', src, 'cellAverages', true);
sVEM2cart = VEM2D(G,0,2,bc_VEM, 'src', src);


%% add Sources
srcCellVol = 0.001607510288066;
r = sqrt(srcCellVol/pi);
srcCells = find( (G.cells.centroids(:,1)-G.cells.centroids(G.cells.tag,1)).^2 ...
               + (G.cells.centroids(:,2)-G.cells.centroids(G.cells.tag,2)).^2 <= r);
QScaled = Q/numel(srcCells);
src = addSource([],srcCells,Qscaled*ones(numel(srcCells),1));

%% Solve Laplace for composite Pebi

sTPFAcomp = incompTPFA(sInit, G, trans, fluid, 'src', src, 'bc', bc_MRST);
sMIMcomp  = solveIncompFlow(sInit, G, S, fluid,'src', src, 'bc', bc_MRST);
sVEM1comp = VEM2D(G,0,1,bc_VEM, 'src', src, 'cellAverages', true);
sVEM2comp = VEM2D(G,0,2,bc_VEM, 'src', src);

%% add Sources
srcCellVol = 9.564213523552989e-04;
r = sqrt(srcCellVol/pi);
srcCells = find( (G.cells.centroids(:,1)-G.cells.centroids(G.cells.tag,1)).^2 ...
               + (G.cells.centroids(:,2)-G.cells.centroids(G.cells.tag,2)).^2 <= r);
QScaled = Q/numel(srcCells);
src = addSource([],srcCells,Qscaled*ones(numel(srcCells),1));

%% Solve Laplace for composite Pebi

sTPFApebi = incompTPFA(sInit, G, trans, fluid, 'src', src, 'bc', bc_MRST);
sMIMpebi  = solveIncompFlow(sInit, G, S, fluid,'src', src, 'bc', bc_MRST);
sVEM1pebi = VEM2D(G,0,1,bc_VEM, 'src', src, 'cellAverages', true);
sVEM2pebi = VEM2D(G,0,2,bc_VEM, 'src', src);

%% Save solutions
Gr = G;
save('../data/pressureDropExactConstSzSrc.mat', 'Gr', ...
     'sTPFAcart', 'sMIMcart', 'sVEM1cart', 'sVEM2cart', ...
     'sTPFAcomp', 'sMIMcomp', 'sVEM1comp', 'sVEM2comp', ...
     'sTPFApebi', 'sMIMpebi', 'sVEM1pebi', 'sVEM2pebi');