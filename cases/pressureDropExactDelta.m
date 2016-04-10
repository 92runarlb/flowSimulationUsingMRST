clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/vem/mat/VEM2D/stable/')
run('../../project-mechanics-fractures/mystartup.m')

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

bc_VEM = VEM_addBC(G, [], boundaryEdges, 'pressure', 0);
bc_MRST = addBC([], boundaryEdges, 'pressure', 0);

            
%% Set fluid and rock properties
gravity reset off 

fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);


%% add Sources
Q = 100;
srcCells = find(G.cells.tag);
src = addSource([],srcCells,Q);

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);
%% Solve Laplace for coarse Cartesian

sTPFAr = incompTPFA(sInit, G, trans, fluid, 'src', src, 'bc', bc_MRST);
sMIMr  = solveIncompFlow(sInit, G, S, fluid,'src', src, 'bc', bc_MRST);
sVEM1r = VEM2D_v3(G,0,1,bc_VEM, 'src', src, 'findCellAverages', true);
sVEM2r = VEM2D_v3(G,0,2,bc_VEM, 'src', src);

%% Save solutions
Gr = G;
save('../data/pressureDropExactDelta.mat', 'Gr', ...
     'sTPFAr', 'sMIMr', 'sVEM1r', 'sVEM2r');