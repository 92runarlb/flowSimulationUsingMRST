clc; clear all; close all

mrstModule add mimetic
mrstModule add streamlines
addpath('../../vem/mat/VEM2D/')

xmax = 1;
ymax = 1;
nx   = 100;
ny   = 100;

G = cartGrid([nx,ny],[xmax,ymax]);

G.nodes.coords = twister(G.nodes.coords);
G = sortEdges(G);
G = computeVEM2DGeometry(G);

%%  Set BC
tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
left = abs(G.faces.centroids(boundaryEdges,1)) < tol;
right= abs(G.faces.centroids(boundaryEdges,1) - xmax) < tol;
             
neuman = ~left& ~right;

bc_MRST = addBC([], boundaryEdges(left), 'pressure', 0);
bc_MRST = addBC(bc_MRST, boundaryEdges(right), 'pressure', 100);
bc_MRST = addBC(bc_MRST, boundaryEdges(neuman), 'flux', 0);

bc_VEM = VEM2D_addBC([], G, boundaryEdges(left), 'pressure', 0);
bc_VEM = VEM2D_addBC(bc_VEM, G, boundaryEdges(right), 'pressure', 100);
bc_VEM = VEM2D_addBC(bc_VEM, G, boundaryEdges(neuman), 'flux', 0);


%% Set fluid and rock properties
gravity reset off 

fluid = initSingleFluid('mu' , 2    , ...
                        'rho', 1);

rock.poro = ones(G.cells.num,1);

rock.perm = ones([G.cells.num,1]);



%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);
%% Solve Laplace

%% TPFA
sTPFA = incompTPFA(sInit, G, trans, fluid, 'bc',bc_MRST);
% seed = (nx:nx-1:nx*ny).';
% SfTPFA = pollock(G, sTPFA, seed,'substeps', 1);
% SbTPFA = pollock(G, sTPFA, seed,'substeps', 1,'reverse', true);

%% MIMETIC
sMIM  = incompMimetic(sInit, G, S, fluid,'bc',bc_MRST);
% seed = (nx:nx-1:nx*ny).';
% SfMIM = pollock(G, sMIM, seed,'substeps', 1);
% SbMIM = pollock(G, sMIM, seed,'substeps', 1,'reverse', true);
%% VEM1
sVEM1 = VEM2D(G,0,1,bc_VEM,'cellAverages',true);
%% VEM2
sVEM2 = VEM2D(G,0,2,bc_VEM);
% [sVEM2] = fluxApprox(G,sMIM, rock,fluid, bc_VEM);
% SfVEM2 = pollock(G, sVEM2, seed,'substeps', 1);
% SbVEM2 = pollock(G, sVEM2, seed,'substeps', 1,'reverse', true);
%% Plotting
%% TPFA
subplot(2,2,1)
hold on
plotCellData(G,sTPFA.pressure,'edgecolor','none');
colormap('jet')
% hf=streamline(SfTPFA);
% hb=streamline(SbTPFA);
% set([hf;hb],'Color','k');

title('TPFA');
colorbar;
%axis equal
%% MIMETIC
subplot(2,2,2)
plotCellData(G,sMIM.pressure,'edgecolor','none');
hold on
% hf=streamline(SfMIM);
% hb=streamline(SbMIM);
% set([hf;hb],'Color','k');
title('Mimetic');
colorbar;
%axis equal
%% VEM1
subplot(2,2,3)
plotCellData(G,sVEM1.cellMoments,'edgecolor','none');
title('VEM 1st order');
colorbar;
%axis equal
%% VEM 2
subplot(2,2,4)
plotCellData(G,sVEM2.cellMoments,'edgecolor','none');
%hf=streamline(SfVEM2);
%hb=streamline(SbVEM2);
%set([hf;hb],'Color','k');
%title('VEM 2nd order');
colorbar;
%axis equal
