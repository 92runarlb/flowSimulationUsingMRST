clc; clear all; close all

mrstModule add mimetic
mrstModule add streamlines
addpath('../../vem/mat/VEM2D/stable/')
addpath('../../pebi/voronoi2D')
run('../../project-mechanics-fractures/mystartup.m')

xmax = 2;
ymax = 1;

nx   = 80;
ny   = 40;

G = cartGrid([nx,ny],[xmax,ymax]);
makeSkew = @(c) c(:,1) + .45*(1-(c(:,1)-1).^2).*(1-c(:,2));
G.nodes.coords(:,1) = makeSkew(G.nodes.coords);
G = computeGeometry(G);

G = sortEdges(G);
G = mrstGridWithFullMappings(G);
G = VEM2D_makeInternalBoundary(G, find(G.faces.tag));
G = computeVEM2DGeometry(G);


%%  Set BC
tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
bot = abs(G.faces.centroids(boundaryEdges,2)) < tol;
top = abs(G.faces.centroids(boundaryEdges,2) - ymax) < tol;
             
neuman = ~bot& ~top;

bc_MRST = addBC([], boundaryEdges(bot), 'pressure', 0);
bc_MRST = addBC(bc_MRST, boundaryEdges(top), 'pressure', 100);
bc_MRST = addBC(bc_MRST, boundaryEdges(neuman), 'flux', 0);

bc_VEM = VEM_addBC(G,[], boundaryEdges(bot), 'pressure', 0);
bc_VEM = VEM_addBC(G,bc_VEM, boundaryEdges(top), 'pressure', 100);
bc_VEM = VEM_addBC(G,bc_VEM, boundaryEdges(neuman), 'flux', 0);
%% Set fluid and rock properties
gravity reset off 

fluid = initSingleFluid('mu' , 2    , ...
                        'rho', 1);

rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);



%% Initialize state
sInit = initState(G, [], 0);
S     = computeMimeticIP(G, rock);
trans = computeTrans(G,rock);
%% Solve Laplace

sTPFA = incompTPFA(sInit, G, trans, fluid, 'bc', bc_MRST);
sMIM  = incompMimetic(sInit, G, S, fluid,'bc', bc_MRST);
sVEM2 = VEM2D_v3(G,0,2,bc_VEM);


plotCellData(G,sTPFA.pressure,'edgecolor','none')
colormap(jet)
colorbar()
figure
plotCellData(G,sMIM.pressure,'edgecolor','none')
colormap(jet)
colorbar()
figure()
plotCellData(G,sVEM2.cellMoments,'edgecolor','none')
colormap(jet)
colorbar()