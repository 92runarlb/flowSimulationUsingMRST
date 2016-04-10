clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/vem/mat/VEM2D/stable/')
run('../../project-mechanics-fractures/mystartup.m')

%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid

fileNames = {'../data/grids/pressureDropGridCart.mat', ...
            '../data/grids/pressureDropGridCompPebi.mat', ...
            '../data/grids/pressureDropGridPebi.mat'};
refName   = '../data/pressureDropExactDelta';
gT = 2;

load(fileNames{gT});
load(refName);



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
src = addSource([],srcCells(1),Q);

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);

%% Solve Laplace
sTPFA = incompTPFA(sInit, G, trans, fluid, 'src', src, 'bc', bc_MRST);
sMIM  = solveIncompFlow(sInit, G, S, fluid,'src', src, 'bc', bc_MRST);
sVEM1 = VEM2D_v3(G,0,1,bc_VEM, 'src', src, 'findCellAverages', true);
sVEM2 = VEM2D_v3(G,0,2,bc_VEM, 'src', src);


tit = {'TPFA', 'Mimetic', 'VEM 1st order', 'VEM 2nd order'};
plotPressDrop(G,sTPFA.pressure, sMIM.pressure, sVEM1.cellMoments, sVEM2.cellMoments, tit);
for i = 1:4
    subplot(2,2,i)
    colorbar;
end


%% TRUST REGION
d = G.cells.diameters(srcCells);
cx = G.cells.centroids(srcCells,1);
cy = G.cells.centroids(srcCells,2);

tC = (G.cells.centroids(:,1)-cx).^2 +(G.cells.centroids(:,2)-cy).^2>d/2;  

%% PLOT DIFFERENCES
% sVEM2.cellMoments(G.cells.tag) = 1;
% sVEM1.cellMoments(G.cells.tag) = 1;
% sMIM.pressure(G.cells.tag) = 1;
% sTPFA.pressure(G.cells.tag) = 1;

VEM2_MIM = (abs(sVEM2.cellMoments-sMIM.pressure))./sVEM2.cellMoments;
VEM2_TPFA = (abs(sVEM2.cellMoments-sTPFA.pressure))./sVEM2.cellMoments;
VEM1_MIM = (abs(sVEM1.cellMoments-sMIM.pressure))./sVEM1.cellMoments;
VEM1_TPFA = (abs(sVEM1.cellMoments-sTPFA.pressure))./sVEM1.cellMoments;

figure;

cScale = [0,min(max([VEM2_MIM,VEM2_TPFA,VEM1_MIM,VEM1_TPFA],[],1))];

% Create colorbar
%hp4 = get(subplot(2,2,4),'Position');
%h = colorbar('Position', [hp4(1)+hp4(3),  hp4(2),  0.03,  hp4(2)+hp4(3)*2.1]);
% cScale = [0,max([VEM2_MIM;VEM2_TPFA;VEM1_MIM;VEM1_TPFA])];
% set(h, 'ylim', cScale)

tit = {'VEM1 - MIM', 'VEM1 - TPFA', 'VEM2 - MIM', 'VEM2-TPFA'};
plotPressDrop(G,VEM1_MIM, VEM1_TPFA,VEM2_MIM, VEM2_TPFA,tit)

for i = 1:4
     subplot(2,2,i)
     caxis manual
     caxis(cScale)
    colorbar()
end





%% PLOT TRUE ERROR
sVEM2rFun = scatteredInterpolant(Gr.cells.centroids(:,1), Gr.cells.centroids(:,2), ...
                 sVEM2r.cellMoments);
sVEM2rInt = sVEM2rFun(G.cells.centroids(:,1), G.cells.centroids(:,2));

sVEM1rFun = scatteredInterpolant(Gr.cells.centroids(:,1), Gr.cells.centroids(:,2), ...
                 sVEM1r.cellMoments);
sVEM1rInt = sVEM1rFun(G.cells.centroids(:,1), G.cells.centroids(:,2));

sTPFArFun = scatteredInterpolant(Gr.cells.centroids(:,1), Gr.cells.centroids(:,2), ...
                 sTPFAr.pressure);
sTPFArInt = sTPFArFun(G.cells.centroids(:,1), G.cells.centroids(:,2));

sMIMrFun = scatteredInterpolant(Gr.cells.centroids(:,1), Gr.cells.centroids(:,2), ...
                 sMIMr.pressure);
sMIMrInt = sMIMrFun(G.cells.centroids(:,1), G.cells.centroids(:,2));


VEM1_err = (abs(sVEM1.cellMoments-sVEM1rInt));
VEM2_err = (abs(sVEM2.cellMoments-sVEM2rInt));
MIM_err = (abs(sMIM.pressure-sMIMrInt));
TPFA_err = (abs(sTPFA.pressure-sTPFArInt));


cScale = max([TPFA_err(tC),MIM_err(tC),VEM1_err(tC),VEM2_err(tC)],[],1);

figure()
tit = {'TPFA', 'Mimetic', 'VEM 1st order', 'VEM 2nd order'};
plotPressDrop(G,TPFA_err, MIM_err,VEM1_err, VEM2_err,tit)


for i = 1:4
     subplot(2,2,i)
      caxis manual
      caxis([0,cScale(i)])
    colorbar()
end

%% TRASH

% switch gT
%     case 1
%         sMIMr = sMIMcart;
%         sTPFAr = sTPFAcart;
%         sVEM1r = sVEM1cart;
%         sVEM2r = sVEM2cart;
%     case 2
%         sMIMr = sMIMcomp;
%         sTPFAr = sTPFAcomp;
%         sVEM1r = sVEM1comp;
%         sVEM2r = sVEM2comp;
%     case 3
%         sMIMr = sMIMpebi;
%         sTPFAr = sTPFApebi;
%         sVEM1r = sVEM1pebi;
%         sVEM2r = sVEM2pebi;
% end


