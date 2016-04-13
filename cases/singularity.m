clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/mat/VEM2D/')
addpath('../../pebiGridding/voronoi2D/')

%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid

xMax = 1; yMax = 1;

C = [xMax/2, yMax/2];
wellLine = {C};                % Set source center

gT = 3;

switch gT
  case 1
    
    nx = 61;
    G = cartGrid([nx,nx],[xMax,yMax]);
    G = computeGeometry(G);
    w1 = wellLine{1};
    D = pdist2(G.cells.centroids, w1);
    [~, I] = min(D, [], 1);
    G.cells.tag = false(G.cells.num,1);
    G.cells.tag(I') = true(size(I'));
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);
    
  case 2 
    
    gridSize = xMax/20;                   % Size of gridcells
    mlqtMax = 2;                            % Set number of reminement levels
    wellGridSize = 0.75/2^mlqtMax;          % Set relative well grid size
    mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                            % Size around wells to be refined 
    G = compositePebiGrid(gridSize, [xMax, yMax], 'wellLines', wellLine, ...
                         'wellGridFactor', wellGridSize, ...
                         'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes);
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);
    
  case 3
      
    gridSize = xMax/15;
    wellGridSize = 0.7/2^2;
    epsilon = gridSize*.7;
    G = pebiGrid(gridSize, [xMax, yMax], 'wellLines', wellLine,   ...
                'wellGridFactor', wellGridSize, 'wellRefinement',true);
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);
    
    otherwise
      
    error('unknown grid case')
    
end                  


%%  Set BC

gD = @(X) log(1./sqrt((X(:,1)-.5).^2 + (X(:,2)-.5).^2));
tol = 1e-6;
boundaryEdges = find(G.faces.neighbors(:,1) == 0 | G.faces.neighbors(:,2) == 0);
bc_VEM = VEM2D_addBC([], G, boundaryEdges, 'pressure', gD);
bc_MRST = addBC([], boundaryEdges, 'pressure', gD(G.faces.centroids(boundaryEdges,:)));
            
%% Set fluid and rock properties
gravity reset off 
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);

%% add Sources
Q = 1;
srcCells = find(G.cells.tag);
src = addSource([],srcCells(1),Q);

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);

%% Solve Laplace
sTPFA = incompTPFA(sInit, G, trans, fluid, 'src', src, 'bc', bc_MRST);
sMIM  = solveIncompFlow(sInit, G, S, fluid,'src', src, 'bc', bc_MRST);
sVEM1 = VEM2D(G,0,1,bc_VEM, 'src', src, 'cellAverages', true);
sVEM2 = VEM2D(G,0,2,bc_VEM, 'src', src);


%%  Plot solutions along radial line

r = .2;     % Trustregion
n = 20;     % ~Number of points

firstQuadTr = G.cells.centroids(:,1) > .5 & G.cells.centroids(:,2) > .5 ...
              & sum(bsxfun(@minus, G.cells.centroids, C).^2,2) > r^2;
dist = abs(G.cells.centroids(:,1) - G.cells.centroids(:,2));
dist = dist(firstQuadTr);
lineCells = find(firstQuadTr);
sortDist = sort(dist);
lineCells = lineCells(ismember(dist,sortDist(1:n)));

XC = G.cells.centroids(lineCells,:);
rC = sqrt(sum(bsxfun(@minus,XC,C).^2,2));
VEM1Vals = sVEM1.cellMoments(lineCells);
VEM2Vals = sVEM2.cellMoments(lineCells);

TPFAVals = sTPFA.pressure(lineCells);
MIMVals = sMIM.pressure(lineCells);

XL = repmat((.5:.001:xMax)',1,2);
rL = sqrt(sum(bsxfun(@minus, XL, C).^2,2));

subplot(1,2,1)
plot(rL,gD(XL),'-.')
hold on
plot(rC, VEM1Vals, '.', rC, VEM2Vals, 'sq', rC, TPFAVals, '+', rC, MIMVals, 'o');
yLMax = max([VEM1Vals; VEM2Vals; TPFAVals; MIMVals]);
yLMin = min([VEM1Vals; VEM2Vals; TPFAVals; MIMVals]);
yLength = 1.1*(yLMax - yLMin); yMid = (yLMin + yLMax)/2;
ylim([yMid - yLength/2, yMid + yLength/2]);
xlim([0,1]);
legend('Exact solution', '1st order VEM', '2nd order VEM', 'TPFA', 'MFD');
xlabel('Distance from center'); ylabel('Pressure')
subplot(1,2,2)
plotGrid(G);
hold on
plot(XL(:,1),XL(:,2),'linewidth', 1);
h1 = plot(XC(:,1), XC(:,2), 'o');
legend(h1, 'Datapoints');
axis equal

% 
% firstQuadTr = G.nodes.coords(:,1) > .5 & G.nodes.coords(:,2) > .5 ...
%               & sum(bsxfun(@minus, G.nodes.coords, C).^2,2) > r^2;
% dist = abs(G.nodes.coords(:,1) - G.nodes.coords(:,2));
% dist = dist(firstQuadTr);
% lineNodes = find(firstQuadTr);
% sortDist = sort(dist);
% lineNodes = lineNodes(ismember(dist,sortDist(1:n)));
% 
% XN = G.nodes.coords(lineNodes,:);
% rN = sqrt(sum(bsxfun(@minus,XN,C).^2,2));
% VEM1Vals = sVEM1.nodeValues(lineNodes);
% VEM2Vals = sVEM2.nodeValues(lineNodes);