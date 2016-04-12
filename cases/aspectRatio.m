clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/mat/VEM2D/')
mult = [1,5,25];
m = numel(mult);
n = 10;
errVec = zeros(m,4);
for i = 1:m
G = cartGrid([n,mult(i)*n],[1,1]);
G = sortEdges(G);
G = computeVEM2DGeometry(G);

f = @(X) zeros(size(X,1),1);
C = -[.2,.2];
gD = @(X) -log(1./(sqrt(sum(bsxfun(@minus, X, C).^2,2)))) + 10;

boundaryEdges = find(any(G.faces.neighbors == 0,2));

bc_VEM = VEM2D_addBC([], G, boundaryEdges, 'pressure', gD);
bc_MRST = addBC([], boundaryEdges, 'pressure',gD(G.faces.centroids(boundaryEdges,:)));

gravity reset off
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.poro = ones(G.cells.num,1);
rock.perm = ones([G.cells.num,1]);

%% Initialize state
sInit = initState(G, [], 0, [0.0,1]);
S     = computeMimeticIP(G, rock, 'Verbose', true);
trans = computeTrans(G,rock);
%% Solve Laplace

sTPFA = incompTPFA(sInit, G, trans, fluid, 'bc', bc_MRST);
sMIM  = solveIncompFlow(sInit, G, S, fluid, 'bc', bc_MRST);
sVEM1 = VEM2D(G,0,1,bc_VEM, 'findCellAverages', true);
sVEM2 = VEM2D(G,0,2,bc_VEM);

faceAvg = polygonInt(G,1:G.cells.num,gD,7)./G.cells.volumes;

NMIM = sqrt(G.cells.num);
NK1 = sqrt(G.nodes.num);
NK2 = sqrt(G.nodes.num + G.faces.num  + G.cells.num);

u = [gD(G.nodes.coords); gD(G.faces.centroids); faceAvg];
errVec(i,:) = [norm(sTPFA.pressure-faceAvg)/NMIM, ...
               norm(sMIM.pressure-faceAvg)/NMIM, ...
               norm(sVEM1.nodeValues-u(1:G.nodes.num))/NK1...
               norm([sVEM2.nodeValues; sVEM2.edgeValues; sVEM2.cellMoments]-u)/NK2];

end

errVec = bsxfun(@rdivide, errVec, errVec(1,:));

x = 1:m;
plot(x,errVec(:,1), x,errVec(:,2),x,errVec(:,3),x,errVec(:,4));

legend('TPFA', 'Mimetic', 'VEM1', 'VEM2')