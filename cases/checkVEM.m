clc; clear; close all;

mrstModule add mimetic
addpath('../../vem/vem/mat/VEM2D/stable/')
addpath('../../vem/vem/mat/')
addpath('../../vem/vem/mat/VEM2D/')
run('../../project-mechanics-fractures/mystartup.m')

f = @(X) zeros(size(X,1),1);
C = -[.2,.2];
gD = @(X) -log(1./(sqrt(sum(bsxfun(@minus, X, C).^2,2))));


nVec = [10,20,40,60];
nIt = numel(nVec);
errVec = zeros(nIt, 3);

for i = 1:nIt
    
    n = nVec(i);
    G = cartGrid([n,n], [1,1]);
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);

    boundaryEdges = find(any(G.faces.neighbors == 0,2));

    bc_VEM = VEM_addBC(G, [], boundaryEdges, 'pressure', gD);

    sVEM1 = VEM2D_v3(G,0,1,bc_VEM);
    sVEM2 = VEM2D_v3(G,0,2,bc_VEM);

	u = [gD(G.nodes.coords); gD(G.faces.centroids); polygonInt_v2(G,1:G.cells.num,gD,7)./G.cells.volumes];

    err1 = sVEM1.nodeValues-u(1:G.nodes.num);
    err2 = [sVEM2.nodeValues; sVEM2.edgeValues; sVEM2.cellMoments] - u;
    
    h = mean(G.cells.diameters); 
    
    errVec(i,:) = [h, h*norm(err1, 'inf'), h*norm(err2, 'inf')];
    
end

loglog(errVec(:,1), errVec(:,2), errVec(:,2), errVec(:,3));
legend('VEM 1st order', 'VEM 2nd order');
p1 = polyfit(log(errVec(:,1)), log(errVec(:,2)),1);
p2 = polyfit(log(errVec(:,1)), log(errVec(:,3)),1);

p1(1)
p2(1)