close all; clear; addpath ../voronoi2D/
% this script performs a simple 2 phase(oil and water) flow simulation in
% Using a composite grid. The case consist of two wells (one injector and
% one producer) and a linear no-flow fault seperating them.

%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid

gT = 1;


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
%% Set simulation parameters
T      = 120*second();    % End time
dT     = T/120;           % Time steps
dTplot = 1:1:T;           % Plot at these time steps

%% Generate grid


                 
%Set internal boundary
G = makeInternalBoundary(G, find(G.faces.tag));
G = computeGeometry(G);

%% Set fluid and rock properties
gravity reset off 

fluid = initSimpleFluid('mu' , [   1,   5]*centi*poise     , ...
                        'rho', [1000, 700]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

rock.poro = ones(G.cells.num,1)*0.15;
rock.perm = ones([G.cells.num,1])*100*milli*darcy;

%% add Sources
srcCells = find(G.cells.tag);
pv = sum(poreVolume(G, rock));
src = addSource([],srcCells(1),0.01*pv,'sat',[1,0]);
src = addSource(src, srcCells(2), -0.01*pv,'sat',[0,1]);

            
%% Solve
state = initState(G, [], 0, [0.0,1]);
trans = computeTrans(G,rock);
state = incompTPFA(state, G, trans, fluid, 'src', src);
% Prepare plotting of saturations
clf;
hold on
plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
axis off equal, view([-120,30]), colormap(flipud(jet))

colorbar; hs = []; ha=[]; zoom(1.3);

% Start the main loop
t  = 0;  plotNo = 1;
while t < T,
   state = implicitTransport(state, G, dT, rock, fluid, 'src', src);

   % Check for inconsistent saturations
   assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);
   
   % Update solution of pressure equation.
   %state = incompMimetic(state, G3D, S, fluid, 'wells', W);
   state = incompTPFA(state, G, trans, fluid, 'src', src);
   % Increase time and continue if we do not want to plot saturations
   t = t + dT;   
   if ( t + dT <= dTplot(plotNo)), continue, end
   delete([hs, ha])
   hs = plotCellData(G, state.s(:,1), find(state.s(:,1) >= 0.0));
   ha = annotation('textbox', [0.1714 0.8214 0.5000 0.1000], 'LineStyle', 'none', ...
                   'String', ['Water saturation at ', ...
                              num2str(convertTo(t,second)), ' s']);
   fig = gcf();
   set(findall(fig,'-property','FontSize'),'FontSize',14) 
   view(0, 90), drawnow, caxis([0 1])
  
   plotNo = plotNo+1;
end