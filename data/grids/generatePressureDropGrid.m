%--------------------------------------------------------------------------
%   Generates grids for pressure drop comparison.
%--------------------------------------------------------------------------

clc; clear; close all;

mrstModule add mimetic
addpath('~/NTNU/5/master/vem/mat/')
addpath('~/NTNU/5/master/vem/mat/VEM2D/stable/')
run('~/NTNU/5/master/project-mechanics-fractures/mystartup.m')



%% Chose grid
% Chose the type of grid.
% gT = 1      Coarse cartesian
% gT = 2      composite pebi
% gT = 3      fully unstructured grid
% gT = 4      Fine Cartesian

gT = 2;


%% Grid parameters
% shared parameters
xmax = 1;                              % Set grid dimentions
ymax = 1;

wellLine = {[xmax/2,ymax/2]};                % Set source center

switch gT
  case 1
    
    nx = sqrt(184);
    G = cartGrid([nx,nx],[xmax,ymax]);
    G = computeGeometry(G);
    w1 = wellLine{1};
    D = pdist2(G.cells.centroids, w1);
    [~, I] = min(D, [], 1);
    G.cells.tag = false(G.cells.num,1);
    G.cells.tag(I') = true(size(I'));
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);
    save('pressureDropGridCart.mat');
    
  case 2 
    
    gridSize = xmax*1/9;                   % Size of gridcells
    mlqtMax = 2;                            % Set number of reminement levels
    wellGridSize = 0.75/2^mlqtMax;          % Set relative well grid size
    mlqtSizes = 2.0*linspace(gridSize,gridSize*wellGridSize,mlqtMax+1)';
                                            % Size around wells to be refined 
    G = compositePebiGrid(gridSize, [xmax, ymax], 'wellLines', wellLine, ...
                         'wellGridFactor', wellGridSize, ...
                         'mlqtMaxLevel', 2, 'mlqtLevelSteps', mlqtSizes);
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);
    save('pressureDropGridCompPebi.mat');
    
  case 3
      
    gridSize = xmax/8;
    wellGridSize = 0.7/2^2;
    epsilon = gridSize*.7;
    G = pebiGrid(gridSize, [xmax, ymax], 'wellLines', wellLine,   ...
                'wellGridFactor', wellGridSize, 'wellRefinement',true);
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);
    save('pressureDropGridPebi.mat');
    
  case 4
    
    nx = 301;
    G = cartGrid([nx,nx],[xmax,ymax]);
    G = computeGeometry(G);
    w1 = wellLine{1};
    D = pdist2(G.cells.centroids, w1);
    [~, I] = min(D, [], 1);
    G.cells.tag = false(G.cells.num,1);
    G.cells.tag(I') = true(size(I'));
    G = sortEdges(G);
    G = computeVEM2DGeometry(G);
    save('pressureDropGridCartFine.mat');
    
    otherwise
      
    error('unknown grid case')
    
end                  

plotGrid(G);