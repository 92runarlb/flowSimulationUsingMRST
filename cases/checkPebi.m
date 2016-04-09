clc; clear all; close all

G = pebiGrid(1/10, [1,1]);
plotGrid(G);
hold on
plot(G.nodes.coords(:,1), G.nodes.coords(:,2),'*')

