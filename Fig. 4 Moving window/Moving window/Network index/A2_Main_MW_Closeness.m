%% Here, we take data from 2002 to 2012 as an exemple to show how to caluculate the network indices
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analyses
clc
clear
close all

%% Import data
load Rhomatrix_cell_limited.mat

%% Calculate the clossness of the lake network
Closeness = [];
for i = 1:length(Rhomatrix_cell_limited)
    Rhomatrix_undirect = Rhomatrix_cell_limited{i,1};
    G = graph(Rhomatrix_undirect);
    Closeness(:,i) = centrality(G, 'closeness', 'Cost', G.Edges.Weight);
end
save Closeness Closeness
