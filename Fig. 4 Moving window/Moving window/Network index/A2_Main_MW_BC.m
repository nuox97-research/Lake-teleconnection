%%  Here, we take data from 2002 to 2012 as an exemple to show how to caluculate the network indices
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analysesclc
clear
close all

%% Import data
load Rhomatrix_cell_limited.mat

%% Calculate the betweenness centrality of the lake network
BC = [];
for i = 1:length(Rhomatrix_cell_limited)
    Rhomatrix_undirect = Rhomatrix_cell_limited{i,1};
    G = graph(Rhomatrix_undirect);
    BC(:,i) = centrality(G, 'betweenness', 'Cost', G.Edges.Weight);
end
save BC BC
