%% Here, we take data from 2002 to 2012 as an exemple to show how to caluculate the network indices
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analyses
clc
clear
close all

%% Import data
load Rhomatrix_cell_limited.mat

%% Calculate the modularity of the lake network
Modularity = [];
for i = 1:length(Rhomatrix_cell_limited)
    Rhomatrix_undirect = Rhomatrix_cell_limited{i,1};

    % python  
    py.importlib.import_module('community');
    py.importlib.import_module('networkx');
    py.importlib.import_module('numpy'); % 导入numpy模块

    Rho_array = py.numpy.array(Rhomatrix_undirect);
    G = py.networkx.from_numpy_array(Rho_array);

    % Louvain
    partition = py.community.best_partition(G, pyargs('weight', 'weight'));

    % calculate the modularity
    modularity = py.community.modularity(partition, G, pyargs('weight', 'weight'));

    Modularity(i,1) = modularity;
end
save Modularity Modularity
