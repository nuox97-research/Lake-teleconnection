%% Here, we take data from 2002 to 2012 as an exemple to show how to caluculate the network indices
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analyses
clc
clear
close all

%% Import data
load Rhomatrix_cell_limited.mat

LakeInfo = readtable('\TSI 02-12\CCM_Geo\LakeGeoInfo.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
N = size(LakeInfo,1);

%% Calculate the motifs index of the lake network
Triangle_values = [];
Triangle_weight_mean = zeros(length(Rhomatrix_cell_limited),1);
Triangle_weight_median = zeros(length(Rhomatrix_cell_limited),1);

for i = 1:length(Rhomatrix_cell_limited)
    Rhomatrix_undirect = Rhomatrix_cell_limited{i,1};

    % python  
    py.importlib.import_module('networkx');
    py.importlib.import_module('numpy'); % 导入numpy模块

    Rho_array = py.numpy.array(Rhomatrix_undirect);
    G = py.networkx.from_numpy_array(Rho_array);

    % triangles
    triangles = py.networkx.triangles(G);
    Triangle_values(:,i) = double(py.array.array('d', triangles.values()));

    % calculate the total weight of triangle motifs
    Gra = graph(Rhomatrix_undirect);
    node_weights_geo = zeros(numnodes(Gra), 1);
    for m = 1:numnodes(Gra)
        neighborNodes = neighbors(Gra, m); % get neighbors of the node i
        for j = 1:length(neighborNodes)
            for k = j+1:length(neighborNodes)
                v = neighborNodes(j);
                w = neighborNodes(k);
                edge_idx_vw = findedge(Gra, v, w);
                if edge_idx_vw > 0
                    % get the weights of three edges 
                    w_uv = Gra.Edges.Weight(findedge(Gra, m, v));
                    w_vw = Gra.Edges.Weight(findedge(Gra, v, w));
                    w_wu = Gra.Edges.Weight(findedge(Gra, w, m));
                    % ​​calculate geometric means and accumulate the results
                    geo_mean = (w_uv * w_vw * w_wu)^(1/3);
                    node_weights_geo(m) = node_weights_geo(m) + geo_mean;
                end
            end
        end
    end
    Triangle_weight_mean(i,1) = mean(node_weights_geo);
    Triangle_weight_median(i,1) = median(node_weights_geo);
end
save Triangle_values Triangle_values
save Triangle_weight_mean Triangle_weight_mean
save Triangle_weight_median Triangle_weight_median


