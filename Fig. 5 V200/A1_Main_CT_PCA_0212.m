clc
clear
close all

%% loading data
CT_raw = readtable('Fullfill_CT_monthly.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
CT_name = CT_raw.Properties.VariableNames;
CT_raw = table2array(CT_raw);

Time = datenum(datetime('2000-01-01'):calmonths(1):datetime('2012-12-01'))';

%% PCA
CT = zscore(deseason(CT_raw,Time),[],1);
    
[coeff, score, latent, tsquared, explained] = pca(CT);

PCA_ts_0212 = score(:, 1:4);
save PCA_ts_0212 PCA_ts_0212

