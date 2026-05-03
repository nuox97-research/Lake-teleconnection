%%  Here, we take data from 2002 to 2012 as an exemple to show how to caluculate the network indices
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analyses
%%% Alongside this workflow for processing 2002-2012 data, 
%%% we directly provide finalized results for the TSI (2016-2024) and FUI datasets in this folder to support downstream simulations
clc
clear
close all

%% load data
load Deg.mat
load BC.mat
load CC.mat
load Closeness.mat
load Triangle_values.mat
load Triangle_weight_mean.mat

%% Degree
Degmean = mean(Deg);
Degstd = std(Deg);
t_critical = tinv(1 - 0.05/2, size(Deg,1)-1);
ci_half_width = t_critical * (Degstd / sqrt(size(Deg,1)));
ci_lower = Degmean - ci_half_width;
ci_upper = Degmean + ci_half_width;

Degstat_0212 = array2table([Degmean;ci_lower;ci_upper]');
Degstat_0212.Properties.VariableNames = {'Degmean','Deglowci','Deghighci'};
save Degstat_0212 Degstat_0212

%% BC
BCmean = mean(BC);
BCstd = std(BC);
t_critical = tinv(1 - 0.05/2, size(BC,1)-1);
ci_half_width = t_critical * (BCstd / sqrt(size(BC,1)));
ci_lower = BCmean - ci_half_width;
ci_upper = BCmean + ci_half_width;

BCstat_0212 = array2table([BCmean;ci_lower;ci_upper]');
BCstat_0212.Properties.VariableNames = {'BCmean','BClowci','BChighci'};
save BCstat_0212 BCstat_0212

%% CC
CCmean = mean(CC);
CCstd = std(CC);
t_critical = tinv(1 - 0.05/2, size(CC,1)-1);
ci_half_width = t_critical * (CCstd / sqrt(size(CC,1)));
ci_lower = CCmean - ci_half_width;
ci_upper = CCmean + ci_half_width;

CCstat_0212 = array2table([CCmean;ci_lower;ci_upper]');
CCstat_0212.Properties.VariableNames = {'CCmean','CClowci','CChighci'};
save CCstat_0212 CCstat_0212

%% Closeness
Closemean = mean(Closeness);
Closestd = std(Closeness);
t_critical = tinv(1 - 0.05/2, size(Closeness,1)-1);
ci_half_width = t_critical * (Closestd / sqrt(size(Closeness,1)));
ci_lower = Closemean - ci_half_width;
ci_upper = Closemean + ci_half_width;

Closenessstat_0212 = array2table([Closemean;ci_lower;ci_upper]');
Closenessstat_0212.Properties.VariableNames = {'Closemean','Closelowci','Closehighci'};
save Closenessstat_0212 Closenessstat_0212

%% Triangle values
Trianvalmean = mean(Triangle_values);
Trianvalstd = std(Triangle_values);
t_critical = tinv(1 - 0.05/2, size(Triangle_values,1)-1);
ci_half_width = t_critical * (Trianvalstd / sqrt(size(Triangle_values,1)));
ci_lower = Trianvalmean - ci_half_width;
ci_upper = Trianvalmean + ci_half_width;

Trianvalstat_0212 = array2table([Trianvalmean;ci_lower;ci_upper]');
Trianvalstat_0212.Properties.VariableNames = {'Trianvalmean','Trianvallowci','Trianvalhighci'};
save Trianvalstat_0212 Trianvalstat_0212

%% Triangle weight
Trianweimean = Triangle_weight_mean;

Trianweistat_0212 = array2table(Trianweimean);
Trianweistat_0212.Properties.VariableNames = {'Trianweimean'};
save Trianweistat_0212 Trianweistat_0212