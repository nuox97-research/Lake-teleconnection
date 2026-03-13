%%  Here, we take data from 2002 to 2012 as an exemple to calculate calculate CT-Climate-TSI relations
%%% The computational procedures remain consistent for TSI (2016-2024) analyses
%%% Alongside this workflow for processing 2002-2012 data, 
%%% we directly provide finalized results for the TSI (2016-2024) datasets in this folder to support downstream simulations
clc
clear
close all

%% Load data
fontsize = 22;
textsize = 15;

Lake_Class_result = readtable('Lake_Climate_CT_Smap_0212.csv','ReadVariableNames',true,'VariableNamingRule','preserve');

%% Classification
Lakecom = Lake_Class_result(:,1:2);
[Lakecom_uni,ia,ic] = unique(Lakecom,"rows","stable");

Classification_full_0212 = [table2array(Lakecom_uni),zeros(size(Lakecom_uni,1),13)];
Classification_simple_0212 = [table2array(Lakecom_uni),zeros(size(Lakecom_uni,1),9)];
for i = 1:size(Lakecom_uni,1)
    i
    Com_spc = Lake_Class_result(ic == i,:);
    class = unique(Com_spc.class);
    if any(class<9)
        class(class>=9) = [];
        Classification_full_0212(i,class+2) = 1;
        Classification_simple_0212(i,class+2) = 1;
    else
        Classification_full_0212(i,class+2) = 1;
        Classification_simple_0212(i,9+2) = 1;
    end
end
Classification_full_0212 = array2table(Classification_full_0212);
Classification_full_0212.Properties.VariableNames = {'lib','target','class1','class2','class3','class4',...
    'class5','class6','class7','class8','class9','class10','class11','class12','class13'};
save Classification_full_0212 Classification_full_0212

Classification_simple_0212 = array2table(Classification_simple_0212);
Classification_simple_0212.Properties.VariableNames = {'lib','target','class1','class2','class3','class4',...
    'class5','class6','class7','class8','class9'};
save Classification_simple_0212 Classification_simple_0212


