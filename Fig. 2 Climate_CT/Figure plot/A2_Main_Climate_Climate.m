clc
clear
close all

%% Load data
fontsize = 22;
textsize = 15;

% 02-12
load ('D:\Lake teleconnection code\Fig. 1 Teleconnection\Lake teleconnection_CCM\Sync_harm_deSeason_0212.mat') 
Sync_harm_0212 = Sync_harm_deSeason_0212(Sync_harm_deSeason_0212.Optlag_rho>0.4 & Sync_harm_deSeason_0212.Dist>3500,:);
Lakecom_0212 = unique(Sync_harm_0212(:,1:2),'rows');

% 16-24
load ('D:\Lake teleconnection code\Fig. 1 Teleconnection\Lake teleconnection_CCM\Sync_harm_deSeason_1624.mat') 
Sync_harm_1624 = Sync_harm_deSeason_1624(Sync_harm_deSeason_1624.Optlag_rho>0.4 & Sync_harm_deSeason_1624.Dist>3500,:);
Lakecom_1624 = unique(Sync_harm_1624(:,1:2),'rows');

% Classification
Classification_full_0212 = struct2array(load('D:\Lake teleconnection code\Fig. 2 Climate_CT\TSI_Climate_CT\Classification_full_0212.mat'));
Classification_full_1624 = struct2array(load('D:\Lake teleconnection code\Fig. 2 Climate_CT\TSI_Climate_CT\Classification_full_1624.mat'));

Classification_full_0212 = Classification_full_0212(ismember(Classification_full_0212(:,1:2),Lakecom_0212),:);
Classification_full_1624 = Classification_full_1624(ismember(Classification_full_1624(:,1:2),Lakecom_1624),:);

Classification_full = [Classification_full_0212;Classification_full_1624];

% filter 
Lake_Class_result_0212 = readtable('D:\Lake teleconnection code\Fig. 2 Climate_CT\TSI_Climate_CT\Lake_Climate_CT_Smap_0212.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
Lake_Class_result_1624 = readtable('D:\Lake teleconnection code\Fig. 2 Climate_CT\TSI_Climate_CT\Lake_Climate_CT_Smap_1624.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
Lake_Class_result_1624.Properties.VariableNames(end) = {'class'};

Lake_Class_result_0212 = Lake_Class_result_0212(ismember(Lake_Class_result_0212(:,1:2),Lakecom_0212),:);
Lake_Class_result_1624 = Lake_Class_result_1624(ismember(Lake_Class_result_1624(:,1:2),Lakecom_1624),:);

% filter classification simple
Classification_simple_0212 = struct2array(load('D:\Lake teleconnection code\Fig. 2 Climate_CT\TSI_Climate_CT\Classification_simple_0212.mat'));
Classification_simple_1624 = struct2array(load('D:\Lake teleconnection code\Fig. 2 Climate_CT\TSI_Climate_CT\Classification_simple_1624.mat'));

Classification_simple_0212 = Classification_simple_0212(ismember(Classification_simple_0212(:,1:2),Lakecom_0212),:);
Classification_simple_1624 = Classification_simple_1624(ismember(Classification_simple_1624(:,1:2),Lakecom_1624),:);

% get not linked pairs
Classification_notlinked_0212 = Classification_simple_0212(Classification_simple_0212.class9 == 1,:);
Classification_notlinked_1624 = Classification_simple_1624(Classification_simple_1624.class9 == 1,:);

%% CF-CF
% eliminate the pairs that are not linked
% for 02-12
CFcom_notlinked_0212 = table2array(Classification_notlinked_0212(:,1:2));
cont = 1;
while cont
    [~,loc] = ismember(CFcom_notlinked_0212,table2array(Lake_Class_result_0212(:,1:2)),'rows');
    loc = loc(loc>0);
    Lake_Class_result_0212(loc,:) = [];
    if ~any(loc)
        cont = 0;
    end
end
Lake_Class_result_0212(Lake_Class_result_0212.class > 8,:) = [];

% for 16-24
CFcom_notlinked_1624 = table2array(Classification_notlinked_1624(:,1:2));
cont = 1;
while cont
    [~,loc] = ismember(CFcom_notlinked_1624,table2array(Lake_Class_result_1624(:,1:2)),'rows');
    loc = loc(loc>0);
    Lake_Class_result_1624(loc,:) = [];
    if ~any(loc)
        cont = 0;
    end
end
Lake_Class_result_1624(Lake_Class_result_1624.class > 8,:) = [];

% to get the CF combination
Lake_Class_result = [Lake_Class_result_0212;Lake_Class_result_1624];
CFcom = table2array(Lake_Class_result(:,3:4));
[uniqueComb, ~, idx] = unique(CFcom, 'rows'); % CFcom(ia,:) = uniqueComb; CFcom = uniqueComb(idx,:
frequency = accumarray(idx, 1);
tab = table(uniqueComb, frequency);

CF_CF = zeros(9,9);
count = 1;
for i = 1:9
    for j = 1:9
        CF_CF(i,j) = frequency(count)/size(CFcom,1)*100;
        count = count+1;
    end
end

% change the position for precipitation and Tmax
CF_CF_change = CF_CF([1,7,3:6,2,8:9],[1,7,3:6,2,8:9]);

%% figure
CFnames = {'T','Tmax','SLH','SM','SSH','SP','P','SR','WS'};

figure
CMP = corrMatPlot_hold_forCF(CF_CF_change,'type','ssq');
Mycolor = flipud(cbrewer2('seq','RdYlBu'));
CMP=CMP.setLabelStr(CFnames);
CMP=CMP.setColorMap(Mycolor);
CMP=CMP.draw();
CMP.setCorrTxt('Color',[0 0 0],'Fontsize',22)
CMP.setXLabel('FontSize',25)
CMP.setYLabel('FontSize',25)
set(gcf,'Position',[600 237 699 562])



