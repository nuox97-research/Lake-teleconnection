clc
clear
close all

%% Load data
fontsize = 24;
textsize = 16;

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

%% CT
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

% to get the CT combination
Lake_Class_result = [Lake_Class_result_0212;Lake_Class_result_1624];
CTcom = table2array(Lake_Class_result(:,5:6));
[uniqueComb, ~, idx] = unique(CTcom, 'rows'); % CFcom(ia,:) = uniqueComb; CFcom = uniqueComb(idx,:
frequency = accumarray(idx, 1);
tab = table(uniqueComb, frequency);

CT_CT = zeros(14,14);
count = 1;
for i = 1:size(tab,1)
    row = tab.uniqueComb(count,1);
    col = tab.uniqueComb(count,2);
    CT_CT(row,col) = tab.frequency(count)/size(CTcom,1)*100;
    count = count+1;
end

% to get the CT percent
CTcontrol = CTcom(:);
CTpercent = tabulate(CTcontrol);

CTnames = {'AO','PDO','PNA','WP','EPNP','ENSO','AMO','TNA','TSA','NAO','EA','CAR','IOD','AAO'};

% sort
[CTsorted,order_id] = sortrows(CTpercent,2,'descend');
CTsorted = CTsorted(:,3);
save order_id order_id

%% figure
% the percent of CT
figure
cmap = addcolor([159:171,104]);

xtips = [];
for i = 1:14
    h = bar(i,CTsorted(i),0.9,'FaceColor',[.75 .75 .75],'EdgeColor','none');
    hold on
    xtips = [xtips,h.XEndPoints];
end
% ylim([0 0.45])
box off
set(gca,'FontName','Arial','Fontsize',fontsize,...
    'XTick',1:14,'XTickLabel',CTnames(order_id),...
    'XTickLabelRotation',55)
% xtips = h.XEndPoints;
ytips = CTsorted+0.01;
labels = string(CTsorted);
for i = 1:length(xtips)
    text(xtips(i),ytips(i),sprintf('%.1f',labels(i)),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',textsize,'FontName','Arial')
end
ylabel('Percentage (%)','FontSize',fontsize)
set(gcf,'Position',[619 280 721 562],'Color',[1 1 1])


