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

%% Frequency
Lake_Class1 = Classification_full(Classification_full.class1 > 0 | Classification_full.class3 > 0,:);
Lake_Class2 = Classification_full(Classification_full.class2 > 0 | Classification_full.class4 > 0,:);
Lake_Class3 = Classification_full(Classification_full.class5 > 0 | Classification_full.class7 > 0,:);
Lake_Class4 = Classification_full(Classification_full.class6 > 0 | Classification_full.class8 > 0,:);

Lake_Class_fre_union = [size(Lake_Class1,1),size(Lake_Class3,1),size(Lake_Class2,1),size(Lake_Class4,1)]./size(Classification_full,1)*100;

figure
h = bar(1:4,Lake_Class_fre_union,0.8,'FaceColor',[.75 .75 .75],'EdgeColor',[.5 .5 .5]);
box off
% ylim([0 60])
set(gca,'FontName','Arial','Fontsize',fontsize,...
    'XTick',1:4,'XTickLabel',{'Class 1','Class 2','Class 3','Class 4'},'XTickLabelRotation',30)
xtips = h.XEndPoints;
ytips = Lake_Class_fre_union+0.01;
labels = string(Lake_Class_fre_union);
for i = 1:length(xtips)
    text(xtips(i),ytips(i),sprintf('%.2f',labels(i)),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',textsize,'FontName','Arial')
end
ylabel('Percentage (%)','FontSize',fontsize)
set(gcf,'color',[1 1 1])
