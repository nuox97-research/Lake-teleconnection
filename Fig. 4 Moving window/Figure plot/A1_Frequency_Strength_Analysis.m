clc 
clear
close all

%% Import data related to TSI
addpath('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Analysis')
load Filetime0212.mat
load Filetime1624.mat

Time_all_TSI = datenum(datetime('2003-01-01'):calmonths(1):datetime('2023-01-01'));

timemid_0212 = datenum(datetime('2004-05-21'):calmonths(6):datetime('2009-11-21'));
timemin_0212 = min(timemid_0212);
timemid_0212 = timemid_0212'-timemin_0212+1;
timemid_1624 = datenum(datetime('2018-04-21'):calmonths(6):datetime('2022-04-21'));
timemid_1624 = timemid_1624'-timemin_0212+1;

timemid_TSI = [timemid_0212;timemid_1624];
Time_all_TSI = Time_all_TSI-timemin_0212;

load Fre_Same_Timefilterfill.mat
load CCMmean_Same_TimefilterFill.mat

Fre_0212 = Fre_Same_Timefilterfill(1:size(Filetime0212,1));
Fre_1624 = Fre_Same_Timefilterfill(1+size(Filetime0212,1):size(Filetime1624,1)+size(Filetime0212,1));

CCMrho_mean_0212 = CCMmean_Same_TimefilterFill(1:size(Filetime0212,1));
CCMrho_mean_1624 = CCMmean_Same_TimefilterFill(1+size(Filetime0212,1):size(Filetime1624,1)+size(Filetime0212,1));

%% Import data related to FUI
% here we simply provide calculated data of FUI 
% but the calculation methodology is consistent with standard TSI procedures
Fre_TimefilterFill = struct2array(load("D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Analysis\Fre_FUI.mat") );
CCMrho_mean_TimefilterFill = struct2array(load("D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Analysis\CCMrho_mean_FUI.mat"));

timemid_start = datetime('2000-08-01')+calmonths(18);
timemid_end = datetime('2015-02-01')+calmonths(18);
timemid_FUI = datenum(timemid_start:calmonths(6):timemid_end);
timemid_FUI = timemid_FUI-timemin_0212;

Time_all = datenum(datetime('2000-01-01'):calmonths(1):datetime('2023-01-01'));
Time_all = Time_all-timemin_0212;

timelabel = {'2000','2002','2004','2006','2008','2010','2012',...
    '2014','2016','2018','2020','2022'};

FUI_breaktime = datenum(datetime('2013-11-01'))-timemin_0212;

%% import data of network indices 
% DC
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Degstat_0212.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Degstat_1624.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Degstat_FUI.mat')

% BC
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\BCstat_0212.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\BCstat_1624.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\BCstat_FUI.mat')

% CC
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\CCstat_0212.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\CCstat_1624.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\CCstat_FUI.mat')

% Closeness
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Closenessstat_0212.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Closenessstat_1624.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Closenessstat_FUI.mat')

% Triangle value
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Trianvalstat_0212.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Trianvalstat_1624.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Trianvalstat_FUI.mat')

% Triangle weight
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Trianweistat_0212.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Trianweistat_1624.mat')
load('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Network index\Trianweistat_FUI.mat')

%% Figure plot parameters
figure
hold on
titlesize = 20;
ticksize = 16;
Allposition = [0.1558,0.8571,0.3347,0.1026;... %1
    0.4894,0.8571,0.3347,0.1026; %2
    0.1558,0.7561,0.3347,0.1026; %3
    0.4894,0.7561,0.3347,0.1026; %4
    0.1558,0.6085,0.3347,0.1026; %5
    0.4894,0.6085,0.3347,0.1026; %6
    0.1558,0.5069,0.3347,0.1026; %7
    0.4894,0.5069,0.3347,0.1026; %8
    0.1558,0.3609,0.3347,0.1026; %9
    0.4894,0.3609,0.3347,0.1026; %10
    0.1558,0.2592,0.3347,0.1026; %11
    0.4894,0.2592,0.3347,0.1026; %12
    0.1558,0.1133,0.3347,0.1026; %13
    0.4894,0.1133,0.3347,0.1026; %14
    0.1558,0.0117,0.3347,0.1026; %15
    0.4894,0.0117,0.3347,0.1026;]; %16

%% TSI Teleconnection 
% TSI Fre
subplot(821,'position',Allposition(1,:))
Xfre_0212 = timemid_0212;
[Y_fit, delta] = fitconint(Xfre_0212,Fre_0212,0.05);
fill([Xfre_0212; flipud(Xfre_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Xfre_0212,Fre_0212,50,addcolor(159),"filled")
plot(Xfre_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

Xfre_1624 = timemid_1624;
[Y_fit, delta] = fitconint(Xfre_1624,Fre_1624,0.05);
fill([Xfre_1624; flipud(Xfre_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(Xfre_1624,Fre_1624,50,addcolor(159),"filled")
plot(Xfre_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])
ylim([min(Fre_0212)*0.95 max(Fre_1624)*1.05])

plot(FUI_breaktime*ones(size(0:0.001:0.9)),0:0.001:0.9,'k-.','LineWidth',1.5)

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize, ...
    'Xtick',Time_all(1:24:end),'XTickLabel','','YTick',0.1:0.02:0.18,'YTickLabel',(0.1:0.02:0.18)*100)
title('Relative frequency','FontSize',titlesize)
% xlabel('')
% ylabel('Percentage (%)');
set(gcf,'color','w') 

% TSI CCMrho
subplot(822,'position',Allposition(2,:))
X_0212 = timemid_0212;
[Y_fit delta] = fitconint(X_0212,CCMrho_mean_0212,0.05);
fill([X_0212; flipud(X_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % 半透明填充
hold on
scatter(X_0212,CCMrho_mean_0212,50,addcolor(159),"filled")
plot(X_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

X_1624 = timemid_1624;
[Y_fit delta] = fitconint(X_1624,CCMrho_mean_1624,0.05);
fill([X_1624; flipud(X_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(X_1624,CCMrho_mean_1624,50,addcolor(159),"filled")
plot(X_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);

xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])
ylim([0.495 0.525])
plot(FUI_breaktime*ones(size(0.1:0.001:0.9)),0.1:0.001:0.9,'k-.','LineWidth',1.5)

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,...
    'Xtick',Time_all(1:24:end),'XTickLabel','', ...
    'Ytick',0.5:0.01:0.53)
title('Teleconnection strength','FontSize',titlesize)
% xlabel('')
% ylabel('CCM {\it\rho}');
ax = gca;
set(ax, 'YAxisLocation', 'right');

%% FUI Teleconnection 
% FUI Fre
subplot(823,'position',Allposition(3,:))
ipT = findchangepts(Fre_TimefilterFill,'Statistic','mean','MaxNumChanges',1);
Xfre_low = 1:(ipT-1);
Xfre_high = ipT:length(Fre_TimefilterFill);
Fre_low = Fre_TimefilterFill(Xfre_low);
Fre_high = Fre_TimefilterFill(Xfre_high);
Time_low = timemid_FUI(Xfre_low)';
Time_high = timemid_FUI(Xfre_high)';

Xfre_low = Time_low;
[Y_fit delta] = fitconint(Xfre_low,Fre_low,0.05);
fill([Xfre_low; flipud(Xfre_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Xfre_low,Fre_low,50,addcolor(163),"filled")
plot(Xfre_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

Xfre_high = Time_high;
[Y_fit delta] = fitconint(Xfre_high,Fre_high,0.05);
fill([Xfre_high; flipud(Xfre_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(Xfre_high,Fre_high,50,addcolor(163),"filled")
plot(Xfre_high, Y_fit,'Color',addcolor(164),'LineWidth',2);

plot(FUI_breaktime*ones(size(0.08:0.001:0.13)),0.08:0.001:0.13,'k-.','LineWidth',1.5)

xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])
ylim([0.08 0.13])
set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,...
    'Xtick',Time_all(1:24:end),'XTickLabel','', ...
    'YTick',0.08:0.02:0.13,'YTickLabel',(0.08:0.02:0.13)*100)   
xlabel('');

% FUI CCMrho
subplot(824,'position',Allposition(4,:))
ipT = findchangepts(CCMrho_mean_TimefilterFill,'Statistic','mean','MaxNumChanges',1);
X_low = 1:(ipT-1);
X_high = ipT:length(CCMrho_mean_TimefilterFill);
CCM_low = CCMrho_mean_TimefilterFill(X_low);
CCM_high = CCMrho_mean_TimefilterFill(X_high);
Time_low = timemid_FUI(X_low)';
Time_high = timemid_FUI(X_high)';

[Y_fit delta] = fitconint(Time_low,CCM_low,0.05);
fill([Time_low; flipud(Time_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Time_low,CCM_low,50,addcolor(163),"filled")
plot(Time_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

[Y_fit delta] = fitconint(Time_high,CCM_high,0.05);
fill([Time_high; flipud(Time_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
scatter(Time_high,CCM_high,50,addcolor(163),"filled")
plot(Time_high, Y_fit,'Color',addcolor(164),'LineWidth',2);
set(gca,'YTick',0.51:0.01:0.54,'YTickLabel',0.51:0.01:0.54)
ylim([0.51 0.535])

plot(FUI_breaktime*ones(size(0.1:0.001:0.6)),0.1:0.001:0.6,'k-.','LineWidth',1.5)

% xlabel('Time','FontSize',25)
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])
set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,...
    'Xtick',Time_all(1:24:end),'XTickLabel','', 'XTickLabelRotation',45)
% ylabel('CCM {\it\rho}','FontSize',25,'Rotation',270);
ax = gca;
set(ax, 'YAxisLocation', 'right');
xlabel('')

%% Network DC
% TSI
subplot(825,'position',Allposition(5,:))
[Y_fit delta] = fitconint(timemid_0212,Degstat_0212.Degmean,0.05);
fill([timemid_0212; flipud(timemid_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(timemid_0212,Degstat_0212.Degmean,50,addcolor(159),"filled")
plot(timemid_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

[Y_fit delta] = fitconint(timemid_1624,Degstat_1624.Degmean,0.05);
fill([timemid_1624; flipud(timemid_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(timemid_1624,Degstat_1624.Degmean,50,addcolor(159),"filled")
plot(timemid_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:200)),0:200,'k-.','LineWidth',1.5)
ylim([min(Degstat_0212.Degmean)*0.95 max(Degstat_1624.Degmean)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize, ...
    'Xtick',Time_all(1:24:end),'XTickLabel','')
title('Degree centrality','FontSize',titlesize)
xlabel('')
% ylabel('Percentage (%)','FontSize',25,'Position',[-1971.1,0.273,-1]);
set(gcf,'color','w') 

% FUI
subplot(827,'position',Allposition(7,:))
ipT = findchangepts(Degstat_FUI.Degmean,'Statistic','mean','MaxNumChanges',1);
X_low = 1:(ipT-1);
X_high = ipT:length(Degstat_FUI.Degmean);
Deg_low = Degstat_FUI.Degmean(X_low);
Deg_high = Degstat_FUI.Degmean(X_high);
Time_low = timemid_FUI(X_low)';
Time_high = timemid_FUI(X_high)';

[Y_fit delta] = fitconint(Time_low,Deg_low,0.05);
fill([Time_low; flipud(Time_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Time_low,Deg_low,50,addcolor(163),"filled")
plot(Time_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

[Y_fit delta] = fitconint(Time_high,Deg_high,0.05);
fill([Time_high; flipud(Time_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
scatter(Time_high,Deg_high,50,addcolor(163),"filled")
plot(Time_high, Y_fit,'Color',addcolor(164),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:200)),0:200,'k-.','LineWidth',1.5)
plot((timemid_FUI(ipT)-30*3)*ones(size(0:200)),0:200,'k--','LineWidth',1.5)
ylim([min(Deg_low)*0.95 max(Deg_high)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,...
    'Xtick',Time_all(1:24:end),'XTickLabel','','Ytick',28:4:36)    
xlabel('');
% ylabel('Percentage (%)','FontSize',25)

%% Network BC
% TSI
subplot(826,'position',Allposition(6,:))
[Y_fit delta] = fitconint(timemid_0212,BCstat_0212.BCmean,0.05);
fill([timemid_0212; flipud(timemid_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(timemid_0212,BCstat_0212.BCmean,50,addcolor(159),"filled")
plot(timemid_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

[Y_fit delta] = fitconint(timemid_1624,BCstat_1624.BCmean,0.05);
fill([timemid_1624; flipud(timemid_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(timemid_1624,BCstat_1624.BCmean,50,addcolor(159),"filled")
plot(timemid_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:500)),0:500,'k-.','LineWidth',1.5)
ylim([min(BCstat_1624.BCmean)*0.95 max(BCstat_0212.BCmean)*1.03])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize, ...
    'Xtick',Time_all(1:24:end),'XTickLabel','')
title('Betweenness centrality','FontSize',titlesize)
xlabel('')
% ylabel('Percentage (%)','FontSize',25,'Position',[-1971.1,0.273,-1]);
ax = gca;
set(ax, 'YAxisLocation', 'right');

% FUI
subplot(828,'position',Allposition(8,:))
ipT = findchangepts(BCstat_FUI.BCmean,'Statistic','mean','MaxNumChanges',1);
X_low = 1:(ipT-1);
X_high = ipT:length(BCstat_FUI.BCmean);
BC_low = BCstat_FUI.BCmean(X_low);
BC_high = BCstat_FUI.BCmean(X_high);
Time_low = timemid_FUI(X_low)';
Time_high = timemid_FUI(X_high)';

[Y_fit delta] = fitconint(Time_low,BC_low,0.05);
fill([Time_low; flipud(Time_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Time_low,BC_low,50,addcolor(163),"filled")
plot(Time_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

[Y_fit delta] = fitconint(Time_high,BC_high,0.05);
fill([Time_high; flipud(Time_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(Time_high,BC_high,50,addcolor(163),"filled")
plot(Time_high, Y_fit,'Color',addcolor(164),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:500)),0:500,'k-.','LineWidth',1.5)
plot((timemid_FUI(ipT)-25*3)*ones(size(0:500)),0:500,'k--','LineWidth',1.5)
ylim([min(BC_low)*0.95 max(BC_high)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,...
    'Xtick',Time_all(1:24:end),'XTickLabel','') 
xlabel('');
% ylabel('Percentage (%)','FontSize',25)
ax = gca;
set(ax, 'YAxisLocation', 'right');

%% Network CC
% FUI
subplot(8,2,11,'position',Allposition(11,:))
ipT = findchangepts(CCstat_FUI.CCmean,'Statistic','mean','MaxNumChanges',1);
X_low = 1:(ipT-1);
X_high = ipT:length(CCstat_FUI.CCmean);
CC_low = CCstat_FUI.CCmean(X_low);
CC_high = CCstat_FUI.CCmean(X_high);
Time_low = timemid_FUI(X_low)';
Time_high = timemid_FUI(X_high)';

[Y_fit delta] = fitconint(Time_low,CC_low,0.05);
fill([Time_low; flipud(Time_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Time_low,CC_low,50,addcolor(163),"filled")
plot(Time_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

[Y_fit delta] = fitconint(Time_high,CC_high,0.05);
fill([Time_high; flipud(Time_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(Time_high,CC_high,50,addcolor(163),"filled")
plot(Time_high, Y_fit,'Color',addcolor(164),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:500)),0:500,'k-.','LineWidth',1.5)
plot((timemid_FUI(ipT)-25*3)*ones(size(0:500)),0:500,'k--','LineWidth',1.5)
ylim([min(CC_low)*0.95 max(CC_high)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,...
    'Xtick',Time_all(1:24:end),'XTickLabel','') 

% TSI
subplot(8,2,9,'position',Allposition(9,:))
[Y_fit delta] = fitconint(timemid_0212,CCstat_0212.CCmean,0.05);
fill([timemid_0212; flipud(timemid_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(timemid_0212,CCstat_0212.CCmean,50,addcolor(159),"filled")
plot(timemid_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

[Y_fit delta] = fitconint(timemid_1624,CCstat_1624.CCmean,0.05);
fill([timemid_1624; flipud(timemid_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(timemid_1624,CCstat_1624.CCmean,50,addcolor(159),"filled")
plot(timemid_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:500)),0:500,'k-.','LineWidth',1.5)
ylim([min(CCstat_0212.CCmean)*0.95 max(CCstat_1624.CCmean)*1.06])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize, ...
    'Xtick',Time_all(1:24:end),'XTickLabel','','YTick',0.15:0.02:0.21,'YTickLabel',0.15:0.02:0.21)
title('Clustering centrality','FontSize',titlesize)
xlabel('')
% ylabel('Percentage (%)','FontSize',25,'Position',[-1971.1,0.273,-1]);   

%% Network Closeness
% FUI
subplot(8,2,12,'position',Allposition(12,:))
ipT = findchangepts(Closenessstat_FUI.Closemean,'Statistic','mean','MaxNumChanges',1);
X_low = 1:(ipT-1);
X_high = ipT:length(Closenessstat_FUI.Closemean);
Closeness_low = Closenessstat_FUI.Closemean(X_low);
Closeness_high = Closenessstat_FUI.Closemean(X_high);
Time_low = timemid_FUI(X_low)';
Time_high = timemid_FUI(X_high)';

[Y_fit delta] = fitconint(Time_low,Closeness_low,0.05);
fill([Time_low; flipud(Time_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Time_low,Closeness_low,50,addcolor(163),"filled")
plot(Time_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

[Y_fit delta] = fitconint(Time_high,Closeness_high,0.05);
fill([Time_high; flipud(Time_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(Time_high,Closeness_high,50,addcolor(163),"filled")
plot(Time_high, Y_fit,'Color',addcolor(164),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:200)),0:200,'k-.','LineWidth',1.5)
plot((timemid_FUI(ipT)-25*3)*ones(size(0:200)),0:200,'k--','LineWidth',1.5)
ylim([min(Closeness_low)*0.95 max(Closeness_high)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,...
    'Xtick',Time_all(1:24:end),'XTickLabel','')
xlabel('');
ax = gca;
set(ax, 'YAxisLocation', 'right');

% TSI
subplot(8,2,10,'position',Allposition(10,:))
[Y_fit delta] = fitconint(timemid_0212,Closenessstat_0212.Closemean,0.05);
fill([timemid_0212; flipud(timemid_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(timemid_0212,Closenessstat_0212.Closemean,50,addcolor(159),"filled")
plot(timemid_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

[Y_fit delta] = fitconint(timemid_1624,Closenessstat_1624.Closemean,0.05);
fill([timemid_1624; flipud(timemid_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(timemid_1624,Closenessstat_1624.Closemean,50,addcolor(159),"filled")
plot(timemid_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:200)),0:200,'k-.','LineWidth',1.5)
ylim([min(Closenessstat_0212.Closemean)*0.95 max(Closenessstat_1624.Closemean)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize, ...
    'Xtick',Time_all(1:24:end),'XTickLabel','')
title('Closeness centrality','FontSize',titlesize)
xlabel('')
ax = gca;
set(ax, 'YAxisLocation', 'right');

%% Network Triangle value
% FUI
subplot(8,2,15,'position',Allposition(15,:))
ipT = findchangepts(Trianvalstat_FUI.Trianvalmean,'Statistic','mean','MaxNumChanges',1);
X_low = 1:(ipT-1);
X_high = ipT:length(Trianvalstat_FUI.Trianvalmean);
Triangleval_low = Trianvalstat_FUI.Trianvalmean(X_low);
Triangleval_high = Trianvalstat_FUI.Trianvalmean(X_high);
Time_low = timemid_FUI(X_low)';
Time_high = timemid_FUI(X_high)';

[Y_fit delta] = fitconint(Time_low,Triangleval_low,0.05);
fill([Time_low; flipud(Time_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Time_low,Triangleval_low,50,addcolor(163),"filled")
plot(Time_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

[Y_fit delta] = fitconint(Time_high,Triangleval_high,0.05);
fill([Time_high; flipud(Time_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(Time_high,Triangleval_high,50,addcolor(163),"filled")
plot(Time_high, Y_fit,'Color',addcolor(164),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:50000)),0:50000,'k-.','LineWidth',1.5)
plot((timemid_FUI(ipT)-30*3)*ones(size(0:50000)),0:50000,'k--','LineWidth',1.5)
ylim([min(Triangleval_low)*0.95 max(Triangleval_high)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,'Xtick',Time_all(1:24:end))  %,...
    %'Xtick',Time_all(1:24:end),'XTickLabel',timelabel) 
    
% TSI
subplot(8,2,13,'position',Allposition(13,:))
[Y_fit delta] = fitconint(timemid_0212,Trianvalstat_0212.Trianvalmean,0.05);
fill([timemid_0212; flipud(timemid_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(timemid_0212,Trianvalstat_0212.Trianvalmean,50,addcolor(159),"filled")
plot(timemid_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

[Y_fit delta] = fitconint(timemid_1624,Trianvalstat_1624.Trianvalmean,0.05);
fill([timemid_1624; flipud(timemid_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(timemid_1624,Trianvalstat_1624.Trianvalmean,50,addcolor(159),"filled")
plot(timemid_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:50000)),0:50000,'k-.','LineWidth',1.5)
ylim([min(Trianvalstat_0212.Trianvalmean)*0.95 max(Trianvalstat_1624.Trianvalmean)*1.08])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize, ...
    'Xtick',Time_all(1:24:end),'XTickLabel','','Ytick',500:500:2500,'YTickLabel',500:500:2500)
title('Triangle number','FontSize',titlesize)

%% Network Triangle weight
% FUI
subplot(8,2,16,'position',Allposition(16,:))
ipT = findchangepts(Trianweistat_FUI.Trianweimean,'Statistic','mean','MaxNumChanges',1);
X_low = 1:(ipT-1);
X_high = ipT:length(Trianweistat_FUI.Trianweimean);
Trianglewei_low = Trianweistat_FUI.Trianweimean(X_low);
Trianglewei_high = Trianweistat_FUI.Trianweimean(X_high);
Time_low = timemid_FUI(X_low)';
Time_high = timemid_FUI(X_high)';

[Y_fit delta] = fitconint(Time_low,Trianglewei_low,0.05);
fill([Time_low; flipud(Time_low)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(Time_low,Trianglewei_low,50,addcolor(163),"filled")
plot(Time_low, Y_fit,'Color',addcolor(164),'LineWidth',2);

[Y_fit delta] = fitconint(Time_high,Trianglewei_high,0.05);
fill([Time_high; flipud(Time_high)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(163), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(Time_high,Trianglewei_high,50,addcolor(163),"filled")
plot(Time_high, Y_fit,'Color',addcolor(164),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:70000)),0:70000,'k-.','LineWidth',1.5)
plot((timemid_FUI(ipT)-30*3)*ones(size(0:70000)),0:70000,'k--','LineWidth',1.5)
ylim([min(Trianglewei_low)*0.95 max(Trianglewei_high)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize,'Xtick',Time_all(1:24:end)) %,...
    %,'XTickLabel',timelabel)   
xlabel('');
ax = gca;
set(ax, 'YAxisLocation', 'right');

% TSI
subplot(8,2,14,'position',Allposition(14,:))
[Y_fit delta] = fitconint(timemid_0212,Trianweistat_0212.Trianweimean,0.05);
fill([timemid_0212; flipud(timemid_0212)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
scatter(timemid_0212,Trianweistat_0212.Trianweimean,50,addcolor(159),"filled")
plot(timemid_0212, Y_fit,'Color',addcolor(160),'LineWidth',2);

[Y_fit delta] = fitconint(timemid_1624,Trianweistat_1624.Trianweimean,0.05);
fill([timemid_1624; flipud(timemid_1624)], [Y_fit - delta; flipud(Y_fit + delta)], ...
    addcolor(159), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
scatter(timemid_1624,Trianweistat_1624.Trianweimean,50,addcolor(159),"filled")
plot(timemid_1624, Y_fit,'Color',addcolor(160),'LineWidth',2);
xlim([min(timemid_FUI)-180 max(timemid_TSI)+180])

plot(FUI_breaktime*ones(size(0:70000)),0:70000,'k-.','LineWidth',1.5)
ylim([min(Trianweistat_0212.Trianweimean)*0.95 max(Trianweistat_1624.Trianweimean)*1.05])

set(gca,'FontName','Arial','FontWeight','normal','FontSize',ticksize, ...
    'Xtick',Time_all(1:24:end),'XTickLabel','')
title('Triangle weight','FontSize',titlesize)
xlabel('')
ax = gca;
set(ax, 'YAxisLocation', 'right');

set(gcf,'Position',[603 50 815 945],'Color',[1 1 1])



