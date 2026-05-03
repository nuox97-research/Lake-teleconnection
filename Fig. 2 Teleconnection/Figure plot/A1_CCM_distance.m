clc 
clear
close all

%% Figure color
Color = [0 44 252;34 98 255;48 154 255;103 203 255;152 238 255;202 255 255;...
    255 254 203;254 238 152;254 205 101;255 152 48;255 99 25;252 43 4]/256;
n = 100;

x = 1:size(Color, 1); 
xi = linspace(1, size(Color, 1), n);
new_colors = zeros(n, 3);  
for i = 1:3  
    new_colors(:, i) = interp1(x, Color(:, i), xi, 'linear');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% CCM  TSI %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Import data
load('D:\Lake teleconnection code\Fig. 1 Teleconnection\Lake teleconnection_CCM\Sync_harm_deSeason_0212.mat')
load('D:\Lake teleconnection code\Fig. 1 Teleconnection\Lake teleconnection_CCM\Sync_harm_deSeason_1624.mat')

Sync_harm = [Sync_harm_deSeason_0212;Sync_harm_deSeason_1624];

%%%% Figure
[CData,~,~,~,~,colorList] = density2C(Sync_harm.Optlag_rho,Sync_harm.Dist,0:0.01:1,0:200:20000,Color);
figure1 = figure('Colormap',new_colors,'Color',[1 1 1],'Position',[362 166 1169 701]);
% main panel
ax1 = axes('Parent',figure1,'Position',[0.2638,0.1741,0.5,0.588]);
hold(ax1,'on')
colormap(colorList)
scatter(Sync_harm.Optlag_rho,Sync_harm.Dist,15,CData,'filled')
plot(0:0.1:1,3500*ones(size(0:0.1:1)),'--','Color',[.5 .5 .5],'LineWidth',2)
set(ax1,'FontName','Arial','FontSize',22,'FontWeight','normal',...
    'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1],'XTickLabel',...
    ['0  ';'0.2';'0.4';'0.6';'0.8';'1  '],...
    'Ytick',[5000 10000 15000 20000]);
xlabel('Distance (km)','FontWeight','normal','FontSize',20,'FontName','Arial');
ylabel({'CCM {\rho}'},'FontWeight','normal','FontSize',20,'FontName','Arial');
xlim([0 1])

% CCM rho distribution
ax4 = axes('Parent',figure1,'Position',[0.2638 0.8327 0.5 0.15]); % left bottom width height
hold(ax4,'on');
mu_lake = mean(log(Sync_harm.Optlag_rho));
sigma_lake = std(log(Sync_harm.Optlag_rho));
x = linspace(0, 1, 1000)'; % 定义绘图范围
pdf_lake = lognpdf(x, mu_lake, sigma_lake);
% 
histogram(Sync_harm.Optlag_rho, 30, 'Normalization', 'pdf', ...
    'FaceColor', Color(6,:), 'FaceAlpha', 0.4, 'EdgeColor', Color(2,:)); % 绘制归一化直方图
% plot(x, pdf_lake,'Color',[.4 .4 .4], 'LineWidth', 2, 'DisplayName', 'Total Distribution'); % 总分布
set(ax4,'TickDir','out','XColor','none','YColor','none','YTickLabel','');
ax4.XLim=ax1.XLim;

% Distance distribution
ax2 = axes('Parent',figure1,'Position',[0.7766 0.1741 0.15 0.588]);
hold(ax2,'on')
gmm = fitgmdist(log(Sync_harm.Dist), 2, 'CovarianceType','full',...
    'Start','plus','Options',statset('MaxIter',1000));
% 混合分布参数
comp_prop = gmm.ComponentProportion;
mu_components = gmm.mu;
sigma_components = squeeze(sqrt(gmm.Sigma));
% 生成PDF曲线
x_dist = linspace(0, 20000, 1000)';
pdf1 = comp_prop(1) * lognpdf(x_dist, mu_components(1), sigma_components(1));
pdf2 = comp_prop(2) * lognpdf(x_dist, mu_components(2), sigma_components(2));
% pdf_total = pdf1 + pdf2;
% 绘图
histogram(Sync_harm.Dist, 50, 'Normalization', 'pdf', 'Orientation', 'horizontal',...
    'FaceColor', Color(6,:), 'FaceAlpha', 0.4, 'EdgeColor', Color(4,:)); % 绘制归一化直方图
plot(pdf1,x_dist,'Color','#467B29', 'LineWidth', 2, 'DisplayName', 'Component 1'); % 分量 1
plot(pdf2,x_dist,'Color','#9E4601', 'LineWidth', 2, 'DisplayName', 'Component 2'); % 分量 2
% plot(x, pdf_total+3e-6, 'Color',[.4 .4 .4], 'LineWidth', 2, 'DisplayName', 'Total Distribution'); % 总分布
set(ax2,'TickDir','out','XColor','none','YColor','none','YTickLabel','');
ax2.YLim=ax1.YLim;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% CCM  FUI %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Import data
load('D:\Lake teleconnection code\Fig. 1 Teleconnection\Lake teleconnection_CCM\Sync_FUI.mat')

%%%% Figure
[CData,~,~,~,~,colorList] = density2C(Sync_FUI.Optlag_rho,Sync_FUI.Dist,0:0.01:1,0:200:20000,Color);
figure2 = figure('Colormap',new_colors,'Color',[1 1 1],'Position',[362 166 1169 701]);
% main panel
ax1 = axes('Parent',figure2,'Position',[0.2638,0.1741,0.5,0.588]);
hold(ax1,'on')
colormap(colorList)
scatter(Sync_FUI.Optlag_rho,Sync_FUI.Dist,15,CData,'filled')
% contourf(XMesh,YMesh,ZMesh,10,'EdgeColor','none')
set(ax1,'FontName','Arial','FontSize',22,'FontWeight','normal',...
    'TickDir','out','XTick',[0 0.2 0.4 0.6 0.8 1],'XTickLabel',...
    ['0  ';'0.2';'0.4';'0.6';'0.8';'1  '],...
    'Ytick',[5000 10000 15000 20000]);
% xlabel('Distance (km)','FontWeight','normal','FontSize',20,'FontName','Arial');
% ylabel({'Synchronous strength','between lake pairs'},'FontWeight','normal','FontSize',20,'FontName','Arial',...
    % 'position',[-4809.929038383438,0.500000476837158,-0.999999999999986]);
xlim([0 1])

% CCM rho distribution
ax4 = axes('Parent',figure2,'Position',[0.2638 0.8327 0.5 0.15]); % left bottom width height
hold(ax4,'on');
mu_lake = mean(log(Sync_FUI.Optlag_rho));
sigma_lake = std(log(Sync_FUI.Optlag_rho));
% 绘制分布
x = linspace(0, 1, 1000)'; % 定义绘图范围
pdf_lake = lognpdf(x, mu_lake, sigma_lake);
% 绘图
histogram(Sync_FUI.Optlag_rho, 50, 'Normalization', 'pdf', ...
    'FaceColor', Color(6,:), 'FaceAlpha', 0.4, 'EdgeColor', Color(4,:)); % 绘制归一化直方图
% plot(x, pdf1, 'Color','#467B29', 'LineWidth', 2, 'DisplayName', 'Component 1'); % 分量 1
% plot(x, pdf2, 'Color','#9E4601', 'LineWidth', 2, 'DisplayName', 'Component 2'); % 分量 2
plot(x, pdf_lake,'Color','#467B29', 'LineWidth', 2, 'DisplayName', 'Total Distribution'); % 总分布
set(ax4,'TickDir','out','XColor','none','YColor','none','YTickLabel','');
ax4.XLim=ax1.XLim;

ax2 = axes('Parent',figure2,'Position',[0.7766 0.1741 0.15 0.588]);
hold(ax2,'on')
gmm = fitgmdist(log(Sync_FUI.Dist), 2, 'CovarianceType','full',...
    'Start','plus','Options',statset('MaxIter',1000));
% 混合分布参数
comp_prop = gmm.ComponentProportion;
mu_components = gmm.mu;
sigma_components = squeeze(sqrt(gmm.Sigma));
% 生成PDF曲线
x_dist = linspace(0, 20000, 1000)';
pdf1 = comp_prop(1) * lognpdf(x_dist, mu_components(1), sigma_components(1));
pdf2 = comp_prop(2) * lognpdf(x_dist, mu_components(2), sigma_components(2));
% pdf_total = pdf1 + pdf2;
% 绘图
histogram(Sync_FUI.Dist, 50, 'Normalization', 'pdf', 'Orientation', 'horizontal',...
    'FaceColor', Color(6,:), 'FaceAlpha', 0.4, 'EdgeColor', Color(4,:)); % 绘制归一化直方图
plot(pdf1,x_dist,'Color','#467B29', 'LineWidth', 2, 'DisplayName', 'Component 1'); % 分量 1
plot(pdf2,x_dist,'Color','#9E4601', 'LineWidth', 2, 'DisplayName', 'Component 2'); % 分量 2
% plot(x, pdf_total+3e-6, 'Color',[.4 .4 .4], 'LineWidth', 2, 'DisplayName', 'Total Distribution'); % 总分布
set(ax2,'TickDir','out','XColor','none','YColor','none','YTickLabel','');
ax2.YLim=ax1.YLim;
