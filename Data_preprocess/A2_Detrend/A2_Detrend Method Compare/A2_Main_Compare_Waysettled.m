clc
clear
close all

%% Import data
Methodcompare = readtable('Methods_Compare.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
X2Y = Methodcompare(:,[1,2,3:2:7]);
Y2X = Methodcompare(:,[1,2,4:2:8]);

X2Y.Properties.VariableNames(1:5) = {'lib', 'tar', 'stl', 'harm', 'mean'};

Y2X.Properties.VariableNames(1:5) = {'lib', 'tar', 'stl', 'harm', 'mean'};

%% X2Y
% stl
CCM_STL_X2Y = zeros(6,6);
for i = 1:size(X2Y,1)
    way_X = X2Y.lib(i);
    way_Y = X2Y.tar(i);
    CCM_STL_X2Y(way_X,way_Y) = CCM_STL_X2Y(way_X,way_Y) + X2Y.stl(i);
end
CCM_STL_X2Y = CCM_STL_X2Y./50;
CMP = corrMatPlot_hold(CCM_STL_X2Y,'type','ssq');
Mycolor = cbrewer2('seq','Blues');
CMP=CMP.setLabelStr({'Simple','Multiple','Varyamp','Varytime','Varyamp\_time','Extra'});
CMP=CMP.setColorMap(Mycolor);
CMP=CMP.draw();
CMP.setCorrTxt('Color',[1 1 1])
set(gcf,'Position',[471 279 871 560])

% harm
CCM_Harm_X2Y = zeros(6,6);
for i = 1:size(X2Y,1)
    way_X = X2Y.lib(i);
    way_Y = X2Y.tar(i);
    CCM_Harm_X2Y(way_X,way_Y) = CCM_Harm_X2Y(way_X,way_Y) + X2Y.harm(i);
end
CCM_Harm_X2Y = CCM_Harm_X2Y./50;
CMP = corrMatPlot_hold(CCM_Harm_X2Y,'type','ssq');
Mycolor = cbrewer2('seq','Blues');
CMP=CMP.setLabelStr({'Simple','Multiple','Varyamp','Varytime','Varyamp\_time','Extra'});
CMP=CMP.setColorMap(Mycolor);
CMP=CMP.draw();
CMP.setCorrTxt('Color',[1 1 1])
set(gcf,'Position',[471 279 871 560])

% mean
CCM_Mean_X2Y = zeros(6,6);
for i = 1:size(X2Y,1)
    way_X = X2Y.lib(i);
    way_Y = X2Y.tar(i);
    CCM_Mean_X2Y(way_X,way_Y) = CCM_Mean_X2Y(way_X,way_Y)+X2Y.mean(i);
end
CCM_Mean_X2Y = CCM_Mean_X2Y./50;
CMP = corrMatPlot_hold(CCM_Mean_X2Y,'type','ssq');
Mycolor = cbrewer2('seq','Blues');
CMP=CMP.setLabelStr({'Simple','Multiple','Varyamp','Varytime','Varyamp\_time','Extra'});
CMP=CMP.setColorMap(Mycolor);
CMP=CMP.draw();
CMP.setCorrTxt('Color',[1 1 1])
set(gcf,'Position',[471 279 871 560])

%% Y2X
% stl
CCM_STL_Y2X = zeros(6,6);
for i = 1:size(Y2X,1)
    way_X = Y2X.lib(i);
    way_Y = Y2X.tar(i);
    CCM_STL_Y2X(way_X,way_Y) = CCM_STL_Y2X(way_X,way_Y) + Y2X.stl(i);
end
CCM_STL_Y2X = CCM_STL_Y2X./50;
CMP = corrMatPlot_hold(CCM_STL_Y2X,'type','ssq');
Mycolor = cbrewer2('seq','Blues');
CMP=CMP.setLabelStr({'Simple','Multiple','Varyamp','Varytime','Varyamp\_time','Extra'});
CMP=CMP.setColorMap(Mycolor);
CMP=CMP.draw();
CMP.setCorrTxt('Color',[0 0 0])
set(gcf,'Position',[471 279 871 560])

% harm
CCM_Harm_Y2X = zeros(6,6);
for i = 1:size(Y2X,1)
    way_X = Y2X.lib(i);
    way_Y = Y2X.tar(i);
    CCM_Harm_Y2X(way_X,way_Y) = CCM_Harm_Y2X(way_X,way_Y)+Y2X.harm(i);
end
CCM_Harm_Y2X = CCM_Harm_Y2X./50;
CMP = corrMatPlot_hold(CCM_Harm_Y2X,'type','ssq');
Mycolor = cbrewer2('seq','Blues');
CMP=CMP.setLabelStr({'Simple','Multiple','Varyamp','Varytime','Varyamp\_time','Extra'});
CMP=CMP.setColorMap(Mycolor);
CMP=CMP.draw();
CMP.setCorrTxt('Color',[0 0 0])
set(gcf,'Position',[471 279 871 560])

% mean
CCM_Mean_Y2X = zeros(6,6);
for i = 1:size(Y2X,1)
    way_X = Y2X.lib(i);
    way_Y = Y2X.tar(i);
    CCM_Mean_Y2X(way_X,way_Y) = CCM_Mean_Y2X(way_X,way_Y)+Y2X.mean(i);
end
CCM_Mean_Y2X = CCM_Mean_Y2X./50;
CMP = corrMatPlot_hold(CCM_Mean_Y2X,'type','ssq');
Mycolor = cbrewer2('seq','Blues');
CMP=CMP.setLabelStr({'Simple','Multiple','Varyamp','Varytime','Varyamp\_time','Extra'});
CMP=CMP.setColorMap(Mycolor);
CMP=CMP.draw();
CMP.setCorrTxt('Color',[0 0 0])
set(gcf,'Position',[471 279 871 560])