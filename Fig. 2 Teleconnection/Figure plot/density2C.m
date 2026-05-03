function [CData,h,XMesh,YMesh,ZMesh,colorList] = density2C(X,Y,XList,YList,colorList)
% 输入：
% X,Y 散点坐标
% XList,YList 用来构造密度曲面网格的序列，其实就是把XLim，YLim分成小份，例如XList=0:0.1:10
% colorList 颜色表mx3数组，用来构造将高度映射到颜色函数的数据表
% 
% 输出：
% CData 各个点对应颜色
% h 各个点对应核密度
% XMesh,YMesh,ZMesh 核密度曲面数据
% colorList 插值后更细密的颜色表
%
[XMesh,YMesh] = meshgrid(XList,YList);
XYi = [XMesh(:) YMesh(:)];
F = ksdensity([X,Y],XYi);
ZMesh = zeros(size(XMesh));
ZMesh(1:length(F)) = F;

h = interp2(XMesh,YMesh,ZMesh,X,Y);
if nargin<5
    colorList=[0.2700    0         0.3300
               0.2700    0.2300    0.5100
               0.1900    0.4100    0.5600
               0.1200    0.5600    0.5500
               0.2100    0.7200    0.4700
               0.5600    0.8400    0.2700
               0.9900    0.9100    0.1300];
end
colorFunc = colorFuncFactory(colorList);
CData = colorFunc((h-min(h))./(max(h)-min(h)));
colorList = colorFunc(linspace(0,1,100)');

function colorFunc = colorFuncFactory(colorList)
    x = (0:size(colorList,1)-1)./(size(colorList,1)-1);
    y1 = colorList(:,1);
    y2 = colorList(:,2);
    y3 = colorList(:,3);
    colorFunc = @(X)[interp1(x,y1,X,'pchip'),interp1(x,y2,X,'pchip'),interp1(x,y3,X,'pchip')];
end
end
