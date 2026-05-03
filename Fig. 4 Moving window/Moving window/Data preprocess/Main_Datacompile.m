%% This m.file is for preparing time-series chunks for moving window analysis
%% Here, we only take data of 2002-2012 as an exemple
clc
clear

%% Import data
TSI = readtable('\Lake teleconnection code\Data_preprocess_new\A2_Detrend\De_Seasonal_FilterFilled.csv','ReadVariableNames',true,'VariableNamingRule','preserve');

load('\Lake teleconnection code\Data_preprocess_new\A1_Datafilling\TS_Time_FILLED.mat')

%% MW = 6 months, Each window = 3y
Time = datetime(TS_Time_FILLED,'ConvertFrom','datenum');

count = 1;
MW_start = 1; 
MW_end = find(Time == Time(MW_start)+calyears(3));

while 1
    TSI_MW = TSI(:,MW_start+3:MW_end+3);
    TSI_MW = [TSI(:,1:3),TSI_MW];
    writetable(TSI_MW,[strcat('TSI_MW',num2str(count)),'_',datestr(Time(MW_start),'yyyymmdd'),'_',datestr(Time(MW_end),'yyyymmdd'),'.csv'])
    Date = array2table(Time(MW_start:MW_end));
    Date.Properties.VariableNames(1) = {'Date'}; 
    writetable(Date,['Date','_','MW',num2str(count),'.csv'])
    MW_start = find(Time == Time(MW_start)+calmonths(6));
    MW_end = find(Time == Time(MW_start)+calyears(3));
    TSI_MW = table();
    if isempty(MW_end)
        break
    end
    count = count+1;
end





