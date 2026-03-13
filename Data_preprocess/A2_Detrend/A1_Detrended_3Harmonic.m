%% This m.file is for removing seasonality of time series 
%% by a harmonic fitting approach from Smith & Boers, 2023

% Reference: Smith T., Boers N. Reliability of vegetation resilience estimates
% depends on biomass density. Nature ecology & evolution 7, 1799-1808.
clc
clear

%% Import data
load TS_FILLTERED_FILLED.mat
load TS_Time_FILLED.mat

TS_Time_FILLED = datenum(TS_Time_FILLED);

%% deSeason
% get time series data
TSI_data = table2array(TS_FILLTERED_FILLED(:,4:end));

% a 3-order harmonic function fit to the data as a deseasoner
harm_freq = 1:3;
Time = TS_Time_FILLED;
Timestamp = Time(1);

order = 3;
xt = (Time-Timestamp)/365.25;
x_rad = xt * 2 * pi;
nr_indep = order*2 + 2;
Indep = zeros(size(TSI_data,2), nr_indep);
Indep(:,1) = ones(size(TSI_data,2),1);
Indep(:,2) = x_rad;

i = 3;
for freq = harm_freq
    cos_freq = cos(x_rad * freq);
    sin_freq = sin(x_rad * freq);
    Indep(:,i) = cos_freq;
    i = i + 1;
    Indep(:,i) = sin_freq;
    i = i + 1;
end

for m = 1:size(TSI_data,1)
    [b,bint,r,rint,stats] = regress(TSI_data(m,:)',Indep);
    fitted = Indep*b;
    Seasonal(m,:) = fitted';
end

De_Seasonal = [TS_FILLTERED_FILLED.Hylak_id,TS_FILLTERED_FILLED.Lat,TS_FILLTERED_FILLED.Lon,TSI_data-Seasonal];
De_Seasonal = array2table(De_Seasonal);
De_Seasonal.Properties.VariableNames(1:end) = TS_FILLTERED_FILLED.Properties.VariableNames;
writetable(De_Seasonal,'De_Seasonal_FilterFilled.csv')

%% deYearSeason
% Rolling mean detrender
YearTrend = zeros(size(TSI_data));
j = 1;
for i = 1:size(TSI_data,1)
    % i
    x = TSI_data(i,:);
    n = length(x);
    xs = zeros(size(x));
    for j = 1:floor(n/2)
        xs(j) = nanmean(x(1:j + floor(n/2)));
    end
    for j = ceil(n/2)+1:n
        xs(j) = nanmean(x(j-floor(n/2)+1:end));
    end
    for j = ceil(n/2)
        xs(j) = nanmean(x(j-floor(n/2)+1 : j+floor(n/2)));
    end
    YearTrend(i,:) = xs;
end
Detrend = TSI_data - YearTrend;

% a 3-order harmonic function fit to the data as a deseasoner
harm_freq = 1:3;
Time = datenum(TS_Time_FILLED);
Timestamp = Time(1);

order = 3;
xt = (Time-Timestamp)/365.25;
x_rad = xt * 2 * pi;
nr_indep = order*2 + 2;
Indep = zeros(size(TSI_data,2), nr_indep);
Indep(:,1) = ones(size(Detrend,2),1);
Indep(:,2) = x_rad;

i = 3;
for freq = harm_freq
    cos_freq = cos(x_rad * freq);
    sin_freq = sin(x_rad * freq);
    Indep(:,i) = cos_freq;
    i = i + 1;
    Indep(:,i) = sin_freq;
    i = i + 1;
end

for m = 1:size(Detrend,1)
    [b,bint,r,rint,stats] = regress(Detrend(m,:)',Indep);
    fitted = Indep*b;
    Seasonal(m,:) = fitted';
end

De_Year_Seasonal = [TS_FILLTERED_FILLED.Hylak_id,TS_FILLTERED_FILLED.Lat,TS_FILLTERED_FILLED.Lon,Detrend-Seasonal];
De_Year_Seasonal = array2table(De_Year_Seasonal);
De_Year_Seasonal.Properties.VariableNames(1:end) = TS_FILLTERED_FILLED.Properties.VariableNames;
writetable(De_Year_Seasonal,'De_Year_Seasonal_FilterFilled.csv')

