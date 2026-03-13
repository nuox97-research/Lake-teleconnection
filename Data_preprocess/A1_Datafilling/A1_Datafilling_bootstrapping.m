%% This m.file is for filling in the data gaps in the time series 
%% using a bootstrapping approach inside a 1-y-long moving window

% This procedure is modified from Gilarranz et al., 2022, PNAS
% Reference: Gilarranz L J, et al., 2022. Regime shifts, trends, and variability 
% of lake productivity at a global scale. PNAS 119, e2116413119

clc
clear

%% Loading data
load TSI_TIMEFILTERED.mat
load TIME_CLEAN.mat

Lake_TSI = TSI_TIMEFILTERED;

Nlakes = size(Lake_TSI,1); % Number of lakes in the dataset
HalfMW = 18; % This goes for a moving window of 37 points --> 370 days (we have 1 datapoint every ~10 days)
BootSize = 1000; % number of randomizations for de bootstrapping

% get time
for i = 1:size(TIME_CLEAN,1)
    timestr = strcat(extractBetween(TIME_CLEAN(i),1,4),'-',extractBetween(TIME_CLEAN(i),5,6),'-',extractBetween(TIME_CLEAN(i),7,8));
    Time(i,1) = datetime(timestr,'InputFormat','uuuu-MM-dd');
end

%% Bootstrapping filling
% create structure to save information about time series
TS_FILLTERED_FILLED = nan(Nlakes,size(Lake_TSI,2)-HalfMW*2-3); % to save filled time series
Hylak_id = Lake_TSI.Hylak_id; % lake ID
Lat = Lake_TSI.Lat; % Latitude
Lon = Lake_TSI.Lon; % Longitude

% to indicate if a lake is discarded or not
DISCARDED = zeros(Nlakes,1);

for f = 1:Nlakes
    % get time series
    TSI_spc = table2array(Lake_TSI(f,4:end));

    % fill the gaps
    FilledTimeseries = nan(length(TSI_spc),1);
    ngaps = sum(isnan(FilledTimeseries));
    ngapsprior = ngaps;
    cont = 0;
    while ngaps > 0
        % bootstrapp data inside a 1yr moving window
        [MeanRSnBT,RandomSampBT,TimeMwbootsBT] = meanmwbootstrap(HalfMW,BootSize,TSI_spc,Time);
        % MeanRSnBT --> this is the bootstrapped time series with the moving window
        % TimeMwbootsBT --> numeric time code that corresponds to each of the values in MeanRSn
        % RandomSampBT=nan(BootSize,Ndays); %store the bootstrapped values (each row is a bootstrapped time series)
        % ALL_RandomSampBT(:,:,f)=RandomSampBT(:,HalfMW+1:length(TSI_spc)-HalfMW);

        % fill the gaps
        Findnan = find(isnan(TSI_spc));
        FilledTimeseries = TSI_spc;
        FilledTimeseries(Findnan) = MeanRSnBT(Findnan);

        % after this step we might still have some gaps in the time series, so we iterate until we have none
        ngaps = sum(isnan(FilledTimeseries(HalfMW+1:size(FilledTimeseries,2)-HalfMW)));
        if ngaps == ngapsprior
            cont = cont+1;
        end
        TSI_spc = FilledTimeseries;
        ngapsprior = ngaps;
        if cont == 5 
            break,
        end
    end
    FilledTimeseries = FilledTimeseries(HalfMW+1:length(TSI_spc)-HalfMW); % discard the first and last half year.
    if sum(isnan(FilledTimeseries)) > 0 %if there is any NaN left in the time series 
        DISCARDED(f,1) = 1; % we discard that lake, continuing the for loop with the next one.
        continue,
    end
    TS_FILLTERED_FILLED(f,:) = FilledTimeseries;
end

% rename the columns
TS_FILLTERED_FILLED = array2table([Hylak_id,Lat,Lon,TS_FILLTERED_FILLED]);
TS_FILLTERED_FILLED.Properties.VariableNames(1:3) = {'Hylak_id','Lat','Lon'};
TS_FILLTERED_FILLED.Properties.VariableNames(4:end) = string(Time(HalfMW+1:length(Time)-HalfMW));

% update the time
TS_Time_FILLED = Time(HalfMW+1:length(Time)-HalfMW);

% discard lakes which still have NaN
TS_FILLTERED_FILLED(logical(DISCARDED),:) = [];

% save time series
save TS_FILLTERED_FILLED TS_FILLTERED_FILLED
save TS_Time_FILLED TS_Time_FILLED % it also works for detrended

