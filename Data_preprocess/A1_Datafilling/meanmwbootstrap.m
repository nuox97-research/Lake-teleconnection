function [MeanRSnBT,RandomSampBT,TimeMwbootsBT] = meanmwbootstrap(HalfMW,BootSize,TimesSeries,Time)
% WE RE-CALCULATE TSI INSIDE A MOVING WINDOW USING BOOTSTRAPPING
  
%INPUT:
  % - HalfMW = 18 --> This goes for a moving window of 37 points --> 370 días | HalfMW = 9 --> moving window half a year | HalfMW = 36 --> 2 years
  % - BootSize: number of randomizations for de bootstrapping
  % - TimesSeries: de-trended time series
  % - Time: time in datetime format

%   meanmwbootstrap(HalfMW,BootSize,cTSI,Time);
%   meanmwbootstrap(HalfMW,BootSize,TimesSeries,Time)

  MWsize = (HalfMW*2)+1; %Size of the moving window
  Ndays = size(TimesSeries,2);
  ConfiVsTime = nan(Ndays,2); %95% confidence intervals around the VARIANCE.
  RandomSampBT = nan(BootSize,Ndays); %store the bootstrapped values
  for d = HalfMW+1:Ndays-HalfMW % 19-338
    Chunk = TimesSeries(d-HalfMW:d+HalfMW);
    NDataChunk = MWsize-nnz(isnan(Chunk));   %number of values in Chunk that doesn't contain nan
    %Bootstrapping
    for r=1:BootSize
      RandChunk = randsample(Chunk,NDataChunk);  %select "NDataChunk" number of elements in Chunk following a uniformly random distribution
      RandomSampBT(r,d) = nanmean(RandChunk);
    end
    %get 95% confidence intervals and store them
    Confi = quantile(RandomSampBT(:,d),[0.0250 0.975]);
    ConfiVsTime(d,:)=Confi;
  end
  
  %averag value (one point per date).
  MeanRSnBT = mean(RandomSampBT)'; 
  TimeMwbootsBT = datenum(Time);
  
%   figure
%   plot(TimeMwbootsBT,MeanRSnBT,'LineWidth',1.5)
%   %white background
%   set(gcf,'color',[1 1 1]);
%   set(gcf,'inverthardcopy','off');
%   %font size
%   set(gca,'fontsize',12,'fontname','Helvetica');
%   set(gca,'fontweight','normal');
%   %axis length
%   %xlim([cdata{1,1}(1,1) cdata{1,1}(size(cdata{1,1},1),1)])
%   %tick length
%   set(gca,'ticklength',0.5*get(gca,'ticklength'))
%   %axis line width
%   %  set(gca,'LineWidth',1.5);
%   %figure size
%   set(gcf,'units','centimeters')
%   set(gcf,'Units','centimeters','Position',[0,0,21,10])   %los dos �ltimos d�gitos, el primero es la anchura, y el segundo la altura  
%   %axis lablels
%   xlabel('Time','fontsize',12,'fontname','Helvetica')
%   ylabel('Bootstrapped TSI','fontsize',12,'fontname','Helvetica')