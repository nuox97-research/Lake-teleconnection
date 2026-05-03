%%  Here, we take data from 2002 to 2012 as an exemple to eliminate the directional influences in CCM calculations
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analyses
%%% Alongside this workflow for processing 2002-2012 data, 
%%% we directly provide finalized results for the TSI (2016-2024) and FUI datasets in this folder to support downstream simulations
clc 
clear

%% Import data
% here we take data from 2002 to 2012 as an exemple
CCM_Harmonic_deSeason_0212 = readtable('CCMresult_df.csv','ReadVariableNames',true,'VariableNamingRule','preserve');

%% Synchronization
Sync_harm_deSeason_0212 = table();
for i = 1:size(CCM_Harmonic_deSeason_0212,1)
    i
    lib = CCM_Harmonic_deSeason_0212.lib(i);
    tar = CCM_Harmonic_deSeason_0212.target(i);
    Twin = CCM_Harmonic_deSeason_0212(CCM_Harmonic_deSeason_0212.lib == tar & CCM_Harmonic_deSeason_0212.target == lib,:);
    if isempty(Twin)
        Sync_harm_deSeason_0212 = [Sync_harm_deSeason_0212;CCM_Harmonic_deSeason_0212(i,:)];
    else
        Subset = [CCM_Harmonic_deSeason_0212(i,:);Twin];
        if isempty(Sync_harm_deSeason_0212)
            N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
            Sync_harm_deSeason_0212 = [Sync_harm_deSeason_0212;Subset(N,:)];
        elseif isempty(find((Sync_harm_deSeason_0212.lib == tar & Sync_harm_deSeason_0212.target == lib) | (Sync_harm_deSeason_0212.lib == lib & Sync_harm_deSeason_0212.target == tar)))
            N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
            Sync_harm_deSeason_0212 = [Sync_harm_deSeason_0212;Subset(N,:)];
        end
    end
end
save Sync_harm_deSeason_0212 Sync_harm_deSeason_0212


