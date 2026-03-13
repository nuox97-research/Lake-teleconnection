%%  Here, we take data from 2002 to 2012 as an exemple to calculate teleconnection frequencies and strength
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analyses
%%% Alongside this workflow for processing 2002-2012 data, 
%%% we directly provide finalized results for the TSI (2016-2024) and FUI datasets in this folder to support downstream simulations
clc 
clear
close all

%% import data
LakeInfo_0212 = readtable('D:\Lake teleconnection code\Fig. 1 Teleconnection\Geoinformation\LakeGeoInfo_0212.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
LakeInfo_1624 = readtable('D:\Lake teleconnection code\Fig. 1 Teleconnection\Geoinformation\LakeGeoInfo_1624.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
Lakesame = intersect(LakeInfo_0212.Hylak_id,LakeInfo_1624.Hylak_id);

%% 2002-2012
% select the 2002-2012 files
% ensure the R file in the "Moving window CCM" folder should be run first
Filelist_0212 = dir('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Moving window CCM\*.csv');

Filetime0212 = [];
MW_0212 = cell(length(Filelist_0212),1);
Fre_0212 = zeros(length(Filelist_0212),1);
CCMrho_mean_0212 = zeros(length(Filelist_0212),1);
for i = 1:length(Filelist_0212)
    i
    Filenamespc = Filelist_0212(i).name;
    Filetimenum = regexp(Filenamespc,'\d+','match');   
    Filetime0212 = [Filetime0212;[Filetimenum{1,1},'-',Filetimenum{1,2}]];
    CCM_spc = readtable([Filelist_0212(i).folder,'\',Filelist_0212(i).name],'ReadVariableNames',true,'VariableNamingRule','preserve');
    CCM_spc(find(~ismember(CCM_spc.lib,Lakesame) | ~ismember(CCM_spc.target,Lakesame)),:) = [];
    CCM_spc = CCM_spc(CCM_spc.Optlag_rho>0.4 & CCM_spc.Dist>3500,:);
    Sync_ccm = table();
    for k = 1:size(CCM_spc,1)
        lib = CCM_spc.lib(k);
        tar = CCM_spc.target(k);
        Twin = CCM_spc(CCM_spc.lib == tar & CCM_spc.target == lib,:);
        if isempty(Twin)
            Sync_ccm = [Sync_ccm;CCM_spc(k,:)];
        else
            Subset = [CCM_spc(k,:);Twin];
            if isempty(Sync_ccm)
                N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
                Sync_ccm = [Sync_ccm;Subset(N,:)];
            elseif isempty(find((Sync_ccm.lib == tar & Sync_ccm.target == lib) | (Sync_ccm.lib == lib & Sync_ccm.target == tar)))
                N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
                Sync_ccm = [Sync_ccm;Subset(N,:)];
            end
        end
    end

    MW_0212{i,1} = Sync_ccm; 
    Fre_0212(i,1) = size(Sync_ccm,1)/(height(Lakesame)*(height(Lakesame)-1)/2);
    CCMrho_mean_0212(i,1) = mean(Sync_ccm.Optlag_rho);
end

%% 2016-2024
% select the 2016-2024 files
% here is just an example, users should modify file paths to their actual locations
Filelist_1624 = dir('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Moving window CCM\*.csv');

Filetime1624 = [];
MW_1624 = cell(length(Filelist_1624),1);
Fre_1624 = zeros(length(Filelist_1624),1);
CCMrho_mean_1624 = zeros(length(Filelist_1624),1);
for i = 1:length(Filelist_1624)
    i
    Filenamespc = Filelist_1624(i).name;
    Filetimenum = regexp(Filenamespc,'\d+','match');   
    Filetime1624 = [Filetime1624;[Filetimenum{1,1},'-',Filetimenum{1,2}]];
    CCM_spc = readtable([Filelist_1624(i).folder,'\',Filelist_1624(i).name],'ReadVariableNames',true,'VariableNamingRule','preserve');
    CCM_spc(find(~ismember(CCM_spc.lib,Lakesame) | ~ismember(CCM_spc.target,Lakesame)),:) = [];
    CCM_spc = CCM_spc(CCM_spc.Optlag_rho>0.4 & CCM_spc.Dist>3500,:);
    Sync_ccm = table();
    for k = 1:size(CCM_spc,1)
        lib = CCM_spc.lib(k);
        tar = CCM_spc.target(k);
        Twin = CCM_spc(CCM_spc.lib == tar & CCM_spc.target == lib,:);
        if isempty(Twin)
            Sync_ccm = [Sync_ccm;CCM_spc(k,:)];
        else
            Subset = [CCM_spc(k,:);Twin];
            if isempty(Sync_ccm)
                N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
                Sync_ccm = [Sync_ccm;Subset(N,:)];
            elseif isempty(find((Sync_ccm.lib == tar & Sync_ccm.target == lib) | (Sync_ccm.lib == lib & Sync_ccm.target == tar)))
                N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
                Sync_ccm = [Sync_ccm;Subset(N,:)];
            end
        end
    end
    MW_1624{i,1} = Sync_ccm;     
    Fre_1624(i,1) = size(Sync_ccm,1)/(height(Lakesame)*(height(Lakesame)-1)/2);
    CCMrho_mean_1624(i,1) = mean(Sync_ccm.Optlag_rho);
end

%%
Fre_Same_Timefilterfill = [Fre_0212;Fre_1624];
CCMmean_Same_TimefilterFill = [CCMrho_mean_0212;CCMrho_mean_1624];

save MW_0212 MW_0212
save MW_1624 MW_1624
save Fre_Same_Timefilterfill Fre_Same_Timefilterfill
save CCMmean_Same_TimefilterFill CCMmean_Same_TimefilterFill

save Filetime0212 Filetime0212
save Filetime1624 Filetime1624
