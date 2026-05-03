%%  Here, we take data from 2002 to 2012 as an exemple to show how to caluculate the network indices
%%% The computational procedures remain consistent for both TSI (2016-2024) and FUI analyses
clc
clear
close all

%% Import data
dinfo = dir('D:\Lake teleconnection code\Fig. 3 Moving window\Moving window\Moving window CCM\CCM*.csv');  
filenum = size(dinfo,1);

% extract lakes with data available in both the 2002–2012 and 2016–2024 periods
LakeInfo_0212 = readtable('D:\Lake teleconnection code\Fig. 1 Teleconnection\Geoinformation\LakeGeoInfo_0212.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
LakeInfo_1624 = readtable('D:\Lake teleconnection code\Fig. 1 Teleconnection\Geoinformation\LakeGeoInfo_1624.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
Lakesame = intersect(LakeInfo_0212.Hylak_id,LakeInfo_1624.Hylak_id);
LakeInfo = LakeInfo_0212(ismember(LakeInfo_0212.Hylak_id,Lakesame),:); % get the lake information

%% Calculate the degree centrality of the lake network
Deg = [];
Rhomatrix_cell_limited = cell(filenum,1);
Dist_cell_limited = cell(filenum,1);
parfor i = 1:filenum 
    filename = [dinfo(i).folder,'\',dinfo(i).name];  
    Table = readtable(filename,'ReadVariableNames',true,'VariableNamingRule','preserve');
    Table = Table(ismember(Table.lib,Lakesame)&ismember(Table.target,Lakesame),:);
    Table = Table(Table.Optlag_rho>0.4 & Table.Dist>3500,:); % only consider lake pairs with teleconnection strength > 0.4 and distances > 3500km

    % get the sync influences  
    Sync_harm = table();
    for m = 1:size(Table,1)
        lib = Table.lib(m);
        tar = Table.target(m);
        Twin = Table(Table.lib == tar & Table.target == lib,:);
        if isempty(Twin)
            Sync_harm = [Sync_harm;Table(m,:)];
        else
            Subset = [Table(m,:);Twin];
            if isempty(Sync_harm)
                N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
                Sync_harm = [Sync_harm;Subset(N,:)];
            elseif isempty(find((Sync_harm.lib == tar & Sync_harm.target == lib) | (Sync_harm.lib == lib & Sync_harm.target == tar)))
                N = find(Subset.Optlag_rho == max(Subset.Optlag_rho),1);
                Sync_harm = [Sync_harm;Subset(N,:)];
            end
        end
    end

    % get the matrix of CCM
    Rhomatrix = zeros(length(Lakesame),length(Lakesame));
    for k = 1:size(Sync_harm,1)
        lib = Sync_harm.lib(k);
        tar = Sync_harm.target(k);
        libnum = find(LakeInfo.Hylak_id == lib);
        tarnum = find(LakeInfo.Hylak_id == tar);
        Rhomatrix(tarnum,libnum) = Sync_harm.Optlag_rho(k);
    end
    Trilow = tril(Rhomatrix);
    Transup = Trilow';
    Trilup = triu(Rhomatrix);
    Translow = Trilup';
    Rhomatrix_undirect = Rhomatrix + Transup + Translow;
    Rhomatrix_cell_limited{i,1} = Rhomatrix_undirect;

    % calculate the degree centrality
    Deg(:,i) = sum(Rhomatrix_undirect,2); % others is the reason of i
end
save Rhomatrix_cell_limited Rhomatrix_cell_limited
save Deg Deg

