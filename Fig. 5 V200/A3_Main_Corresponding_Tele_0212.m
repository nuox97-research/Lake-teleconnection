clc
clear

%%
LakeInfo = readtable('LakeGeoInfo.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
Sync_CCM = struct2array(load('Sync_harm_deSeason.mat'));
Sync_CCM_strong = Sync_CCM(Sync_CCM.Optlag_rho>0.4 & Sync_CCM.Dist>3500,:);

Positive_points = readtable('Positive_points_0212.csv','ReadVariableNames',true,'VariableNamingRule','preserve');
Positive_points = Positive_points.Hylak_id;
Negative_points = table2array(LakeInfo(~ismember(LakeInfo.Hylak_id,Positive_points),'Hylak_id'));

%%
for i = 1:height(Sync_CCM_strong)
    lib = Sync_CCM_strong.lib(i);
    tar = Sync_CCM_strong.target(i);

    libloc = find(LakeInfo.Hylak_id == lib);
    tarloc = find(LakeInfo.Hylak_id == tar);

    OLON = LakeInfo.repr_poi_1(tarloc);
    OLAT = LakeInfo.repr_point(tarloc);
    DLON = LakeInfo.repr_poi_1(libloc);
    DLAT = LakeInfo.repr_point(libloc);

    if ismember(lib,Positive_points) & ismember(tar,Positive_points)
        ODLOC_0212(i,:) = [lib,tar,OLON,OLAT,DLON,DLAT,1];
    elseif ismember(lib,Negative_points) & ismember(tar,Negative_points)
        ODLOC_0212(i,:) = [lib,tar,OLON,OLAT,DLON,DLAT,-1];
    else 
        ODLOC_0212(i,:) = [lib,tar,OLON,OLAT,DLON,DLAT,0];
    end
end
ODLOC_0212 = array2table(ODLOC_0212);
ODLOC_0212.Properties.VariableNames = {'lib','target','OLON','OLAT','DLON','DLAT','class'};

tabulate(ODLOC_0212.class);

writetable(ODLOC_0212,'ODLOC_0212.csv')