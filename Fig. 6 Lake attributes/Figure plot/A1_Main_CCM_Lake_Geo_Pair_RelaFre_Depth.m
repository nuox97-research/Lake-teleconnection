%% We take lake depth as an exemple
clc
clear
close all

%% Figure parameter
ax1position = [0.1898 0.17 0.6 0.6];
ax2position = [0.1898 0.82 0.6 0.15];
ax3position = [0.8271 0.17 0.15 0.6];
ticksize = 18;
labelsize = 26;
legendsize = 20;

%% Loading data
% CCM
load('D:\Lake teleconnection code\Fig. 1 Teleconnection\Lake teleconnection_CCM\Sync_harm_deSeason_0212.mat')
load('D:\Lake teleconnection code\Fig. 1 Teleconnection\Lake teleconnection_CCM\Sync_harm_deSeason_1624.mat')

Sync_harm_0212 = Sync_harm_deSeason_0212;
Sync_harm_1624 = Sync_harm_deSeason_1624;

Sync_harm = [Sync_harm_0212;Sync_harm_1624];
Sync_harm = Sync_harm(Sync_harm.Optlag_rho>0.4 & Sync_harm.Dist>3500,:);

related_lake = unique([Sync_harm.lib;Sync_harm.target]);

% Lake informations
LakeInfo = readtable('LakeGeoInfo.xlsx','ReadVariableNames',true,'VariableNamingRule','preserve');
LakeInfo = LakeInfo(ismember(LakeInfo.Hylak_id,related_lake),:);

%% For Lake Depth
Depth_lib = zeros(size(LakeInfo,1),1);
Depth_target = zeros(size(LakeInfo,1),1);
parfor i = 1:size(Sync_harm,1)
    lib = Sync_harm.lib(i);
    target = Sync_harm.target(i);

    libloc = find(LakeInfo.Hylak_id == lib);
    targetloc = find(LakeInfo.Hylak_id == target);

    if ~isempty(libloc) & ~isempty(targetloc)
        Depth_lib(i,1) = LakeInfo.Depth_avg(libloc);
        Depth_target(i,1) = LakeInfo.Depth_avg(targetloc);
    else
        Depth_lib(i,1) = NaN;
        Depth_target(i,1) = NaN;
    end

end
Sync_harm.libDepth = Depth_lib;
Sync_harm.tarDepth = Depth_target;
Sync_harm(isnan(Sync_harm.libDepth)|isnan(Sync_harm.tarDepth),:) = [];

Sync_harm.meanDepth = mean([Sync_harm.libDepth,Sync_harm.tarDepth],2);

%% All possible pairs
n_lakes = height(LakeInfo);
all_pairs = nchoosek(LakeInfo.Hylak_id, 2);
all_possible_pairs = zeros(size(all_pairs, 1), 3);
parfor i = 1:size(all_pairs, 1)
    lib = all_pairs(i, 1);
    tar = all_pairs(i, 2);

    libloc = find(LakeInfo.Hylak_id == lib);
    tarloc = find(LakeInfo.Hylak_id == tar);

    all_possible_pairs(i, :) = [lib, tar, ...
        mean([LakeInfo.Depth_avg(libloc),LakeInfo.Depth_avg(tarloc)])];
end

%% 
n_observed_pairs = height(Sync_harm);
n_total_pairs = size(all_possible_pairs, 1);

n_bins = 25;
Depth_all = all_possible_pairs(:, 3);
Depth_sync = Sync_harm.meanDepth;

% the bin edge
[lakenum, bin_edges] = histcounts(log10([Depth_all; Depth_sync]), n_bins+1);
bin_edges(find(lakenum<100)+1) = [];
bin_edges = 10.^bin_edges;
bin_centers = sqrt(bin_edges(1:end-1) .* bin_edges(2:end));

% Count the expected and observed pairs in each bin
expected_counts = zeros(length(bin_edges)-1, 1);
observed_counts = zeros(length(bin_edges)-1, 1);
observed_strength_mean = zeros(length(bin_edges)-1, 1);
observed_strength_std = zeros(length(bin_edges)-1, 1);

for bin_idx = 1:length(bin_edges)-1
    lower_bound = bin_edges(bin_idx);
    upper_bound = bin_edges(bin_idx+1);
    
    % number of all possible pairs in this bin (expected)
    in_bin_all = (Depth_all >= lower_bound) & (Depth_all < upper_bound);
    expected_counts(bin_idx) = (sum(in_bin_all)/ n_total_pairs) * n_observed_pairs;
    
    % number of observed pairs in this bin (observed)
    in_bin_sync = (Depth_sync >= lower_bound) & (Depth_sync < upper_bound);
    observed_counts(bin_idx) = sum(in_bin_sync);
    
    % calculate the mean teleconnection strength
    if expected_counts(bin_idx) > 1 && observed_counts(bin_idx) > 0
        sync_in_bin = Sync_harm.Optlag_rho(in_bin_sync);
        observed_strength_mean(bin_idx) = mean(sync_in_bin);
        observed_strength_std(bin_idx) = std(sync_in_bin);
    else
        observed_strength_mean(bin_idx) = NaN;
        observed_strength_std(bin_idx) = NaN;        
    end
end

% calculate the frequency (observed/expected)
observed_frequency = (observed_counts ./ expected_counts);
observed_frequency(isnan(observed_frequency) | isinf(observed_frequency)) = 0;

%% Confidence interval (bootstraping)
n_bootstrap = 1000;  % Bootstrap time
alpha = 0.05; 
bootstrap_frequency = zeros(length(bin_centers), n_bootstrap);
bootstrap_strength = zeros(length(bin_centers), n_bootstrap);

% bootstrapping
fprintf('Bootstrapping (n=%d):\n', n_bootstrap);
progress_step = ceil(n_bootstrap / 100);
rng(123);
for b = 1:n_bootstrap
    if mod(b, progress_step) == 0
        fprintf('Finished %d/%d (%.0f%%)\n', b, n_bootstrap, b/n_bootstrap*100);
    end
    
    % 1. Resample observed teleconnection pairs with replacement
    bootstrap_indices = randsample(n_observed_pairs, n_observed_pairs, true);
    bootstrap_sample = Sync_harm(bootstrap_indices, :);
    
    % 2. Resample all teleconnection pairs with replacement
    bootstrap_expected_indices = randsample(n_total_pairs, n_total_pairs, true);
    bootstrap_all = Depth_all(bootstrap_expected_indices);
    
    % 3. Calculate the statistics for each bin.
    for bin_idx = 1:length(bin_centers)
        lower_bound = bin_edges(bin_idx);
        upper_bound = bin_edges(bin_idx+1);
        
        % observed
        in_bin_obs = (bootstrap_sample.meanDepth >= lower_bound) & ...
                     (bootstrap_sample.meanDepth < upper_bound);
        obs_count = sum(in_bin_obs);
        
        % expected
        in_bin_exp = (bootstrap_all >= lower_bound) & ...
                     (bootstrap_all < upper_bound);
        exp_count = sum(in_bin_exp) / n_total_pairs * n_observed_pairs;
        
        % frequency
        if exp_count > 1
            bootstrap_frequency(bin_idx, b) = obs_count / exp_count;
        else
            bootstrap_frequency(bin_idx, b) = 0;
        end
        
        % mean teleconnection strength
        if obs_count > 0
            strength_vals = bootstrap_sample.Optlag_rho(in_bin_obs);
            bootstrap_strength(bin_idx, b) = mean(strength_vals);
        else
            bootstrap_strength(bin_idx, b) = NaN;
        end
    end
end

% Calculate the confidence interval for frequency and teleconnection strength
frequency_ci_lower = zeros(length(bin_centers), 1);
frequency_ci_upper = zeros(length(bin_centers), 1);

strength_ci_lower = zeros(length(bin_centers), 1);
strength_ci_upper = zeros(length(bin_centers), 1);

% percentile
lower_percentile = 5;
upper_percentile = 95;

for bin_idx = 1:length(bin_centers)
    % confidence interval for frequency
    freq_vals = bootstrap_frequency(bin_idx, :);
    freq_vals = freq_vals(~isnan(freq_vals) & ~isinf(freq_vals));
    
    if ~isempty(freq_vals)
        frequency_ci_lower(bin_idx) = prctile(freq_vals, lower_percentile);
        frequency_ci_upper(bin_idx) = prctile(freq_vals, upper_percentile);
    end
    
    % confidence interval for teleconnection strength
    strength_vals = bootstrap_strength(bin_idx, :);
    strength_vals = strength_vals(~isnan(strength_vals));
    
    if ~isempty(strength_vals)
        strength_ci_lower(bin_idx) = prctile(strength_vals, lower_percentile);
        strength_ci_upper(bin_idx) = prctile(strength_vals, upper_percentile);
    end
end

%% Figure
valididx = (observed_frequency>0);

figure1 = figure('Color',[1 1 1],'Position',[649 342 629 498]);
axes1 = axes('Parent',figure1,...
    'Position',[0.133921565168676 0.149759031716121 0.732180717039437 0.759176711255766]);
hold(axes1,'on');

yyaxis left;
fill([bin_centers(valididx) fliplr(bin_centers(valididx))],...
     [frequency_ci_lower(valididx)' fliplr(frequency_ci_upper(valididx)')],...
     [0 158 181]/255, 'FaceAlpha',0.3, 'EdgeColor','none'); % 蓝色半透明填充
plot(bin_centers(valididx), observed_frequency(valididx),'-','LineWidth',2,'MarkerSize',25,'Color','#009EB5');
yline(1,'--','Color',[.5 .5 .5],'LineWidth',1.5)
set(gca,'Fontname','Arial','Fontsize', ticksize, ...
    'Ycolor','#009EB5','XScale','log','XTick',[1 10 100 1000 1e4 1e5 1e6])
xlim([0.25 400])
ylim([0.2 1.5])
xlabel('Mean lake depth of lake pairs (m)','FontSize',labelsize)
ylabel('Relative occurrence frequency','FontSize',labelsize)

yyaxis right;
pa = patch([bin_centers(valididx) fliplr(bin_centers(valididx))], ...
    [strength_ci_lower(valididx)' fliplr(strength_ci_upper(valididx)')],'r');
pa.FaceColor = 'k'; 
pa.FaceAlpha = 0.3;  
pa.LineStyle = 'none';  
hold on
plot(bin_centers(valididx), observed_strength_mean(valididx),'-','LineWidth',2,'Color','k'); % '#FF0000'
ylabel('Mean teleconnection strength','FontSize',labelsize); 
set(gca,'Ycolor','k') %#FF0000
ylim([0 1])


