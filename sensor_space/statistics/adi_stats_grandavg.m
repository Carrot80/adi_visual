function [] = adi_stats(path2data)

load([path2data filesep 'grandavg_like.mat'])
load([path2data filesep 'grandavg_dislike.mat'])

load([path2data filesep 'avg_subjects_dislike.mat'])
load([path2data filesep 'avg_subjects_like.mat'])

cfg=[];
grandavg = ft_timelockgrandaverage(cfg, avg_like, avg_dislike);

cfg = [];
cfg.method  = 'within';
cfg.keepindividual = 'yes';
grandavg_like_group = ft_timelockgrandaverage(cfg, avg_subjects_like{:});
grandavg_dislike_group = ft_timelockgrandaverage(cfg, avg_subjects_dislike{:});


%% global field power:

cfg = [];
cfg.method = 'power';
[gmf_like] = ft_globalmeanfield(cfg, avg_like);
[gmf_dislike] = ft_globalmeanfield(cfg, avg_dislike);

figure('Renderer', 'painters', 'Position', [300 300 800 200])
plot(gmf_like.time, gmf_like.avg, 'r', 'LineWidth',1)
axis tight
hold on
plot(gmf_dislike.time, gmf_dislike.avg, 'b', 'LineWidth',1)
legend({'like'; 'dislike'})
title('global field power ')
legend('boxoff') 
savefig([path2data 'global_field_power_avg.fig' ]);

%% global field amplitude
cfg = [];
cfg.method = 'amplitude'; %'amplitude', power
[gfa_like] = ft_globalmeanfield(cfg, avg_like);
[gfa_dislike] = ft_globalmeanfield(cfg, avg_dislike);
[gfa] = ft_globalmeanfield(cfg, grandavg);

% plot global field amplitude
figure('Renderer', 'painters', 'Position', [300 300 800 200])
plot(gfa_like.time, gfa_like.avg, 'r', 'LineWidth',1)
axis tight
hold on
plot(gfa_dislike.time, gfa_dislike.avg, 'b', 'LineWidth',1)
% hold on
legend({'like'; 'dislike'})
% title('global field amplitude ')
set(gcf,'color','w');
box off
legend({'like';  'dislike'}, 'boxoff')
legend('boxoff') 
xlabel('Time [s]')
ylabel('Global Field Amplitude [t]')
savefig([path2data 'global_field_amplitude_avg.fig' ]);

figure('Renderer', 'painters', 'Position', [300 300 800 200])
plot(gfa.time, gfa.avg, 'k', 'LineWidth',1)
axis tight
xlabel('Time [s]')
ylabel('Global Field Amplitude [t]')
title('global field amplitude like and dislike together')
set(gcf,'color','w');
box off
savefig([path2data 'global_field_amplitude_both_conditions.fig' ]);
 
%% topoplot für Publikation (like, dislike zusammen)

% baseline:
cfg = [];
cfg.xlim  = [-0.5 -0.03];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
cfg.layout = '4D248.lay';
figure; ft_topoplotER(cfg, grandavg);
colorbar; 
set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );
savefig([path2data 'Topoplot_baseline.fig' ]);

% comp 1:
cfg = [];
cfg.xlim  = [0.03 0.15];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
cfg.layout = '4D248.lay';
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );
savefig([path2data 'Topoplot_comp1.fig' ]);

% comp2:
cfg = [];
cfg.xlim  = [0.15 0.56];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
cfg.layout = '4D248.lay';
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );
savefig([path2data 'Topoplot_comp2.fig' ]);

% comp 3:
cfg = [];
cfg.xlim  = [0.56 0.73];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
cfg.layout = '4D248.lay';
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );
savefig([path2data 'Topoplot_comp3.fig' ]);

%%
% calculate t-statistic at each time point:
% dependent samples ttest laut tutorial 

cfg = [];
cfg.showlabels  = 'yes';
cfg.layout      = '4D248_helmet.mat';
figure; ft_multiplotER(cfg, grandavg_like_group, grandavg_dislike_group)

% define the parameters for the statistical comparison
cfg = [];
% cfg.channel     = 'all';
cfg.latency     = [0.03 0.15];%[0.09 0.09]; %[0.56 0.74];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'max';%'FDR';
cfg.numrandomization = 10000;

Nsub = numel(avg_subjects_dislike);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
stat_max_distribution_comp1 = ft_timelockstatistics(cfg,avg_subjects_like{:} ,avg_subjects_dislike{:});   % don't forget the {:}!
cfg.latency     = [0.15 0.56];
stat_max_distribution_comp2 = ft_timelockstatistics(cfg,avg_subjects_like{:} ,avg_subjects_dislike{:});   % don't forget the {:}!
save([path2data 'stat_max_distribution_comp1.mat' ], 'stat_max_distribution_comp1');
save([path2data 'stat_max_distribution_comp2.mat' ], 'stat_max_distribution_comp2');

cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectlike_vs_dislike = ft_math(cfg, avg_like, avg_dislike);

cfg=[];
cfg.latency = [0.03 0.15];
raweffectlike_vs_dislike_comp1 = ft_selectdata(cfg, raweffectlike_vs_dislike);
cfg.latency = [0.15 0.56];
raweffectlike_vs_dislike_comp2 = ft_selectdata(cfg, raweffectlike_vs_dislike);
cfg.latency = [0.56 0.73];
raweffectlike_vs_dislike_comp3 = ft_selectdata(cfg, raweffectlike_vs_dislike);

cfg = [];
% cfg.style     = 'blank';
cfg.layout    = '4D248.lay';
cfg.highlight = 'on';
cfg.zlim=[-0.15 0.15]
cfg.highlightchannel = find(stat_max_distribution_comp1.mask);
% cfg.comment   = 'yes';
cfg.colorbar = 'yes';
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp1); set(gcf,'color','w');
savefig([path2data 'stat_max_distribution_comp1.fig' ]);

cfg.highlightchannel = find(stat_max_distribution_comp2.mask);
cfg.zlim=[-0.15 0.15]
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp2); set(gcf,'color','w');% 
savefig([path2data 'stat_max_distribution_comp2.fig' ]);
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp3) ;set(gcf,'color','w');% 
savefig([path2data 'stat_max_distribution_comp3.fig' ]);



cfg = [];
% cfg.style     = 'blank';
cfg.layout    = '4D248.lay';
set(gcf,'color','w');
cfg.highlight = 'on';
cfg.highlightchannel = find(stat_max_distribution.mask);
cfg.comment   = 'yes';
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp1) % 
title('analytic without mcc')
savefig([path2data 'stats_analytic_without_mcc.fig' ]);
cfg.highlightchannel = find(stat_uncorrected.mask);
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp1) % 
savefig([path2data 'stats_analytic_FDR.fig' ]);



%% cluster permutation test:

cfg_neighb        = [];
cfg_neighb.method = 'distance';
cfg_neighb.template = 'bti248_neighb.mat';
cfg_neighb.layout = '4D248.lay';
neighbours        = ft_prepare_neighbours(cfg_neighb);

cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'depsamplesT'; % use the independent samples T-statistic as a measure to
                               % evaluate the effect at the sample level
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
                               % permutation distribution.
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
                               % required for a selected sample to be included
                               % in the clustering algorithm (default=0).
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.05;               % alpha level of the permutation test
cfg.numrandomization = 5000;      % number of draws from the permutation distribution
cfg.correcttail = 'alpha';
subj = size(grandavg_like_group.individual,1);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
     
cfg.neighbours    = neighbours;                             
% cfg.channel       = {'MEG'};     
cfg.latency       = [0.03 0.15];     
     % time interval over which the experimental

% cfg.avgovertime = 'yes';                                 % conditions must be compared (in seconds)   
cfg.avgovertime = 'no'
[stat_comp1] = ft_timelockstatistics(cfg, grandavg_like_group, grandavg_dislike_group);
cfg.latency       = [0.15 0.56]; 
[stat_comp2] = ft_timelockstatistics(cfg, grandavg_like_group, grandavg_dislike_group);
cfg.latency       = [0.56 0.74]; 
[stat_comp3] = ft_timelockstatistics(cfg, grandavg_like_group, grandavg_dislike_group);

%%

mask = find(sum(stat_comp1.mask, 2))

% stat.posclusters(1)
% stat.negclusters(1)

% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
% if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.
stat_comp3.cfg.alpha = 0.05; % sonst zu konservativ

% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos_cluster_pvals = [stat_comp1.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat_comp1.cfg.alpha);
pos = ismember(stat_comp1.posclusterslabelmat, pos_signif_clust);


% and now for the negative clusters...
neg_cluster_pvals = [stat_comp1.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat_comp1.cfg.alpha);
neg = ismember(stat_comp1.negclusterslabelmat, neg_signif_clust);

% pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
% neg = stat.negclusterslabelmat == 1;

timestep = 0.005; % timestep between time windows for each subplot (in seconds)
sampling_rate = 1017.25; % Data has a temporal resolution of 300 Hz
sample_count = length(stat_comp1.time);
% number of temporal samples in the statistics object
j = [0.03:timestep:0.15]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count+timestep*sampling_rate]; % t

[i1,i2] = match_str(raweffectlike_vs_dislike_comp1.label, stat_comp1.label);

for k = 1:25
   subplot(5,5,k);
   cfg = [];
   cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
%    cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel is in a to-be-plotted cluster, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).

   % Next, check which channels are in the clusters over the
   % entire time interval of interest.
   pos_int = zeros(numel(raweffectlike_vs_dislike_comp1.label),1);
   neg_int = zeros(numel(raweffectlike_vs_dislike_comp1.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

   cfg.highlight   = 'on';
   % Get the index of the to-be-highlighted channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.comment     = 'xlim';
   cfg.commentpos  = 'title';
   cfg.layout      = '4D248.lay';
   cfg.interactive = 'no';
   cfg.figure      = 'gca'; % plots in the current axes, here in a subplot
   ft_topoplotER(cfg, raweffectlike_vs_dislike_comp1);
end
%% figure comp1:


posclusterslabelmat = stat_comp1.posclusterslabelmat;
posclusterslabelmat(find(posclusterslabelmat~=1))=0;
indx_cluster1 = find(posclusterslabelmat==1);
pos_cluster = find(sum(posclusterslabelmat,2));

negclusterslabelmat = stat_comp1.negclusterslabelmat;
negclusterslabelmat(find(negclusterslabelmat~=1))=0;
indx_cluster1 = find(negclusterslabelmat==1);
neg_cluster = find(sum(negclusterslabelmat,2));

cfg = [];
cfg.xlim=[0.03 0.15]; 
% cfg.style     = 'blank';
cfg.zlim=[-0.12 0.12];
cfg.highlight = 'on';
cfg.highlightchannel = [pos_cluster;neg_cluster];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
savefig([path2data 'cluster_permutation_test_comp1.fig' ]);

%% comp 2 figure:

mask = find(sum(stat_comp2.mask, 2))

stat_comp2.posclusters(1)
stat_comp2.negclusters(1)

% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
% if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.
stat_comp2.cfg.alpha = 0.05; % sonst zu konservativ

% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos_cluster_pvals = [stat_comp2.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat_comp2.cfg.alpha);
pos = ismember(stat_comp2.posclusterslabelmat, pos_signif_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat_comp2.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat_comp2.cfg.alpha);
neg = ismember(stat_comp2.negclusterslabelmat, neg_signif_clust);

% pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
% neg = stat.negclusterslabelmat == 1;

timestep = 0.005; % timestep between time windows for each subplot (in seconds)
sampling_rate = 1017.25; % Data has a temporal resolution of 300 Hz
sample_count = length(stat_comp2.time);
% number of temporal samples in the statistics object
j = [0.15:timestep:0.56]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count+timestep*sampling_rate]; % t

[i1,i2] = match_str(raweffectlike_vs_dislike_comp2.label, stat_comp2.label);

for k = 1:25
   subplot(5,5,k);
   cfg = [];
   cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
%    cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel is in a to-be-plotted cluster, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).

   % Next, check which channels are in the clusters over the
   % entire time interval of interest.
   pos_int = zeros(numel(raweffectlike_vs_dislike_comp2.label),1);
   neg_int = zeros(numel(raweffectlike_vs_dislike_comp2.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

   cfg.highlight   = 'on';
   % Get the index of the to-be-highlighted channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.comment     = 'xlim';
   cfg.commentpos  = 'title';
   cfg.layout      = '4D248.lay';
   cfg.interactive = 'no';
   cfg.figure      = 'gca'; % plots in the current axes, here in a subplot
   ft_topoplotER(cfg, raweffectlike_vs_dislike_comp2);
end





%% 
cfg = [];
cfg.xlim=[0.03 0.15]; 
% cfg.style     = 'blank';
% cfg.zlim=[-0.12 0.12];

cfg.highlight = 'on';
cfg.highlightchannel = [mask];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
savefig([path2data 'cluster_permutation_test_comp1.fig' ]);



%% figure comp2 
posclusterslabelmat = stat_comp2.posclusterslabelmat;
indx_cluster4 = find(stat_comp2.posclusterslabelmat==4);
posclusterslabelmat(indx_cluster4)=0;
indx_cluster3 = find(stat_comp2.posclusterslabelmat==3);
posclusterslabelmat(indx_cluster3)=0;
indx_cluster2 = find(stat_comp2.posclusterslabelmat==2);
posclusterslabelmat(indx_cluster2)=0;
pos_cluster = find(sum(posclusterslabelmat,2));

negclusterslabelmat = stat_comp2.negclusterslabelmat;
negclusterslabelmat(find(negclusterslabelmat~=1))=0;
indx_cluster1 = find(negclusterslabelmat==1);
neg_cluster = find(sum(negclusterslabelmat,2));

cfg = [];
cfg.xlim=[0.15 0.56]; 
% cfg.style     = 'blank';
cfg.zlim=[-0.12 0.12];
cfg.highlight = 'on';
cfg.highlightchannel = [pos_cluster; neg_cluster];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
% savefig([path2data 'cluster_permutation_test_comp2.fig' ]);

sign_chans_cluster_perm_comp2=stat_comp2.label(unique([pos_cluster; neg_cluster]));

%% comp 3: n.s.

%%

stat_comp2.cfg.alpha = 0.025; % da drei tests gerechnet werden (1 pro componente)
pos_cluster_pvals = [stat_comp2.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat_comp2.cfg.alpha);
pos = ismember(stat_comp2.posclusterslabelmat, pos_signif_clust);
neg_cluster_pvals = [stat_comp2.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat_comp2.cfg.alpha);
neg = ismember(stat_comp2.negclusterslabelmat, neg_signif_clust);


cfg = [];
cfg.xlim=[0.15 0.56]; 
% cfg.xlim=[0.56 0.74]; 
% cfg.style     = 'blank';
cfg.zlim=[-0.12 0.12];
pos_int = zeros(numel(raweffectlike_vs_dislike_comp1.label),1);
pos_int(i1) = all(pos(i2, 1), 2);  
neg_int = zeros(numel(raweffectlike_vs_dislike.label),1);
neg_int(i1) = all(neg(i2, 1), 2); 
cfg.highlight = 'on';
cfg.highlightchannel = find(pos_int);
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
savefig([path2data 'cluster_permutation_test_comp2.fig' ]);

%% comp 3:

stat_comp3.cfg.alpha = 0.025; % da drei tests gerechnet werden (1 pro componente)
pos_cluster_pvals = [stat_comp3.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat_comp3.cfg.alpha);
pos = ismember(stat_comp3.posclusterslabelmat, pos_signif_clust);
neg_cluster_pvals = [stat_comp3.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat_comp3.cfg.alpha);
neg = ismember(stat_comp3.negclusterslabelmat, neg_signif_clust);


cfg = [];
cfg.xlim=[0.56 0.74]; 
% cfg.style     = 'blank';
cfg.zlim=[-0.12 0.12];
pos_int = zeros(numel(raweffectlike_vs_dislike_comp1.label),1);
pos_int(i1) = all(pos(i2, 1), 2);  
neg_int = zeros(numel(raweffectlike_vs_dislike.label),1);
neg_int(i1) = all(neg(i2, 1), 2); 
% cfg.highlight = 'on';
cfg.highlightchannel = find(pos_int);
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
savefig([path2data 'cluster_permutation_test_comp3.fig' ]);












%%
cfg = [];
cfg.style     = 'blank';
cfg.layout    = '4D248.lay';
cfg.highlight = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.colorbar           = 'yes';
cfg.highlightchannel = find(stat.mask);
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp1) % 
title('analytic without mcc')
savefig([path2data 'stats_analytic_without_mcc.fig' ]);
cfg.highlightchannel = find(stat_FDR.mask);
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp1) % 
savefig([path2data 'stats_analytic_FDR.fig' ]);

end