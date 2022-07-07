function [] = adi_stats_grandavg_clearly_rated(path2data)

load([path2data filesep 'grandavg_like.mat'])
load([path2data filesep 'grandavg_dislike.mat'])

load([path2data filesep 'avg_dislike_subjects.mat'])
load([path2data filesep 'avg_like_subjects.mat'])

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
plot(gmf_like.time, gmf_like.avg, 'r', 'LineWidth',1))
axis tight
hold on
plot(gmf_dislike.time, gmf_dislike.avg, 'b', 'LineWidth',1))
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
cfg.xlim  = [0.56 0.74];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
cfg.layout = '4D248.lay';
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );
savefig([path2data 'Topoplot_comp3.fig' ]);

%%
% define the parameters for the statistical comparison
cfg = [];
cfg.channel     = 'all';
cfg.latency     = [0.03 0.15]; %[0.56 0.74];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';%'FDR';

Nsub = 28;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,avg_subjects_like{:} ,avg_subjects_dislike{:});   % don't forget the {:}!

cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectlike_vs_dislike = ft_math(cfg, avg_like, avg_dislike);

cfg=[];
cfg.latency = [0.03 0.15];
raweffectlike_vs_dislike_comp1 = ft_selectdata(cfg, raweffectlike_vs_dislike)

cfg = [];
% cfg.style     = 'blank';
cfg.layout    = '4D248.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'yes';
figure; ft_topoplotER(cfg, raweffectlike_vs_dislike_comp1) % 
title('analytic without mcc')
savefig([path2data 'stats_analytic_without_mcc.fig' ]);

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
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 5000;      % number of draws from the permutation distribution

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
% cfg.latency       = [0.15 0.56];      % time interval over which the experimental

cfg.avgovertime = 'yes';                                 % conditions must be compared (in seconds)    
[stat] = ft_timelockstatistics(cfg, grandavg_like_group, grandavg_dislike_group);

% stat.posclusters(1)
% stat.negclusters(1)

% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
% if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.
stat.cfg.alpha = 0.025; % da drei tests gerechnet werden (1 pro componente)

% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);


% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);

% pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
% neg = stat.negclusterslabelmat == 1;

% timestep = 0.01; % timestep between time windows for each subplot (in seconds)
% sampling_rate = 1017.25; % Data has a temporal resolution of 300 Hz
% sample_count = length(stat.time);
% number of temporal samples in the statistics object
% j = [0.03:timestep:0.15]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*sampling_rate:sample_count+timestep*sampling_rate]; % t

[i1,i2] = match_str(raweffectlike_vs_dislike_comp1.label, stat.label);

% plot avg:
 
cfg = [];
cfg.xlim=[0.03 0.15]; 
% cfg.zlim=[-0.15 0.15];
pos_int = zeros(numel(raweffectlike_vs_dislike_comp1.label),1);
pos_int(i1) = all(pos(i2, 1), 2);  
% neg_int = zeros(numel(raweffectlike_vs_dislike.label),1);
% neg_int(i1) = all(neg(i2, 1), 2); 
cfg.highlight = 'on';
cfg.highlightchannel = find(pos_int);
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
savefig([path2data 'cluster_permutation_test.fig' ]);

end