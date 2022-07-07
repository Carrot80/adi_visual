
function [] = adi_ERF_freqanalysis(subjectpath, path2data, condition)

switch condition
    case 'like'
    trialcount_like = [];
    avg_subjects_like = struct([]);
    for ii = [2:length(subjectpath)] % wichtig, adi_04 rausnehmen, da andere Kanalsortierung
        %% like
        filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2data 'Neu_Like*.mat']);
        all_runs_subject = struct([]);
        for kk = 1:length(filename)
            load ([filename(kk).folder filesep filename(kk).name])
            cfg = [];   
            cfg.demean = 'yes';
            cfg.baselinewindow  = [-0.5 -0.030];
            data = ft_preprocessing(cfg, cleanMEG_interp);

            %% z-transformation:
             [data] = subfun_ztransform(data);

             % realignment per subject:    
             if kk==1
                 gradfile = data.grad;
             else
                 cfg=[];
                 cfg.template = gradfile;
                 cfg.inwardshift = 1;
                 load([subjectpath(ii).folder filesep subjectpath(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\vol\vol.mat'])
                 cfg.headmodel = vol;
                 data = ft_megrealign(cfg, data);
             end
            

             %% 
             all_runs_subject{length(all_runs_subject)+1} = data;
             trialcount_like(length(trialcount_like)+1) = numel(data.trial);
             clear cleanMEG_interp 

        end
            clear vol gradfile
            close all

        %% appenddata 
        
        count_runs = length(all_runs_subject);
        cfg = [];
        like_all_runs = ft_appenddata(cfg, all_runs_subject{:});
        
        %% dislike:
        all_runs_subject_dislike = [];
        filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2data 'Neu_Dislike*.mat']);
        all_runs_subject = struct([]);
        for kk = 1:length(filename)
            load ([filename(kk).folder filesep filename(kk).name])
            cfg = [];   
            cfg.demean = 'yes';
            cfg.baselinewindow  = [-0.5 -0.030];
            data = ft_preprocessing(cfg, cleanMEG_interp);

            %% z-transformation:
             [data] = subfun_ztransform(data);

             % realignment per subject:    
             if kk==1
                 gradfile = data.grad;
             else
                 cfg=[];
                 cfg.template = gradfile;
                 cfg.inwardshift = 1;
                 load([subjectpath(ii).folder filesep subjectpath(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\vol\vol.mat'])
                 cfg.headmodel = vol;
                 data = ft_megrealign(cfg, data);
             end
            

             %% 
             all_runs_subject_dislike{length(all_runs_subject_dislike)+1} = data;
             trialcount_dislike(length(trialcount_like)+1) = numel(data.trial);
             clear cleanMEG_interp 

        end
            clear vol gradfile
            close all

        %% appenddata 
        
        count_runs = length(all_runs_subject);
        cfg = [];
        dislike_all_runs = ft_appenddata(cfg, all_runs_subject_dislike{:});
        
        cfg=[]; 
        like_dislike_all_runs = ft_appenddata (cfg, like_all_runs, dislike_all_runs)
        path2save = [subjectpath(ii).folder filesep subjectpath(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\freqanalysis\toi_bl_0.01_0.5\'];
        if ~exist(path2save, 'dir')
            mkdir(path2save)
        end
        
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = 'MEG';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:2:46;                         % analysis 2 to 30 Hz in steps of 2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
        cfg.toi          = -0.5:0.01:0.5;    
%         cfg.keeptrials = 'yes'; % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
        TFRhann_like = ft_freqanalysis(cfg, like_all_runs);
        TFRhann_dislike = ft_freqanalysis(cfg, dislike_all_runs);
        TFRhann_like_dislike = ft_freqanalysis(cfg, like_dislike_all_runs);
        save([path2save 'TFRhann_dislike.mat'], 'TFRhann_dislike')
        save([path2save 'TFRhann_like.mat'], 'TFRhann_like')
        save([path2save 'TFRhann_like_dislike.mat'], 'TFRhann_like_dislike')
    
        
%         cfg = [];
%         cfg.baseline     = [-0.5 0];
%         cfg.baselinetype = 'absolute';
% %         cfg.zlim         = [-2.5e-27 2.5e-27];
%         cfg.showlabels   = 'yes';
%         cfg.layout       = '4D248_helmet.mat';
%         figure
%         ft_multiplotTFR(cfg, TFRhann_like_dislike);
%         savefig([path2save 'TFRhann_like_dislike.fig'])
        
        cfg = [];
        cfg.zlim = 'maxabs';
        cfg.xlim = [-0.5 1];
        cfg.baseline = [-0.5 0];
        % cfg.zlim = [-0.3 0.3];
        cfg.maskstyle    = 'saturation';
        cfg.baselinetype = 'relchange'; %'relchange'; % absolute
        cfg.colormap = 'jet';
        cfg.layout       = '4d248_helmet.mat';
        % cfg.channel      = TRF_wave_all.label(ind_left) ;% 'A132';
        cfg.interactive  = 'yes';
        cfg.layout       = '4d248_helmet.mat';
        cfg.title = 'TFRhann_like_dislike';
        ft_singleplotTFR(cfg, TFRhann_like_dislike); 
        savefig([path2save 'TFRhann_like_dislike.fig'])
        %% 

%         cfg = [];
%         cfg.channel    = 'MEG';
%         cfg.method     = 'wavelet';
%         cfg.width      = 7; % 7 cycles better than 5 cycles
%         cfg.output     = 'pow';
%         cfg.foi        = 1:2:46;
%         cfg.toi        = -0.5:0.05:1.5;
%         cfg.keeptrials = 'yes';
%         TRF_wavelet_like_7_cycles = ft_freqanalysis(cfg, like_all_runs);
%         TRF_wavelet_dislike_7_cycles = ft_freqanalysis(cfg, dislike_all_runs);
%         TRF_wave_like_dislike_7_cycles = ft_freqanalysis(cfg, like_dislike_all_runs);
%         
% 
%         cfg = [];
%         cfg.zlim = 'maxabs';
%         cfg.xlim = [-0.5 1];
%         cfg.baseline = [-0.5 0];
%         cfg.zlim = [-0.3 0.3];
%         cfg.maskstyle    = 'saturation';
%         cfg.baselinetype = 'relchange'; %'relchange'; % absolute
%         cfg.colormap = 'jet';
%         cfg.layout       = '4d248_helmet.mat';
% %         cfg.channel      = TRF_wave_all.label(ind_left) ;% 'A132';
%         cfg.interactive  = 'yes';
%         cfg.layout       = '4d248_helmet.mat';
%         cfg.title = 'TRF_wave_like_dislike_7_cycles';
%         ft_singleplotTFR(cfg, TRF_wave_like_dislike_7_cycles); 
%         savefig([path2save 'TRF_wave_like_dislike.fig'])
%         
%         save([path2save 'TRF_wavelet_like_7_cycles.mat'], 'TRF_wavelet_like_7_cycles')
%         save([path2save 'TRF_wavelet_dislike_7_cycles.mat'], 'TRF_wavelet_dislike_7_cycles')
%         save([path2save 'TRF_wave_like_dislike_7_cycles.mat'], 'TRF_wave_like_dislike_7_cycles')
    
        
        
%         cfg = [];
%         cfg.baseline     = [-0.5 0];
%         cfg.baselinetype = 'absolute';
% %         cfg.zlim         = [-2e-25 2e-25];
%         cfg.showlabels   = 'yes';
%         cfg.layout       = '4D248_helmet.mat';
%         cfg.colorbar     = 'yes';
%         figure
%         ft_multiplotTFR(cfg, TRF_wavelet_like)
%         ft_multiplotTFR(cfg, TRF_wavelet_dislike)
%         
%         
%         cfg = [];
%         cfg.latency          = 'all';
%         cfg.frequency        = 'all';
%         cfg.method           = 'montecarlo';
%         cfg.statistic        = 'ft_statfun_depsamplesT';
%         cfg.correctm         = 'cluster';
%         cfg.clusteralpha     = 0.05;
%         cfg.clusterstatistic = 'maxsum';
%         cfg.minnbchan        = 2;
%         cfg.tail             = 0;
%         cfg.clustertail      = 0;
%         cfg.alpha            = 0.025;
%         cfg.numrandomization = 500;
%         % prepare_neighbours determines what sensors may form clusters
%         cfg_neighb.method    = 'distance';
%         cfg_neighb.layout = '4D248.lay';
%         cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, TRF_wavelet_like);
% 
%         design = zeros(1,size(TRF_wavelet_dislike.powspctrm,1) + size(TRF_wavelet_like.powspctrm,1));
%         design(1,1:size(TRF_wavelet_dislike.powspctrm,1)) = 1;
%         design(1,(size(TRF_wavelet_dislike.powspctrm,1)+1):(size(TRF_wavelet_dislike.powspctrm,1)+...
%         size(TRF_wavelet_like.powspctrm,1))) = 2;
% 
%         cfg.design           = design;
%         cfg.ivar             = 1;
% 
%         [stat] = ft_freqstatistics(cfg, TRF_wavelet_dislike, TRF_wavelet_like);


    end

    






end















% %%
% 
% load([path2data 'avg_subjects_dislike.mat' ])
% load([path2data 'avg_subjects_like.mat' ])
% 
% 
% cfg=[];
% cfg.keepindividual = 'yes';
% avg_subjects_like = ft_timelockgrandaverage(cfg, avg_subjects_like{:})
% avg_subjects_dislike = ft_timelockgrandaverage(cfg, avg_subjects_dislike{:})
% 
% cfg = [];
% cfg.toilim = [0 .5];
% activation = ft_redefinetrial(cfg, avg_subjects_like);
% activation.grad = template.mean_grad;
% 
% cfg = [];
% cfg.toilim = [-0.5 0];
% baseline = ft_redefinetrial(cfg, avg_subjects_like);
% baseline.time = activation.time;
% baseline.grad = template.mean_grad;
% 
% cfg_neighb = [];
% cfg_neighb.planarmethod = 'sincos';
% % cfg_neighb.layout = '4D248.lay';
% % prepare_neighbours determines with what sensors the planar gradient is computed
% cfg_neighb.method    = 'distance';
% cfg_neighb.neighbours       = ft_prepare_neighbours(cfg_neighb, activation);
% activation_planar = ft_megplanar(cfg_neighb, activation);
% baseline_planar   = ft_megplanar(cfg_neighb, baseline);
% 
% cfg = [];
% cfg.output = 'pow';
% cfg.channel = 'all';
% cfg.method = 'mtmconvol';
% cfg.taper = 'hanning';
% cfg.foi = 17;
% cfg.toi = [0:0.05:0.5];
% cfg.t_ftimwin = 4./cfg.foi; %7 cycles
% cfg.keeptrials = 'yes';
% 
% freq_hanning_activation = ft_freqanalysis(cfg, activation);
% freq_hanning_baseline  = ft_freqanalysis(cfg, baseline);
% 
% cfg = [];
% cfg.channel          = {'all'};
% cfg.latency          = [0 0.5];
% cfg.method           = 'montecarlo';
% cfg.frequency        = 17;
% cfg.statistic        = 'ft_statfun_actvsblT'; %'ft_statfun_depsamplesT'; % ''; %
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 2;
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 1000;
% % prepare_neighbours determines what sensors may form clusters
% cfg_neighb.method    = 'distance';
% cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freq_hanning_activation);
% 
% ntrials = size(freq_hanning_activation.powspctrm,1);
% design  = zeros(2,2*ntrials);
% design(1,1:ntrials) = 1;
% design(1,ntrials+1:2*ntrials) = 2;
% design(2,1:ntrials) = [1:ntrials];
% design(2,ntrials+1:2*ntrials) = [1:ntrials];
% 
% cfg.design   = design;
% cfg.ivar     = 1;
% cfg.uvar     = 2;
% 
% [stat] = ft_freqstatistics(cfg, freq_hanning_activation, freq_hanning_baseline);
% 
% cfg = [];
% cfg.alpha  = 0.025;
% cfg.parameter = 'stat';
% % cfg.zlim   = [-4 4];
% cfg.layout = '4D248_helmet.mat';
% ft_clusterplot(cfg, stat);  % funktioniert







end


function [output_data_zscore] = subfun_ztransform(input_data)

baseline_zero_samples = nearest(input_data.time{1,1},-0.03);
baseline_dur_samples = nearest(input_data.time{1,1}, -0.5);

output_data_zscore = input_data;

for pp=1:length(input_data)
    output_data_zscore(pp).trial = [];
    output_data_zscore(pp).trial = cell(1, length(input_data(pp).trial));
    for kk = 1:length(input_data(pp).trial)
        output_data_zscore(pp).trial{kk} = zeros(size(input_data(pp).trial{kk},1),size(input_data(pp).trial{kk},2));
        for oo = 1:size(input_data(pp).trial{kk},1)   
            M = mean(input_data(pp).trial{kk}(oo, baseline_dur_samples:baseline_zero_samples));
            STD = std(input_data(pp).trial{kk}(oo,baseline_dur_samples:baseline_zero_samples));
            output_data_zscore(pp).trial{kk}(oo,:) = (input_data(pp).trial{kk}(oo,:)- M)/STD;
            clearvars M STD
        end
    end
end

clear input_data


end
