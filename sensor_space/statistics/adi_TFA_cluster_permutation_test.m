function [] = adi_ERF_cluster_permutation_test(subjectpath, path2data)


TRF_wave_dislike = struct([]);
TRF_wave_like = struct([]);

%   for ii = 2:length(subjectpath) 
% 
%         load ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2data filesep 'TRF_wavelet_dislike_7_cycles.mat'])
%         load ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2data filesep 'TRF_wavelet_like_7_cycles.mat'])
% 
%         TRF_wave_dislike{ii-1}=TRF_wavelet_dislike_7_cycles;
%         TRF_wave_like{ii-1}=TRF_wavelet_like_7_cycles;
% 
%         clear TRF_wavelet_dislike_7_cycles TRF_wavelet_like_5_cycles
%   end

    for ii = 2:length(subjectpath) 

        load ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2data filesep 'TFRhann_like.mat'])
        load ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2data filesep 'TFRhann_dislike.mat'])

        TRF_wave_dislike{ii-1}=TFRhann_dislike;
        TRF_wave_like{ii-1}=TFRhann_like;

        clear TFRhann_like TFRhann_dislike
  end
  
 cfg =[];
 cfg.keepindividual = 'yes';
 TRF_wave_dislike_all = ft_freqgrandaverage(cfg, TRF_wave_dislike{:});
 TRF_wave_like_all = ft_freqgrandaverage(cfg, TRF_wave_like{:});
 TRF_wave_like_dislike_all = ft_freqgrandaverage(cfg, TRF_wave_like{:}, TRF_wave_dislike{:});

 channelselection = {'A45'; 'A46'; 'A47'; 'A48'; 'A49'; 'A50'; 'A51'; 'A52'; 'A53';....
    'A70'; 'A71'; 'A72'; 'A73'; 'A74'; 'A75'; 'A76'; 'A77'; 'A78'; 'A79'; 'A80';  ...
    'A99'; 'A100'; 'A101'; 'A102'; 'A103'; 'A104'; 'A105'; 'A106'; 'A107'; 'A108'; 'A109'; 'A110'; 'A111'; ...
    'A131'; 'A132'; 'A133'; 'A134'; 'A135'; 'A136'; 'A137'; 'A138'; 'A139'; 'A140'; 'A141'; 'A142'; 'A143'; ...
    'A159'; 'A160'; 'A161'; 'A162'; 'A163'; 'A164'; 'A165'; 'A166'; 'A167'; 'A168'; 'A169'; 'A170'; ...
    'A180'; 'A181'; 'A182'; 'A183'; 'A184'; 'A185'; 'A186'; 'A187'; 'A188'; 'A189'; 'A190'; 'A191'; 'A192'; ...
     'A214'; 'A215'; 'A216'; 'A217'; 'A218'; 'A219'; 'A220'; 'A221'; 'A222'; 'A223'; 'A224'; 'A225'; 'A226'; ...
     'A234';  'A235';  'A236';  'A237';  'A238';  'A239';  'A240';  'A241';  'A242';  'A243'; 'A244'};
cfg = [];
cfg.zlim = 'maxabs';
% cfg.xlim = [-0.5 1];
cfg.baseline = [-0.5 0];
cfg.zlim = [-0.3 0.3];
cfg.maskstyle    = 'saturation';
cfg.baselinetype = 'relchange'; %'relchange'; % absolute
cfg.colormap = 'jet';
cfg.layout       = '4d248_helmet.mat';
cfg.channel      =  channelselection;
cfg.interactive  = 'yes';
cfg.layout       = '4d248_helmet.mat';
cfg.title = 'dislike';
ft_singleplotTFR(cfg, TRF_wave_dislike_all); 
set(gcf,'Position',[100 100 400 300])
cfg.title = 'like';
ft_singleplotTFR(cfg, TRF_wave_like_all); 
set(gcf,'Position',[100 100 400 300])
cfg.title = 'like and dislike';
ft_singleplotTFR(cfg, TRF_wave_like_dislike_all); 
set(gcf,'Position',[100 100 400 300])

cfg=[];
cfg.layout = '4D248.lay';
ft_layoutplot(cfg, TRF_wave_like)

figure;
scatter(lay.pos(:,1), lay.pos(:,2))
hold on
scatter(lay.pos(:,1), lay.pos(:,2), lay.label)


%% 
cfg = [];
cfg.channel          = channelselection%{'all'};
% cfg.latency          = [0.03 0.15];
cfg.method           = 'montecarlo';
% cfg.frequency        = 17;
cfg.statistic        = 'ft_statfun_depsamplesT'; %'ft_statfun_depsamplesT'; % ''; %
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'distance';
cfg_neighb.template      = 'bti248_neighb.mat';
gradfile = load('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\template_gradfile\template.mat');
TRF_wave_dislike_all.grad = gradfile.template.mean_grad;
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, TRF_wave_dislike_all);

ntrials = size(TRF_wave_like_all.powspctrm,1);
design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials) = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

% design = zeros(1,size(freqFIC_planar_cmb.powspctrm,1) + size(freqFC_planar_cmb.powspctrm,1));
% design(1,1:size(freqFIC_planar_cmb.powspctrm,1)) = 1;
% design(1,(size(freqFIC_planar_cmb.powspctrm,1)+1):(size(freqFIC_planar_cmb.powspctrm,1)+...
% size(freqFC_planar_cmb.powspctrm,1))) = 2;

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;

[stat] = ft_freqstatistics(cfg, TRF_wave_like_all, TRF_wave_dislike_all);

% plotting:
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
% cfg.zlim   = [-4 4];
cfg.layout = '4d248_helmet.mat';
ft_clusterplot(cfg, stat);








end