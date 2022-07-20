function [sign_chans_cluster_perm] = adi_mvpa_select_channels(path2data)

%svm
chans_svm = load([path2data 'svm_sign_channels_0.03_0.15ms.mat']);
chans_svm_comp1 = chans_svm.sign_channels; 
chans_svm_peak = load([path2data 'svm_sign_channels_0.09_0.1ms.mat']);
svm_chan_indices_all = unique([chans_svm_comp1.index chans_svm_peak.sign_channels.index]);
save([path2data 'svm_chan_indices_all.mat'], 'svm_chan_indices_all')
% lda

chans_lda = load([path2data 'lda_sign_channels_0.03_0.15ms.mat']);
chans_lda_comp1 = chans_lda.sign_channels; 
chans_lda_peak = load([path2data 'lda_sign_channels_0.09_0.1ms.mat']);
lda_chan_indices_all = unique([chans_lda_comp1.index chans_lda_peak.sign_channels.index]);
save([path2data 'lda_chan_indices_all.mat'], 'lda_chan_indices_all')

load(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\cluster_permutation_test\stats_comp1_nonaveraged.mat'])
pos_cluster = find(sum(stat_comp1.posclusterslabelmat, 2));
mask_pos_neg = find(sum(stat_comp1.mask, 2));
sign_chans_cluster_perm = stat_comp1.label(mask_pos_neg);
% negative mit eingeschlossen, da sonst zu konservativ, auch wenn p-wert
% von 0.025 nicht erreicht wird; intuitiv würde ich sagen, dass sonst
% informative Kanäle fehlen

%% component 2:
load(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\cluster_permutation_test\stats_comp2_nonaveraged.mat'])
pos_cluster = find(sum(stat_comp2.posclusterslabelmat, 2));
mask_pos_neg = find(sum(stat_comp2.mask, 2));
sign_chans_cluster_perm = stat_comp2.label(mask_pos_neg);














end