function [] = adi_permutation_compute_p_values(subjectdir, config)
% Gespräch mit Stefan: keine p-Werte mitteln
% classifier nicht mitteln
% compute_mean_accuracy_per_subj(subjectdir, config)
% compute_max_accuracy_per_subj(subjectdir, config)
% compute_mean_accuracy_across_subj(subjectdir, config)

%% sort mean accuracy across subjects without pvalues

% sort_mean_accuracy_across_subj_without_pvalues(subjectdir, config)


%% compute mean accuracy with p-values 
acc_svm_all_subj = [];

for ii = 1:length(subjectdir)
    
    if config.timerange(1,1) == 0.09 && exist(['E:\Arbeit\adidas\data_analysis\visual_stimuli\single_subjects\' subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x1_trainfold_only_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_sample_orig.mat'])
            load(['E:\Arbeit\adidas\data_analysis\visual_stimuli\single_subjects\' subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x1_trainfold_only_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_sample_orig.mat'])
    elseif config.timerange(1,1) == 0.09 && ~exist(['E:\Arbeit\adidas\data_analysis\visual_stimuli\single_subjects\' subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x1_trainfold_only_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_sample_orig.mat'])
           load(['E:\Arbeit\adidas\data_analysis\visual_stimuli\single_subjects\' subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x1_trainfold_only_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) '_ms.mat'], 'perf')
    else
        load(['E:\Arbeit\adidas\data_analysis\visual_stimuli\single_subjects\' subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x1_trainfold_only_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) '_ms.mat'], 'perf')
    end
    accuracy_lda = mean(perf.lda.accuracy,1);
    accuracy_svm = mean(perf.svm.accuracy,1);

    acc_lda_all_subj(ii,:)=accuracy_lda; 
    acc_svm_all_subj(ii,:)=accuracy_svm; 

    balldesign =config.balldesign;
    
    if config.timerange(1,1) == 0.09 && exist([config.path2subjects subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\mvpa\searchlight\permutation_distributation' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_svm_sample_orig.mat'])
        perf_perm_svm = load([config.path2subjects subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\mvpa\searchlight\permutation_distributation' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_svm_sample_orig.mat'])
        
    else
        perf_perm_svm = load([config.path2subjects subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\mvpa\searchlight' '\permutation_distribution_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_svm.mat']);
    end
%     perf_perm_lda = load(['E:\Arbeit\adidas\data_analysis\visual_stimuli\single_subjects\' subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight' '\permutation_distribution_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) '.mat']);
    perf_perm_lda = load([config.path2subjects subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\mvpa\searchlight' '\permutation_distribution_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms.mat']);

    %   svm:     
    perm_svm_balldesigns = perf_perm_svm.perf_perm.svm.(balldesign{1}).accuracy; 
    for kk = 2:length(balldesign)
        perm_svm_balldesigns = (perm_svm_balldesigns+perf_perm_svm.perf_perm.svm.(balldesign{kk}).accuracy)./2;
    end
    
    prc_95_svm = prctile(perm_svm_balldesigns, 95);
    sign_chans_svm_subj = find(accuracy_svm >=prc_95_svm);
    sign_chans_svm_all_subj.([subjectdir(ii).name]).channels = sign_chans_svm_subj;
    
    switch ii
        case 1
            sign_chans_svm = find(accuracy_svm >=prc_95_svm);
        otherwise 
            indx_sign_chans_svm = find(accuracy_svm >=prc_95_svm);
            sign_chans_svm = cat(2,sign_chans_svm, indx_sign_chans_svm);
    end
    
    clear accuracy_svm prc_95_svm indx_sign_chans_svm

    % lda:
    
    perm_lda_balldesigns = perf_perm_lda.perf_perm.lda.(balldesign{1}).accuracy; 
    for kk = 2:length(balldesign)
        perm_lda_balldesigns = (perm_lda_balldesigns+perf_perm_lda.perf_perm.lda.(balldesign{kk}).accuracy)./2;
    end
    
    prc_95_lda = prctile(perm_lda_balldesigns, 95);
    sign_chans_lda_subj = find(accuracy_lda >=prc_95_lda);
    sign_chans_lda_all_subj.([subjectdir(ii).name]).channels = sign_chans_lda_subj;
    
    switch ii
        case 1
            sign_chans_lda = find(accuracy_lda >=prc_95_lda);
        otherwise 
            indx_sign_chans_lda = find(accuracy_lda >=prc_95_lda);
            sign_chans_lda = cat(2,sign_chans_lda, indx_sign_chans_lda);
    end
end

mean_accuracy_svm = mean(acc_svm_all_subj,1);
acc_svm = [];
acc_svm.accuracy = mean(acc_svm_all_subj,1)'
acc_svm.label = perf.features;
acc_svm.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
ft_topoplotER(cfg, acc_svm); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\svm_mean_accuracy_all_subj' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms.fig']);

mean_accuracy_lda = mean(acc_lda_all_subj,1);
acc_lda = [];
acc_lda.accuracy = mean(acc_lda_all_subj,1)'
acc_lda.label = perf.features;
acc_lda.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
ft_topoplotER(cfg, acc_lda); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\lda_mean_accuracy_all_subj' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms.fig']);




sign_chans_count_svm = tabulate(sign_chans_svm);  
stat = [];
stat.mvpa.perf = sign_chans_count_svm(:,2)./29*100;
stat.accuracy = sign_chans_count_svm(:,2)./29*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(unique(sign_chans_svm));
cfg.title = ['svm - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\svm_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms.fig']);
sign_channels = [];
sign_channels.channels = perf.features(unique(sign_chans_svm));
sign_channels.index = unique(sign_chans_svm);
save(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\svm_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms.mat'], 'sign_channels')

sign_chans_count_lda = tabulate(sign_chans_lda);  
stat.mvpa.perf = sign_chans_count_lda(:,2)./29*100;
stat.accuracy = sign_chans_count_lda(:,2)./29*100;
cfg.highlightchannel = perf.features(unique(sign_chans_lda));
cfg.title = ['lda - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\lda_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms.fig']);
sign_channels = [];
sign_channels.channels = perf.features(unique(sign_chans_lda));
sign_channels.index = unique(sign_chans_lda);
save(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\lda_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms.mat'], 'sign_channels')

%% plot channels that are significant in at least 2 subjects
indx_one_svm=find(sign_chans_count_svm(:,2)==1);
sign_channelnumbers_svm = unique(sign_chans_svm);
[~, location_index_svm] = ismember(indx_one_svm,sign_channelnumbers_svm)
sign_channelnumbers_svm(location_index_svm)=[];

stat = [];
stat.mvpa.perf = sign_chans_count_svm(:,2)./29*100;
stat.accuracy = sign_chans_count_svm(:,2)./29*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(sign_channelnumbers_svm);
cfg.title = ['svm - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=2 sig. channels in subjects']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\svm_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_2_subj.fig']);


indx_one_lda=find(sign_chans_count_lda(:,2)==1);
sign_channelnumbers_lda = unique(sign_chans_lda);
[~, location_index_lda] = ismember(indx_one_lda,sign_channelnumbers_lda)
sign_channelnumbers_lda(location_index_lda)=[];
stat.mvpa.perf = sign_chans_count_lda(:,2)./29*100;
stat.accuracy = sign_chans_count_lda(:,2)./29*100;
cfg.highlightchannel = perf.features(sign_channelnumbers_lda);
cfg.title = ['lda - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=2 sig. channels in subjects'];
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\lda_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_2_subj.fig']);

%% plot channels that are significant in at least 3 subjects
% svm:
indx_sign_xsubj_svm=find(sign_chans_count_svm(:,2)>2); % mind x subjects sollen signifikanten kanal haben
stat = [];
stat.mvpa.perf = sign_chans_count_svm(:,2)./29*100;
stat.accuracy = sign_chans_count_svm(:,2)./29*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(indx_sign_xsubj_svm);
cfg.title = ['svm - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=3 sig. channels in subjects']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\svm_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_3_subj.fig']);

%lda:
indx_sign_xsubj_lda=find(sign_chans_count_lda(:,2)>2); % mind x subjects sollen signifikanten kanal haben
stat = [];
stat.mvpa.perf = sign_chans_count_lda(:,2)./29*100;
stat.accuracy = sign_chans_count_lda(:,2)./29*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(indx_sign_xsubj_lda);
cfg.title = ['lda - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=3 sig. channels in subjects']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\lda_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_3_subj.fig']);
%% plot channels that are significant in at least 4 subjects
% svm:
indx_sign_xsubj_svm=find(sign_chans_count_svm(:,2)>3); % mind x subjects sollen signifikanten kanal haben
stat = [];
stat.mvpa.perf = sign_chans_count_svm(:,2)./29*100;
stat.accuracy = sign_chans_count_svm(:,2)./29*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(indx_sign_xsubj_svm);
cfg.title = ['svm - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=4 sig. channels in subjects']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\svm_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_4_subj.fig']);
indx_sign_xsubj_svm=find(sign_chans_count_svm(:,2)>4); % mind x subjects sollen signifikanten kanal haben
cfg.highlightchannel = perf.features(indx_sign_xsubj_svm);
cfg.title = ['svm - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=5 sig. channels in subjects']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\svm_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_5_subj.fig']);



%lda:
indx_sign_xsubj_lda=find(sign_chans_count_lda(:,2)>3); % mind x subjects sollen signifikanten kanal haben
stat = [];
stat.mvpa.perf = sign_chans_count_lda(:,2)./29*100;
stat.accuracy = sign_chans_count_lda(:,2)./29*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(indx_sign_xsubj_lda);
cfg.title = ['lda - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=4 sig. channels in subjects']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\lda_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_4_subj.fig']);
indx_sign_xsubj_lda=find(sign_chans_count_lda(:,2)>4);
cfg.highlightchannel = perf.features(indx_sign_xsubj_lda);
cfg.title = ['lda - ' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) ' >=5 sig. channels in subjects']; 
ft_topoplotER(cfg, stat); 
set(gcf,'color','w');
savefig(['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\searchlight\lda_sign_channels_' num2str(config.timerange(1)) '_' num2str(config.timerange(2)) 'ms_at_least_5_subj.fig']);

% für Publikation: mind 3-4 Probanden sollten signifikanten Kanal haben

%% Stefan: Accuracy-Werte in eine Rangfolge bringen und schauen, ob es "Abbruchkante" gibt: war nicht wirklich der Fall

acc_all_classifiers_all_subj = mean(mean_acc_all_classifiers_all_subj);

[sort_acc, chan_location] = sort(acc_all_classifiers_all_subj, 'descend')

figure;
scatter(1:248, sort_acc); % ersten 47 Kanäle nehmen

stat = [];
stat.mvpa.perf = acc_all_classifiers_all_subj'*100;
stat.accuracy = acc_all_classifiers_all_subj'*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(chan_location(1:47));
cfg.title = ['significant channels in subjects']; 
set(gcf,'color','w');
figure; ft_topoplotER(cfg, stat); 


%%



 % 2. calculation of p value with FDR-correction
p = nan(1, numel(perm_timerange));
for tt = perm_timerange
    [indx, acc_taller_orig] = find(acc_perm_mean_balldesigns(:,tt) >  acc_mean_balldesigns_orig_data(tt));
    p(tt-(perm_timerange(1)-1)) = (sum(acc_taller_orig)+1)./(1000+1);  % p=(b+1)/(m+1) => see Dima et al(2018), Human Brain Mapping
end
p_ = nan(1, size(acc_perm_mean_balldesigns,2));
p_(perm_timerange) = p;
sign_p = p_<0.05;

p_all_subjects(ii,:) = p(1,:);

% confidence intervals:
for kk = 1:length(acc_perf_balldesigns)
   CI_orig(kk).lda = ci(acc_perf_balldesigns(:,kk));
end
for kk = 1:length(CI_orig)
    ci_orig_lower(kk) = CI_orig(kk).lda(1);
    ci_orig_upper(kk) = CI_orig(kk).lda(2);
end

for kk = 1:size(acc_perm_mean_balldesigns,2)
   CI_perm(kk).lda = ci(acc_perm_mean_balldesigns(:,kk));
end

for kk = 1:length(CI_perm)
    ci_perm_lower(kk) = CI_perm(kk).lda(1);
    ci_perm_upper(kk) = CI_perm(kk).lda(2);
end
    
figure; hold on;
plot(time, acc_mean_balldesigns_orig_data)
plot(time, ci_orig_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_orig_upper, 'color', [0.5 0.5 0.5])
plot(time, nanmean(acc_perm_mean_balldesigns))
plot(time, sign_p*0.4, 'r*'  )
ylim([0.1 1])
title('accuracy uncorrected')
xlabel('time')
ylabel('accuracy')


transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_upper,fliplr(ci_orig_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

filled_perm = [ci_perm_lower,fliplr(ci_perm_upper)];
fillhandle = fill(xpoints,filled_perm,'r');%plot the data

transparency_color = 0.5;
transparency_edge = 0.5;
edge = 'r';
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

plot(time, nanmean(acc_perm_mean_balldesigns), 'r')
plot(time, ci_perm_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_perm_upper, 'color', [0.5 0.5 0.5])
    

savefig([subjectdir(ii).folder filesep subjectdir(ii).name filesep config.subfolder '\accuracy_orig_data_uncorrected.fig'])
close 
% FDR-correction:
p(2,:) = perm_timerange;
[sort_p, indx_p] = sort(p(1,:)); % 
rank_p = 1:sum(numel(sort_p));%182;
p_adjusted = nan(1,numel(rank_p));
p_adjusted(end) = sort_p(end);  % größter p-Wert kommt ans Ende 
%     
for tt = fliplr(rank_p(1:end-1))
    temp = sort_p(tt)*(numel(sort_p)/rank_p(tt));
    if temp < p_adjusted(rank_p(tt+1))
        p_adjusted(tt) = temp;
    elseif temp > p_adjusted(rank_p(tt+1))
        p_adjusted(tt) = p_adjusted(tt+1);
    else
        p_adjusted(tt) = temp;
    end
    clear temp
end
    
[~, old_order] = sort(indx_p);
p_adjusted_ = p_adjusted(old_order);

p_adjusted_ = [nan(1, perm_timerange(1)-1) p_adjusted_  nan(1,385-numel(nan(1, perm_timerange(1)-1))-numel(p_adjusted_) )];
indx_p_sign = find(p_adjusted_<0.05);
p_sign = nan(1, 385);
p_sign(indx_p_sign) = 1;

p_adjusted_all_subj(ii, :) = p_adjusted;


 % plot significance as stars:
figure; hold on;
plot(time, acc_mean_balldesigns_orig_data)
plot(time, ci_orig_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_orig_upper, 'color', [0.5 0.5 0.5])
plot(time, nanmean(acc_perm_mean_balldesigns), 'r')
plot(time, ci_perm_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_perm_upper, 'color', [0.5 0.5 0.5])
plot(time, p_sign*0.73, 'r*'  )
transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_upper,fliplr(ci_orig_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

ylim([0.1 1])
title('accuracy FDR corrected')
xlabel('time')
ylabel('accuracy')
savefig([subjectdir(ii).folder filesep subjectdir(ii).name filesep config.subfolder '\accuracy_orig_data_fdr_ci.fig'])
close  
clear acc_mean_balldesigns_orig_data acc_perf_balldesigns acc_perm_mean_balldesigns perf perf_perm p_adjusted CI_orig ci_perm_lower ci_perm_upper CI_perm ci_orig_lower ci_orig_upper p


% figure;
% plot(time, nanmean(acc_mean_balldesigns_orig_data_all_subjects(:,:))) 
% hold on
% plot(time, nanmean(acc_perm_mean_balldesigns_all_subjects))
% for kk = 1:length(acc_mean_balldesigns_orig_data_all_subjects)
%    CI_orig_all_subjects(kk).lda = ci(acc_mean_balldesigns_orig_data_all_subjects(:,kk));
% end
% for kk = 1:length(CI_orig_all_subjects)
%     ci_orig_all_subj_lower(kk) = CI_orig_all_subjects(kk).lda(1);
%     ci_orig_all_subj_upper(kk) = CI_orig_all_subjects(kk).lda(2);
% end


transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_all_subj_upper,fliplr(ci_orig_all_subj_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color



for kk = 1:size(acc_perm_mean_balldesigns_all_subjects, 2)
   CI_perm_all_subjects(kk).lda = ci(acc_perm_mean_balldesigns_all_subjects(:,kk));
end
for kk = 1:size(acc_perm_mean_balldesigns_all_subjects, 2)
    ci_perm_all_subj_lower(kk) = CI_perm_all_subjects(kk).lda(1);
    ci_perm_all_subj_upper(kk) = CI_perm_all_subjects(kk).lda(2);
end
transparency_color = 0.3;
transparency_edge = 0.3;
filled_perm = [ci_perm_all_subj_lower,fliplr(ci_perm_all_subj_upper)];
fillhandle = fill(xpoints,filled_perm,'r');%plot the data
edge = 'r';
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

% find corrected p_values per subject to figure out the proportion of subjects with significant samples:
p_sign_all_subj = zeros(size(p_adjusted_all_subj, 1), size(p_adjusted_all_subj, 2));
for kk = 1:size(p_adjusted_all_subj, 1)
    indx_p_sign_all_subj = find(p_adjusted_all_subj(kk,:)<0.05);
    p_sign_all_subj(kk,indx_p_sign_all_subj) = 1;
end

perf.accuracy_sign = zeros(1, size(time,2));
 perm_timerange = nearest(time, timerange(1)):1:nearest(time, timerange(2));
for tt = perm_timerange
    prc_ = prctile(acc_perm_mean_balldesigns_all_subjects(:,tt), 95);
    if acc_mean_balldesigns_orig_data_all_subjects(1,tt) > prc_
       perf.accuracy_sign(tt) = 1;
    else
       perf.accuracy_sign(tt) = 0;
    end
    clear prc_
end

% find uncorrected p_values per subject to figure out the proportion of subjects with significant samples:
p_sign_all_subj = zeros(size(p_all_subjects, 1), size(p_all_subjects, 2));
for kk = 1:size(p_all_subjects, 1)
    indx_p_sign_all_subj = find(p_all_subjects(kk,:)<0.05);
    p_sign_all_subj(kk,indx_p_sign_all_subj) = 1;
end
sum(p_sign_all_subj)
temp = nan(1, size(acc_perm_mean_balldesigns_all_subjects,2));
temp(perm_timerange)=sum(p_sign_all_subj);

 % 2. calculation of p value over subjects with FDR-correction
p_all_subj = nan(1, numel(perm_timerange));
for tt = perm_timerange
    [indx_all_subj, acc_taller_orig_all_subj] = find(acc_perm_mean_balldesigns_all_subjects(:,tt) >  nanmean(acc_mean_balldesigns_orig_data_all_subjects(:,tt)));
    p_all_subj(tt-(perm_timerange(1)-1)) = (sum(acc_taller_orig_all_subj)+1)./(1000+1);  % p=(b+1)/(m+1) => see Dima et al(2018), Human Brain Mapping
end
p_all_subj_ = nan(1, size(acc_perm_mean_balldesigns_all_subjects,2));
p_all_subj_(perm_timerange) = p_all_subj;
sign_p_all_subj = p_all_subj_<=0.05;

plot(time, sign_p_all_subj*0.45, 'r*')
ylim([0.38 0.68])
ylabel('accuracy')
xlabel('time')
box off

temp = nan(28, size(acc_perm_mean_balldesigns_all_subjects,2));
temp(:, perm_timerange) = p_all_subjects;
p_comp_1 = temp(:, nearest(time, 0.03): nearest(time, 0.15));

test = zeros(size(p_comp_1,1), size(p_comp_1,2));
for kk = 1:size(p_comp_1,1)
    ind = find(p_comp_1(kk,:)<=0.05);
    test(kk,ind) = 1;
end
sum(test,2)
% ci:
% saved in 'E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\significance_testing_mvpa'

end

function [] = mkfigures(subjectdir, config)

data_balldesign_all_subj = nan(length(subjectdir), 385);
data_perm_balldesign_all_subj = nan(1000, 385);

for ii = 1:length(subjectdir)

    if exist('file', 'var')
        load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\performance_real_data.mat'], 'perf')

        data_balldesign_all_subj(ii,:) = perf.(balldesign{1}).lda.accuracy;

        perf_perm.(balldesign{1}) = load([config.path2subjects subjectdir(ii).name filesep ...
                config.subfolder '\performance_permutation_null_distribution_' balldesign{1} '.mat'], ['perm_' balldesign{1}]);
    else
        load([subjectdir(ii).folder filesep subjectdir(ii).name filesep config.subfolder '\performance_real_data.mat'], 'perf')

        data_balldesign_all_subj(ii,:) = perf.(config.balldesign{1}).lda.accuracy;

        perf_perm.(config.balldesign{1}) = load([subjectdir(ii).folder filesep subjectdir(ii).name filesep config.subfolder '\performance_permutation_null_distribution_' config.balldesign{1} '.mat'], ['perm_' config.balldesign{1}]);
        
    end
    
switch ii
    case 1
        data_perm_balldesign_all_subj =  perf_perm.(config.balldesign{1}).(['perm_' config.balldesign{1}]).accuracy;   
    otherwise
        data_perm_balldesign_all_subj = (data_perm_balldesign_all_subj+perf_perm.(config.balldesign{1}).(['perm_' config.balldesign{1}]).accuracy)./2;
end
end
fsample = 256.001;
fsample_sec = 1/fsample;
time = -0.5:fsample_sec:1;

% confidence intervals:
for kk = 1:length(data_balldesign_all_subj)
   CI_orig(kk).lda = ci(data_balldesign_all_subj(:,kk));
end
for kk = 1:length(data_balldesign_all_subj)
    ci_orig_lower(kk) = CI_orig(kk).lda(1);
    ci_orig_upper(kk) = CI_orig(kk).lda(2);
end

for kk = 1:size(data_perm_balldesign_all_subj,2)
   CI_perm(kk).lda = ci(data_perm_balldesign_all_subj(:,kk));
end

for kk = 1:length(CI_perm)
    ci_perm_lower(kk) = CI_perm(kk).lda(1);
    ci_perm_upper(kk) = CI_perm(kk).lda(2);
end

figure; hold on;
plot(time, nanmean(data_balldesign_all_subj))
plot(time, nanmean(data_perm_balldesign_all_subj))

transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_upper,fliplr(ci_orig_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

filled_perm = [ci_perm_lower,fliplr(ci_perm_upper)];
fillhandle = fill(xpoints,filled_perm,'r');%plot the data

transparency_color = 0.5;
transparency_edge = 0.5;
edge = 'r';
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

perm_timerange = nearest(time, config.timerange(1)):1:nearest(time, config.timerange(2));

p_all_subj = nan(1, numel(perm_timerange));
for tt = perm_timerange
    [indx_all_subj, acc_taller_orig_all_subj] = find(data_perm_balldesign_all_subj(:,tt) >  nanmean(data_balldesign_all_subj(:,tt)));
    p_all_subj(tt-(perm_timerange(1)-1)) = (sum(acc_taller_orig_all_subj)+1)./(1000+1);  % p=(b+1)/(m+1) => see Dima et al(2018), Human Brain Mapping
end
p_all_subj_ = nan(1, size(data_perm_balldesign_all_subj,2));
p_all_subj_(perm_timerange) = p_all_subj;
sign_p_all_subj = p_all_subj_<=0.05;

plot(time, sign_p_all_subj*0.45, 'r*')
ylim([0.3 0.85])
ylabel('accuracy')
xlabel('time')
box off
set(gcf, 'Position', [434 664 762 303])
gcf

% save
 
end
 
function [] = compute_max_accuracy_per_subj(subjectdir, config)

    for ii = 1:length(subjectdir)

        load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x2_trainfold_only_0.03_0.15_ms.mat'], 'perf')
        load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\permutation_distribution_0.03_0.15.mat'])

        max_accuracy_lda = max(perf.lda.accuracy);
        max_accuracy_svm = max(perf.svm.accuracy);
        max_accuracy_logreg = max(perf.logreg.accuracy);
        max_acc_all_classifiers = (max_accuracy_lda+max_accuracy_svm+max_accuracy_logreg)./3;

        for kk = 1:length(balldesign)
            permutation_all_classifiers.(balldesign{kk}) = (perf_perm.lda.(balldesign{kk}).accuracy+perf_perm.svm.(balldesign{kk}).accuracy+perf_perm.logreg.(balldesign{kk}).accuracy)./3;
        end

        for kk = 1:length(balldesign)
            perm_balldesigns(kk,:,:) = permutation_all_classifiers.(balldesign{kk});
        end

        max_perm_balldesign = squeeze(max(perm_balldesigns));
        prc_95_max = prctile(max_perm_balldesign, 95);

        sign_chans_max = find(max_acc_all_classifiers >=prc_95_max);
        stat = [];
        % stat.time = perf.lda.time;
        stat.mvpa.perf = max_acc_all_classifiers';
        stat.accuracy = max_acc_all_classifiers';
        stat.label = perf.features;
        stat.dimord = 'chan';
        cfg              = [];
        cfg.parameter    = 'accuracy';
        cfg.layout       = '4D248_helmet.mat';            
        %    cfg.xlim         = [0, 0];
        cfg.colorbar     = 'yes';
        cfg.highlight = 'on';
        cfg.highlightchannel = perf.features(sign_chans_max);
        cfg.title = ['max significant channels ' subjectdir(ii).name]; 
        figure; ft_topoplotER(cfg, stat); 
        savefig([config.path2subjects subjectdir(ii).name filesep config.subfolder 'sign_sensors_max.fig'] )
    end

end

function [] = compute_mean_accuracy_per_subj(subjectdir, config)
    for ii = 1:length(subjectdir)

        load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x2_trainfold_only_0.03_0.15_ms.mat'], 'perf')
        load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\permutation_distribution_0.03_0.15.mat'])

        % Mittelung aller Bälle und Classifier:
        balldesign = {'gbf'; 'gbs'; 'gbv'; 'ggf'; 'ggs'; 'ggv'; 'rwf'; 'rws'; 'rwv'}; % Reihenfolge der abgespeicherten Balldesigns

        mean_accuracy_lda = mean(perf.lda.accuracy,1);
        mean_accuracy_svm = mean(perf.svm.accuracy,1);
        mean_accuracy_logreg = mean(perf.logreg.accuracy,1);
      
        % mittlere Accuracy über alle bälle nehmen?
        % permutation:
        % lda:
        perm_lda_balldesigns = perf_perm.lda.(balldesign{1}).accuracy; 
        for kk = 2:length(balldesign)
            perm_lda_balldesigns = (perm_lda_balldesigns+perf_perm.lda.(balldesign{kk}).accuracy)./2;
        end
        
        perm_svm_balldesigns = perf_perm.svm.(balldesign{1}).accuracy; 
        for kk = 2:length(balldesign)
            perm_svm_balldesigns = (perm_svm_balldesigns+perf_perm.svm.(balldesign{kk}).accuracy)./2;
        end
        
        perm_logreg_balldesigns = perf_perm.logreg.(balldesign{1}).accuracy; 
        for kk = 2:length(balldesign)
            perm_logreg_balldesigns = (perm_logreg_balldesigns+perf_perm.logreg.(balldesign{kk}).accuracy)./2;
        end
    % 	
    % lda:
        prc_95_lda = prctile(perm_lda_balldesigns, 95);
        sign_chans_lda = find(mean_accuracy_lda >=prc_95_lda);
        
        stat = [];
        % stat.time = perf.lda.time;
        stat.mvpa.perf = mean_accuracy_lda'*100;
        stat.accuracy = mean_accuracy_lda'*100;
        stat.label = perf.features;
        stat.dimord = 'chan';
        cfg              = [];
        cfg.parameter    = 'accuracy';
        cfg.layout       = '4D248_helmet.mat';            
        %    cfg.xlim         = [0, 0];
        cfg.colorbar     = 'yes';
        cfg.highlight = 'on';
        cfg.highlightchannel = perf.features(sign_chans_lda);
        cfg.title = ['significant channels lda ' subjectdir(ii).name]; 
        % figure; 
        ft_topoplotER(cfg, stat); 
        savefig([config.path2subjects subjectdir(ii).name filesep config.subfolder 'sign_sensors_lda.fig'] )
        
         % svm:
        prc_95_svm = prctile(perm_svm_balldesigns, 95);
        sign_chans_svm = find(mean_accuracy_svm >=prc_95_svm);
        
        stat = [];
        % stat.time = perf.lda.time;
        stat.mvpa.perf = mean_accuracy_svm'*100;
        stat.accuracy = mean_accuracy_svm'*100;
        stat.label = perf.features;
        stat.dimord = 'chan';
        cfg              = [];
        cfg.parameter    = 'accuracy';
        cfg.layout       = '4D248_helmet.mat';            
        %    cfg.xlim         = [0, 0];
        cfg.colorbar     = 'yes';
        cfg.highlight = 'on';
        cfg.highlightchannel = perf.features(sign_chans_svm);
        cfg.title = ['significant channels svm ' subjectdir(ii).name]; 
        % figure; 
        ft_topoplotER(cfg, stat); 
        savefig([config.path2subjects subjectdir(ii).name filesep config.subfolder 'sign_sensors_svm.fig'] )
        
        
        % logreg:
        prc_95_logreg = prctile(perm_logreg_balldesigns, 95);
        sign_chans_logreg = find(mean_accuracy_logreg >=prc_95_logreg);
     
        stat = [];
        % stat.time = perf.lda.time;
        stat.mvpa.perf = mean_accuracy_logreg'*100;
        stat.accuracy = mean_accuracy_logreg'*100;
        stat.label = perf.features;
        stat.dimord = 'chan';
        cfg              = [];
        cfg.parameter    = 'accuracy';
        cfg.layout       = '4D248_helmet.mat';            
        %    cfg.xlim         = [0, 0];
        cfg.colorbar     = 'yes';
        cfg.highlight = 'on';
        cfg.highlightchannel = perf.features(sign_chans_logreg);
        cfg.title = ['significant channels logreg ' subjectdir(ii).name]; 
        % figure; 
        ft_topoplotER(cfg, stat); 
        savefig([config.path2subjects subjectdir(ii).name filesep config.subfolder 'sign_sensors_logreg.fig'] )

    end
end

function [] = compute_mean_accuracy_across_subj(subjectdir, config)
    for ii = 1:length(subjectdir)

        load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x2_trainfold_only_0.03_0.15_ms.mat'], 'perf')
        mean_accuracy_lda = mean(perf.lda.accuracy,1);
        mean_accuracy_svm = mean(perf.svm.accuracy,1);
        mean_accuracy_logreg = mean(perf.logreg.accuracy,1);
        mean_acc_all_classifiers = (mean_accuracy_lda+mean_accuracy_svm+mean_accuracy_logreg)./3;
        mean_acc_all_classifiers_all(ii,:)=mean_acc_all_classifiers; 

        load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\permutation_distribution_0.03_0.15.mat'])

        for kk = 1:length(balldesign)
            permutation_all_classifiers.(balldesign{kk}) = (perf_perm.lda.(balldesign{kk}).accuracy+perf_perm.svm.(balldesign{kk}).accuracy+perf_perm.logreg.(balldesign{kk}).accuracy)./3
        end

        perm_mean_classifiers_balldesigns = permutation_all_classifiers.(balldesign{1}); 
        for kk = 2:length(balldesign)
            perm_mean_classifiers_balldesigns = (perm_mean_classifiers_balldesigns+permutation_all_classifiers.(balldesign{kk}))./2;
        end

        switch ii
            case 1
                perm_mean_all_subjects = perm_mean_classifiers_balldesigns;
            otherwise
                perm_mean_all_subjects = (perm_mean_all_subjects+perm_mean_classifiers_balldesigns)./2;
        end
        clear perm_mean_classifiers_balldesigns
    end

    mean_subjects_mean_acc_all_classifiers_all = mean(mean_acc_all_classifiers_all);
    prc_95_mean_subj = prctile(perm_mean_all_subjects, 95);

    sign_chans_mean_all_subj = find(mean_subjects_mean_acc_all_classifiers_all >=prc_95_mean_subj);

    stat = [];
    stat.mvpa.perf = mean_subjects_mean_acc_all_classifiers_all'*100;
    stat.accuracy = mean_subjects_mean_acc_all_classifiers_all'*100;
    stat.label = perf.features;
    stat.dimord = 'chan';
    cfg              = [];
    cfg.parameter    = 'accuracy';
    cfg.layout       = '4D248_helmet.mat';            
    %    cfg.xlim         = [0, 0];
    cfg.colorbar     = 'yes';
    cfg.highlight = 'on';
    cfg.highlightchannel = perf.features(sign_chans_mean_all_subj);
    cfg.title = ['significant channels mean all subjects']; 
    figure; ft_topoplotER(cfg, stat); 
end

function [] = sort_mean_accuracy_across_subj_without_pvalues(subjectdir, config)

% nach Gespräch mit Stefan: accuracy-Werte in eine Rangfolge bringen und
% schauen, ob es abbruchkante gibt ähnlich wie bei PCA
for ii = 1:length(subjectdir)
 load([config.path2subjects subjectdir(ii).name filesep config.subfolder '\perf_pseudotrials_equalnum_x2_trainfold_only_0.09_0.1ms.mat'], 'perf')
        mean_accuracy_lda = mean(perf.lda.accuracy,1);
        mean_accuracy_svm = mean(perf.svm.accuracy,1);
        mean_accuracy_logreg = mean(perf.logreg.accuracy,1);
        mean_acc_all_classifiers = (mean_accuracy_lda+mean_accuracy_svm+mean_accuracy_logreg)./3;
        mean_acc_all_classifiers_all(ii,:)=mean_acc_all_classifiers; 
end

mean_acc_all_classifiers_all_subj = mean(mean_acc_all_classifiers_all,1);

[sort_acc, chan_location] = sort(mean_acc_all_classifiers_all_subj, 'descend');

figure;
scatter(1:248, sort_acc); 
find(sort_acc>0.55)

stat = [];
stat.mvpa.perf = mean_acc_all_classifiers_all_subj'*100;
stat.accuracy = mean_acc_all_classifiers_all_subj'*100;
stat.label = perf.features;
stat.dimord = 'chan';
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
cfg.highlight = 'on';
cfg.highlightchannel = perf.features(chan_location(1:54));
cfg.title = ['significant channels in subjects']; 
set(gcf,'color','w');
figure; ft_topoplotER(cfg, stat); 

end
