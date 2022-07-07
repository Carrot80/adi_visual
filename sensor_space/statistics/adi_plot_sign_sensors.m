function adi_plot_sign_sensors(path2data, path2stats)
% plot significant sensors in component 1 and 2, see Dima et al., 2012, figure 2B
% noch nicht fertig, unklar, ob überhaupt sinnvoll für Publikation
load([path2data filesep 'avg_subjects_dislike.mat'])
load([path2data filesep 'avg_subjects_like.mat'])

load([path2data filesep 'grandavg_like.mat'])
load([path2data filesep 'grandavg_dislike.mat'])

load([path2stats 'stat_max_distribution_comp1.mat'])
% load([path2stats 'stats_comp2.mat'])

ind_sign_labels_comp1 = find(stat_max_distribution_comp1.mask)
sign_labels_comp1 = stat_max_distribution_comp1.label(ind_sign_labels_comp1);

load('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\template_gradfile\template.mat')

ind_sign_labels_comp1_occ_right = ind_sign_labels_comp1(1:4);
scatter3(mean_grad.chanpos(1:248,1), mean_grad.chanpos(1:248,2), mean_grad.chanpos(1:248,3))
hold on
scatter3(mean_grad.chanpos(ind_sign_labels_comp1_occ_right,1), mean_grad.chanpos(ind_sign_labels_comp1_occ_right,2), mean_grad.chanpos(ind_sign_labels_comp1_occ_right,3), 'r')
ind_sign_labels_comp1_occ_left = ind_sign_labels_comp1(5:8);
hold off
ind_sign_labels_comp1_occ_left = ind_sign_labels_comp1(5:8);
scatter3(mean_grad.chanpos(1:248,1), mean_grad.chanpos(1:248,2), mean_grad.chanpos(1:248,3))
hold on
scatter3(mean_grad.chanpos(ind_sign_labels_comp1_occ_left,1), mean_grad.chanpos(ind_sign_labels_comp1_occ_left,2), mean_grad.chanpos(ind_sign_labels_comp1_occ_left,3), 'r')



cfg = [];
cfg.method  = 'within';
cfg.keepindividual = 'yes';
grandavg_like_group = ft_timelockgrandaverage(cfg, avg_subjects_like{:});
grandavg_dislike_group = ft_timelockgrandaverage(cfg, avg_subjects_dislike{:});

cfg = [];
cfg.channel = stat_max_distribution_comp1.label(ind_sign_labels_comp1);
cfg.showlegend = 'yes';
figure; ft_singleplotER(cfg, grandavg_dislike_group, grandavg_like_group)
legend('boxoff')
set(gcf,'color','w');

cfg = [];
cfg.channel = stat_max_distribution_comp1.label(ind_sign_labels_comp1_occ_right);
cfg.showlegend = 'yes';
cfg.ylim = [-0.5 0.5];
figure; ft_singleplotER(cfg, grandavg_dislike_group, grandavg_like_group)
legend('boxoff')
set(gcf,'color','w');
title('Right Occipital')

cfg = [];
cfg.channel = stat_max_distribution_comp1.label(ind_sign_labels_comp1_occ_left);
%cfg.showlegend = 'yes';
cfg.ylim = [-0.5 0.5];
figure; ft_singleplotER(cfg, grandavg_dislike_group, grandavg_like_group)
%legend('boxoff')
set(gcf,'color','w');
title('Left Temporo-Parietal')

%% global field amplitude
avg_like_sel = avg_like;
avg_like_sel.avg = avg_like.avg(ind_sign_labels_comp1_occ_right,:);
avg_dislike_sel = avg_dislike;
avg_dislike_sel.avg = avg_dislike.avg(ind_sign_labels_comp1_occ_right,:);

cfg = [];
cfg.method = 'amplitude'; %'amplitude', power
[gfa_like_amplitude] = ft_globalmeanfield(cfg, avg_like_sel);
[gfa_dislike_amplitude] = ft_globalmeanfield(cfg, avg_dislike_sel);

figure; 
plot(gfa_like_amplitude.time, gfa_like_amplitude.avg, 'r')
hold on;
plot(gfa_dislike_amplitude.time, gfa_dislike_amplitude.avg, 'b')

cfg = [];
cfg.method = 'power'; %'amplitude', power
[gfa_like_power] = ft_globalmeanfield(cfg, avg_like_sel);
[gfa_dislike_power] = ft_globalmeanfield(cfg, avg_dislike_sel);

figure; 
plot(gfa_like_power.time, gfa_like_power.avg, 'r')
hold on;
plot(gfa_dislike_power.time, gfa_dislike_power.avg, 'b')

%% amplitude

like_selected_sensors = grandavg_like_group.individual(:,ind_sign_labels_comp1_occ_right,:);
dislike_selected_sensors = grandavg_dislike_group.individual(:,ind_sign_labels_comp1_occ_right,:);


for kk=1:29
    rms_like_subjects(kk,:,:)=rms(like_selected_sensors(kk,:,:))
    rms_dislike_subjects(kk,:,:)=rms(dislike_selected_sensors(kk,:,:))
end








%%



like_selected_sensors = squeeze(mean(like_selected_sensors,1));
dislike_selected_sensors = squeeze(mean(dislike_selected_sensors,1));



figure; 
hold on;
plot(grandavg_like_group.time, mean(rms_like_subjects,1), 'r', 'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(grandavg_like_group.time, mean(rms_like_subjects,1)+SEM_like, 'color', [0.5 0.5 0.5]) % unklar, ob SEM noch mal durch 2 zu teilen ist
plot(grandavg_like_group.time, mean(rms_like_subjects,1)-SEM_like, 'color', [0.5 0.5 0.5])





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

SEM_dislike = std(rms_dislike_subjects,[],1)/sqrt(size(rms_dislike_subjects,1));
SEM_like = std(rms_like_subjects,[],1)/sqrt(size(rms_like_subjects,1));

figure; 
errorbar(grandavg_like_group.time, mean(rms_like_subjects,1), std(rms_like_subjects,[],1)/sqrt(size(rms_like_subjects,1)));
hold on;
errorbar(grandavg_dislike_group.time, mean(rms_dislike_subjects,1), std(rms_dislike_subjects,[],1)/sqrt(size(rms_dislike_subjects,1)));
plot(grandavg_like_group.time, mean(rms_like_subjects,1), 'b', 'MarkerFaceColor','k','MarkerEdgeColor','k')
plot(grandavg_like_group.time, mean(rms_dislike_subjects,1), 'k', 'MarkerFaceColor','b','MarkerEdgeColor','k')


% plot(grandavg_like_group.time, mean(rms_dislike_subjects,1), 'b', 'MarkerFaceColor','b','MarkerEdgeColor','k')


transparency_color = 0.1;
edge = 'r';
transparency_edge = 0.1;
filled = [mean(rms_like_subjects,1)+SEM_like,fliplr(mean(rms_like_subjects,1)-SEM_like)];
xpoints = [grandavg_like_group.time, fliplr(grandavg_like_group.time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color


end