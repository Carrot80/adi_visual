function [grandavg, trl_count] = adi_grandavg_group_realignment_per_subject(subjectpath, path2inputfile, condition)


% load template grad-file:

switch condition
    case 'like'
 
    trialcount_like = [];
    avg_subjects_like = struct([]);
    for ii = 2:length(subjectpath) % wichtig, adi_04 rausnehmen, da andere Kanalsortierung
        %% like
        filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2inputfile 'Neu_Like*.mat']);
        all_runs_subject = struct([]);
        for kk = 1:length(filename)
            load ([filename(kk).folder filesep filename(kk).name])
%             [trials] = kh_trial2dat(cleanMEG_interp.trial);
%             if any(any(any(isnan(trials))))
%                 error(['NaNs in ' filename(kk).name ' ' subjectpath(ii).name]);
%             end
            clear trials
            cfg = [];
            cfg.lpfilter = 'yes';
            cfg.lpfreq   = 45;     
            cfg.demean = 'yes';
            cfg.baselinewindow  = [-0.5 -0.030];
            data = ft_preprocessing(cfg, cleanMEG_interp);

            %% select data

            cfg = [];
            cfg.latency = [-0.5 1];
            data = ft_selectdata(cfg, data);

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
            

             %% avg per run
             cfg = [];
             avg = ft_timelockanalysis(cfg, data);
             all_runs_subject{length(all_runs_subject)+1} = avg;
             trialcount_like(length(trialcount_like)+1) = numel(data.trial);
             clear cleanMEG_interp avg

        end
            clear vol gradfile
            close all

        %% grandavg like 
        count_runs = length(all_runs_subject);
        cfg = [];
        cfg.method  = 'within';
        switch count_runs
            case 3
                avg_like = ft_timelockgrandaverage(cfg, all_runs_subject{1}, all_runs_subject{2}, all_runs_subject{3});
            case 2
                avg_like = ft_timelockgrandaverage(cfg, all_runs_subject{1}, all_runs_subject{2});
            case 1
                avg_like = all_runs_subject{1};
        end
%         if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save], 'dir')
%             mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save])
%         end       
%         figure
%         plot(avg_like.time, avg_like.avg)
%         title('grandavg like')
%         axis tight
        avg_subjects_like{length(avg_subjects_like)+1} = avg_like;
        clear avg_like
    end

        %% grandavg like  
        cfg = [];
        cfg.method  = 'across';
        avg_like = ft_timelockgrandaverage(cfg, avg_subjects_like{:});
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\grandavg_like.mat', 'avg_like');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\avg_subjects_like.mat', 'avg_subjects_like');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\trialcount_like.mat', 'trialcount_like');
      
        figure
        plot(avg_like.time, avg_like.avg)
        title('grandavg like')
        axis tight
        savefig (['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\grandavg_like.fig' ]);
                    


       case 'dislike'       
            trialcount_dislike = [];
            avg_subjects_dislike = struct([]);
            for ii = 2:length(subjectpath) % wichtig, adi_04 rausnehmen, da andere Kanalsortierung
               
                filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2inputfile 'Neu_D*.mat']);
                all_runs_subject = struct([]);
                for kk = 1:length(filename)
                    load ([filename(kk).folder filesep filename(kk).name])
                    [trials] = kh_trial2dat(cleanMEG_interp.trial);
                    if any(any(any(isnan(trials))))
                        error(['NaNs in ' filename(kk).name ' ' subjectpath(ii).name]);
                    end
                    clear trials
                    cfg = [];
                    cfg.lpfilter = 'yes';
                    cfg.lpfreq   = 45;     
                    cfg.demean = 'yes';
                    cfg.baselinewindow  = [-0.5 -0.030];
                    data = ft_preprocessing(cfg, cleanMEG_interp);

                    %% select data

                    cfg = [];
                    cfg.latency = [-0.5 1];
                    data = ft_selectdata(cfg, data);

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

                     
                     %% avg per run
                     cfg = [];
                     avg = ft_timelockanalysis(cfg, data);
                     all_runs_subject{length(all_runs_subject)+1} = avg;
                     trialcount_dislike(length(trialcount_dislike)+1) = numel(data.trial);
                     clear cleanMEG_interp avg

                end
                    clear vol 
                    close all

                %% grandavg dislike 
                count_runs = length(all_runs_subject);
                cfg = [];
                cfg.method  = 'within';
                switch count_runs
                    case 6
                        avg_dislike = ft_timelockgrandaverage(cfg, all_runs_subject{1}, all_runs_subject{2}, all_runs_subject{3}, all_runs_subject{4}, all_runs_subject{5},  all_runs_subject{6});
                    case 5
                        avg_dislike = ft_timelockgrandaverage(cfg, all_runs_subject{1}, all_runs_subject{2}, all_runs_subject{3}, all_runs_subject{4}, all_runs_subject{5} );
                    case 4
                        avg_dislike = ft_timelockgrandaverage(cfg, all_runs_subject{1}, all_runs_subject{2}, all_runs_subject{3}, all_runs_subject{4});
                    case 3
                        avg_dislike = ft_timelockgrandaverage(cfg, all_runs_subject{1}, all_runs_subject{2}, all_runs_subject{3});
                    case 2
                        avg_dislike = ft_timelockgrandaverage(cfg, all_runs_subject{1}, all_runs_subject{2});
                    case 1
                        avg_dislike = all_runs_subject{1};
                end
        %         if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save], 'dir')
        %             mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save])
        %         end       
        %         figure
        %         plot(avg_like.time, avg_like.avg)
        %         title('grandavg like')
        %         axis tight
                avg_subjects_dislike{length(avg_subjects_dislike)+1} = avg_dislike;
                clear avg_dislike
            end

        %% grandavg dislike  
        cfg = [];
        cfg.method  = 'across';
        avg_dislike = ft_timelockgrandaverage(cfg, avg_subjects_dislike{:});
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\dontcare_und_dislike_zusammen\grandavg_dislike.mat', 'avg_dislike');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\dontcare_und_dislike_zusammen\avg_subjects_dislike.mat', 'avg_subjects_dislike');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\dontcare_und_dislike_zusammen\\trialcount_dislike.mat', 'trialcount_dislike');
        
        figure
        plot(avg_dislike.time, avg_dislike.avg)
        title('grandavg dislike')
        axis tight
        savefig (['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_all_trials\realignment_per_subject\dontcare_und_dislike_zusammen\grandavg_dislike.fig' ]);
                    
      end                

        
%% global field power: 
cfg = [];
cfg.method = 'power';
[gmf_like] = ft_globalmeanfield(cfg, avg_like);
[gmf_dislike] = ft_globalmeanfield(cfg, avg_dislike);
if exist('avg_dontcare', 'var')
    [gmf_dontcare] = ft_globalmeanfield(cfg, avg_dontcare);
    figure('Renderer', 'painters', 'Position', [300 300 800 200])
    plot(gmf_like.time, gmf_like.avg, 'r')
    axis tight
    hold on
    plot(gmf_dislike.time, gmf_dislike.avg, 'b')
    hold on
    plot(gmf_dontcare.time, gmf_dontcare.avg, 'k:')
    legend({'like'; 'dislike'; 'dontcare'})
    legend('boxoff') 
    title('global field power ')

    savefig([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'global_field_power_avg.fig' ]);
    clear avg_dontcare
else
    figure('Renderer', 'painters', 'Position', [300 300 800 200])
    plot(gmf_like.time, gmf_like.avg)
    axis tight
    hold on
    plot(gmf_dislike.time, gmf_dislike.avg)
    legend({'like'; 'dislike'})
    title('global field power ')
    legend('boxoff') 
    savefig([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'global_field_power_avg.fig' ]);

end


%% grandavg: 
          
cfg = [];
cfg.method  = 'within';
grandavg_like_group = ft_timelockgrandaverage(cfg, grandavg_like.avg);
grandavg_dislike_group = ft_timelockgrandaverage(cfg, grandavg_dislike.avg);

ind=[];
for kk = 1:numel(grandavg_dontcare)
    ind(kk) = isempty(grandavg_dontcare(kk).avg);
end

grandavg_dontcare(find(ind)) = [];
grandavg_dontcare_group = ft_timelockgrandaverage(cfg, grandavg_dontcare.avg);     

cfg = [];
cfg.method  = 'within';
grandavg = ft_timelockgrandaverage(cfg, avg_like, avg_dislike); 

figure('Renderer', 'painters', 'Position', [300 300 700 300]); plot(grandavg.time, grandavg.avg)

save (['W:\neurochirurgie\science\adidas\Kirsten\data_analysis\visual_stimuli\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_like.mat' ], 'grandavg_like_group');
save (['W:\neurochirurgie\science\adidas\Kirsten\data_analysis\visual_stimuli\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_dislike.mat' ], 'grandavg_dislike_group');
save (['W:\neurochirurgie\science\adidas\Kirsten\data_analysis\visual_stimuli\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_dontcare.mat' ], 'grandavg_dontcare_group');

%% Erstellung Abbildung für Publikation:

%% global field amplitude
cfg = [];
cfg.method = 'power'; %'amplitude', power
[gmf_like] = ft_globalmeanfield(cfg, grandavg_like_group);
[gmf_dislike] = ft_globalmeanfield(cfg, grandavg_dislike_group);
[gmf] = ft_globalmeanfield(cfg, grandavg);
[gmf_dontcare] = ft_globalmeanfield(cfg, grandavg_dontcare_group);

% plot global field amplitude

figure('Renderer', 'painters', 'Position', [300 300 800 200])
plot(gmf_like.time, gmf_like.avg, 'r', 'LineWidth',1)
axis tight
hold on
plot(gmf_dislike.time, gmf_dislike.avg, 'b', 'LineWidth',1)
% hold on
% plot(gmf_dontcare.time, gmf_dontcare.avg, 'k:', 'LineWidth',1)
legend({'like'; 'dislike'})
% title('global field amplitude ')
set(gcf,'color','w');
box off
legend({'like';  'dislike'}, 'boxoff')
legend('boxoff') 
xlabel('Time [s]')
ylabel('Global Field Amplitude [t]')

%% global field power
cfg = [];
cfg.method = 'power'; %'amplitude', power
[gmf_like] = ft_globalmeanfield(cfg, grandavg_like_group);
[gmf_dislike] = ft_globalmeanfield(cfg, grandavg_dislike_group);
[gmf] = ft_globalmeanfield(cfg, grandavg);
[gmf_dontcare] = ft_globalmeanfield(cfg, grandavg_dontcare_group);

% plot global field power
figure('Renderer', 'painters', 'Position', [300 300 800 200])
plot(gmf_like.time, gmf_like.avg, 'r', 'LineWidth',1)
axis tight
hold on
plot(gmf_dislike.time, gmf_dislike.avg, 'b', 'LineWidth',1)
hold on
plot(gmf_dontcare.time, gmf_dontcare.avg, 'k:', 'LineWidth',1)
legend({'like'; 'dislike'; 'dontcare'})
% title('global field power ')
set(gcf,'color','w');
box off
legend({'like';  'dislike'; 'dontcare'}, 'boxoff')
legend('boxoff') 
xlabel('Time [s]')
ylabel('Global Field Power')

%%

figure('Renderer', 'painters', 'Position', [300 300 800 200])
plot(gmf.time, gmf.avg, 'k', 'LineWidth',1.2)
hold on
plot(gmf.time, gmf_like.avg, 'r--', 'LineWidth',1.2)
hold on
plot(gmf.time, gmf_dislike.avg, 'b:', 'LineWidth',1.2)
axis tight
set(gcf,'color','w');
box off
legend({'all trials like and dislike'; 'like';  'dislike'}, 'boxoff')
legend('boxoff') 
    
%%
[rms_like] = rms(avg_like.avg);
[rms_dislike] = rms(avg_dislike.avg);
[rms_dontcare] = rms(grandavg_dontcare_group.avg);
[rms_avg] = rms(grandavg.avg);

figure('Renderer', 'painters', 'Position', [300 300 800 200])
plot(gmf_like.time, rms_like, 'r', 'LineWidth',1.2)
hold on
plot(gmf_like.time, rms_dislike, 'b', 'LineWidth',1.2)
% hold on
% plot(gmf.time, rms_dontcare, 'k:', 'LineWidth',1.2)
axis tight
set(gcf,'color','w');
box off
legend({'like';  'dislike'}, 'boxoff')
legend('boxoff') 
xlabel('Time [s]')
ylabel('Global Mean Amplitude') % GMA ist das gleiche wie RMS

%%
figure
plot(gmf.time, rms_avg, 'k', 'LineWidth',1.2)
hold on
plot(gmf.time, rms_like, 'r--', 'LineWidth',1.2)
hold on
plot(gmf.time, rms_dislike, 'b:', 'LineWidth',1.2)
axis tight
set(gcf,'color','w');
box off
legend({'all trials like and dislike'; 'like';  'dislike'}, 'boxoff')
legend('boxoff') 

%evtl. wie in language paper rms mit Signifikanztests berechnen

figure
plot(gmf.time, rms_avg, 'k', 'LineWidth',1.2)
hold on
plot(gmf.time, rms_like, 'r--', 'LineWidth',1.2)
hold on
plot(gmf.time, rms_dislike, 'b:', 'LineWidth',1.2)
axis tight
set(gcf,'color','w');
box off
legend({'all trials like and dislike'; 'like';  'dislike'}, 'Location','northeast');
legend('boxoff') 
xlabel('Time [s]')
ylabel('Global Mean Amplitude')
f = gcf;
f.Position(3:4) = [800 300]; 

%% topoplot für Publikation (like, dislike zusammen)


% baseline:
cfg = [];
cfg.xlim  = [-0.5 -0.03];
cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
figure; ft_topoplotER(cfg, grandavg);
colorbar; 
set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );

% comp 1:
cfg = [];
cfg.xlim  = [0.03 0.15];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );

% comp 2:

cfg = [];
cfg.xlim  = [0.15 0.56];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );

% comp 3:
cfg = [];
cfg.xlim  = [0.56 0.74];
% cfg.zlim = [-0.4 0.4];
cfg.comment = 'no';
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); 
set(gcf, 'Position', [885   539   363   288] );

%%


cfg = [];
cfg.xlim  = [0.0 0.1];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.1 0.2];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.2 0.3];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.3 0.4];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.4 0.5];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.5 0.6];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.6 0.7];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );


%% components condition like
close all
% baseline:
cfg = [];
cfg.xlim  = [-0.5 -0.03];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group);colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 1:
cfg = [];
cfg.xlim  = [0.03 0.15];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 2:

cfg = [];
cfg.xlim  = [0.15 0.56];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 3:
cfg = [];
cfg.xlim  = [0.56 0.74];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

%% components condition dislike
close all
% baseline:
cfg = [];
cfg.xlim  = [-0.5 -0.03];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group);colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 1:
cfg = [];
cfg.xlim  = [0.03 0.15];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 2:

cfg = [];
cfg.xlim  = [0.15 0.56];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 3:
cfg = [];
cfg.xlim  = [0.56 0.74];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

%%

figure
plot(grandavg_like_group.time, grandavg_like_group.avg)
axis tight
title('grandavg like')
figure
plot(grandavg_dislike_group.time, grandavg_dislike_group.avg)
axis tight
title('grandavg dislike')
        
savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_like.fig' ]);
        
%%
% comp1:

cfg = [];
cfg.channel   = 'all';
cfg.latency   = [.03 .15];
cfg.parameter = 'avg';
GA_like_comp1         = ft_timelockgrandaverage(cfg,grandavg_like.avg);
GA_dislike_comp1       = ft_timelockgrandaverage(cfg,grandavg_dislike.avg);

cfg = [];
cfg.showlabels  = 'yes';
cfg.layout      = '4D248.lay';
figure; ft_multiplotER(cfg,GA_like_comp1, GA_dislike_comp1)

% define the parameters for the statistical comparison
cfg = [];
cfg.channel     = 'all';
cfg.latency     = [0.03 0.15]; %[0.56 0.74];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'FDR';

Nsub = 29;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike.avg);   % don't forget the {:}!

cfg = [];
% cfg.style     = 'blank';
cfg.layout    = '4D248.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'yes';
figure; ft_topoplotER(cfg, GA_like_comp1) % 
title('analytic FDR')

%%
cfg = [];
cfg.channel     = 'all';
cfg.latency     = [0.03 0.15];
% cfg.latency     = [0.15 0.56];
% cfg.latency     = [0.56 0.74];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;

Nsub = 30;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike.avg); 

cfg = [];
cfg.style     = 'blank';
cfg.layout    = '4D248.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_like_comp1)
title('Nonparametric: significant without multiple comparison correction')
   

%% grandavg dislike und dontcare zusammengenommen:
cfg = [];
cfg.method  = 'within';
grandavg_dislike_dontcare_group = ft_timelockgrandaverage(cfg, grandavg_dislike_group, grandavg_dontcare_group)

cfg = [];
cfg.channel     = 'all';
% cfg.latency     = [0.56 0.74];
cfg.latency     = [0.03 0.16];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'FDR';

Nsub = 30;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike_dontcare_group);   % don't forget the {:}!



%%
        
        
chan = 52;
time = [0.3 0.7];

% find the time points for the effect of interest in the grand average data
timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));

% select the individual subject data from the time points and calculate the mean
for isub = 1:10
    values_FIC(isub)  = mean(allsubjFIC{isub}.avg(chan,timesel_FIC));
    values_FC(isub)  = mean(allsubjFC{isub}.avg(chan,timesel_FC));
end

% plot to see the effect in each subject
M = [values_FC',values_FIC'];
figure; plot(M','o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
        'subj7', 'subj8', 'subj9', 'subj10'}, 'location','EastOutside');        
        
        
        
        
        
        
        
        
%%

         
         %% grandavg like and dislike together:
         
        cfg = [];
        cfg.method  = 'within';
        grandavg_group = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group);
        figure
        plot(grandavg_group.time, grandavg_group.avg)
        axis tight
        title('grandavg like, dislike')
        cfg = [];
        cfg.method = 'power';
        [gmf] = ft_globalmeanfield(cfg, grandavg_group);
        figure
        plot(gmf.time, gmf.avg,  'Color',[0.17, 0.17, 0.17], 'LineWidth',1)
        hold on
        
        cfg = [];
        cfg.method  = 'within';
        grandavg_group_alltrls = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group,grandavg_dontcare_group);
        cfg = [];
        cfg.method = 'power';
        [gmf_alltrls] = ft_globalmeanfield(cfg, grandavg_group_alltrls);
        
        figure
        plot(gmf_alltrls.time, gmf_alltrls.avg,  'k', 'LineWidth',1)
        set(gcf,'color','w');
        box off

        
        
 %% cluster permutation test:
 
cfg = [];
cfg.method  = 'within';
cfg.keepindividual = 'yes';
grandavg_like_group = ft_timelockgrandaverage(cfg, avg_subjects_like{:});

cfg = [];
cfg.method  = 'within';
cfg.keepindividual = 'yes';
grandavg_dislike_group = ft_timelockgrandaverage(cfg, avg_subjects_dislike{:});

grandavg_like_group.grad = template.mean_grad;
cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, grandavg_like_group);

cfg = [];
avglike = ft_timelockanalysis(cfg, grandavg_like_group);
avgdislike = ft_timelockanalysis(cfg, grandavg_dislike_group);
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectlike_vs_dislike = ft_math(cfg, avglike, avgdislike);

%% %% cluster permutation test:
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
%cfg.neighbours = neighbours;   % see below
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
     
cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
                                 % which other sensors it can form clusters
% cfg.channel       = {'MEG'};     % cell-array with selected channel labels
% cfg.latency       = [0.03 0.15];      % time interval over which the experimental
cfg.latency       = [0.15 0.56];      % time interval over which the experimental

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

[i1,i2] = match_str(raweffectlike_vs_dislike.label, stat.label);

% % plot
% for k = 1:length(j)-1;
%    subplot(3,round(length(j)/3),k);
%    cfg = [];
%    cfg.xlim=[j(k) j(k+1)];
%    pos_int = zeros(numel(raweffectlike_vs_dislike.label),1);
%    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
%    cfg.highlight = 'on';
%    cfg.highlightchannel = find(pos_int);
%    cfg.comment = 'xlim';
%    cfg.commentpos = 'title';
%    cfg.layout = '4D248.lay';
%    ft_topoplotER(cfg, raweffectlike_vs_dislike);
% end



% plot avg:
 
cfg = [];
cfg.xlim=[0.03 0.15]; 
cfg.zlim=[-0.15 0.15];
pos_int = zeros(numel(raweffectlike_vs_dislike.label),1);
pos_int(i1) = all(pos(i2, 1), 2);  
% neg_int = zeros(numel(raweffectlike_vs_dislike.label),1);
% neg_int(i1) = all(neg(i2, 1), 2); 
cfg.highlight = 'on';
cfg.highlightchannel = find(pos_int);
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
clear pos pos_int

%%

for kk = 1:length(raweffectlike_vs_dislike.avg)
    [stat(kk).h,stat(kk).pvalue,stat(kk).ci,stat(kk).stats] = ttest(raweffectlike_vs_dislike.avg(:,kk), 0, 0.05) % H0: mean = 0, alpha 0.05
    t_stat(kk) = stat(kk).stats.tstat;
end

figure;
plot(raweffectlike_vs_dislike.time, raweffectlike_vs_dislike.avg)

figure;
plot(raweffectlike_vs_dislike.time, mean(abs(grandavg_like_group.avg)), 'r')
hold on
plot(raweffectlike_vs_dislike.time, mean(abs(grandavg_dislike_group.avg)), 'b')

figure;
plot(raweffectlike_vs_dislike.time, mean(abs(avglike.avg)), 'r')
hold on
plot(raweffectlike_vs_dislike.time, mean(abs(avgdislike.avg)), 'b')

figure;
plot(raweffectlike_vs_dislike.time, t_stat)
hold on
plot(raweffectlike_vs_dislike.time, zeros(1, length(raweffectlike_vs_dislike.time)), '-')
hold on
plot(raweffectlike_vs_dislike.time, 1.96*ones(1, length(raweffectlike_vs_dislike.time)), '-')
hold on
plot(raweffectlike_vs_dislike.time, -1.96*ones(1, length(raweffectlike_vs_dislike.time)), '-')

%%
%%
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'depsamplesT'; % use the independent samples T-statistic as a measure to
                               % evaluate the effect at the sample level
cfg.correctm = 'no';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
                               % permutation distribution.
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
cfg.avgoverchan = 'no';                               % required for a selected sample to be included
                               % in the clustering algorithm (default=0).
%cfg.neighbours = neighbours;   % see below
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 5000;      % number of draws from the permutation distribution

subj = size(grandavg_dislike_group.individual,1);
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
     
cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
                                 % which other sensors it can form clusters
% cfg.channel       = {'MEG'};     % cell-array with selected channel labels
cfg.latency       = [0.08 0.08];      % time interval over which the experimental
cfg.avgovertime = 'yes';                                 % conditions must be compared (in seconds)    
% [stat] = ft_timelockstatistics(cfg, grandavg_like_group_trls, grandavg_dislike_group_trls);
[stat] = ft_timelockstatistics(cfg, avg_subjects_like{:}, avg_subjects_dislike{:});

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









