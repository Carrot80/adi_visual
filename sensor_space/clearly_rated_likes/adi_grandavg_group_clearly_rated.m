function [grandavg, trl_count] = adi_grandavg_group_clearly_rated(subjectpath, grandavg, path2inputfile, condition);


response_tbl = readtable(['E:\Arbeit\adidas\data_analysis\balldesign_ratings.xlsx']);
likes = response_tbl.eindeutigeLikes;
dislikes = response_tbl.eindeutigeDislikes;


switch condition
    case 'like'
        counter = 1;
        avg_subjects = struct([]);
        for ii = 2:length(subjectpath) % wichtig, adi_04 rausnehmen, da andere Kanalsortierung
            %% like
            filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2inputfile 'Neu_Like*.mat']);
            all_runs_subject = struct([]);
            for kk = 1:length(filename)
                load ([filename(kk).folder filesep filename(kk).name])
                
                [likes_clearly_rated] = find_unambiguous_trials(cleanMEG_interp, subjectpath(ii).name, response_tbl, likes)
                
                % realignment per subject:
%                 if kk==1
%                     
%                     
%                 end
                if numel(likes_clearly_rated.trial) > 0

                    [trials] = kh_trial2dat(likes_clearly_rated.trial);
            
                    if any(any(any(isnan(trials))))
                        error(['NaNs in ' filename(kk).name ' ' subjectpath(ii).name]);
                    end
                    clear trials
                    cfg = [];
                    cfg.lpfilter = 'yes';
                    cfg.lpfreq   = 45;     
                    cfg.demean = 'yes';
                    cfg.baselinewindow  = [-0.5 -0.030];
                    data = ft_preprocessing(cfg, likes_clearly_rated);

                    %% select data

                    cfg = [];
                    cfg.latency = [-0.5 1];
                    data = ft_selectdata(cfg, data);

                    %% z-transformation:
                     [data] = subfun_ztransform(data);

                     %% avg per run
                     cfg = [];
                     avg = ft_timelockanalysis(cfg, data);
                     all_runs_subject{length(all_runs_subject)+1} = avg;
                     trialcount_like(counter) = numel(data.trial);
                     clear cleanMEG_interp data avg
                     counter = counter+1;
                end

                end
                    clear vol 
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
            if ~isempty(all_runs_subject)
                avg_subjects_like{length(avg_subjects)+1} = avg_like;
                clear avg_like
            end
      
            
        end

        %% grandavg like  
        cfg = [];
        cfg.method  = 'across';
        avg_like = ft_timelockgrandaverage(cfg, avg_subjects_like{:});
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\grandavg_like.mat', 'avg_like');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\avg_like_subjects.mat', 'avg_subjects_like');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\trialcount_like.mat', 'trialcount_like');
      
        figure
        plot(avg_like.time, avg_like.avg)
        title('grandavg like')
        axis tight
        savefig (['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\grandavg_like.fig' ]);



    case 'dislike'
        
        counter = 1;
        avg_subjects = struct([]);
        for ii = 2:length(subjectpath) % wichtig, adi_04 rausnehmen, da andere Kanalsortierung
            %% dislike
            filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2inputfile 'Neu_Dislike*.mat']);
            all_runs_subject = struct([]);
            for kk = 1:length(filename)
                load ([filename(kk).folder filesep filename(kk).name])
                
                [dislikes_clearly_rated] = find_unambiguous_trials(cleanMEG_interp, subjectpath(ii).name, response_tbl, dislikes)
                
                % realignment per subject:
%                 if kk==1
%                     
%                     
%                 end
                if numel(dislikes_clearly_rated.trial) > 0

                    [trials] = kh_trial2dat(dislikes_clearly_rated.trial);
            
                    if any(any(any(isnan(trials))))
                        error(['NaNs in ' filename(kk).name ' ' subjectpath(ii).name]);
                    end
                    clear trials
                    cfg = [];
                    cfg.lpfilter = 'yes';
                    cfg.lpfreq   = 45;     
                    cfg.demean = 'yes';
                    cfg.baselinewindow  = [-0.5 -0.030];
                    data = ft_preprocessing(cfg, dislikes_clearly_rated);

                    %% select data

                    cfg = [];
                    cfg.latency = [-0.5 1];
                    data = ft_selectdata(cfg, data);

                    %% z-transformation:
                     [data] = subfun_ztransform(data);

                     %% avg per run
                     cfg = [];
                     avg = ft_timelockanalysis(cfg, data);
                     all_runs_subject{length(all_runs_subject)+1} = avg;
                     trialcount_dislike(counter) = numel(data.trial);
                     clear cleanMEG_interp data avg
                     counter = counter+1;
                end

                end
                    clear vol 
                    close all

            %% grandavg dislike 
            count_runs = length(all_runs_subject);
            cfg = [];
            cfg.method  = 'within';
            switch count_runs
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
            if ~isempty(all_runs_subject)
                avg_subjects_dislike{length(avg_subjects)+1} = avg_dislike;
                clear avg_dislike
            end
      
            
        end
        
        cfg = [];
        cfg.method  = 'across';
        avg_dislike = ft_timelockgrandaverage(cfg, avg_subjects{:});
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\grandavg_dislike.mat', 'avg_dislike');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\avg_dislike_subjects.mat', 'avg_subjects_dislike');
        save ('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\trialcount_dislike.mat', 'trialcount_dislike');
       
        figure
        plot(avg_dislike.time, avg_dislike.avg)
        title('grandavg dislike')
        axis tight
        savefig (['E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\grandavg_dislike.fig' ]);

end



end

function [cleanMEG_interp] = find_unambiguous_trials(cleanMEG_interp, subject, response_tbl, condition)

response_tbl.Subject
ind_subj = strcmp(response_tbl.Subject, subject);
ind_balldesign = response_tbl.balldesign(ind_subj);
ind_keeptrials = condition(ind_subj);
ind_balldesign(find(ind_keeptrials))=[];

delete_trials = zeros(1, length(cleanMEG_interp.trial));

for pp=1:length(cleanMEG_interp.trial)
    try
        balldesign_trial=cleanMEG_interp.trialinfo.balldesign_short{1,pp}{1,1};
    catch
        
    end
%     if any(strcmp(ind_balldesign, balldesign_trial))
%         delete_trials(pp) = 1;
%     end 
    if any(strcmp(balldesign_trial, ind_balldesign))
        delete_trials(pp) = 1;
    end 
end

cleanMEG_interp.trial(find(delete_trials)) = [];
cleanMEG_interp.time(find(delete_trials)) = [];
cleanMEG_interp.trialinfo.balldesign_short(find(delete_trials)) = [];
cleanMEG_interp.trialinfo.balldesign(find(delete_trials)) = [];
cleanMEG_interp.trialinfo.response_label(find(delete_trials)) = [];
cleanMEG_interp.trialinfo.response(find(delete_trials)) = [];
cleanMEG_interp.trialinfo.triggerlabel(find(delete_trials)) = [];
cleanMEG_interp.trialinfo.triggerchannel(find(delete_trials)) = [];
cleanMEG_interp.trialinfo.responsechannel(find(delete_trials)) = [];

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

