function [] = adi_leave_out_exemplar_singlsubj(subjectdir, config)

%%
% created: 21.01.2020
% modified: 26.2.2020 (built in channelselection)
% modified 1.4.2022 ==> realignent of head positions to common sensor
% position

for ii = 1:length(subjectdir)
    
    if 2 == exist([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\session.mat'], 'file')
        
        load ([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\session.mat'])
        
    else
    
        dir_data = dir([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\' '*.mat']);
        ind_ = [];
        for kk = 1:length(dir_data)
             ind_(kk) = contains(dir_data(kk).name, 'Dontcare');
        end
        
        indx_dontcare = find(ind_);
        if ~isempty(indx_dontcare)
            dir_data(indx_dontcare) = [];
        end
        for kk = 1:length(dir_data)
               
            load ([dir_data(kk).folder filesep dir_data(kk).name], 'cleanMEG_interp')
            for pp = 1:length(cleanMEG_interp.trial)
                cleanMEG_interp.trialinfo.run{pp} = dir_data(kk).name(end-4);
            end
            if 1 == isfield(cleanMEG_interp, 'additional_cleaning')
                cleanMEG_interp = rmfield(cleanMEG_interp,  'additional_cleaning');
            end
            [data_bpfreq_res_sel] = adi_bpfilter(cleanMEG_interp, 'bp1_45Hz');         
            clear cleanMEG_interp 
        

        %% realignment per subject:
        
            if endsWith(dir_data(kk).name, '_1.mat') || kk==1
                    gradfile = data_bpfreq_res_sel.grad;
                    data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'ChannelFlag_Bst');
                    if isfield(data_bpfreq_res_sel, 'additional_cleaning')
                        data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'additional_cleaning');
                    end
                else
                    cfg=[];
                    cfg.template = gradfile;
                    cfg.inwardshift = 1;
                    load([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\vol\vol.mat'])
                    cfg.headmodel = vol;
                    data_bpfreq_res_sel_ra = ft_megrealign(cfg, data_bpfreq_res_sel);
                    data_bpfreq_res_sel_ra = rmfield(data_bpfreq_res_sel_ra, 'sampleinfo');
                    data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'ChannelFlag_Bst');
                    if isfield(data_bpfreq_res_sel, 'additional_cleaning')
                        data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'additional_cleaning');
                    end
                    close all
                    fn_ra = fieldnames(data_bpfreq_res_sel_ra);
                    fn_orig = fieldnames(data_bpfreq_res_sel);
                    diff_fn = setdiff(fn_orig, fn_ra);
                    data_bpfreq_res_sel_ra.trialinfo = data_bpfreq_res_sel.trialinfo;
                    data_bpfreq_res_sel_ra.dimord = data_bpfreq_res_sel.dimord;  
                    data_bpfreq_res_sel_ra=orderfields(data_bpfreq_res_sel_ra);
                    data_bpfreq_res_sel=orderfields(data_bpfreq_res_sel);
            end   
            if endsWith(dir_data(kk).name, '_1.mat') || kk==1
                   data_realigned(kk) = data_bpfreq_res_sel;
               else
                   data_realigned(kk) = data_bpfreq_res_sel_ra;
            end
                
                clear cleanMEG_interp data_bpfreq_res_sel  data_bpfreq_res_sel_ra
        end
          
        %%
        session = data_realigned(1);
        session.response_label = data_realigned(1).trialinfo.response_label;
        session.run = data_realigned(1).trialinfo.run;
        session.balldesign = data_realigned(1).trialinfo.balldesign_short;
        session = rmfield(session, 'trialinfo');

        for kk = 2:length(data_realigned)
            session.trial = cat(2, session.trial, data_realigned(kk).trial);
            session.time = cat(2, session.time, data_realigned(kk).time);
%             session.cfg = cat(2, session.cfg, data_realigned(kk).cfg);
            session.response_label = cat(2, session.response_label, data_realigned(kk).trialinfo.response_label);
            session.balldesign = cat(2, session.balldesign, data_realigned(kk).trialinfo.balldesign_short);
            session.run = cat(2, session.run, data_realigned(kk).trialinfo.run);
        end
        clear data_realigned
        for kk = 1:length(session.response_label)
            switch session.response_label{kk}
                case 'like' % 'Volley'
                    session.labels(kk) = 1;
                case 'Neu_Like' % 'Volley'
                    session.labels(kk) = 1;
                case 'dislike' % 'Space'
                    session.labels(kk) = 2;
                case 'Neu_Dislike' % 'Space'
                    session.labels(kk) = 2;
            end
        end

      %% z-transform trials  

        [session] = adi_ztrans_sensorspace(session);
    
    end
 
     % select time interval:
     if ~isempty(config.channelselection)
        cfg = [];
        cfg.latency = config.time_interval;
        cfg.channel = config.channelselection;
        session = ft_selectdata(cfg, session);
     end
     
    % bei Proband 5 und 28 sind zwei Antworten mit 0 gelabelt. Diese trials
    % müssen entfernt werden
    indx_wrong_label = find(session.labels == 0);
    if ~isempty(indx_wrong_label)
        session.trial(indx_wrong_label)=[];
        session.time(indx_wrong_label)=[];
        session.labels(indx_wrong_label)=[];
        session.balldesign(indx_wrong_label)=[];
        session.run(indx_wrong_label)=[];
    end

    like_dislike_ratings = [];
    for kk = 1:length(session.balldesign) 
        balldesign_array(kk)= session.balldesign{kk};
    end
    
    [count_balldesign, balls] = histcounts(categorical(balldesign_array), categorical(unique(balldesign_array)));
    
    for kk = 1:length(balls)
        like_dislike_ratings.(balls{kk}) = session.labels(find(strcmp(balldesign_array, balls{kk})));
    end
    session.like_dislike_ratings = like_dislike_ratings; 
    
    clear balldesign_array
    
%     session = rmfield(session, 'trial');
%     if isfield(config, 'channelselection')
%         chan_ind = zeros(1, size(config.channelselection,1));
%         for jj = 1:length(config.channelselection)
%             chan_ind(jj) = find(strcmp(config.channelselection{jj}, session.label));
%         end
%         data = data(:, chan_ind, :);
%     end

    %% MVPA  
    mvpa_leave_out_balldesign(session, config, subjectdir(ii));
   
    
    %%
    clear session data
      
   
   end




end

function [] = mvpa_leave_out_balldesign(session, config, subjectdir)

%% check if toolbox is in path:

if isempty(strfind(path,'E:\Arbeit\Matlab\Bibliotheken\MVPA-Light-master'))
    addpath(genpath('E:\Arbeit\Matlab\Bibliotheken\MVPA-Light-master\'));
end
data_trials = kh_trial2dat(session.trial);

%% Get default hyperparameters for the logreg and lda classifier

param_lda = mv_get_hyperparameter('lda');
param_svm = mv_get_hyperparameter('svm');

    %% built train and test data based on CV folds
cf_lda = cell(length(config.balldesign),length(session.time{1}));
cf_svm = cell(length(config.balldesign),length(session.time{1}));

%% Crossvalidation folds

CV = adi_crossval_leaveExemplarOut(session.balldesign, session.labels);

%%
indx_time_1 = nearest(session.time{1,1}, config.time_interval(1)); 
indx_time_2 = nearest(session.time{1,1}, config.time_interval(2));  

n_CV = length(CV);

%% data_bootstrap_trainfold verwenden für die Klassifikation, so dass nur labels getauscht werden und die Daten nicht verändert werden
   
for kk = 1:n_CV 
      
    % labels des testset werden beibehalten
    test_fold = data_trials(CV(kk).testset,:,:);
    clabel_test_fold = session.labels(CV(kk).testset);
   
    train_fold = data_trials;
    train_fold(CV(kk).testset,:,:) = [];  
    clabel_train_fold = session.labels;
    clabel_train_fold(CV(kk).testset) = [];
    
%% MVPA 
       
       
%%  bootstrap trainfold    
    trainfold_like = train_fold(clabel_train_fold == 1,:,:);
    trainfold_dislike = train_fold(clabel_train_fold == 2,:,:);
    max_numel_trials = max([size(trainfold_like,1), size(trainfold_dislike,1)]);
    trainfold = [];
    cfg = [];
    cfg.mode = 'bootstrap';
    cfg.averages = 4;
    cfg.repetitions = max_numel_trials;
    bootstrap_trainfold_like = fte_subaverage(cfg, trainfold_like);
    bootstrap_trainfold_dislike = fte_subaverage(cfg, trainfold_dislike);
    trainfold_like = [];
    trainfold_dislike = [];
    data_bootstrap_trainfold = cat(1, bootstrap_trainfold_like, bootstrap_trainfold_dislike);
    clabel_train_fold = cat(2, ones(1, size(bootstrap_trainfold_like,1)), 2*ones(1, size(bootstrap_trainfold_dislike,1)));
    bootstrap_trainfold_like = [];
    bootstrap_trainfold_dislike = [];
    
%%  lda for single time points  
  
    
   %%  
    size_random_dataset = 1000;
    clabel_train_fold_rand = zeros(size_random_dataset, length(clabel_train_fold));
    
    for pp = 1:size_random_dataset
        clabel_train_fold_rand(pp,:) = clabel_train_fold(randperm(length(clabel_train_fold)));
    end

    disp(['creating null distribution for balldesign ' config.balldesign{kk}])
    
    
    %% lda
    accuracy_lda = nan(size_random_dataset, size(session.time{1},2));
%     auc = nan(size_random_dataset, size(session.time{1},2));
%     confusion_percentage = cell(size_random_dataset, size(session.time{1},2));
    confusion_lda = cell(size_random_dataset, size(session.time{1},2));
    f1score_lda = cell(size_random_dataset, size(session.time{1},2));
 
    for tt = indx_time_1:indx_time_2
        data_bootstrap_trainfold_tt = data_bootstrap_trainfold(:,:,tt);
        test_fold_tt = test_fold(:,:,tt);
        parfor pp = 1:size(clabel_train_fold_rand,1)
            cf_perm_lda = train_lda(param_lda, data_bootstrap_trainfold_tt, clabel_train_fold_rand(pp,:));
            [predlabel_perm, dval_perm] = test_lda(cf_perm_lda, test_fold_tt); 

            % Calculate AUC and Accuracy and other performance measures:
%             auc(pp, tt) = mv_calculate_performance('auc', 'dval', dval_perm, clabel_test_fold);
            accuracy_lda(pp, tt) = mv_calculate_performance('acc', 'dval', dval_perm, clabel_test_fold);
%             confusion_percentage{pp, tt} = mv_calculate_performance('confusion', 'clabel', predlabel_perm', clabel_test_fold);
            confusion_lda{pp, tt} = confusionmat(clabel_test_fold, predlabel_perm');
            f1score_lda{pp, tt} = kh_mv_calculate_performance_alt('f1', 'clabel', predlabel_perm', clabel_test_fold);
        end     
    end
    
    perf_perm.lda.(config.balldesign{kk}).accuracy = accuracy_lda;   
    perf_perm.lda.(config.balldesign{kk}).confusion = confusion_lda;
    perf_perm.lda.(config.balldesign{kk}).f1score = f1score_lda;
    perf_perm.lda.(config.balldesign{kk}).number_of_trials.testset = numel(clabel_test_fold);  
    perf_perm.lda.(config.balldesign{kk}).number_of_trials.trainingsset = numel(clabel_train_fold);
    perf_perm.lda.(config.balldesign{kk}).time_interval = [num2str(config.time_interval(1)) '_' num2str(config.time_interval(2))];
 
  
     %% svm
    accuracy_svm = nan(size_random_dataset, size(session.time{1},2));
%     auc = nan(size_random_dataset, size(session.time{1},2));
    confusion_svm = cell(size_random_dataset, size(session.time{1},2));
    f1score_svm = cell(size_random_dataset, size(session.time{1},2));
 
    for tt = indx_time_1:indx_time_2
        data_bootstrap_trainfold_tt = data_bootstrap_trainfold(:,:,tt);
        test_fold_tt = test_fold(:,:,tt);
        parfor pp = 1:size(clabel_train_fold_rand,1)
            cf_perm_svm = train_svm(param_svm, data_bootstrap_trainfold_tt, clabel_train_fold_rand(pp,:));
            [predlabel_perm, dval_perm] = test_svm(cf_perm_svm, test_fold_tt); 

            % Calculate AUC and Accuracy and other performance measures:
%             auc(pp, tt) = mv_calculate_performance('auc', 'dval', dval_perm, clabel_test_fold);
            accuracy_svm(pp, tt) = mv_calculate_performance('acc', 'dval', dval_perm, clabel_test_fold);
%             confusion_percentage{pp, tt} = mv_calculate_performance('confusion', 'clabel', predlabel_perm', clabel_test_fold);
            confusion_svm{pp, tt} = confusionmat(clabel_test_fold, predlabel_perm');
            f1score_svm{pp, tt} = kh_mv_calculate_performance_alt('f1', 'clabel', predlabel_perm', clabel_test_fold);
        end     
    end
    
    perf_perm.svm.(config.balldesign{kk}).accuracy = accuracy_svm;   
    perf_perm.svm.(config.balldesign{kk}).confusion = confusion_svm;
    perf_perm.svm.(config.balldesign{kk}).f1score = f1score_svm;
    perf_perm.svm.(config.balldesign{kk}).number_of_trials.testset = numel(clabel_test_fold);  
    perf_perm.svm.(config.balldesign{kk}).number_of_trials.trainingsset = numel(clabel_train_fold);
    perf_perm.svm.(config.balldesign{kk}).time_interval = [num2str(config.time_interval(1)) '_' num2str(config.time_interval(2))];

  
    clear test_fold train_fold clabel_train_fold clabel_test_fold
    
end

perf_perm.CV = CV;
perf_perm.config = config;
perf_perm.config.methods = 'perf_without_pca_bootstrap_4pseudotrials_per_condition_in_trainfold';
perf_perm.config.comment = 'da testset trials mit sowohl like als auch dislike-Antworten enthalten kann, wurden nur trials aus dem trainfold zur Bildung von Pseudotrials gemittelt. Um unbalanced data zu balancieren, wurde oversampling der kleineren Bedingung (Like/dislike) durchgeführt';
perf_perm.ratings = session.like_dislike_ratings;
   
if ~exist([config.path2save.mainfolder subjectdir.name filesep config.path2save.subfolder], 'dir')
    mkdir([config.path2save.mainfolder subjectdir.name filesep config.path2save.subfolder])
end
save([config.path2save.mainfolder subjectdir.name filesep config.path2save.subfolder 'performance_permutation_null_distribution_cluster_perm_' num2str(config.time_interval(1)) '_' num2str(config.time_interval(2)) 'ms.mat'], 'perf_perm')


  
end

function [CV] = adi_crossval_leaveExemplarOut(balldesign_short, labels)

for kk=1:length(balldesign_short) 
    balldesign_array(kk) = balldesign_short{kk};
end

[count_balldesign, balls] = histcounts(categorical(balldesign_array), categorical(unique(balldesign_array)));

for kk = 1:length(balls)
    index_design.(balls{kk}) = find(strcmp(balldesign_array, balls(kk)));
end


for kk = 1:length(balls)
    CV(kk).testset = index_design.(balls{kk});
    CV(kk).trainingsset = 1:length(balldesign_short);
    CV(kk).trainingsset(CV(kk).testset) = [];
    CV(kk).design = balls{kk};
end

for kk = 1:length(balls)
    CV(kk).labels_trainingsset = labels;
    CV(kk).labels_trainingsset(CV(kk).testset) = [];
    CV(kk).labels_testset = labels(index_design.(balls{kk}));
    CV(kk).balldesign_trainingsset = balldesign_short;
    CV(kk).balldesign_trainingsset(CV(kk).testset) = [];
    CV(kk).trialnumer_likes_trainingsset = numel(find(CV(kk).labels_trainingsset==1));
    CV(kk).trialnumer_dislikes_trainingsset = numel(find(CV(kk).labels_trainingsset==2));
    CV(kk).ratio_trialnumer_likes_dislikes_trainingsset = CV(kk).trialnumer_likes_trainingsset/CV(kk).trialnumer_dislikes_trainingsset ;
    
end

end


function [data_bpfreq_res_sel] = adi_bpfilter(filename, bpname)


switch bpname
    case 'bp1_95Hz'
       data_bpfreq_res_sel =  filename;
        return
    case 'bp1_45Hz'
        bpfreq = [1 45];
    case'1_5_45Hz'
        bpfreq = [1.5 45];
    case 'bp2_45Hz'
        bpfreq = [2 45];
    case 'bp3_45Hz'
        bpfreq = [3 45];
    case 'delta'
        bpfreq = 4;
    case 'theta'
        bpfreq = [4 8];
    case 'alpha'
        bpfreq = [8 13];
    case 'beta'
        bpfreq = [13 25];
    case 'low_gamma'
        bpfreq = [25 45];
    case 'high_gamma'
        bpfreq = [55 90];
    case 'bp10-45Hz'
        bpfreq = [10 45];
end

cfg = [];
cfg.keeptrials = 'yes';
cfg.vartrllength = 2;

cfg.trials  = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') %|| 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

try
    [data_bpfreq] = ft_preprocessing(cfg, filename); 
    [warnMsg, warnID] = lastwarn;
    if ~isempty(warnMsg)
       warnMsg
    end
catch
    
    warnMsg
end

cfg =[];
cfg.resamplefs = 256;
% cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);

for k = 1:length(data_bpfreq_res_sel.trial)
    data_bpfreq_res_sel.grad.label(249:end) = [];
    data_bpfreq_res_sel.grad.chanori(249:end, :) = [];
    data_bpfreq_res_sel.grad.chanpos(249:end, :) = [];
    data_bpfreq_res_sel.grad.tra(249:end, :) = [];
    data_bpfreq_res_sel.label(249:end) = [];
end

fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_res_sel);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_sel_res);

for k=1:length(diff_fieldnames)
    data_bpfreq_res_sel.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
end
% fn_filename{end+1}='cfg';
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, fn_filename);
clearvars filename data_bpfreq data_bpfreq_res 

end



