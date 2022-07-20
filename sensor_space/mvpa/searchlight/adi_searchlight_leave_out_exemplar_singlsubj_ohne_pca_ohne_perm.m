function [] = adi_leave_out_exemplar_singlsubj(mainpath, subjectdir, filename, balldesign, time_interval, neighbours)


for ii = 1:length(subjectdir)
    
    
    if ~exist([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\'], 'dir')
        mkdir(([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\']))
    end
    if 2 == exist([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\session_fsample_orig.mat'], 'file')
        
        load ([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\session_fsample_orig.mat'])
        
    else
    
        dir_data = dir([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\' '*.mat']);
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

                [data_bpfreq_res_sel] = adi_bpfilter(cleanMEG_interp, 'bp1_45Hz');

                if endsWith(dir_data(kk).name, '_1.mat') || kk==1
                    gradfile = data_bpfreq_res_sel.grad;
                    data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'ChannelFlag_Bst');
                    data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'additional_cleaning');
                else
                    cfg=[];
                    cfg.template = gradfile;
                    cfg.inwardshift = 1;
                    load([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\vol\vol.mat'])
                    cfg.headmodel = vol;
                    data_bpfreq_res_sel_ra = ft_megrealign(cfg, data_bpfreq_res_sel);
%                     data_bpfreq_res_sel_ra = rmfield(data_bpfreq_res_sel_ra, 'sampleinfo');
                    data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'ChannelFlag_Bst');
                    data_bpfreq_res_sel = rmfield(data_bpfreq_res_sel, 'additional_cleaning');
                    fn_ra = fieldnames(data_bpfreq_res_sel_ra);
                    fn_orig = fieldnames(data_bpfreq_res_sel);
                    diff_fn = setdiff(fn_orig, fn_ra);
                    data_bpfreq_res_sel_ra.trialinfo = data_bpfreq_res_sel.trialinfo;
                    data_bpfreq_res_sel_ra.dimord = data_bpfreq_res_sel.dimord;  
                    data_bpfreq_res_sel_ra=orderfields(data_bpfreq_res_sel_ra);
                    data_bpfreq_res_sel=orderfields(data_bpfreq_res_sel);
                end   
                if endsWith(dir_data(kk).name, '_1.mat') || kk==1
                   data(kk) = data_bpfreq_res_sel;
                else
                   data(kk) = data_bpfreq_res_sel_ra;
                end
                
                clear cleanMEG_interp data_bpfreq_res_sel  data_bpfreq_res_sel_ra
               
        end
         
        close all
        
        session = data(1);
        session.response_label = data(1).trialinfo.response_label;
        session.run = data(1).trialinfo.run;
        session.balldesign = data(1).trialinfo.balldesign_short;
        session = rmfield(session, 'trialinfo');

        for kk = 2:length(data)
            session.trial = cat(2, session.trial, data(kk).trial);
            session.time = cat(2, session.time, data(kk).time);
%             session.cfg = cat(2, session.cfg, data(kk).cfg);
            session.response_label = cat(2, session.response_label, data(kk).trialinfo.response_label);
            session.balldesign = cat(2, session.balldesign, data(kk).trialinfo.balldesign_short);
            session.run = cat(2, session.run, data(kk).trialinfo.run);
        end
        
        clear data
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
%         save ([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realignment_per_subject\session.mat'])
%    if ~exist(['D:\Arbeit\adidas\data_analysis\single_subjects\' subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realignment_per_subject\'], 'dir')
%         mkdir(['D:\Arbeit\adidas\data_analysis\single_subjects\' subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realignment_per_subject\'])
%    end
        save (['C:\Users\herfurkn1\Documents\adidas\data_analysis\single_subjects\' subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\mvpa\session_fsample_orig.mat'], 'session')
        save ([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\session_fsample_orig.mat'], 'session')
    end
  
    %% MVPA  
    

        
%         if ~exist([mainpath subjectdir(ii).name filesep filename  '_' num2str(comp(pp, 1)) '_' num2str(comp(pp, 2)) '_ms_lda.fig'], 'file')
        
            [perf] = mvpa_leave_out_balldesign(session, balldesign, time_interval, neighbours);

            if ~exist([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\'], 'dir')
                mkdir([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\'])
            end
        %     perf.config.methods = 'bootstrap 2*trialnumber';
        %     perf.config.methods = 'bootstrap_pseudotrials';

            perf.number_of_trials = numel(session.trial);
            perf.config.methods = 'perf_without_pca_bootstrap_4pseudotrials_per_condition_trainfold_only';
            perf.config.comment = ' da alle trials mit den uneindeutigen Antworten zur Klassifikation genutzt worden, wurden nur trials aus dem trainfold zur Bildung von Pseudotrials gemittelt; es wurde trialanzahl für like und dislike angeglichen ';
            perf.features = session.label;
            save([mainpath subjectdir(ii).name filesep filename 'perf_pseudotrials_equalnum_x1_trainfold_only_' num2str(time_interval(1)) '_' num2str(time_interval(2)) 'ms_sample_orig.mat'], 'perf')
            save(['C:\Users\herfurkn1\Documents\adidas\data_analysis\single_subjects\' subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\mvpa\searchlight\' 'perf_pseudotrials_equalnum_x1_trainfold_only_' num2str(time_interval(1)) '_' num2str(time_interval(2)) 'ms_sample_orig.mat'], 'perf')
            if ~exist([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\figures\'], 'dir')
                mkdir([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\figures\'])
            end

           %% mean accuracy
           for kk = 1:length(perf.lda.accuracy)
               perf.lda.mean_accuracy(kk) = mean(perf.lda.accuracy(:,kk));
               perf.svm.mean_accuracy(kk,1) = mean(perf.svm.accuracy(:,kk));
               perf.logreg.mean_accuracy(kk,1) = mean(perf.logreg.accuracy(:,kk));
           end
           
           %% mean accuracy across lda, svm, logreg:
%                       filename_figures = [mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\searchlight\realigned_data_per_subject\figures\'];
%            if ~exist(filename_figures, 'dir')
%                 mkdir(filename_figures)
%            end
            
                       
           %% figure result mean all classifiers:
           
           for ss=1:size(perf.lda,2) 
               mean_accuracy_classifiers = (mean(perf.lda(ss).accuracy)+mean(perf.svm(ss).accuracy)+mean(perf.logreg(ss).accuracy))./3;
               stat = [];  
               stat.time = str2double(perf.lda(ss).time(1:4)); 
               stat.mvpa.perf = mean_accuracy_classifiers;
               stat.accuracy = mean_accuracy_classifiers;
               stat.label = session.label;
               stat.dimord = 'chan';

               cfg              = [];
               cfg.parameter    = 'accuracy';
               cfg.layout       = '4D248_helmet.mat';            
            %    cfg.xlim         = [0, 0];
               cfg.colorbar     = 'yes';
               cfg.title = ['Searchlight mean all classifiers' perf.lda(ss).time]; 
               figure; ft_topoplotER(cfg, stat); 
               savefig([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\figures\perf_pseudotrials_equalnum_x1_trainfold_only_' perf.lda(ss).time '_ms_mean_all_classifiers.fig'])
           end
          close all
                    
           %% figure result lda:
           for ss=1:size(perf.lda,2) 
               stat = [];  
               stat.time = str2double(perf.lda(ss).time(1:4)); 
               stat.mvpa.perf = mean(perf.lda(ss).accuracy);
               stat.accuracy = mean(perf.lda(ss).accuracy);
               stat.label = session.label;
               stat.dimord = 'chan';

               cfg              = [];
               cfg.parameter    = 'accuracy';
               cfg.layout       = '4D248_helmet.mat';            
            %    cfg.xlim         = [0, 0];
               cfg.colorbar     = 'yes';
%                cfg.title = ['Searchlight LDA' perf.lda(ss).time]; 
               figure; ft_topoplotER(cfg, stat); 
               savefig([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\figures\perf_pseudotrials_equalnum_x1_trainfold_only_' perf.lda(ss).time '_ms_lda.fig'])
           end
           close all
           %% figure result svm:

           for ss=1:size(perf.svm,2) 
               stat = [];  
               stat.time = str2double(perf.svm(ss).time(1:4)); 
               stat.mvpa.perf = mean(perf.svm(ss).accuracy);
               stat.accuracy = mean(perf.svm(ss).accuracy);
               stat.label = session.label;
               stat.dimord = 'chan';

               cfg              = [];
               cfg.parameter    = 'accuracy';
               cfg.layout       = '4D248_helmet.mat';            
            %    cfg.xlim         = [0, 0];
               cfg.colorbar     = 'yes';
               cfg.title = ['Searchlight svm' perf.svm(ss).time]; 
               figure; ft_topoplotER(cfg, stat); 
               savefig([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\figures\perf_pseudotrials_equalnum_x1_trainfold_only_' perf.svm(ss).time '_ms_svm.fig'])
           end
           close all

            %% figure result logreg:

           for ss=1:size(perf.logreg,2) 
               stat = [];  
               stat.time = str2double(perf.logreg(ss).time(1:4)); 
               stat.mvpa.perf = mean(perf.logreg(ss).accuracy);
               stat.accuracy = mean(perf.logreg(ss).accuracy);
               stat.label = session.label;
               stat.dimord = 'chan';

               cfg              = [];
               cfg.parameter    = 'accuracy';
               cfg.layout       = '4D248_helmet.mat';            
            %    cfg.xlim         = [0, 0];
               cfg.colorbar     = 'yes';
               cfg.title = ['Searchlight logreg' perf.logreg(ss).time]; 
               figure; ft_topoplotER(cfg, stat); 
               savefig([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\searchlight\figures\perf_pseudotrials_equalnum_x1_trainfold_only_' perf.svm(ss).time '_ms_logreg.fig'])
           end
           close all
 
  
            %% Konfidenzintervall
        %    for kk = 1:length(perf.lda.accuracy)
        %        perf.lda.ci_acc(kk) = ci(perf.lda.accuracy(:,kk));
        %        perf.svm.ci_acc(kk) = ci(perf.svm.accuracy(:,kk));
        %        perf.logreg.ci_acc(kk) = ci(perf.logreg.accuracy(:,kk));
        %    end

                 
%         end

 

    end
        clear session
    

end




function [perf] = mvpa_leave_out_balldesign(session, balldesign, time_interval, neighbours)

%% Get default hyperparameters for the logreg and lda classifier
param_logreg = mv_get_hyperparameter('logreg');

param_lda = mv_get_hyperparameter('lda');

param_svm = mv_get_hyperparameter('svm');

    %% built train and test data based on CV folds
cf_logreg = cell(length(balldesign), length(session.label));
cf_lda = cell(length(balldesign), length(session.label));
cf_svm = cell(length(balldesign), length(session.label));


%% Crossvalidation folds

CV = adi_crossval_leaveExemplarOut(session.balldesign);
data_trials = kh_trial2dat(session.trial);
%%

for kk = 1:length(CV)
    
    test_fold = data_trials(CV(kk).testset,:,:);
    clabel_test_fold = session.labels(CV(kk).testset);
    % hier überprüfen, ob testset eindeutig einer condition zugeordnet ist:
    % evtl. morgen weitermachen
%     if ismember
%         size(unique(clabel_test_fold),2) > 1
%         C = categorical(clabel_test_fold,[1 2 0],{'like','dislike','dontcare'});
%         if str
%         [group_count, b] = histcounts(C);
%         exist()
%     end
%     
    
    train_fold = data_trials;
    train_fold(CV(kk).testset,:,:) = [];  
    clabel_train_fold = session.labels;
    clabel_train_fold(CV(kk).testset) = [];
    
    %% MVPA 
    
        %% PCA (z-transformation noch einbauen)
%         [pca_tt.coeff,pca_tt.score,pca_tt.latent,pca_tt.tsquared,pca_tt.explained] = pca(train_fold(:,:,tt), 'rows', 'all');
        
%         sum_explained = 0;
%         idx = 0;
%         while sum_explained < 95 
%             idx = idx + 1;
%             sum_explained = sum_explained + pca_tt.explained(idx);
%         end

%         train_pca = train_fold(:,:,tt)* pca_tt.coeff(:,1:idx);
%         test_pca = test_fold(:,:,tt)* pca_tt.coeff(:,1:idx);
        
%%  bootstrap trainfold    
    [trials_trainfold] = kh_trial2dat(train_fold);
 
    trials_trainfold_like = [];
    trials_trainfold_like.trial = trials_trainfold(clabel_train_fold == 1);
    trials_trainfold_like.time = session.time(1:length(trials_trainfold_like.trial));
    trials_trainfold_like.fsample = session.fsample;
    trials_trainfold_like.dimord = 'rpt_chan_time';
    trials_trainfold_like.istimelock = 1;
    trials_trainfold_like.label = session.label;
    
    trials_trainfold_dislike = [];
    trials_trainfold_dislike.trial = trials_trainfold(clabel_train_fold == 2);
    trials_trainfold_dislike.time = session.time(1:length(trials_trainfold_dislike.trial));
    trials_trainfold_dislike.fsample = session.fsample;
    trials_trainfold_dislike.dimord = 'rpt_chan_time';
    trials_trainfold_dislike.istimelock = 1;
    trials_trainfold_dislike.label = session.label;
    
    max_numel_trials = max([numel(trials_trainfold_like.trial), numel(trials_trainfold_dislike.trial)]);
    
    cfg = [];
    cfg.mode = 'bootstrap';
    cfg.averages = 4;
    cfg.repetitions = max_numel_trials*1;
    bootstrap_trainfold_like = fte_subaverage(cfg, trials_trainfold_like);
    bootstrap_trainfold_dislike = fte_subaverage(cfg, trials_trainfold_dislike);
    
    data_bootstrap_trainfold_like = kh_trial2dat(bootstrap_trainfold_like.trial);
    data_bootstrap_trainfold_dislike = kh_trial2dat(bootstrap_trainfold_dislike.trial);
    
    data_bootstrap_trainfold = cat(1, data_bootstrap_trainfold_like, data_bootstrap_trainfold_dislike);
    clabel_train_fold = cat(2, ones(1, size(data_bootstrap_trainfold_like,1)), 2*ones(1, size(data_bootstrap_trainfold_dislike,1)));
    
%     clabel_train_fold(find(clabel_train_fold == 1))

    %% bootstrap testfold: geht nur bei eindeutigen Antworten, deshalb erst einmal weglassen:
%     [trls_testfold] = kh_trial2dat(test_fold);
%  
%     trials_testfold = [];
%     trials_testfold.trial = trls_testfold;
%     trials_testfold.time = session.time(1:length(trials_testfold.trial));
%     trials_testfold.fsample = session.fsample;
%     trials_testfold.dimord = 'rpt_chan_time';
%     trials_testfold.istimelock = 1;
%     trials_testfold.label = session.label;
%     
%     cfg = [];
%     cfg.mode = 'bootstrap';
%     cfg.averages = 4;
%     cfg.repetitions = numel(trials_testfold.trial);
%     bootstrap_testfold = fte_subaverage(cfg, trials_testfold);
%     data_bootstrap_testfold = kh_trial2dat(bootstrap_testfold.trial);
%%

% We want to classify on time window of components 
% time_idx = find(session.time{1} >= comp(1)  &  session.time{1} <= comp(2));
% nachbarschaft 5 cm laut Stefan bzw. noch besser template von Fieldtrip verwenden:
% Zeitkomponente wird gemittelt:

%% Zeitintervalle verkleinert:
% time = 0.03:0.02:0.15;
% time=[0.07 0.1];

for ss = 1:numel(time_interval)-1

time_idx = find(session.time{1} >= time_interval(ss)  &  session.time{1} <= time_interval(ss+1));    
X_train = mean(data_bootstrap_trainfold(:,:,time_idx),3);
X_test = mean(test_fold(:,:,time_idx),3);

 %% LDA - Loop across features
for ff = 1:length(session.label)
    
    current_sensor = session.label(ff);
    current_neighb = neighbours(str2num(current_sensor{1}(2:end))).neighblabel;
    for pp = 1:length(current_neighb)
        neighb_index(pp) = find(strcmp(session.label, current_neighb(pp)));
    end
    nb = [ff neighb_index];

    fprintf('Classifying using feature %d with neighbours %s\n', ff, mat2str(setdiff(nb,ff)))
    
    Xfeat_train = reshape(X_train(:,nb,:), size(X_train,1), []);
    Xfeat_test = reshape(X_test(:,nb,:), size(X_test,1), []);
    
    ind_zero = find(clabel_test_fold==0);
    if ~isempty(ind_zero)
        clabel_test_fold(ind_zero) = 2;
    end
    %% LDA: 
    % Store the current state of the random number generator
    rng_state = rng;
    cf_lda{kk, ff} = train_lda(param_lda, Xfeat_train, clabel_train_fold);
    [predlabel, dval] = test_lda(cf_lda{kk, ff}, Xfeat_test); 
    % To calculate classification accuracy, compare the predicted labels to
    % the true labels and take the mean
    % Calculate AUC and Accuracy
    perf.lda(ss).auc(kk, ff) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
    perf.lda(ss).accuracy(kk, ff) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    perf.lda(ss).time = [num2str(time_interval(ss)) '_' num2str(time_interval(ss+1))];
    
    perf.lda(ss).confusion{kk, ff} = confusionmat(clabel_test_fold, predlabel', 'order', [1, 2]);
    
    if 1==size(perf.lda(ss).confusion{kk, ff},1) && isequal(clabel_test_fold, predlabel') 
        confusion = [];
        confusion(1,1) = perf.lda(ss).confusion{kk, ff};
        confusion(2,2) = 0;
        confusion(1,2)= 0;
        confusion(2,1)= 0;
        perf.lda(ss).confusion{kk, ff} = confusion; 
    end
    TP=perf.lda(ss).confusion{kk, ff}(1,1); 
    TN=perf.lda(ss).confusion{kk, ff}(2,2);
    FP=perf.lda(ss).confusion{kk, ff}(2,1);
    FN=perf.lda(ss).confusion{kk, ff}(1,2);
    F1_Score = 2*TP./(2*TP+FP+FN);
    perf.lda(ss).f1score{kk, ff} =  F1_Score;
    clear F1_Score TP FP FN TN
    
    %% SVM:
    % Store the current state of the random number generator
    rng_state = rng;
    cf_svm{kk, ff} = train_svm(param_svm, Xfeat_train, clabel_train_fold);
    [predlabel, dval] = test_svm(cf_svm{kk,ff}, Xfeat_test);
    % To calculate classification accuracy, compare the predicted labels to
    % the true labels and take the mean

    % Calculate AUC and Accuracy
    perf.svm(ss).auc(kk, ff) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
    perf.svm(ss).accuracy(kk, ff) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    perf.svm(ss).time = [num2str(time_interval(ss)) '_' num2str(time_interval(ss+1))];
    perf.svm(ss).confusion{kk, ff} = confusionmat(clabel_test_fold, predlabel', 'order', [1, 2]);
    TP=perf.svm(ss).confusion{kk, ff}(1,1);
    if 1==size(perf.svm(ss).confusion{kk, ff},1) && isequal(clabel_test_fold, predlabel') 
        confusion = [];
        confusion(1,1) = perf.svm(ss).confusion{kk, ff};
        confusion(2,2) = 0;
        confusion(1,2)= 0;
        confusion(2,1)= 0;
        perf.svm(ss).confusion{kk, ff} = confusion; 
    end
    FP=perf.svm(ss).confusion{kk, ff}(2,1);
    FN=perf.svm(ss).confusion{kk, ff}(1,2);
    F1_Score = 2*TP./(2*TP+FP+FN);
    perf.svm(ss).f1score{kk, ff} =  F1_Score;
    perf.svm(ss).time = [num2str(time_interval(ss)) '_' num2str(time_interval(ss+1))];
    clear F1_Score TP FP FN TN
    
    %% logreg
    % Store the current state of the random number generator
    rng_state = rng;
    cf_logreg{kk, ff} = train_logreg(param_logreg, Xfeat_train, clabel_train_fold);
    [predlabel, dval, prob] = test_logreg(cf_logreg{kk, ff}, Xfeat_test);
        
    % To calculate classification accuracy, compare the predicted labels to
    % the true labels and take the mean

    % Calculate AUC and Accuracy
    perf.logreg(ss).auc(kk, ff) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
    perf.logreg(ss).accuracy(kk, ff) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    perf.logreg(ss).confusion{kk, ff} = confusionmat(clabel_test_fold, predlabel', 'order', [1, 2]);
    TP=perf.logreg(ss).confusion{kk, ff}(1,1);
    if 1==size(perf.logreg(ss).confusion{kk, ff},1) && isequal(clabel_test_fold, predlabel') 
        confusion = [];
        confusion(1,1) = perf.logreg(ss).confusion{kk, ff};
        confusion(2,2) = 0;
        confusion(1,2)= 0;
        confusion(2,1)= 0;
        perf.logreg(ss).confusion{kk, ff} = confusion; 
    end
    FP=perf.logreg(ss).confusion{kk, ff}(2,1);
    FN=perf.logreg(ss).confusion{kk, ff}(1,2);
    F1_Score = 2*TP./(2*TP+FP+FN);
    perf.logreg(ss).f1score{kk, ff} =  F1_Score;   
    perf.logreg(ss).time = [num2str(time_interval(ss)) '_' num2str(time_interval(ss+1))];
    clear F1_Score TP FP FN TN

end

end
    clear Xfeat_test Xfeat_train clabel_train_fold clabel_test_fold neighb_index nb
    
end

 end


function [CV] = adi_crossval_leaveExemplarOut(balldesign_short)

for kk=1:length(balldesign_short) 
    balldesign_array(kk)=balldesign_short{kk};
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
    cfg.demean = 'yes';
    cfg.baselinewindow  = [-0.5 -0.030];
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

% cfg =[];
% cfg.resamplefs = 256;
% % cfg.detrend = 'no';
% [data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq);

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
fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_res_sel);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_data_bpfreq_sel_res, fn_filename);

if ~isempty(diff_fieldnames)
    for k=1:length(diff_fieldnames)
%         data_bpfreq_res_sel.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
        filename.(diff_fieldnames{k}) = data_bpfreq_res_sel.(diff_fieldnames{k});
    end   
    
end
fn_filename = fieldnames(filename);
% fn_filename{end+1}='cfg';

data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, fn_filename);

if ~isfield(data_bpfreq_res_sel, 'additional_cleaning')
        data_bpfreq_res_sel = setfield(data_bpfreq_res_sel, 'additional_cleaning', 'no');
end

right_fieldorder = {'label'; 'time';'trial'; 'ChannelFlag_Bst'; 'additional_cleaning'; 'cfg'; 'dimord'; 'fsample'; 'grad'; 'trialinfo'; 'sampleinfo'};
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, right_fieldorder);

clearvars filename data_bpfreq data_bpfreq_res 

end

%%
%  perf_perm = creating_zero_distr(clabel_train_fold, balldesign{kk}, session, time_idx)
   

%     %% creating null distribution:
% 
%     size_random_dataset = 1000;
%     clabel_train_fold_rand = zeros(size_random_dataset, length(clabel_train_fold));
% 
%     for pp = 1:size_random_dataset
%         clabel_train_fold_rand(pp,:) = clabel_train_fold(randperm(length(clabel_train_fold)));
%     end
% 
%     disp(['creating null distribution for balldesign ' balldesign{kk}])
% 
%     accuracy = nan(size_random_dataset, size(session.time{1},2));
% %     auc = nan(size_random_dataset, size(session.time{1},2));
%     confusion_percentage = cell(size_random_dataset, size(session.time{1},2));
%     confusion = cell(size_random_dataset, size(session.time{1},2));
%     f1score = cell(size_random_dataset, size(session.time{1},2));
% 
%     
%     data_bootstrap_trainfold_tt = mean(data_bootstrap_trainfold(:,:,time_idx),3);
%     test_fold_tt = mean(test_fold(:,:,time_idx),3);
% %         parfor pp = 1:size(clabel_train_fold_rand,1)
% %   lda:
%     parfor pp = 1:size(clabel_train_fold_rand,1)
% %     for pp = 1:size(clabel_train_fold_rand,1)
%         cf_perm_lda = train_lda(param_lda, data_bootstrap_trainfold_tt, clabel_train_fold_rand(pp,:));
%         [predlabel_perm_lda, dval_perm_lda] = test_lda(cf_perm_lda, test_fold_tt); 
% 
%         % Calculate AUC and Accuracy and other performance measures:
%         accuracy_lda(pp, ss) = mv_calculate_performance('acc', 'dval', dval_perm_lda, clabel_test_fold);
% 
%     end
%     
% %     perf_perm.(config.balldesign{kk}).auc = auc;
%     perf_perm.lda.(balldesign{kk})(ss).mean_accuracy = mean(accuracy_lda);    
%     perf_perm.lda.(balldesign{kk})(ss).number_of_trials.testset = numel(clabel_test_fold);  
%     perf_perm.lda.(balldesign{kk})(ss).number_of_trials.trainingsset = numel(clabel_train_fold);
%     
%     %svm:
%     parfor pp = 1:size(clabel_train_fold_rand,1)
% %     for pp = 1:size(clabel_train_fold_rand,1)
%         cf_perm_svm = train_lda(param_lda, data_bootstrap_trainfold_tt, clabel_train_fold_rand(pp,:));
%         [predlabel_perm_svm, dval_perm_svm] = test_lda(cf_perm_svm, test_fold_tt); 
% 
%         % Calculate AUC and Accuracy and other performance measures:
%         accuracy_svm(pp, ss) = mv_calculate_performance('acc', 'dval', dval_perm_svm, clabel_test_fold);
% 
%     end
%  
% 
%     perf_perm.svm.(balldesign{kk})(ss).mean_accuracy = mean(accuracy_svm);    
%     perf_perm.svm.(balldesign{kk})(ss).number_of_trials.testset = numel(clabel_test_fold);  
%     perf_perm.svm.(balldesign{kk})(ss).number_of_trials.trainingsset = numel(clabel_train_fold);
%     
%         %logreg:
%     parfor pp = 1:size(clabel_train_fold_rand,1)
% %     for pp = 1:size(clabel_train_fold_rand,1)
%         cf_perm_logreg = train_lda(param_lda, data_bootstrap_trainfold_tt, clabel_train_fold_rand(pp,:));
%         [predlabel_perm_logreg, dval_perm_logreg] = test_logreg(cf_perm_logreg, test_fold_tt); 
% 
%         % Calculate AUC and Accuracy and other performance measures:
%         accuracy_logreg(pp, ss) = mv_calculate_performance('acc', 'dval', dval_perm_logreg, clabel_test_fold);
% 
%     end
%     
% %     perf_perm.(config.balldesign{kk}).auc = auc;
%     perf_perm.logreg.(balldesign{kk})(ss).mean_accuracy = mean(accuracy_logreg);    
%     perf_perm.logreg.(balldesign{kk})(ss).number_of_trials.testset = numel(clabel_test_fold);  
%     perf_perm.logreg.(balldesign{kk})(ss).number_of_trials.trainingsset = numel(clabel_train_fold);