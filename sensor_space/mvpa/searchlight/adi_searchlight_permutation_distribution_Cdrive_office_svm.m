function [] = adi_leave_out_exemplar_singlsubj(mainpath, subjectdir, filename, balldesign, time_interval, neighbours)

 path2Cdrive = 'C:\Users\herfurkn1\Documents\adidas\data_analysis\single_subjects\';
 pathname = 'MEG_analysis\realigned_data_per_subject\mvpa\searchlight\';
% for ii = 3:length(subjectdir)
%     if ~exist([path2Cdrive subjectdir(ii).name filesep pathname], 'dir')
%         mkdir([path2Cdrive subjectdir(ii).name filesep pathname])
%     end
%     session_source = [mainpath subjectdir(ii).name '\MEG_analysis\noisereduced\1_95Hz\mvpa\realigned_data_per_subject\session.mat'];
%     session_destination = [path2Cdrive subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\' ]; 
%     copyfile (session_source, session_destination)
% 
% end

for ii = [2:length(subjectdir)]
     
   	load ([path2Cdrive subjectdir(ii).name filesep 'MEG_analysis\realigned_data_per_subject\session.mat'])
        
    %% MVPA  
    
    [perf_perm] = mvpa_leave_out_balldesign(session, balldesign, time_interval, neighbours, subjectdir(ii).name);


%     perf.config.methods = 'bootstrap 2*trialnumber';
%     perf.config.methods = 'bootstrap_pseudotrials';

    perf_perm.number_of_trials = numel(session.trial);
    perf_perm.config.methods = 'perf_without_pca_bootstrap_4pseudotrials_per_condition_trainfold_only';
    perf_perm.config.comment = ' da alle trials mit den uneindeutigen Antworten zur Klassifikation genutzt worden, wurden nur trials aus dem trainfold zur Bildung von Pseudotrials gemittelt; es wurde trialanzahl f�r like und dislike angeglichen ';
    perf_perm.features = session.label;


    %% Konfidenzintervall
%    for kk = 1:length(perf.lda.accuracy)
%        perf.lda.ci_acc(kk) = ci(perf.lda.accuracy(:,kk));
%        perf.svm.ci_acc(kk) = ci(perf.svm.accuracy(:,kk));
%        perf.logreg.ci_acc(kk) = ci(perf.logreg.accuracy(:,kk));
%    end


 save([path2Cdrive subjectdir(ii).name filesep pathname 'permutation_distribution_' num2str(time_interval(1)) '_' num2str(time_interval(2)) 'ms_svm.mat'], 'perf_perm')


    clear session
    

end


end

function [perf_perm] = mvpa_leave_out_balldesign(session, balldesign, time_interval, neighbours, subjname)

%% Get default hyperparameters for the logreg and lda classifier
param_logreg = mv_get_hyperparameter('logreg');

param_lda = mv_get_hyperparameter('lda');

param_svm = mv_get_hyperparameter('svm');

    %% built train and test data based on CV folds
% cf_logreg = cell(length(balldesign), length(session.label));
% cf_lda = cell(length(balldesign), length(session.label));
cf_svm = cell(length(balldesign), length(session.label));


%% Crossvalidation folds

CV = adi_crossval_leaveExemplarOut(session.balldesign);
data_trials = kh_trial2dat(session.trial);
%%

% % switch subjname
% %     
% %     case 'nl_adi_21'
%         indx = 2:length(CV);
%     otherwise
%         indx = 1:length(CV);
% end
indx = 1:length(CV);
for kk = indx
    disp(['creating null distribution for balldesign ' balldesign{kk}])
    test_fold = data_trials(CV(kk).testset,:,:);
    clabel_test_fold = session.labels(CV(kk).testset);
    % hier �berpr�fen, ob testset eindeutig einer condition zugeordnet ist:
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
    cfg.repetitions = max_numel_trials*2;
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

%% Zeitintervalle verkleinern:

% time = 0.03:0.04:0.15;
% time=[0.03 0.15]; % Besprch Stefan 

for ss = 1:numel(time_interval)-1

time_idx = find(session.time{1} >= time_interval(ss)  &  session.time{1} <= time_interval(ss+1));    
X_train = mean(data_bootstrap_trainfold(:,:,time_idx),3);
X_test = mean(test_fold(:,:,time_idx),3);


% perf_perm.lda.(balldesign{kk})(ss).accuracy = zeros(1000, 248);
perf_perm.svm.(balldesign{kk})(ss).accuracy = zeros(1000, 248);
% perf_perm.logreg.(balldesign{kk})(ss).accuracy = zeros(1000, 248);

% perf_perm.lda.(balldesign{kk})(ss).number_of_trials.testset = numel(clabel_test_fold);  
% perf_perm.lda.(balldesign{kk})(ss).number_of_trials.trainingsset = numel(clabel_train_fold);
% perf_perm.lda.(balldesign{kk})(ss).time = [num2str(time_interval(ss)) '_' num2str(time_interval(ss+1))];
 
perf_perm.svm.(balldesign{kk})(ss).number_of_trials.testset = numel(clabel_test_fold);  
perf_perm.svm.(balldesign{kk})(ss).number_of_trials.trainingsset = numel(clabel_train_fold);
perf_perm.svm.(balldesign{kk})(ss).time = [num2str(time_interval(ss)) '_' num2str(time_interval(ss+1))];

% perf_perm.logreg.(balldesign{kk})(ss).number_of_trials.testset = numel(clabel_test_fold);  
% perf_perm.logreg.(balldesign{kk})(ss).number_of_trials.trainingsset = numel(clabel_train_fold);
% perf_perm.logreg.(balldesign{kk})(ss).time = [num2str(time_interval(ss)) '_' num2str(time_interval(ss+1))];

tic
 %% SVM - Loop across features
parfor ff = 1:length(session.label)

    current_sensor = session.label(ff);
    current_neighb = neighbours(str2num(current_sensor{1}(2:end))).neighblabel;
    neighb_index = zeros(1,length(current_neighb));
    numel_current_neighbours = 1:length(current_neighb);
    for pp = numel_current_neighbours
        neighb_index(pp) = find(strcmp(session.label, current_neighb(pp)));
    end
%     clear pp
    nb = [ff neighb_index];

    fprintf('Classifying using feature %d with neighbours %s\n', ff, mat2str(setdiff(nb,ff)))
 
    Xfeat_train = reshape(X_train(:,nb,:), size(X_train,1), []);
    Xfeat_test = reshape(X_test(:,nb,:), size(X_test,1), []);
    
 
    %% creating null distribution:

    size_random_dataset = 1000;
    clabel_train_fold_rand = zeros(size_random_dataset, length(clabel_train_fold));

    for pp = 1:size_random_dataset
        clabel_train_fold_rand(pp,:) = clabel_train_fold(randperm(length(clabel_train_fold)));
    end

   

    %     auc = nan(size_random_dataset, size(session.time{1},2));
    confusion_percentage = cell(size_random_dataset, size(session.time{1},2));
    confusion = cell(size_random_dataset, size(session.time{1},2));
    f1score = cell(size_random_dataset, size(session.time{1},2));
%     accuracy_lda = zeros(size(clabel_train_fold_rand,1),size(time_interval,1));
    
    data_bootstrap_trainfold_tt = mean(data_bootstrap_trainfold(:,:,time_idx),3);
    test_fold_tt = mean(test_fold(:,:,time_idx),3);

% parpool('local',4)

%     for pp = 1:size(clabel_train_fold_rand,1)
% numworkers = [1 2 3 4];
% t_local = zeros(size(numworkers));
% for w = 1:numel(numworkers)
% tic
%     parfor (pp = 1:size(clabel_train_fold_rand,1), numworkers(w))
%     parfor pp = 1:size(clabel_train_fold_rand,1)
%         cf_perm_lda = train_lda(param_lda, data_bootstrap_trainfold_tt, clabel_train_fold_rand(pp,:));
%         [predlabel_perm_lda, dval_perm_lda] = test_lda(cf_perm_lda, test_fold_tt); 
% 
%         % Calculate AUC and Accuracy and other performance measures:
%         accuracy_lda(pp, ss) = mv_calculate_performance('acc', 'dval', dval_perm_lda, clabel_test_fold);
% 
%     end
%     
%     t_local(w) = toc;   

% end
%     perf_perm.(config.balldesign{kk}).auc = auc;
%     perf_perm.lda.(balldesign{kk})(ss).accuracy(:,ff) = accuracy_lda;    
%      clear accuracy_lda
    
    %% svm:
    
     accuracy_svm = zeros(size(clabel_train_fold_rand,1),size(time_interval,1));
   
%     parfor pp = 1:size(clabel_train_fold_rand,1)
    for pp = 1:size(clabel_train_fold_rand,1)
        cf_perm_svm = train_svm(param_svm, data_bootstrap_trainfold_tt, clabel_train_fold_rand(pp,:));
        [predlabel_perm_svm, dval_perm_svm] = test_svm(cf_perm_svm, test_fold_tt); 

        % Calculate AUC and Accuracy and other performance measures:
        accuracy_svm(pp, ss) = mv_calculate_performance('acc', 'dval', dval_perm_svm, clabel_test_fold);

    end
    perf_perm_svm(ff) =   accuracy_svm;
%     perf_perm.svm.(balldesign{kk})(ss).accuracy(:,ff) = accuracy_svm;    
%     clear accuracy_svm
    
%    %% logreg:

%  accuracy_logreg = zeros(size(clabel_train_fold_rand,1),size(time_interval,1));
%    
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
% % %     perf_perm.(config.balldesign{kk}).auc = auc;
%     perf_perm.logreg.(balldesign{kk})(ss).accuracy(:,ff) = accuracy_logreg;    
%     clear accuracy_logreg
end
toc
% perf_perm.svm.(balldesign{kk})

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

if ~isfield(data_bpfreq_res_sel, 'additional_cleaning')
        data_bpfreq_res_sel = setfield(data_bpfreq_res_sel, 'additional_cleaning', 'no');
end

right_fieldorder = {'label'; 'time';'trial'; 'ChannelFlag_Bst'; 'additional_cleaning'; 'cfg'; 'dimord'; 'fsample'; 'grad'; 'trialinfo'};
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, right_fieldorder);

clearvars filename data_bpfreq data_bpfreq_res 

end



