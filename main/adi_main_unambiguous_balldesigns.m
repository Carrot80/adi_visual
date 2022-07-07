%% adi main tests: hier neue Auswertung ausprobiert



%% avg realigned mit eindeutigen Bällen: leider wird cluster permutation test nicht signifikant, im Gegensatz zu allen verwendeten Trials, deshalb hier nicht weiter gemacht
clear
path2subj = 'E:\Arbeit\adidas\data_analysis\visual_stimuli\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
path2inputfile = 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\';
path2save = [];
delete_balldesign_subj = [];
grandavg = []; 
% [grandavg, trl_count] = adi_grandavg_group_clearly_rated(subject_list, grandavg, path2inputfile, 'like');
[grandavg, trl_count] = adi_grandavg_group_clearly_rated(subject_list, grandavg, path2inputfile, 'dislike');

%% statistics grand avg

path2data = 'E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\sensor_space\grandavg_clearly_rated_trials\not_realigned\'
adi_stats_grandavg_clearly_rated(path2data)