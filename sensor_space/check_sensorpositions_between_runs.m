 
function [] = check_runs(subject_list, path2folder)

% Coil Positions => funktioniert nicht, da in hs files immer gleiche
% Postion aus Digitalisierung vorhanden ist
% for ii = 2:length(subject_list)
% 
%     dir_hs_files = dir([subject_list(ii).folder filesep subject_list(ii).name filesep path2folder filesep 'hs_file_*']);
%     
%     ind=[];
%     for kk = 1:length(dir_hs_files)
%         ind(kk) = endsWith(dir_hs_files(kk).name,'.fig');
%     end
%     if any(ind)
%         dir_hs_files(find(ind))=[];
%     end
%     
%     for kk=1:length(dir_hs_files)
%            shape = ft_read_headshape([dir_hs_files(kk).folder filesep dir_hs_files(kk).name], 'format', '4d_hs');
%            runs(kk).pos = shape.fid.pos;
%     end
%     euclid_dist = [];
%     for ll = 1:5 % 5 coils
%         p = [5.22	6.99	1.3];
%         q = [5.25	6.97	1.33]; % immer gleiche Position 
%         coils_run1 = runs(1).pos;
%         coils_run2 = runs(2).pos;
%         for ii=1:size(q,1)
%             euclid_dist(ii).run(ll) = sqrt( (coils_run2(:,1)-coils_run1(:,1)).^2 + (coils_run2(:,2)-coils_run1(:,2)).^2 + (coils_run2(:,3)-coils_run1(:,3)).^2);
%         end
%     end
% end
%  

 
for ii = 6:length(subject_list)

    dir_hs_files = dir([subject_list(ii).folder filesep subject_list(ii).name filesep path2folder filesep 'hs_file_*']);
    
    ind=[];
    for kk = 1:length(dir_hs_files)
        ind(kk) = endsWith(dir_hs_files(kk).name,'.fig');
    end
    if any(ind)
        dir_hs_files(find(ind))=[];
    end
    
    figure; hold on;
    for kk=1:length(dir_hs_files)
        
        if exist([dir_hs_files(kk).folder filesep 'Neu_Like500_' dir_hs_files(kk).name(end) '.mat'], 'file')
            shape = ft_read_headshape([dir_hs_files(kk).folder filesep dir_hs_files(kk).name], 'format', '4d_hs');
            load([dir_hs_files(kk).folder filesep 'Neu_Like500_' dir_hs_files(kk).name(end) '.mat'], 'cleanMEG')
            ft_plot_headshape(shape)
            hold on 
            plot3(cleanMEG.grad.chanpos(1:248,1), cleanMEG.grad.chanpos(1:248,2), cleanMEG.grad.chanpos(1:248,3), 'o')      

         % print(gcf, '-append','-dpsc2', ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\headpositions.ps']);
            clear shape cleanMEG
        end
        
    end
    legend('run1', 'run2', 'run3')
    ax = gca;
    set(gca, 'CameraPosition', [-0.242 -2.599 0.105]);  
    title([subject_list(ii).name 'all_runs' ])
    savefig ([dir_hs_files(kk).folder filesep 'Sensorpositions_all_runs.fig']);  
    hold off
%     close
end
 
% 
%  
% 
% for ii=1:length(subject_list)
%         filename = dir([subject_list(i).folder filesep subject_list(i).name filesep  path2inputfile 'Neu_Like*.mat']);
% 
%         for ke=1:length(filename)
%             load ([filename(k).folder filesep filename(k).name])
%             grad_all_runs(kk) = cleanMEG_interp.grad;
%             % delete balldesigns which were not clearly rated:
%             delete_trials = zeros(1, length(cleanMEG_interp.trial));
%             for p=1:length(cleanMEG_interp.trial)
%                 trial=cleanMEG_interp.trialinfo.balldesign_short{1,p}{1,1};
%                 % check if trial should be deleted:
%                 if 1==delete_run.(subjectpath(i).name).(['run' filename(k).name(end-4)]).(trial)
%                    delete_trials(p) = 1;
%                 end
%             end
%             
%             
%             
%         end
%         
%         
% end