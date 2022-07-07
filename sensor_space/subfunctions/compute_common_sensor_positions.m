function [grad_common] = compute_common_sensor_positions(subject_list)
% compute_common_sensor_positions reads in all grad-fields and appends them
% in a struct array; appendstruct is a fieldtrip-subfunction; grad-files
% are averaged; code was taken from ft_megrealign
% to do: je nachdem, ob noch weitere runs/probanden rausgeworfen werden 


grad_all_subj=struct([]);

    for ii = 1:length(subject_list)
        %% like
        filename = dir([subject_list(ii).folder filesep subject_list(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\02_interpolated\Neu_Like*.mat']);
        for kk = 1:length(filename)
            load ([filename(kk).folder filesep filename(kk).name])
            grad_all_subj = appendstruct(grad_all_subj, cleanMEG_interp.grad); 
            clear cleanMEG_interp
         end     
        
    end
    
mean_grad = ft_average_sens(grad_all_subj);
template.mean_grad = mean_grad;
template.comment = 'without nl_adi_04 has different order of_sensor positions';
save('E:\Arbeit\adidas\data_analysis\visual_stimuli\group_analysis\template_gradfile\template.mat');





end




