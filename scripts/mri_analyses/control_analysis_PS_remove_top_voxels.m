% Read in existing pattern matrices created by 
% `RSA_btwn_runs_exclude_outlier_trials.m`, 
% NaN out top N voxels, and write out trial pattern correlation matrices
%
% Halle R. Dimsdale-Zucker

% setup 
initialize_ABCDCon

roi_dirs = {'ashs_left','ashs_right'};
num_vox = 5; % this is the number of voxels to be removed from PS
rois = {'brCA1_body', 'brCA2_3_DG_body'};

% pick removal method
% options: 'mean', 'sd', 'cond_diff'
REMOVE_METHOD_FLAG = 'cond_diff';
if strcmp(REMOVE_METHOD_FLAG, 'mean')
    remove_str = '_by_mean';
elseif strcmp(REMOVE_METHOD_FLAG, 'sd')
    remove_str = '_by_SD';
elseif strcmp(REMOVE_METHOD_FLAG, 'cond_diff')
    remove_str = '_by_condition_diff';
else
    fprintf('Voxel removal method does not match possible file types.')
    break;
end


%  loop across subject
for isub=1:length(subjects)

    fprintf('\n----Working on s%s----\n',num2str(subjects(isub),'%03d'))
    subj_timer = tic;

    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.analMRIDir = analMRIDir;
    b.modelDir = [analMRIDir,'multivariate_sanityCheck'];
        
    % Run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);

    % read in current subjects' trial labels file so can identify bad
    % trials, etc.
    trial_lbls_fname = fullfile(b.modelDir, b.curSubj, [b.curSubj,'_trial_pattern_ids_all_runs.mat']);
    % this will load in a structure called `ids`
    if(exist(trial_lbls_fname, 'file'))
        load(trial_lbls_fname);
        
        for idir = 1:length(roi_dirs)
            for iroi=1:length(rois)
                b.cur_ROI_dir = fullfile(analMRIDir, b.curSubj,'ROIs', roi_dirs{idir});
                cur_roi = rois{iroi};
                cur_hemi = roi_dirs{idir}(strfind(roi_dirs{idir},'_')+1:end);

                % read in existing pattern matrix
                % this gets created by `RSA_btwn_runs_exclude_outlier_trials.m`
                pattern_mtx_fname = fullfile(b.cur_ROI_dir, sprintf('%s_pattern_mtx_no_outlier_trials_all_runs.mat', cur_roi));
                if(exist(pattern_mtx_fname, 'file'))
                    load(pattern_mtx_fname)

                    % read in voxels to exclude  
                    % this gets created by `control_analysis_drop_voxels.R`
                    exclude_voxels_fname = fullfile(b.analMRIDir, b.curSubj, sprintf('%s_%s_top_%d_voxels%s.csv', cur_hemi, cur_roi, num_vox, remove_str));                        

                    
                    if(exist(exclude_voxels_fname, 'file'))
                        exclude_voxels = csvread(exclude_voxels_fname);
                        
                        % this is based on the approach in `RSA_btwn_runs_exclude_outlier_trials.m`
                        % because reading in matrices where bad trials and voxels have
                        % already been removed, no need to do this again
                        % remove rows for voxels to drop
                        pattern_truncated = pattern_all_runs;
                        fname_corr_out = cellstr(fullfile(b.cur_ROI_dir,sprintf('%s_pattern_corr_no_outlier_trials_%02d_truncated_voxels%s_all_runs', cur_roi, num_vox, remove_str)));
                        pattern_truncated(exclude_voxels, :) = [];

                        % take correlation across all voxels and all trials
                        % and save this out too 
                        pattern_corr = corr(pattern_truncated);
                        save(fname_corr_out{:},'pattern_corr')
                    else
                        sprintf('\nExclude voxels file %s does not exist.', exclude_voxels_fname)
                        continue;
                    end %exist(exclude_voxels_fname                     
                else
                    sprintf('\nPattern matrix file %s does not exist.', pattern_mtx_fname)
                    continue;
                end % exist(pattern_mtx_fname
                clear pattern_all_runs pattern_corr
            end %iroi
        end %idir
    else
        sprintf('\n Labels file for subject %s does not exist. Skipping.', b.curSubj)
        continue;
    end %exist(trial_lbls_fname
end %isub

