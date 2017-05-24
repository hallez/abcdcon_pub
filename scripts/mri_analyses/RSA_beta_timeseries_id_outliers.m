% Author: Halle R. Dimsdale-Zucker
%% Setup
initialize_ABCDCon

% variables
num_subjs = length(subjects);

%% Determine outlier cutoffs
% These come from visual inspection of the graphs created in `beta_timeseries_graphs.m`
outlier_zmin = 0.7;
outlier_zmax = 0.85;

%% Loop across subjects
for isub=1:num_subjs
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    fprintf('------Working on subject %s------\n',b.curSubj)
    
    % run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);
    
    % set up variables that may change between subjects
    num_runs = length(b.runs);
    
    b.cur_subj_fpath = fullfile(analMRIDir,...
        'multivariate_sanityCheck',...
        b.curSubj);
    
    %% read in that subjects' summary file (these are written out in `RSA_beta_timeseries_graphs.m`)
    cur_sub_betas_fname = [b.cur_subj_fpath,filesep,sprintf('%s_beta_timeseries_summary.mat',b.curSubj)];
    
    if exist(cur_sub_betas_fname, 'file')==2
        % this will create a structure variable
        % cur_sub_betas.runs
        cur_sub_betas = load(cur_sub_betas_fname);
    else
        fprintf('No betas summary file for %s.\n', b.curSubj)
    end
    
    for irun=1:num_runs
        
        cur_run_id = b.runs{irun};
        cur_run_num = str2double(strtok(cur_run_id, 'run'));
        
        %% identify bad timepoints
        cur_run_outliers = find(outlier_zmin > cur_sub_betas.runs.mean_abs_zbetas_across_voxels(:,:,irun) | cur_sub_betas.runs.mean_abs_zbetas_across_voxels(:,:,irun) > outlier_zmax);
        
        if isempty(cur_run_outliers)
            cur_run_outliers = 0; 
        end

        %% save out badtimepoints to .txt file
        fname_bad_timepoints = sprintf('%s_bad_timepoints_run%02d.dat', b.curSubj,cur_run_num);
        csvwrite([b.cur_subj_fpath, filesep, fname_bad_timepoints], cur_run_outliers)
        
        % clear any variables that may change between loop iterations
        clear cur_run_outliers 
    
    end %irun

    %% clear subject-specific variables
    clear cur_sub_betas_fname cur_sub_betas
    
end %isub =
