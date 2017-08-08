% Author: Halle R. Dimsdale-Zucker
%% Setup
initialize_ABCDCon

% setup lengths of variables
num_subjs = length(subjects);

% setup flags
SAVE_OUT_FLAG = 1;
GRAPH_FLAG = 0;
GRAPH_GROUP_FLAG = 1;
ROIS_FLAG = 1;

% roi variables
roi_dirs = {'ashs_left','ashs_right'};
roi_names = {'brCA2_3_DG_body', 'brCA1_body',...
    'brERC','brsubiculum',...
    'brwhole_hippo'};
num_rois = length(roi_names);

% plotting setup
plots_dir = fullfile(dropbox_dir, 'writeups', 'figures');

% let's time to see how long this takes
all_subjects_timer = tic;

%% Loop across subjects
for isub = 1:num_subjs
    % let's also time the length for each subject 
    cur_subj_timer = tic;
    
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    fprintf('------Working on subject %s------\n',b.curSubj)
    
    % make plotting output directory
    if ~isdir(fullfile(plots_dir, b.curSubj))
        mkdir(fullfile(plots_dir, b.curSubj))
    end
    
    % setup variables within "b" that may change between subjects (ie,
    % re-set to defaults in case they've been overwritten from last
    % subject)
    b.beta_matrix_size = [144,152,52];
    
    % run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);
    
    % setup variables that may change between subjects
    num_runs = length(b.runs);
    num_trials = numRecogTrials/num_runs;
    
    num_voxels = (b.beta_matrix_size(1) * b.beta_matrix_size(2) * b.beta_matrix_size(3));
    
    b.cur_subj_fpath = fullfile(analMRIDir,...
        'multivariate_sanityCheck',...
        b.curSubj);

    for irun = 1:num_runs
        b.cur_run = b.runs{irun};
        b.cur_run_fpath = fullfile(b.cur_subj_fpath,...
                b.cur_run);
        cur_run_betas = nan([b.beta_matrix_size,num_trials]);
        
        for itrial = 1:num_trials
            b.cur_trial = ['trial_', num2str(itrial, '%03d')];
            b.cur_trial_fpath = fullfile(b.cur_run_fpath,...
                b.cur_trial);
            b.cur_beta_fpath = fullfile(b.cur_trial_fpath, 'beta_0001.img');
            
            if exist(b.cur_beta_fpath, 'file')
                b.cur_beta = spm_read_vols(spm_vol(b.cur_beta_fpath));
                cur_run_betas(:,:,:,itrial) = b.cur_beta;
                
                % in case the beta changes in size between trials, clear it
                % just to be safe
                clear b.cur_beta
            end %if exist

        end %itrial
        
        if ROIS_FLAG
            %% look at variability across ROIs that might be of interest
            for idir=1:length(roi_dirs)
                b.cur_ROI_dir = [analMRIDir, filesep, b.curSubj, filesep,...
                    'ROIs', filesep, roi_dirs{idir}];

                for iroi=1:num_rois
                    % loop across br*.nii files 
                    % Read in the .nii file
                    b.cur_roi = spm_read_vols(spm_vol([b.cur_ROI_dir,filesep,roi_names{iroi},'.nii'])); 
                    cur_roi_name = [roi_dirs{idir}, '_', roi_names{iroi}];

                    % check to make sure the current ROI
                    % actually has voxels in it ("1" values)
                    % this comes up w/ anterior/posterior ROIs
                    if sum(unique(b.cur_roi))==1
                        % loop across trials to ensure pulling betas for each
                        % trial
                        for itrial = 1:num_trials
                            itrial_cur_run_betas = cur_run_betas(:,:,:,itrial);
                            rois.(cur_roi_name)(:,itrial,irun) = itrial_cur_run_betas(b.cur_roi > 0);
                        end %itrial

                        %% compute summary stats
                        % this is based on the section for how to do this for
                        % all betas (ie, the non-ROI-specific section below)
                        % take mean across trials, but preserve voxels
                        rois.(['mean_',cur_roi_name])(:,irun) = nanmean(rois.(cur_roi_name)(:,:,irun),2);
                        % nanstd(data, flag, dimension)
                        rois.(['sd_',cur_roi_name])(:,irun) = nanstd(rois.(cur_roi_name)(:,:,irun),1,2);
                        x_minus_mean = bsxfun(@minus,rois.(cur_roi_name)(:,:,irun),rois.(['mean_',cur_roi_name])(:,irun));
                        z_betas = bsxfun(@rdivide, x_minus_mean, rois.(['sd_',cur_roi_name])(:,irun));
                        abs_zbetas = abs(z_betas);
                        rois.(['mean_abs_z_',cur_roi_name])(:,irun) = nanmean(abs_zbetas,2);
                        
                        clear x_minus_mean z_betas abs_zbetas
                    end %sum(unique
                end %iroi
            end %idir
        end %if ROIS_FLAG
                
        % clear up variables that get re-used
        clear cur_run_betas cur_run_betas_reshaped x_minus_mean cur_run_zbetas cur_run_abs_zbetas cur_run_abs_zbetas_reshaped        
    end %irun
    
    %% aggregate subjects to put into a group matrix
    % collapse across run dimension (ie, take mean across levels of runs)    
    if ROIS_FLAG
        for idir = 1:length(roi_dirs)
           for iroi = 1:num_rois
              cur_roi_name = [roi_dirs{idir}, '_', roi_names{iroi}];
              % again, take mean across runs dimension
              all_subj_rois(isub).(['mean_',cur_roi_name])(:,1) = mean(rois.(['mean_',cur_roi_name]),2); 
              all_subj_rois(isub).(['mean_abs_z_',cur_roi_name])(:,1) = mean(rois.(['mean_abs_z_',cur_roi_name]),2);
              all_subj_rois(isub).ids(:,1) = {b.curSubj}; 
           end
        end
    end

    %% graph variability across time (ie, across trials)
    if GRAPH_FLAG
        for iroi = 1:num_rois
            % plot as a histogram
            fig = figure('Name', sprintf('%s %s: Histogram of beta means, by run',b.curSubj, roi_names{iroi}));
            figname = fullfile(plots_dir, b.curSubj, sprintf('%s_histogram_by_run_%s', b.curSubj, roi_names{iroi}));
            for irun=1:num_runs
               subplot(num_runs,1,irun)
               histogram(all_subj_rois(isub).(['mean_',cur_roi_name]))
               title(sprintf('%s %s: Histogram of beta means, Run %d',b.curSubj, roi_names{iroi}, irun))
            end
            add_subplot_axis_labels(fig,'Mean beta values', 'Count') %NB - this is a custom function
            saveas(fig,figname,'eps') 
            clear fig

            % TODO: save out a single file for each subject w/ these figures 
            % close figures so don't overload
            close all
        end %iroi
    end %GRAPH_FLAG
    
    % clear subject-level variable before moving onto next subject
    % this will also help deal w/ the issue that s001 has a different
    % matrix size than any of the other subjects
    if ROIS_FLAG
        clear rois
    end
    
    tEnd_subj =  toc(cur_subj_timer);
    fprintf('Processing %s took %d minutes and %f seconds.\n',b.curSubj,floor(tEnd_subj/60),rem(tEnd_subj,60))
end %isub

keyboard
if(ROIS_FLAG)
    % save out mean betas for all subjects
    save(fullfile(plots_dir, 'allSubj_allROI_means.mat'), 'all_subj_rois')
end
keyboard

tEnd_all_subjects = toc(all_subjects_timer);
fprintf('\nRunning all those subjects took %d minutes and %f seconds. \n',floor(tEnd_all_subjects/60),rem(tEnd_all_subjects,60))

%% graph values across subjects
xlim_buffer = 0.75;
xlim_zscore_buffer = 0.01;
if GRAPH_GROUP_FLAG
    if ROIS_FLAG
        % graph ROIs by subject
        % TODO: make this into a function 
        fig = figure('Name','CA23DG body: Mean abs_z_beta values by subject');
        row_counter = 0;
        for isub=1:num_subjs
            row_counter = row_counter + 1;
            % plot the left hemisphere in the left column
            subplot(num_subjs,2,row_counter) 
            histogram(all_subj_rois(isub).mean_abs_z_ashs_left_brCA2_3_DG_body(:,1))
            title(sprintf('%s: left',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_left_brCA2_3_DG)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_left_brCA2_3_DG)) + xlim_zscore_buffer])

            % plot the right hemisphere in the right column
            row_counter = row_counter + 1;
            subplot(num_subjs,2,row_counter)
            histogram(all_subj_rois(isub).mean_abs_z_ashs_right_brCA2_3_DG_body(:,1))
            title(sprintf('%s: right',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_right_brCA2_3_DG)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_right_brCA2_3_DG)) + xlim_zscore_buffer])
        end
        add_subplot_axis_labels(fig,'Mean absolute value z-scored beta','Count')
        saveas(fig,fullfile(plots_dir, sprintf('all_subj_histogram_CA23DG_body')),'eps') 
        clear fig

        fig = figure('Name','CA1 body: Mean abs_z_beta values by subject');
        row_counter = 0;
        for isub=1:num_subjs
            row_counter = row_counter + 1;
            % plot the left hemisphere in the left column
            subplot(num_subjs,2,row_counter) 
            histogram(all_subj_rois(isub).mean_abs_z_ashs_left_brCA1_body(:,1))
            title(sprintf('%s: left',all_subj_rois(isub).ids{:}))

            % plot the right hemisphere in the right column
            row_counter = row_counter + 1;
            subplot(num_subjs,2,row_counter)
            histogram(all_subj_rois(isub).mean_abs_z_ashs_right_brCA1_body(:,1))
            title(sprintf('%s: right',all_subj_rois(isub).ids{:}))

        end
        add_subplot_axis_labels(fig,'Mean absolute value z-scored beta','Count')
        saveas(fig,fullfile(plots_dir, sprintf('all_subj_histogram_CA1_body')),'eps') 
        clear fig
    end %ROIS_FLAG
end