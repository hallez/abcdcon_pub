% Author: Halle R. Dimsdale-Zucker
%% Setup
initialize_ABCDCon

% setup lengths of variables
num_subjs = length(subjects);

% setup flags
SAVE_OUT_FLAG = 1;
GRAPH_FLAG = 0;
GRAPH_GROUP_FLAG = 1;
ROIS_FLAG = 0;

% roi variables
roi_dirs = {'ashs_left','ashs_right'};
roi_names = {'br35','br35_36','br36',...
    'brCA1','brCA2_3_DG',...
    'brERC','brsubiculum',...
    'brwhole_hippo'};
num_rois = length(roi_names);

% let's time to see how long this takes
all_subjects_timer = tic;

%% Loop across subjects
for isub = 1:num_subjs
    % let's also time the length for each subject 
    cur_subj_timer = tic;
    
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    fprintf('------Working on subject %s------\n',b.curSubj)
    
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
        
        %% compute summary stats 

        % compute means
        runs.mean_betas_across_trials(:,:,:,irun) = nanmean(cur_run_betas,4);
        cur_run_betas_reshaped = reshape(cur_run_betas, [num_voxels, num_trials]);
        runs.mean_betas_across_voxels(1,:,irun) = nanmean(cur_run_betas_reshaped);
        
        % compute SDs
        % nanstd(data, flag, dimension)
        runs.sd_betas(:,:,:,irun) = nanstd(cur_run_betas,1,4);
        
        % compute SEM
        % standard deviation / sqrt(length (ie, in this case, the number of
        % trials because that's the dimension the SD is taken across))
        runs.sem_betas(:,:,:,irun) = runs.sd_betas(:,:,:,irun) / sqrt(size(cur_run_betas,4));
        
        % compute zscores, absolute(zscores), and mean(absolute(zscore))
        % since have NAN values, easiest to manually compute
        % bsxfun comes from: 
        % http://stackoverflow.com/questions/13402637/how-to-subtract-each-item-of-a-matrix-from-each-coressponding-row-of-another-mat
        x_minus_mean = (bsxfun(@minus,cur_run_betas, runs.mean_betas_across_trials(:,:,:,irun)));
        cur_run_zbetas = bsxfun(@rdivide, x_minus_mean, runs.sd_betas(:,:,:,irun));
        cur_run_abs_zbetas = abs(cur_run_zbetas);
        % reshape so have a num_voxels by num_trials matrix
        cur_run_abs_zbetas_reshaped = reshape(cur_run_abs_zbetas, [num_voxels,num_trials]);
        % take mean across voxel dimension
        runs.mean_abs_zbetas_across_voxels(1,:,irun) = nanmean(cur_run_abs_zbetas_reshaped);
        
        % clear up variables that get re-used
        clear cur_run_betas cur_run_betas_reshaped x_minus_mean cur_run_zbetas cur_run_abs_zbetas cur_run_abs_zbetas_reshaped        
    end %irun
    
    %% aggregate subjects to put into a group matrix
    % collapse across run dimension (ie, take mean across levels of runs)
    % should now have a trials x subjects matrix
    all_subj.mean_betas_across_voxels(:,isub) = mean(runs.mean_betas_across_voxels(:,:,:),3);
    all_subj.mean_abs_zbetas_across_voxels(:,isub) = mean(runs.mean_abs_zbetas_across_voxels(:,:,:),3); 
    all_subj.ids(:,isub) = {b.curSubj};
    
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
        %% plot means
        fig = figure('Name',sprintf('%s: Beta means across trials, by run',b.curSubj));
        for irun=1:num_runs
           subplot(num_runs,1,irun)
           bar(runs.mean_betas(:,:,:,irun))
        end
        % NB: this is an HRZ custom function
        add_subplot_axis_labels(fig,'Trial Number','Mean beta value')
        % clear `fig` variable before move onto next plot
        clear fig 

        % plot as a histogram
        fig = figure('Name', sprintf('%s: Histogram of beta means, by run',b.curSubj));
        for irun=1:num_runs
           subplot(num_runs,1,irun)
           histogram(runs.mean_betas(:,:,:,irun))
        end
        add_subplot_axis_labels(fig,'Mean beta values', 'Count')
        clear fig
        
        %% plot SDs
        fig = figure('Name',sprintf('%s: Beta means and SDs across trials, by run',b.curSubj));
        for irun=1:num_runs
           subplot(num_runs,1,irun)
           hold on
           bar(runs.mean_betas(:,:,:,irun))
           errorbar(runs.mean_betas(:,:,:,irun), runs.sd_betas(:,:,:,irun),'r.')
        end
        hold off
        add_subplot_axis_labels(fig,'Trial Number','Mean beta value')
        clear fig
        
        %% plot SEMs
        fig = figure('Name',sprintf('%s: Beta means and SEMs across trials, by run',b.curSubj));
        for irun=1:num_runs
           subplot(num_runs,1,irun)
           hold on
           bar(runs.mean_betas(:,:,:,irun))
           errorbar(runs.mean_betas(:,:,:,irun), runs.sem_betas(:,:,:,irun),'r.')
        end
        hold off
        add_subplot_axis_labels(fig,'Trial Number','Mean beta value')
        clear fig
        
        %% plot zscores
        % plot histogram
        fig = figure('Name',sprintf('%s: Histogram of abs(z_beta) means, by run',b.curSubj));
        for irun=1:num_runs
           subplot(num_runs,1,irun)
           histogram(runs.mean_abs_zbetas_across_voxels(:,:,irun))
        end
        add_subplot_axis_labels(fig,'Mean absolute value z-scored beta value','Count')
        clear fig
        
        % plot barplot
        fig = figure('Name',sprintf('%s: abs(z_beta) means, by run',b.curSubj));
        for irun=1:num_runs
           subplot(num_runs,1,irun)
           bar(runs.mean_abs_zbetas_across_voxels(:,:,irun))
        end
        add_subplot_axis_labels(fig,'Trial number','Mean absolute value z-scored beta value')
        clear fig
        
        % TODO: save out a single file for each subject w/ these figures 
        % close figures so don't overload
        close all
    end %GRAPH_FLAG
   
    if SAVE_OUT_FLAG
        % save out run-level summary info for the current subject
        fname_out = sprintf('%s_beta_timeseries_summary',b.curSubj);
        save([b.cur_subj_fpath,filesep,fname_out],'runs')

        if ROIS_FLAG
            % save out ROI information
            roi_fname_out = sprintf('run%03d_beta_summary_stats_by_rois',irun);
            save([b.cur_subj_fpath,filesep,roi_fname_out],'rois');
        end
    end 
    
    % clear subject-level variable before moving onto next subject
    % this will also help deal w/ the issue that s001 has a different
    % matrix size than any of the other subjects
    clear runs 
    if ROIS_FLAG
        clear rois
    end
    
    tEnd_subj =  toc(cur_subj_timer);
    fprintf('Processing %s took %d minutes and %f seconds.\n',b.curSubj,floor(tEnd_subj/60),rem(tEnd_subj,60))
end %isub

tEnd_all_subjects = toc(all_subjects_timer);
fprintf('\nRunning all those subjects took %d minutes and %f seconds. \n',floor(tEnd_all_subjects/60),rem(tEnd_all_subjects,60))


%% compute summaries across subjects 
all_subj.grand_mean = nanmean(all_subj.mean_betas_across_voxels,2);
all_subj.grand_mean_abs_zbetas = nanmean(all_subj.mean_abs_zbetas_across_voxels,2);

%% graph values across subjects
xlim_buffer = 0.75;
xlim_zscore_buffer = 0.01;
if GRAPH_GROUP_FLAG
    fig = figure('Name','Mean beta values by subject');
    for isub=1:num_subjs
       subplot(num_subjs,1,isub) 
       histogram(all_subj.mean_betas_across_voxels(:,isub))
       title(sprintf('%s',all_subj.ids{isub}))
       % figure out the xlimits across all subjects so consistent
       % min(min()) or max(max())
       % add some buffer so not butting data up to edge of graph
       xlim([min(min(all_subj.mean_betas_across_voxels)) - xlim_buffer,...
           max(max(all_subj.mean_betas_across_voxels)) + xlim_buffer])
    end
    add_subplot_axis_labels(fig,'Mean beta value','Count')
    clear fig
    
    fig = figure('Name','Mean abs_z_beta values by subject');
    for isub=1:num_subjs
       subplot(num_subjs,1,isub) 
       histogram(all_subj.mean_abs_zbetas_across_voxels(:,isub))
       title(sprintf('%s',all_subj.ids{isub}))
       xlim([min(min(all_subj.mean_abs_zbetas_across_voxels)) - xlim_zscore_buffer,...
           max(max(all_subj.mean_abs_zbetas_across_voxels)) + xlim_zscore_buffer])
    end
    add_subplot_axis_labels(fig,'Mean absolute value z-scored beta','Count')
    clear fig
    
    fig = figure('Name',sprintf('Grand means across %d subjects: Betas',num_subjs));
    histogram(all_subj.grand_mean)
    add_subplot_axis_labels(fig,'Grand mean beta value','Count')
    clear fig
    
    fig = figure('Name',sprintf('Grand means across %d subjects: abs_zBetas',num_subjs));
    histogram(all_subj.grand_mean_abs_zbetas)
    add_subplot_axis_labels(fig,'Grand mean absolute value z-scored beta value','Count')
    clear fig
    
    if ROIS_FLAG
        %% graph ROIs by subject
        % TODO: make this into a function 
        fig = figure('Name','Whole_hippo: Mean abs_z_beta values by subject');
        row_counter = 0;
        for isub=1:num_subjs
            row_counter = row_counter + 1;
            % plot the left hemisphere in the left column
            subplot(num_subjs,2,row_counter) 
            histogram(all_subj_rois(isub).mean_abs_z_ashs_left_brwhole_hippo(:,1))
            title(sprintf('%s: left',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_left_brwhole_hippo)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_left_brwhole_hippo)) + xlim_zscore_buffer])

            % plot the right hemisphere in the right column
            row_counter = row_counter + 1;
            subplot(num_subjs,2,row_counter)
            histogram(all_subj_rois(isub).mean_abs_z_ashs_right_brwhole_hippo(:,1))
            title(sprintf('%s: right',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_right_brwhole_hippo)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_right_brwhole_hippo)) + xlim_zscore_buffer])
        end
        add_subplot_axis_labels(fig,'Mean absolute value z-scored beta','Count')
        clear fig

        fig = figure('Name','CA23DG: Mean abs_z_beta values by subject');
        row_counter = 0;
        for isub=1:num_subjs
            row_counter = row_counter + 1;
            % plot the left hemisphere in the left column
            subplot(num_subjs,2,row_counter) 
            histogram(all_subj_rois(isub).mean_abs_z_ashs_left_brCA2_3_DG(:,1))
            title(sprintf('%s: left',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_left_brCA2_3_DG)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_left_brCA2_3_DG)) + xlim_zscore_buffer])

            % plot the right hemisphere in the right column
            row_counter = row_counter + 1;
            subplot(num_subjs,2,row_counter)
            histogram(all_subj_rois(isub).mean_abs_z_ashs_right_brCA2_3_DG(:,1))
            title(sprintf('%s: right',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_right_brCA2_3_DG)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_right_brCA2_3_DG)) + xlim_zscore_buffer])
        end
        add_subplot_axis_labels(fig,'Mean absolute value z-scored beta','Count')
        clear fig

        fig = figure('Name','CA1: Mean abs_z_beta values by subject');
        row_counter = 0;
        for isub=1:num_subjs
            row_counter = row_counter + 1;
            % plot the left hemisphere in the left column
            subplot(num_subjs,2,row_counter) 
            histogram(all_subj_rois(isub).mean_abs_z_ashs_left_brCA1(:,1))
            title(sprintf('%s: left',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_left_brCA1)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_left_brCA1)) + xlim_zscore_buffer])

            % plot the right hemisphere in the right column
            row_counter = row_counter + 1;
            subplot(num_subjs,2,row_counter)
            histogram(all_subj_rois(isub).mean_abs_z_ashs_right_brCA1(:,1))
            title(sprintf('%s: right',all_subj_rois(isub).ids{:}))
%             % set xlimits so consistent
%             xlim([min(min(all_subj_rois.mean_abs_z_ashs_right_brCA1)) - xlim_zscore_buffer,...
%                max(max(all_subj_rois.mean_abs_z_ashs_right_brCA1)) + xlim_zscore_buffer])
        end
        add_subplot_axis_labels(fig,'Mean absolute value z-scored beta','Count')
        clear fig
    end %ROIS_FLAG
end