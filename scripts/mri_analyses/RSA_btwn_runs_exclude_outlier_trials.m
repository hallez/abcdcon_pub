% Loop across ROIs, extract (mask) single-trial betas within these ROIs
% unless the current trial has been identified as an outlier,
% and write out voxel x trial pattern matrices and trial pattern
% correlation matrices
%
% Halle R. Dimsdale-Zucker

%% Setup some basics
initialize_ABCDCon

roi_dirs = {'ashs_left','ashs_right'};

%% Loop across subjects 
for isub=1:length(subjects)

    fprintf('\n----Working on s%s----\n',num2str(subjects(isub),'%03d'))
    subj_timer = tic;

    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.analMRIDir = analMRIDir;
    b.modelDir = [analMRIDir,'multivariate_sanityCheck'];
    % if last subject overwrote b.runs, revert back 
    runbase = 'run';
    runnums = 1:4;
    nruns = length(runnums);
    b.runs = cell(1,nruns);
    for irun=1:nruns
        b.runs(irun) = cellstr([runbase,num2str(runnums(irun))]);
    end %for irun=
    
    % read in current subjects' trial labels file so can identify bad
    % trials, etc.
    trial_lbls_fname = fullfile(b.modelDir, b.curSubj, [b.curSubj,'_trial_pattern_ids_all_runs.mat']);
    % this will load in a structure called `ids`
    load(trial_lbls_fname);
    
    % Run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);
    
    %% Loop across ROI directories
    for idir=1:length(roi_dirs)
    
        fprintf('Working on ROI dir:%s.\n', roi_dirs{idir})
        
        % loop across ROI directories
        b.cur_ROI_dir = [analMRIDir, filesep, b.curSubj,filesep,'ROIs', filesep, roi_dirs{idir}];
        b.all_ROIs = dir(b.cur_ROI_dir);
        b.all_ROIs = {b.all_ROIs.name};
        % make sure only getting binarized ROIs
        idx = regexp(b.all_ROIs,'^b*');
        idx = ~cellfun(@isempty,idx);
        b.all_ROIs = {b.all_ROIs{idx}};
        % make sure only getting ROIs (this is critical if re-running the
        % script since you'll have br*.mat files
        idx2 = regexp(b.all_ROIs,'\>.nii');
        idx2 = ~cellfun(@isempty,idx2);
        b.all_ROIs = {b.all_ROIs{idx2}};
        
        % remove other ROIs we're never going to analyze to speed up the script
        rois_to_include = 1:length(b.all_ROIs);
        for iroi = 1:length(b.all_ROIs)
            % remove "garbage" ROIs that ASHS produces
            if contains_str(b.all_ROIs{iroi}, 'zeros') ||  contains_str(b.all_ROIs{iroi}, 'Clear_Label') || contains_str(b.all_ROIs{iroi}, 'MISC') 
                fprintf('Contains one of zeros, Clear_Label, MISC. ROI is %s.\n', b.all_ROIs{iroi})
                rois_to_include(iroi) = 0;
            end
        end
        
        % make rois_to_include into a logical index
        rois_to_include = logical(rois_to_include);
        % subset so just get ROIs we want
        b.all_ROIs = {b.all_ROIs{rois_to_include}};
        
        %% Loop across ROIs w/in the current ROI directory 
        for iroi=1:length(b.all_ROIs) 
            
            fprintf('Working on ROI %d of %d ROI: %s\n',iroi,length(b.all_ROIs),b.all_ROIs{iroi})
            b.cur_roi=spm_vol([b.cur_ROI_dir,filesep,b.all_ROIs{iroi}]);
            b.cur_roi.img = spm_read_vols(b.cur_roi);
            [path, roi_name, ext] = fileparts(b.cur_roi.fname);
            
            %% Loop across runs
            % initialize a column counter 
            % this will be used to make sure that the current run's pattern
            % gets put into the correct place in the across run matrix
            col_counter = 0;
            
            % initialize a matrix of NaNs where can put pattern
            % matrices from each run
            % WILL ENCOUNTER PROBLEMS WHEN RUNS HAVE DIFFERENT NUMBERS
            % OF ROWS (ie, if needed to delete a row due to NaN values)
            num_trials = 63;
            pattern_all_runs = nan(size(find(b.cur_roi.img > 0),1),num_trials*length(b.runs));
            pattern_ids_all_runs = cell(1,num_trials*length(b.runs));
            z_pattern_all_runs = nan(size(find(b.cur_roi.img > 0),1),num_trials*length(b.runs));
            trial_means_all_runs = nan(1,num_trials*length(b.runs));
                        
            for irun=1:length(b.runs)

                b.curRun = b.runs{irun};
                fprintf('Working on: %s\n',b.curRun);
                
                %% Loop across trials
                curRunDir = [b.modelDir,filesep,b.curSubj,filesep,b.curRun];                
                curDirs = dir(curRunDir);
                trialCounter = 0;

                % figure out how many trials were in the current run (should always
                % be 63, but this is coded for generalizability) 
                for j = 1:size(curDirs,1)
                    s=curDirs(j).name;
                    if strfind(s,'trial_')
                        trialCounter = trialCounter + 1;
                    end %if strfind
                end %j=

                % check to ensure this is a reasonable amount of trials
                if trialCounter <= (numRecogTrials/nruns)
                    numCurTrials = 1:trialCounter;
                else 
                    error('Implausible number of trials for current run: %s.\n',b.curRun)
                end %if trialCounter

                % initialize a variable the size of the current ROI
                % (rows) x number of trials
                cur_pattern = nan(size(find(b.cur_roi.img > 0),1),trialCounter);
                % initialize a variable for the trial IDs 
                cur_pattern_trial_ids = cell(1,trialCounter);

                for itrial=1:trialCounter

                    b.curTrial = ['trial_',num2str(numCurTrials(itrial),'%03d')];

                    % read in the beta_0001.img for the current trial and
                    % run
                    b.curTrial_dir = [b.modelDir,filesep,b.curSubj,filesep,b.curRun,filesep,b.curTrial];
                    b.cur_beta=spm_vol([b.curTrial_dir,filesep,'beta_0001.img']);
                    b.cur_beta.img = spm_read_vols(b.cur_beta);

                    % apply current ROI mask to current beta image,
                    % save results into a voxel x trial pattern matrix
                    cur_pattern(:,itrial) = b.cur_beta.img(b.cur_roi.img > 0);

                    % figure out trial-specific label
                    % NB: could now do this from the IDs file 
                    [TOKEN, REMAIN] = strtok(cellstr(b.cur_beta.descrip),'R');
                    [TOKEN, REMAIN] = strtok(REMAIN,'*');
                    cur_pattern_trial_ids(:,itrial) = TOKEN; 

                end %itrial=  

                if size(cur_pattern,1) > 0
                    % put the current pattern (and ids) into a variable
                    % across runs
                    pattern_all_runs(:,col_counter + 1:col_counter + num_trials) = cur_pattern;
                    pattern_ids_all_runs(1,col_counter + 1:col_counter + num_trials) = cur_pattern_trial_ids;
                    
                    % z-score w/in an roi, w/in a run
                    % this is essentially a way of mean-centering (b/c
                    % z_mean = 0)
                    % 0). this is one way of removing the overall univariate
                    % effect from the patterns
                    % we also have to deal w/ what happens when there are
                    % NaN values
                    % code based on: http://www.mathworks.com/matlabcentral/answers/105736-zscore-of-an-array-with-nan-s
                    % bsxfun based on: http://stackoverflow.com/questions/11027457/z-score-with-nan-values-in-matlab-vectorized
                    if any(isnan(cur_pattern(:)))
                        xmu=nanmean(cur_pattern);
                        xsigma=nanstd(cur_pattern);
                        x_minus_mean = bsxfun(@minus, cur_pattern, xmu);
                        z_cur_pattern = bsxfun(@rdivide, x_minus_mean, xsigma); 
                    else
                        z_cur_pattern = zscore(cur_pattern);
                    end
                    
                    z_pattern_all_runs(:, col_counter + 1:col_counter + num_trials) = z_cur_pattern;
                    
                    % take the mean across all voxels w/in a given trial
                    % we will use this when we want to remove the
                    % univariate effects from the mixed models
                    trial_means = nanmean(cur_pattern);
                    trial_means_all_runs(1, col_counter + 1:col_counter + num_trials) = trial_means;
                end %if size(pattern,1) > 0
                
                col_counter = col_counter + num_trials; 
                
                % clear variables that change for each run loop
                clear cur_pattern cur_pattern_trial_ids z_cur_pattern

            end %irun=
           
            % save out pattern matrix
            % deal w/ when have no patterns
            % will happen for zeros ROIs
            if size(pattern_all_runs,1) > 0
                % first, let's get rid of rows w/ bad voxels and then NaN
                % out trials that need to be removed (either on the basis
                % of behavior or being a "bad" beta)
                % we'll be naughty and just keep overwriting
                % `pattern_all_runs`
                %
                % NB: this elimination of bad voxels (rows) and bad trials 
                % (columns) is done in a 2-step process so don't end up in
                % a situation where have NaN values in rows AND columns
                % because then when go to compute correlations, this would
                % result in a matrix of NaNs
                %
                % based on: http://www.mathworks.com/matlabcentral/answers/68510-remove-rows-or-cols-whose-elements-are-all-nan
                pattern_all_runs = pattern_all_runs(all(~isnan(pattern_all_runs),2),:);
                z_pattern_all_runs = z_pattern_all_runs(all(~isnan(z_pattern_all_runs),2),:);
                
                % now, NaN out bad trials 
                % we do NOT want to eliminate these columns b/c then this
                % will get tricky when want to match up with trial IDs
                pattern_all_runs(:,logical(ids.mem_resp_scored_numeric == 99)) = NaN;
                z_pattern_all_runs(:,logical(ids.mem_resp_scored_numeric == 99)) = NaN;
                
                % NOW, we can finally save some stuff out
                fname_pattern_out = cellstr([b.cur_ROI_dir,filesep,roi_name,'_pattern_mtx_no_outlier_trials_all_runs']);
                save(fname_pattern_out{:},'pattern_all_runs')
                
                fname_ids_out = cellstr([b.cur_ROI_dir,filesep,roi_name,'_pattern_mtx_ids_no_outlier_trials_all_runs']);
                save(fname_ids_out{:},'pattern_ids_all_runs')
                
                % take correlations of the z-scored patterns
                % this should help remove the univariate effects
                z_pattern_corr = corr(z_pattern_all_runs);
                fname_z_corr_out = cellstr([b.cur_ROI_dir,filesep,roi_name,'_z_pattern_corr_no_outlier_trials_all_runs']);
                save(fname_z_corr_out{:}, 'z_pattern_corr')
                                
                % take correlation across all voxels and all trials
                % and save this out too 
                pattern_corr = corr(pattern_all_runs);
                fname_corr_out = cellstr([b.cur_ROI_dir,filesep,roi_name,'_pattern_corr_no_outlier_trials_all_runs']);
                save(fname_corr_out{:},'pattern_corr')
                
                % save out the nan-means for all trials across all runs
                fname_trial_means_out = cellstr([b.cur_ROI_dir, filesep, roi_name, '_trial_nan_means_all_runs']);
                save(fname_trial_means_out{:},'trial_means_all_runs')
            end

            
        end %iroi
        
    end %idir=
    
    tEnd_subj =  toc(subj_timer);
    fprintf('Extracting patterns for %s took %d minutes and %f seconds.\n',b.curSubj,floor(tEnd_subj/60),rem(tEnd_subj,60))
end %isub

