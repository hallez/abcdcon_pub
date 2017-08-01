% Loop across ROIs, extract (mask) FIR betas within these ROIs
%
% Halle R. Dimsdale-Zucker

%% Setup some basics
initialize_ABCDCon

roi_dirs = {'ashs_left','ashs_right'};
b.all_ROIs = {'brCA1_body', 'brCA2_3_DG_body', 'brwhole_hippo.nii'};
FIR_order = 10; % this is a parameter that's set in `control_analysis_first_level_FIR.m`

%% Loop across subjects 
for isub=1:length(subjects)

    fprintf('\n----Working on s%s----\n',num2str(subjects(isub),'%03d'))
    subj_timer = tic;

    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.analMRIDir = analMRIDir;
    b.modelDir = [analMRIDir,'univariate_FIR_byMemory'];
    % if last subject overwrote b.runs, revert back 
    runbase = 'run';
    runnums = 1:4;
    nruns = length(runnums);
    b.runs = cell(1,nruns);
    for irun=1:nruns
        b.runs(irun) = cellstr([runbase,num2str(runnums(irun))]);
    end %for irun=
    
    % Run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);
    
    %% Loop across ROI directories
    for idir=1:length(roi_dirs)
    
        fprintf('Working on ROI dir:%s.\n', roi_dirs{idir})
        b.cur_ROI_dir = fullfile(b.analMRIDir, b.curSubj, 'ROIs', roi_dirs{idir});
        
        %% Loop across ROIs w/in the current ROI directory 
        for iroi=1:length(b.all_ROIs) 
            
            fprintf('Working on ROI %d of %d ROI: %s\n',iroi,length(b.all_ROIs),b.all_ROIs{iroi})
            b.cur_roi=spm_vol(fullfile(b.cur_ROI_dir, [b.all_ROIs{iroi} '.nii']));
            b.cur_roi.img = spm_read_vols(b.cur_roi);
            [path, roi_name, ext] = fileparts(b.cur_roi.fname);
            
            %% Loop across runs            
            for irun=1:length(b.runs)

                b.curRun = b.runs{irun};
                fprintf('Working on: %s\n',b.curRun);
                                
                % read in current regressor file
                load(fullfile(b.modelDir, b.curSubj, sprintf('byMemory_model_%s_%s_regs.mat', b.curSubj, b.curRun)))
                num_regs = length(names); % `names` is read in from the .mat file
                num_betas = 19; %num_regs * FIR_order; % there are technically also betas for motion and spike regressors, but we don't care about those for now
                
                %% Loop across betas
                curRunDir = [b.modelDir,filesep,b.curSubj,filesep,b.curRun];                
                curDirs = dir(curRunDir);

                % initialize a variable the size of the current ROI
                % (rows) x number of betas of interest
                cur_betas = nan(size(find(b.cur_roi.img > 0),1), num_betas);
                % initialize a variable for the trial IDs 
                cur_reg_ids = cell(1,num_betas);

                for ibeta=1:num_betas
                    % read in the current beta
                    b.cur_beta=spm_vol(fullfile(b.modelDir, b.curSubj, sprintf('beta_%04d.img', ibeta)));
                    b.cur_beta.img = spm_read_vols(b.cur_beta);

                    % apply current ROI mask to current beta image,
                    % save results into a voxel x trial pattern matrix
                    cur_betas(:,ibeta) = b.cur_beta.img(b.cur_roi.img > 0);

                    % figure out trial-specific label
                    % NB: could now do this from the IDs file 
                    [TOKEN, TRIAL_LBL] = strtok(cellstr(b.cur_beta.descrip),'R'); % ADAPT THIS SO WORKS FOR ALL REGRESSOR NAMES
                    [REG_NAME, REMAIN] = strtok(TRIAL_LBL,'*');
                    [TOKEN, REG_NUM] = strtok(TRIAL_LBL,'(');
                    cur_reg_ids(:,ibeta) = strcat(REG_NAME, '_', REG_NUM);

                end %ibeta=  
                keyboard
                
                if size(cur_betas,1) > 0
                    % put the current pattern (and ids) into a variable
                    % across runs
                    % FIGURE OUT HOW TO ADAPT IDEA OF `COL_COUNTER`
                    pattern_all_runs(:,col_counter + 1:col_counter + num_trials) = cur_betas;
                    pattern_ids_all_runs(1,col_counter + 1:col_counter + num_trials) = cur_reg_ids;
                    
                    % take the mean across all voxels w/in a given beta
                    trial_means = nanmean(cur_betas);
                    trial_means_all_runs(1, col_counter + 1:col_counter + num_trials) = trial_means;
                end %if size(pattern,1) > 0
                
                col_counter = col_counter + num_trials; 
                
                % clear variables that change for each run loop
                clear cur_pattern cur_reg_ids z_cur_pattern

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

