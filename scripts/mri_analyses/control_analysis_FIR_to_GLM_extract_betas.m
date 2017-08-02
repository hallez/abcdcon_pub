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
            
            % determine how many betas the current subject has
            all_cur_betas = dir(fullfile(b.modelDir, b.curSubj,'beta*.img'));
            num_betas = length(all_cur_betas);

            % initialize a variable the size of the current ROI
            % (rows) x number of betas of interest
            cur_betas = nan(size(find(b.cur_roi.img > 0),1), num_betas);
            % initialize a variable for the trial IDs 
            beta_ids = cell(1,num_betas);

            for ibeta=1:num_betas
                % read in the current beta
                b.cur_beta=spm_vol(fullfile(b.modelDir, b.curSubj, sprintf('beta_%04d.img', ibeta)));
                b.cur_beta.img = spm_read_vols(b.cur_beta);

                % apply current ROI mask to current beta image,
                % save results into a voxel x trial pattern matrix
                cur_betas(:,ibeta) = b.cur_beta.img(b.cur_roi.img > 0);

                % figure out trial-specific label
                pat = '^spm_spm:beta \(\d{4}\) \- Sn\(\d{1}\) (?<trial>\w+)';
                n = regexp(b.cur_beta.descrip, pat, 'names');
                pat2 = '^spm_spm:beta \(\d{4}\) \- Sn\((?<runnum>\d+)';
                n2 = regexp(b.cur_beta.descrip, pat2, 'names');

                beta_ids{:,ibeta} = sprintf('%s_beta%02d_run%s',n.trial, ibeta, n2.runnum);
            end %ibeta= 

            if size(cur_betas,1) > 0
                % take the mean across all voxels w/in a given beta
                beta_means_all_runs = nanmean(cur_betas);
            end %if size(cur_betas,1) > 0

            % save out 
            if size(betas_all_runs,1) > 0
                % save out the nan-means for all trials across all runs
                fname_beta_means_out = cellstr([b.cur_ROI_dir, filesep, roi_name, '_FIR_beta_nan_means_all_runs']);
                save(fname_beta_means_out{:},'beta_means_all_runs')
                
                fname_ids_out = cellstr([b.cur_ROI_dir,filesep,roi_name,'_FIR_beta_ids_all_runs']);
                save(fname_ids_out{:},'beta_ids')
            end

            
        end %iroi
        
    end %idir=
    
    tEnd_subj =  toc(subj_timer);
    fprintf('Extracting betas for %s took %d minutes and %f seconds.\n',b.curSubj,floor(tEnd_subj/60),rem(tEnd_subj,60))
end %isub

