function [] = control_analysis_RSA_single_trial_models_batch_FIR()

% Halle R. Dimsdale-Zucker
% Loop across subjects, runs, and trials to estimate single trial betas for
% multivariate (pattern similarity) analyses. 


%% Specify variables not already in initialize script

initialize_ABCDCon
curdir = pwd;
outLabel = 'single_trial_sanityCheck_model_FIR';

%%
% let's time to see how long this takes
all_models_timer = tic;

for isub=1:length(subjects)
    
    % let's also time the length for each subject 
    cur_subj_timer = tic;
    
    fprintf('\n----Working on s%s----',num2str(subjects(isub),'%03d'))

    % Define variables for individual subjects
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.dataDir = [analMRIDir,b.curSubj,'/'];
    b.analMRIDir = analMRIDir;
    b.modelDir = [analMRIDir,'multivariate_FIR_sanityCheck'];

    % Run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);
    
    % Loop across runs
    for irun=1:length(b.runs)
        b.curRun = b.runs{irun};
        
        % Loop across trials
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
                
        for itrial=1:trialCounter
            
            b.curTrial = ['trial_',num2str(numCurTrials(itrial),'%03d')];
            
            % specify matlabbatch variable with subject-specific inputs
            matlabbatch = batch_job(b);

            % save matlabbatch variable for posterity
            outName = strcat(b.dataDir,outLabel,'_',date);
            save(outName, 'matlabbatch');

            % run matlabbatch job
            cd(b.dataDir);
            try
                spm_jobman('initcfg')
                spm('defaults', 'FMRI');
                spm_jobman('serial', matlabbatch);
            catch
                cd(curdir);
                continue;
            end %try
            
        end %itrial=  
        
    end %irun=
    
    tEnd_subj =  toc(cur_subj_timer);
    fprintf('Running the models for %s took %d minutes and %f seconds.\n',b.curSubj,floor(tEnd_subj/60),rem(tEnd_subj,60))
    
    cd(curdir);
end %isub

tEnd_all_models = toc(all_models_timer);
fprintf('\nRunning all those models took %d minutes and %f seconds. Phew, glad that''s over!\n',floor(tEnd_all_models/60),rem(tEnd_all_models,60))
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
function [matlabbatch]=batch_job(b)
%This function generates the matlabbatch variable: can be copied in
%directly from the batch job output, then modify lines as necessary to
%generalize the paths, etc, using b variables

    %% Specfiy files
    fprintf('Figuring out files for subject %s.\n',b.curSubj)
    matlabbatch{1}.cfg_basicio.file_fplist.dir = {[b.dataDir,b.curRun,filesep]};
    matlabbatch{1}.cfg_basicio.file_fplist.filter = '^rf';
    matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{2}.cfg_basicio.file_fplist.dir = {[b.modelDir,filesep,b.curSubj,filesep,b.curRun,filesep,b.curTrial,filesep]};
    matlabbatch{2}.cfg_basicio.file_fplist.filter = 'regs.mat';
    matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{3}.cfg_basicio.file_fplist.dir = {[b.dataDir,b.curRun,filesep]};
    matlabbatch{3}.cfg_basicio.file_fplist.filter = 'spike_regs_rp.txt';
    matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPList';
    
    %% Run model
    matlabbatch{4}.spm.stats.fmri_spec.dir = {[b.modelDir,filesep,b.curSubj,filesep,b.curRun,filesep,b.curTrial,filesep]};
    matlabbatch{4}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{4}.spm.stats.fmri_spec.timing.RT = 2.01;
    matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep;
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tname = 'Scans';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).sname = 'File Selector (Batch Mode): Selected Files (^rf)';
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1).src_output = substruct('.','files');
    matlabbatch{4}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1) = cfg_dep;
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tname = 'Multiple conditions';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).sname = 'File Selector (Batch Mode): Selected Files (regs.mat)';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi(1).src_output = substruct('.','files');
    matlabbatch{4}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1) = cfg_dep;
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tname = 'Multiple regressors';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).sname = 'File Selector (Batch Mode): Selected Files (spike_regs_rp.txt)';
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg(1).src_output = substruct('.','files');
    matlabbatch{4}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{4}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{4}.spm.stats.fmri_spec.bases.fir.length = 20;
    matlabbatch{4}.spm.stats.fmri_spec.bases.fir.order = 10;
    matlabbatch{4}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{4}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{4}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{4}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% Estimate model
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
    matlabbatch{5}.spm.stats.fmri_est.method.Classical = 1;
    
end %batch_job
