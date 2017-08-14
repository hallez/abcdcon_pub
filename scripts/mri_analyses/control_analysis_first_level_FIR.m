function [] = control_analysis_first_level_FIR()

% Halle R. Dimsdale-Zucker


%% Specify variables not already in initialize script

initialize_ABCDCon
curdir = pwd;

%%%MODEL DEFINITIONS 
allModels(1).name = 'byMemory_model';
allModels(1).type = 'univariate_FIR';

allModels(2).name = 'byMemoryFHitMiss_model';
allModels(2).type = 'univariate_FIR';

% enable selecting to analyze different models
modelSelect = 0;
for imodel=1:size(allModels,2)
    curModel = allModels(imodel).name;
    b.curModel = curModel; % this is needed to feed into SPM batch
    curModelBase = strtok(curModel,'_');
    curModelType = allModels(imodel).type;
    outLabel = sprintf('firstlevel_%s', curModelBase);
    modelSelect = input(['Do you want to analyze ' curModel '? (Y=1,N=0): ']);
    if modelSelect
        break;
    end %if modelSelect
end %imodel=
fprintf('Selected to analyze %s\n',curModel)

outputDir = [analMRIDir curModelType '_' curModelBase filesep];
if ~isdir(outputDir)
    fprintf('Creating model directory: %s\n',curModel)
    mkdir(outputDir)
end %if isdir(

%%
for isub=1:length(subjects)

    fprintf('\n----Working on s%s----',num2str(subjects(isub),'%03d'))

    % Define variables for individual subjects
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.dataDir = [analMRIDir,b.curSubj,'/'];
    b.analMRIDir = analMRIDir;
    b.analysisName = sprintf('%s_model_', curModelBase);
    b.modelDir = [analMRIDir, sprintf('%s_%s', curModelType, curModelBase)];

    % Run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);

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

    cd(curdir);

end %isub

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
function [matlabbatch]=batch_job(b)
%This function generates the matlabbatch variable: can be copied in
%directly from the batch job output, then modify lines as necessary to
%generalize the paths, etc, using b variables

    %% Specfiy files
    fprintf('Figuring out files for subject %s.\n',b.curSubj)
    matlabbatch{1}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run1/']};
    matlabbatch{1}.cfg_basicio.file_fplist.filter = '^rf';
    matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{2}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run2/']};
    matlabbatch{2}.cfg_basicio.file_fplist.filter = '^rf';
    matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{3}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run3/']};
    matlabbatch{3}.cfg_basicio.file_fplist.filter = '^rf';
    matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{4}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run4/']};
    matlabbatch{4}.cfg_basicio.file_fplist.filter = '^rf';
    matlabbatch{4}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{5}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run1/']};
    matlabbatch{5}.cfg_basicio.file_fplist.filter = 'spike_regs_rp.txt';
    matlabbatch{5}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{6}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run2/']};
    matlabbatch{6}.cfg_basicio.file_fplist.filter = 'spike_regs_rp.txt';
    matlabbatch{6}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{7}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run3/']};
    matlabbatch{7}.cfg_basicio.file_fplist.filter = 'spike_regs_rp.txt';
    matlabbatch{7}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{8}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run4/']};
    matlabbatch{8}.cfg_basicio.file_fplist.filter = 'spike_regs_rp.txt';
    matlabbatch{8}.cfg_basicio.file_fplist.rec = 'FPListRec';

    %% First level modelling
    fprintf('Starting first level for subject %s.\n',b.curSubj)
    matlabbatch{9}.spm.stats.fmri_spec.dir = {[b.modelDir,'/',b.curSubj]};
    matlabbatch{9}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{9}.spm.stats.fmri_spec.timing.RT = 2.01;
    matlabbatch{9}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{9}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).tname = 'Scans';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).sname = 'File Selector (Batch Mode): Selected Files (^rf)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).scans(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi = {fullfile(b.modelDir,b.curSubj, sprintf('%s_%s_run1_regs.mat',b.curModel, b.curSubj))};
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).tname = 'Multiple regressors';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).sname = 'File Selector (Batch Mode): Selected Files (spike_regs_rp.txt)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi_reg(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).hpf = 128;
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).tname = 'Scans';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).sname = 'File Selector (Batch Mode): Selected Files (^rf)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).scans(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi = {fullfile(b.modelDir,b.curSubj, sprintf('%s_%s_run2_regs.mat',b.curModel, b.curSubj))};
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).tname = 'Multiple regressors';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).sname = 'File Selector (Batch Mode): Selected Files (spike_regs_rp.txt)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi_reg(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).hpf = 128;
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).tname = 'Scans';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).sname = 'File Selector (Batch Mode): Selected Files (^rf)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).scans(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi = {fullfile(b.modelDir,b.curSubj, sprintf('%s_%s_run3_regs.mat',b.curModel, b.curSubj))};
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).tname = 'Multiple regressors';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).sname = 'File Selector (Batch Mode): Selected Files (spike_regs_rp.txt)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi_reg(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).hpf = 128;
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).tname = 'Scans';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).sname = 'File Selector (Batch Mode): Selected Files (^rf)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).scans(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi = {fullfile(b.modelDir,b.curSubj, sprintf('%s_%s_run4_regs.mat',b.curModel, b.curSubj))};
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).regress = struct('name', {}, 'val', {});
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1) = cfg_dep;
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).tname = 'Multiple regressors';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).sname = 'File Selector (Batch Mode): Selected Files (spike_regs_rp.txt)';
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1});
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi_reg(1).src_output = substruct('.','files');
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).hpf = 128;
    matlabbatch{9}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{9}.spm.stats.fmri_spec.bases.fir.length = 20;
    matlabbatch{9}.spm.stats.fmri_spec.bases.fir.order = 10;
    matlabbatch{9}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{9}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{9}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{9}.spm.stats.fmri_spec.cvi = 'AR(1)';

    %% Estimate
    fprintf('Starting estimation for subject %s.\n',b.curSubj)
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{10}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
    matlabbatch{10}.spm.stats.fmri_est.method.Classical = 1;
end %function [matlabbatch]