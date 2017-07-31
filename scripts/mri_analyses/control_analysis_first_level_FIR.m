function [] = control_analysis_first_level_FIR()

% Halle R. Dimsdale-Zucker


%% Specify variables not already in initialize script

initialize_ABCDCon
curdir = pwd;
outLabel = 'firstlevel_byMemory';
b.replace_cons = 0;

%%
for isub=1:length(subjects)

    fprintf('\n----Working on s%s----',num2str(subjects(isub),'%03d'))

    % Define variables for individual subjects
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.dataDir = [analMRIDir,b.curSubj,'/'];
    b.analMRIDir = analMRIDir;
    b.analysisName = 'byMemory_model_';
    b.modelDir = [analMRIDir,'univariate_FIR_byMemory'];
    b.cons = '_cons.mat';

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
    matlabbatch{9}.spm.stats.fmri_spec.sess(1).multi = {[b.modelDir,'/',b.curSubj,'/byMemory_model_',b.curSubj,'_run1_regs.mat']};
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
    matlabbatch{9}.spm.stats.fmri_spec.sess(2).multi = {[b.modelDir,'/',b.curSubj,'/byMemory_model_',b.curSubj,'_run2_regs.mat']};
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
    matlabbatch{9}.spm.stats.fmri_spec.sess(3).multi = {[b.modelDir,'/',b.curSubj,'/byMemory_model_',b.curSubj,'_run3_regs.mat']};
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
    matlabbatch{9}.spm.stats.fmri_spec.sess(4).multi = {[b.modelDir,'/',b.curSubj,'/byMemory_model_',b.curSubj,'_run4_regs.mat']};
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


%%%%CONTRAST FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = batch_contrasts(b)
%This function evaluates a set of previously-defined t-contrasts. Adapted
%from fmri_contrasts_b (Ken Roberts, Duke Univ, CCN)

%load SPM
load(strcat(b.modelDir,filesep,b.curSubj,filesep,'SPM.mat'));

%create xCon structures for each con
% what is the range of this batch of contrasts?
if ~isempty(SPM.xCon) && b.replace_cons~=1
    num_existing_cons = length(SPM.xCon); %can modify this to overwrite
else
    num_existing_cons = 0;
    SPM.xCon = struct( 'name',{{'init'}}, 'STAT', [1], 'c', [1], 'X0', [1], ...
        'iX0', {{'init'}}, 'X1o', [1], 'eidf', [], 'Vcon', [],  'Vspm', [] );
end;

%load in contrast vectors and names
load(strcat(b.modelDir,filesep,b.curSubj,filesep,b.analysisName,b.curSubj,b.cons));
num_cons_tocalc = length(contrast_vectors);
sessionmeans = zeros(1,length(b.runs)); % this right-pads to account for number of runs for run means to ensure that contrast vector is same length as design matrix. 

%update xCon structures
for iCon = 1:num_cons_tocalc
    SPM.xCon(iCon + num_existing_cons) = spm_FcUtil('Set', contrastNames{iCon}, ...
        'T', 'c', transpose([contrast_vectors{iCon} sessionmeans]), SPM.xX.xKXs);
end;

%evaluate contrasts
Ci = num_existing_cons+1:(num_existing_cons + num_cons_tocalc);
spm_contrasts(SPM, Ci);

end
