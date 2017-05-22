function [] = new_segment_write_deformations_batch()

% Halle R. Dimsdale-Zucker

%% Specify variables not already in initialize script

initialize_ABCDCon
curdir = pwd;
outLabel = 'new_segment'; 

for isub=1:length(subjects)
    
    fprintf('\n----Working on s%s----',num2str(subjects(isub),'%03d'))
    
    % Define variables for individual subjects
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.dataDir = [analMRIDir,b.curSubj,'/'];
    b.analMRIDir = analMRIDir;
    b.analysisName = 'sanityCheck_model_';
    b.modelDir = [analMRIDir,'univariate_sanityCheck'];
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
    matlabbatch{1}.cfg_basicio.file_fplist.dir = {[b.dataDir,'mprage/']};
    matlabbatch{1}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{2}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run1/']};
    matlabbatch{2}.cfg_basicio.file_fplist.filter = '^f.*';
    matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{3}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run1/']};
    matlabbatch{3}.cfg_basicio.file_fplist.filter = '^meanf.*';
    matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{4}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run2/']};
    matlabbatch{4}.cfg_basicio.file_fplist.filter = '^f.*';
    matlabbatch{4}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{5}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run2/']};
    matlabbatch{5}.cfg_basicio.file_fplist.filter = '^meanf.*';
    matlabbatch{5}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{6}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run3/']};
    matlabbatch{6}.cfg_basicio.file_fplist.filter = '^f.*';
    matlabbatch{6}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{7}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run3/']};
    matlabbatch{7}.cfg_basicio.file_fplist.filter = '^meanf.*';
    matlabbatch{7}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{8}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run4/']};
    matlabbatch{8}.cfg_basicio.file_fplist.filter = '^f.*';
    matlabbatch{8}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{9}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run4/']};
    matlabbatch{9}.cfg_basicio.file_fplist.filter = '^meanf.*';
    matlabbatch{9}.cfg_basicio.file_fplist.rec = 'FPListRec';
    matlabbatch{10}.cfg_basicio.file_fplist.dir = {[b.modelDir,filesep,b.curSubj,filesep]};
    matlabbatch{10}.cfg_basicio.file_fplist.filter = '^con.*';
    matlabbatch{10}.cfg_basicio.file_fplist.rec = 'FPListRec';
    
    %% New segment
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1) = cfg_dep;
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).tname = 'Volumes';
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).sname = 'File Selector (Batch Mode): Selected Files (^s)';
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{11}.spm.tools.preproc8.channel.vols(1).src_output = substruct('.','files');
    matlabbatch{11}.spm.tools.preproc8.channel.biasreg = 0.0001;
    matlabbatch{11}.spm.tools.preproc8.channel.biasfwhm = 60;
    matlabbatch{11}.spm.tools.preproc8.channel.write = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(1).tpm = {'/Applications/spm8/toolbox/Seg/TPM.nii,1'};
    matlabbatch{11}.spm.tools.preproc8.tissue(1).ngaus = 2;
    matlabbatch{11}.spm.tools.preproc8.tissue(1).native = [1 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(1).warped = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(2).tpm = {'/Applications/spm8/toolbox/Seg/TPM.nii,2'};
    matlabbatch{11}.spm.tools.preproc8.tissue(2).ngaus = 2;
    matlabbatch{11}.spm.tools.preproc8.tissue(2).native = [1 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(2).warped = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(3).tpm = {'/Applications/spm8/toolbox/Seg/TPM.nii,3'};
    matlabbatch{11}.spm.tools.preproc8.tissue(3).ngaus = 2;
    matlabbatch{11}.spm.tools.preproc8.tissue(3).native = [1 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(3).warped = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(4).tpm = {'/Applications/spm8/toolbox/Seg/TPM.nii,4'};
    matlabbatch{11}.spm.tools.preproc8.tissue(4).ngaus = 3;
    matlabbatch{11}.spm.tools.preproc8.tissue(4).native = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(4).warped = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(5).tpm = {'/Applications/spm8/toolbox/Seg/TPM.nii,5'};
    matlabbatch{11}.spm.tools.preproc8.tissue(5).ngaus = 4;
    matlabbatch{11}.spm.tools.preproc8.tissue(5).native = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(5).warped = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(6).tpm = {'/Applications/spm8/toolbox/Seg/TPM.nii,6'};
    matlabbatch{11}.spm.tools.preproc8.tissue(6).ngaus = 2;
    matlabbatch{11}.spm.tools.preproc8.tissue(6).native = [0 0];
    matlabbatch{11}.spm.tools.preproc8.tissue(6).warped = [0 0];
    matlabbatch{11}.spm.tools.preproc8.warp.mrf = 0;
    matlabbatch{11}.spm.tools.preproc8.warp.reg = 4;
    matlabbatch{11}.spm.tools.preproc8.warp.affreg = 'mni';
    matlabbatch{11}.spm.tools.preproc8.warp.samp = 3;
    matlabbatch{11}.spm.tools.preproc8.warp.write = [0 1];
    
    %% Deformations
    matlabbatch{12}.spm.util.defs.comp{1}.def(1) = cfg_dep;
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).tname = 'Deformation Field';
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).sname = 'New Segment: Forward Deformations';
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).src_exbranch = substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.comp{1}.def(1).src_output = substruct('.','fordef', '()',{':'});
    matlabbatch{12}.spm.util.defs.ofname = '';
    matlabbatch{12}.spm.util.defs.fnames(1) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(1).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(1).sname = 'File Selector (Batch Mode): Selected Files (^s*.nii)';
    matlabbatch{12}.spm.util.defs.fnames(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(1).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(2) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(2).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(2).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(2).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(2).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(2).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(2).sname = 'File Selector (Batch Mode): Selected Files (^f.*)';
    matlabbatch{12}.spm.util.defs.fnames(2).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(2).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(3) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(3).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(3).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(3).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(3).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(3).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(3).sname = 'File Selector (Batch Mode): Selected Files (^meanf.*)';
    matlabbatch{12}.spm.util.defs.fnames(3).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(3).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(4) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(4).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(4).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(4).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(4).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(4).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(4).sname = 'File Selector (Batch Mode): Selected Files (^f.*)';
    matlabbatch{12}.spm.util.defs.fnames(4).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(4).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(5) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(5).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(5).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(5).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(5).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(5).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(5).sname = 'File Selector (Batch Mode): Selected Files (^meanf.*)';
    matlabbatch{12}.spm.util.defs.fnames(5).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(5).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(6) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(6).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(6).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(6).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(6).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(6).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(6).sname = 'File Selector (Batch Mode): Selected Files (^f.*)';
    matlabbatch{12}.spm.util.defs.fnames(6).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(6).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(7) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(7).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(7).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(7).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(7).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(7).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(7).sname = 'File Selector (Batch Mode): Selected Files (^meanf.*)';
    matlabbatch{12}.spm.util.defs.fnames(7).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(7).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(8) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(8).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(8).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(8).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(8).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(8).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(8).sname = 'File Selector (Batch Mode): Selected Files (^f.*)';
    matlabbatch{12}.spm.util.defs.fnames(8).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(8).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(9) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(9).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(9).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(9).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(9).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(9).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(9).sname = 'File Selector (Batch Mode): Selected Files (^meanf.*)';
    matlabbatch{12}.spm.util.defs.fnames(9).src_exbranch = substruct('.','val', '{}',{9}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(9).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(10) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(10).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(10).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.util.defs.fnames(10).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.util.defs.fnames(10).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(10).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(10).sname = 'File Selector (Batch Mode): Selected Files (^con.*)';
    matlabbatch{12}.spm.util.defs.fnames(10).src_exbranch = substruct('.','val', '{}',{10}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(10).src_output = substruct('.','files');
    matlabbatch{12}.spm.util.defs.fnames(11) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(11).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(11).tgt_spec{1}(1).name = 'filter';
    matlabbatch{12}.spm.util.defs.fnames(11).tgt_spec{1}(1).value = 'image';
    matlabbatch{12}.spm.util.defs.fnames(11).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(11).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(11).sname = 'New Segment: c1 Images';
    matlabbatch{12}.spm.util.defs.fnames(11).src_exbranch = substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(11).src_output = substruct('.','tiss', '()',{1}, '.','c', '()',{':'});
    matlabbatch{12}.spm.util.defs.fnames(12) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(12).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(12).tgt_spec{1}(1).name = 'filter';
    matlabbatch{12}.spm.util.defs.fnames(12).tgt_spec{1}(1).value = 'image';
    matlabbatch{12}.spm.util.defs.fnames(12).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(12).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(12).sname = 'New Segment: c2 Images';
    matlabbatch{12}.spm.util.defs.fnames(12).src_exbranch = substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(12).src_output = substruct('.','tiss', '()',{2}, '.','c', '()',{':'});
    matlabbatch{12}.spm.util.defs.fnames(13) = cfg_dep;
    matlabbatch{12}.spm.util.defs.fnames(13).tname = 'Apply to';
    matlabbatch{12}.spm.util.defs.fnames(13).tgt_spec{1}(1).name = 'filter';
    matlabbatch{12}.spm.util.defs.fnames(13).tgt_spec{1}(1).value = 'image';
    matlabbatch{12}.spm.util.defs.fnames(13).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.util.defs.fnames(13).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.util.defs.fnames(13).sname = 'New Segment: c3 Images';
    matlabbatch{12}.spm.util.defs.fnames(13).src_exbranch = substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{12}.spm.util.defs.fnames(13).src_output = substruct('.','tiss', '()',{3}, '.','c', '()',{':'});
    matlabbatch{12}.spm.util.defs.savedir.savesrc = 1;
    matlabbatch{12}.spm.util.defs.interp = 1;

    
end