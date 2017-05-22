function [] = preprocess()

% Runs basic preprocessing steps 
%
% Halle R. Dimsdale-Zucker

%% Specify variables not already in initialize script

initialize_ABCDCon
curdir = pwd;
outLabel = 'preproc_moco_quickcoreg';

for isub=1:length(subjects)

    fprintf('\n----Working on s%s----',num2str(subjects(isub),'%03d'))

    % Define variables for individual subjects
    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.dataDir = [analMRIDir,b.curSubj,'/'];

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
    matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{2}.cfg_basicio.file_fplist.dir = {[b.dataDir, 't2_19/']};
    matlabbatch{2}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{3}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run1/']};
    matlabbatch{3}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{4}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run2/']};
    matlabbatch{4}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{4}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{5}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run3/']};
    matlabbatch{5}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{5}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{6}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run4/']};
    matlabbatch{6}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{6}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{7}.cfg_basicio.file_fplist.dir = {[b.dataDir,'room_localizer/']};
    matlabbatch{7}.cfg_basicio.file_fplist.filter = '^f';
    matlabbatch{7}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{8}.cfg_basicio.file_fplist.dir = {[b.dataDir,'fieldmap1/']};
    matlabbatch{8}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{8}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{9}.cfg_basicio.file_fplist.dir = {[b.dataDir,'fieldmap2/']};
    matlabbatch{9}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{9}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{10}.cfg_basicio.file_fplist.dir = {[b.dataDir,'ep_seg_partial/']};
    matlabbatch{10}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{10}.cfg_basicio.file_fplist.rec = 'FPList';
    matlabbatch{11}.cfg_basicio.file_fplist.dir = {[b.dataDir,'ep_seg_wholebrain/']};
    matlabbatch{11}.cfg_basicio.file_fplist.filter = '^s';
    matlabbatch{11}.cfg_basicio.file_fplist.rec = 'FPList';

    %% Motion correction (realignment + reslice)
    fprintf('Doing motion correction for subject %s.\n',b.curSubj)
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep;
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).tname = 'Session';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
    matlabbatch{12}.spm.spatial.realign.estwrite.data{1}(1).src_output = substruct('.','files');
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep;
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).tname = 'Session';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1});
    matlabbatch{12}.spm.spatial.realign.estwrite.data{2}(1).src_output = substruct('.','files');
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1) = cfg_dep;
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).tname = 'Session';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1});
    matlabbatch{12}.spm.spatial.realign.estwrite.data{3}(1).src_output = substruct('.','files');
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1) = cfg_dep;
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).tname = 'Session';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1});
    matlabbatch{12}.spm.spatial.realign.estwrite.data{4}(1).src_output = substruct('.','files');
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1) = cfg_dep;
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).tname = 'Session';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).sname = 'File Selector (Batch Mode): Selected Files (^f)';
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1});
    matlabbatch{12}.spm.spatial.realign.estwrite.data{5}(1).src_output = substruct('.','files');
    matlabbatch{12}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{12}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{12}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{12}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{12}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{12}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{12}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{12}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{12}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{12}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{12}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{12}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

    %% Coreg
    fprintf('Doing coreg with MPRAGE to allow for displaying results for subject %s.\n',b.curSubj)
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep;
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).tname = 'Reference Image';
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).sname = 'Realign: Estimate & Reslice: Mean Image';
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).src_exbranch = substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{13}.spm.spatial.coreg.estwrite.ref(1).src_output = substruct('.','rmean');
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1) = cfg_dep;
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).tname = 'Source Image';
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).sname = 'File Selector (Batch Mode): Selected Files (^s)';
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{13}.spm.spatial.coreg.estwrite.source(1).src_output = substruct('.','files');
    matlabbatch{13}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{13}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{13}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{13}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{13}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{13}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{13}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{13}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{13}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

    fprintf('Doing coreg with T2 to allow for displaying results for subject %s.\n',b.curSubj)
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep;
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).tname = 'Reference Image';
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).name = 'filter';
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).value = 'image';
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).sname = 'Realign: Estimate & Reslice: Mean Image';
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).src_exbranch = substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{14}.spm.spatial.coreg.estwrite.ref(1).src_output = substruct('.','rmean');
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1) = cfg_dep;
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).tname = 'Source Image';
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).sname = 'File Selector (Batch Mode): Selected Files (^s)';
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
    matlabbatch{14}.spm.spatial.coreg.estwrite.source(1).src_output = substruct('.','files');
    matlabbatch{14}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{14}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{14}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{14}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{14}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{14}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{14}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{14}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{14}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

end
