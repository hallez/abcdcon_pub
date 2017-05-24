function []=RSA_reslice_t2_and_ROIs_batch

% Apply parameters needed to coregister and reslice T2 into mean EPI space
% to all ROIs.
%
% Halle R. Dimsdale-Zucker

initialize_ABCDCon
curdir=pwd; %this is a SPM-specific convention
outLabel = 'coreg_reslice_ROIs'; 

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions 

function [matlabbatch]=batch_job(b)

%% Specify files
matlabbatch{1}.cfg_basicio.file_fplist.dir = {[b.dataDir,'run1/']};
matlabbatch{1}.cfg_basicio.file_fplist.filter = '^meanf';
matlabbatch{1}.cfg_basicio.file_fplist.rec = 'FPList';
matlabbatch{2}.cfg_basicio.file_fplist.dir = {[b.dataDir,'t2_19/']};
matlabbatch{2}.cfg_basicio.file_fplist.filter = '^s';
matlabbatch{2}.cfg_basicio.file_fplist.rec = 'FPList';
matlabbatch{3}.cfg_basicio.file_fplist.dir = {[b.dataDir,'ROIs/']};
% based on http://stackoverflow.com/questions/2116328/regexp-matching-string-not-starting-with-my 
% and http://www.regextester.com/15
% skip over ROIs that start with 'r' or 'br' (ie, don't re-reslice already
% resliced OR already binarized and resliced files
matlabbatch{3}.cfg_basicio.file_fplist.filter = '^((?![br]).|(?![r]).)*$';
matlabbatch{3}.cfg_basicio.file_fplist.rec = 'FPListRec';

%% Coregister and reslice
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep;
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tname = 'Reference Image';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).sname = 'File Selector (Batch Mode): Selected Files (^meanf)';
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.coreg.estwrite.ref(1).src_output = substruct('.','files');
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1) = cfg_dep;
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tname = 'Source Image';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).sname = 'File Selector (Batch Mode): Selected Files (^s)';
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.coreg.estwrite.source(1).src_output = substruct('.','files');
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1) = cfg_dep;
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tname = 'Other Images';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).sname = 'File Selector (Batch Mode): Selected Files (^((?![br]).|(?![r]).)*$)';
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.coreg.estwrite.other(1).src_output = substruct('.','files');
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{4}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{4}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

end