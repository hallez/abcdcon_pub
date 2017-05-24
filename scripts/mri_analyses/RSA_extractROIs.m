function [] = extractROIs()

    % Code to read in traced ROIs (manual or ASHS) and extract ROIs (based on
    % values in labels file) and save out as individual files. 
    %
    % Halle R. Dimsdale-Zucker

    initialize_ABCDCon

    %% Set flags
    % //TODO loop across all options instead of using flags 
    
    % use ASHS segmentation (could set other flags for other tracings,
    % e.g., manual, FreeSurfer, etc.)
    ASHS_FLAG = 1;
    
    % use raw ASHS segmentation (otherwise, use coreg)
    RAW_FLAG = 1; 
    
    % use 1.9 mm T2 (otherwise, use 1.5 mm)
    T2_19_FLAG = 1;
    
    % switch between L and R hemi (relevant only for ASHS)
    % 1=left hemi, 0=right hemi
    LEFT_FLAG = 0;

    %% Loop across subjects
    for isub=1:length(subjects)

        fprintf('\n----Working on s%s----\n',num2str(subjects(isub),'%03d'))

        b.curSubj = ['s' num2str(subjects(isub),'%03d')];
        
        % Run exceptions for subject-specific naming conventions
        b = run_exceptions_ABCDCon(b);

        % Setup base directory for ROI output
        b.ROI_base_dir = [analMRIDir, b.curSubj,filesep,'ROIs'];
        if ~exist(b.ROI_base_dir,'dir');
           fprintf('Making base ROIs directory.\n')
           mkdir(b.ROI_base_dir);
        end % if ~exist

        if ASHS_FLAG
            % deal w/ hemisphere of current interest
            if LEFT_FLAG 
                cur_hemi = 'left';
            else
                cur_hemi = 'right';
            end
            
            % define where split ROIs will be written out to
            b.ROI_dir = [b.ROI_base_dir,filesep,'ashs_',cur_hemi];

            % define where to-be split file lives 
            if RAW_FLAG
                b.dataDir = [analMRIDir,b.curSubj,filesep,'ashs',filesep,'raw',filesep,'final'];
            else 
                b.dataDir = [analMRIDir,b.curSubj,filesep,'ashs',filesep,'coreg',filesep,'final'];
            end

            % define name of the to-be split file 
            if T2_19_FLAG
                b.nii_file = [b.dataDir,filesep,[b.curSubj,'_19_no_coreg_',cur_hemi,'_lfseg_corr_usegray.nii']];
            else
                b.nii_file = [b.dataDir,filesep,[b.curSubj,'_15_no_coreg_',cur_hemi,'_lfseg_corr_usegray.nii']];
            end

            % read in the labels file 
            b.lbls=tdfread([analMRIDir,filesep,'snaplabels_forMatlab.txt'],',');

        else
            % define where split ROIs will be written out to
            b.ROI_dir = [b.ROI_base_dir,filesep,'manual'];

            % Define which T2 was traced and needs to be split 
            if T2_19_FLAG
                b.dataDir = [config.directories.tracings,filesep,b.curSubj,filesep,'t2_19'];
            else
                b.dataDir = [config.directories.tracings,filesep,b.curSubj,filesep,'t2_15'];
            end

            % define name of the to-be split file
            b.nii_file=[b.dataDir,filesep,b.curSubj,'_t2_traced.nii'];

            % read in the labels file
            b.lbls=tdfread([config.directories.tracings,filesep,'HighResMTLLabels_forMatlab.txt'],',');

        end %if ASHS_FLAG

        % Actually run the extraction
        extract(b);

    end %for isub=
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
function []=extract(b)
    
    if ~exist(b.ROI_dir,'dir');
       fprintf('Making current ROI directory.\n')
       mkdir(b.ROI_dir);
    end % if ~exist
    
    % Unzip the .nii.gz file
    % Should create a .nii file with the same name in the pwd
    if exist(b.nii_file,'file')==0
        fprintf('Unzipping %s\n',b.nii_file)
        gunzip([b.nii_file,'.gz']);
    end %if ~exist

    % Read in the .nii file
    b.all_ROIs=spm_vol(b.nii_file);
    b.all_ROIs.img = spm_read_vols(b.all_ROIs);

    % Get info about the labels
    % 'stable' will return them in the order they occur in the file; for
    % ITK-SNAP output, not an issue b/c ordered in descending numerical
    % order (but just a precaution)
    b.lbl_ids=unique(b.lbls.IDX,'stable');
    b.nlbls=length(b.lbl_ids);
    
    % Do a sanity check that the image contains all of the values from the ROI
    % labels
    b.lbl_idx=ismember(b.lbl_ids,unique(b.all_ROIs.img,'stable'));

    for ilbl=1:b.nlbls
        
        % skip over labels that don't exist
        % //TODO make variable names clearer 
        if sum(strcmp(cellstr(b.lbls.LABEL(ilbl,:)),cellstr(b.lbls.LABEL(~b.lbl_idx,:))))>0
            fprintf('Skipping %s. Does not exist in current file.\n',(b.lbls.LABEL(ilbl,:)))
        else
            % set the current label (name, index value)
            cur_lbl_name=b.lbls.LABEL(ilbl,:);
            cur_lbl_idx=b.lbls.IDX(ilbl);

            fprintf('Working on %s.\n',cur_lbl_name);

            %deal with whitespace in cur_lbl_name
            for i=1:length(cur_lbl_name)
                if cur_lbl_name(i)==' '
                    cur_lbl_name=cur_lbl_name(1:i-1);
                    break;
                end %if cur_lbl_name(i)==' '
            end %i=1:length(cur_lbl_name)

            % extract the voxels from b.all_ROIs that match the current label's index
            % value
            cur_ROI=b.all_ROIs.img==cur_lbl_idx;

            % save out new ROI but make sure don't overwrite if already exists
            % If need to check: 
            % imagesc(cur_ROI_nii(:,:,i) where i iterates through the slices in the
            % z dimension
            if exist([b.ROI_dir,filesep,cur_lbl_name,'.nii'],'file')==0
                fname = [b.ROI_dir,filesep,cur_lbl_name,'.nii'];
                tmp_b_all_rois = b.all_ROIs;
                tmp_b_all_rois.fname = fname;
                spm_write_vol(tmp_b_all_rois,cur_ROI);
            end %if exist
        end %sum(strcmp(
        
    end %ilbl=

end
