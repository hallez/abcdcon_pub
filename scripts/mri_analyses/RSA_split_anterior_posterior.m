% Halle R. Dimsdale-Zucker
initialize_ABCDCon

% set up directories, variables, etc. that are shared for all subjects
b.analyzed_mri_dir = analMRIDir;
b.raw_behav_dir = rawBehavDir;
b.rois_to_split = {'CA1', 'CA3', 'CA2_3_DG', 'CA3_DG', 'subiculum', 'whole_hippo'};
b.roi_dirs = {'ashs_left', 'ashs_right'};

for isub = 1:length(subjects)
    b.curSubj = sprintf('s%03d',subjects(isub));
    fprintf('\n----Working on %s----\n',b.curSubj)

    % set up variables that may get overwritten between
    % subejcts when check for exceptions
    b.beta_matrix_size = [144 152 52];
    
    % check for exceptions
    b = run_exceptions_ABCDCon(b);
    
    % set up subject-specific paths
    
    % read in transitions file
    transitions_fname = fullfile(b.raw_behav_dir, filesep, b.curSubj, sprintf('%s_hc_transitions.yml',b.curSubj));
    
    if exist(transitions_fname, 'file')
        transitions = ReadYaml(transitions_fname);
        
        b.left_head = transitions.left.head_slice;
        b.right_head = transitions.right.head_slice;
        b.left_body_tail = transitions.left.body_tail_transition;
        b.right_body_tail = transitions.right.body_tail_transition;
    else 
       warning('No transitions file for %s. Moving onto next subject', b.curSubj) 
       break;
    end
    
    % loop across rois
    for idir = 1:length(b.roi_dirs)
        b.cur_roi_dir = b.roi_dirs{idir};
        % this is a bit of a hack since `strtok` won't look for a delimiter
        % **string**
        b.cur_hemi = b.cur_roi_dir(strfind(b.cur_roi_dir,'_')+1:end);

        for iroi = 1:length(b.rois_to_split)
            b.cur_roi = b.rois_to_split{iroi};
            b.cur_roi_fpath = fullfile(b.analyzed_mri_dir, filesep, b.curSubj, filesep, 'ROIs', b.cur_roi_dir);
            b.cur_roi_fname = fullfile(b.cur_roi_fpath, filesep, sprintf('%s.nii', b.cur_roi));

            if exist(b.cur_roi_fname, 'file')
                % read in existing ROI
                Q = spm_vol(b.cur_roi_fname);
                QV = spm_read_vols(Q);

                if strcmp(b.cur_hemi, 'left')
                    % split into head, body, body + tail, and tail sections
                    QV_head = zeros(size(QV));
                    QV_head(:,:,1:b.left_head) = QV(:,:,1:b.left_head);

                    QV_body = zeros(size(QV));
                    QV_body(:,:,b.left_head+1:b.left_body_tail) = QV(:,:,b.left_head+1:b.left_body_tail);

                    QV_tail = zeros(size(QV));
                    QV_tail(:,:,b.left_body_tail+1:end) = QV(:,:,b.left_body_tail+1:end);
                    
                    QV_body_tail = zeros(size(QV));
                    QV_body_tail(:,:,b.left_head+1:end) = QV(:,:,b.left_head+1:end);
                elseif strcmp(b.cur_hemi, 'right')
                    QV_head = zeros(size(QV));
                    QV_head(:,:,1:b.right_head) = QV(:,:,1:b.right_head);

                    QV_body = zeros(size(QV));
                    QV_body(:,:,b.right_head+1:b.right_body_tail) = QV(:,:,b.right_head+1:b.right_body_tail);

                    QV_tail = zeros(size(QV));
                    QV_tail(:,:,b.right_body_tail+1:end) = QV(:,:,b.right_body_tail+1:end);
                    
                    QV_body_tail = zeros(size(QV));
                    QV_body_tail(:,:,b.right_head+1:end) = QV(:,:,b.right_head+1:end);
                end %if strcmp(b.cur_hemi

                % write out new ROIs
                % use 'Q' to get the structure in the way
                % spm_write_vol likes
                Q_head = Q;
                Q_body = Q;
                Q_tail = Q;
                Q_body_tail = Q;

                Q_head.fname = fullfile(b.cur_roi_fpath, filesep, sprintf('%s_head.nii', b.cur_roi));
                Q_body.fname = fullfile(b.cur_roi_fpath, filesep, sprintf('%s_body.nii', b.cur_roi));
                Q_tail.fname = fullfile(b.cur_roi_fpath, filesep, sprintf('%s_tail.nii', b.cur_roi));
                Q_body_tail.fname = fullfile(b.cur_roi_fpath, filesep, sprintf('%s_body_tail.nii', b.cur_roi));

                spm_write_vol(Q_head, QV_head);
                spm_write_vol(Q_body, QV_body);
                spm_write_vol(Q_tail, QV_tail);
                spm_write_vol(Q_body_tail, QV_body_tail);

                % clean up before moving onto next ROI
                clear Q QV Q_head QV_head Q_body QV_body Q_tail QV_tail Q_body_tail QV_body_tail
            end % if exist(b.cur_roi_fname
        end %iroi

    end % idir


end %isub

