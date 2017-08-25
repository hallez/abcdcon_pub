initialize_ABCDCon
cur_dir = pwd;

ana_dir = analMRIDir; 
output_dir = [analMRIDir, filesep, 'univariate_sanityCheck'];

cons = {'con_0005.img' 'con_0010.img' 'con_0015.img'};
conlabels = {'allRHits_vs_FamHitsANDMiss' 'brownRHitsxFHits_Miss' 'grayRHitsxFHits_Miss'};

outLabel = 'body_ROIs';

ROIs = {'rCA1_body.nii','rCA2_3_DG_body.nii','rwhole_hippo.nii'}; 
roi_dirs = {'ashs_left','ashs_right'};

% Initialize counter
% this determines the row to write to
% row 1 is reserved for the headers
cnt = 2;

for iroi_dir = 1:length(roi_dirs);
    
    cur_roi_dir = roi_dirs{iroi_dir}; 
    fprintf('----\nCurrent directory: %s\n----', cur_roi_dir)
    
    for curROI = 1:length(ROIs);
        
        fprintf('\n\tStarting %s\n',ROIs{curROI})
        [a ROIname c] = fileparts(ROIs{curROI});
        output(1,1) = {'subject'};
        output(1,2) = {'roi_file'};
        output(1,3) = {'contrast'};
        output(1,4) = {'activity'};
        output(1,5) = {'hemi'};

        for  isub = 1:length(subjects)
            
            curSub = ['s' num2str(subjects(isub),'%03d')];
            fprintf('\tWorking on %s\n',curSub)

            % Get subject-specific ROI mask filename
            maskDir = [ana_dir, filesep, curSub, filesep, 'ROIs', filesep, cur_roi_dir]; 
            curMask = [maskDir, filesep, ROIs{curROI}];

            if exist(curMask,'file')

                % Check if mask is zipped; if so, unzip it
                unzip_flag = 0;
                if strfind(curMask,'.gz')
                   gunzip(curMask);
                   curMask = strrep(curMask,'.gz','');
                   unzip_flag = 1;
                end

                % Read in mask
                [roi_y roi_xyz] = spm_read_vols(spm_vol(curMask));

                for curCon = 1:length(cons)

                    % Get contrast values within the mask
                    cur_con_fpath = [ana_dir, filesep, 'univariate_sanityCheck',...
                        filesep, curSub, filesep, cons{curCon}];

                    if exist(cur_con_fpath,'file')

                        [con_y con_xyz] = spm_read_vols(spm_vol(cur_con_fpath));

                        % ensure dimensions are the same between the ROI
                        % and the contrast
                        if size(con_y) == size(roi_y) 
                            % find values in contrast for current ROI
                            t = con_y(roi_y>0);
                            t = t(~isnan(t));

                            % Create output matrix
                            output(cnt,1) = {curSub};
                            output(cnt,2) = {ROIname};
                            output(cnt,3) = {conlabels{curCon}};
                            output(cnt,4) = {mean(t)};
                            output(cnt,5) = {cur_roi_dir};

                            %Increment counter and reset for next subject
                            clear t con_y con_xyz
                                    cnt = cnt + 1;
                        else
                            warning('Dimensions differ for contrast %s and roi %s',...
                                cons{curCon},ROIs{curROI})
                        end %if con_y == roi_y
                        
                    else
                        warning('curCon file does not exist, skipping: %s.\n',cur_con_fpath);
                    end %exist(cur_con_fpath
                end %for curCon

                % Re-zip mask if necessary
                if unzip_flag
                    gzip(curMask);
                    delete(curMask);
                end

                %Reset for next subject
                clear roi_y roi_xyz
                
            else
                warning('curROI file does not exist, skipping: %s.\n',curMask);
            end %exist(curMask
        end %for curSub
    end %for curROI
end %iroi_dir

% Write out results
cd(output_dir)
cell2csv(strcat('cons_',outLabel,'.csv'),output); %using excel on mac

cd(cur_dir)



