% This script generates regressors for running single trial event models.
%
% Halle R. Dimsdale-Zucker

initialize_ABCDCon

%% MODEL INFO
curModel = 'sanityCheck_model';
curModelBase = strtok(curModel,'_');
curModelType = 'multivariate';

outputDir = [analMRIDir curModelType '_' curModelBase filesep];
if ~isdir(outputDir)
    fprintf('Creating model directory: %s\n',curModel)
    mkdir(outputDir)
end %if isdir(

%% OPTIONS
write_reg = 1; %do you want to write out regressor mats?
include_dur = 0; %do you want to include durations? 

%% END OF MODIFICATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(strcat(outputDir,'model_info'));

for isub = 1:length(subjects)
    fprintf('\n----Working on subject s%s----\n',num2str(subjects(isub),'%03d'))
    
    b.curSubj = ['s',num2str(subjects(isub),'%03d')]; 

    if ~isdir([outputDir,b.curSubj])
        mkdir(strcat(outputDir,b.curSubj));
    end
    
    % read in onset file information 
    datafile = [rawBehavDir,b.curSubj,filesep,'onsets_file.csv'];
    conditions = tdfread(datafile,','); 
    
    % check that have information for current model
    if ~isfield(conditions,curModel) 
        error('There is no column for the current model (%s) in onsets file (%s).\n',curModel,datafile);
    end %if 
    excludereg = 99;
    
    % replace values that may change between subjects
    b.runs = {'run1','run2','run3','run4' };
    
    b = run_exceptions_ABCDCon(b);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%Generate regressor files
    for irun = 1:length(b.runs)
        
        cur_run = b.runs{irun};
        cur_run_num = str2double(strtok(cur_run, 'run'));
        
        % create run-level directory
        runDir = [outputDir,filesep,b.curSubj,filesep,cur_run];
        if ~isdir(runDir)
                fprintf('Creating %s directory.\n',cur_run)
                mkdir(runDir)
        end %if isdir(
            
        % generate variable to represent where current functional run's
        % data is stored
        curfuncrun = [analMRIDir,b.curSubj,filesep,cur_run,filesep];
        
        % create a logical index for the current run (using
        % conditions.BlockID)
        curRunIdx = conditions.BlockID==cur_run_num;
        
        % subset information w/in conditions struct for the current run
        % using this logical index
        % NB: a more standard Matlab way to do this would be to have a
        % matrix where columns represent each of these variables (would
        % then also need a key to know which columns represent which
        % variables)
        curConditions.trialType = cellstr(conditions.trialTypeKey(curRunIdx,:));
        curConditions.Onset = conditions.Onset(curRunIdx);
        curConditions.ObjectID = conditions.ObjectID(curRunIdx);
        curConditions.VideoID = cellstr(conditions.VideoID(curRunIdx,:));
        curConditions.objRec_trialNum = conditions.objRec_trialNum(curRunIdx);
        curConditions.sanityCheck_model = conditions.sanityCheck_model(curRunIdx);
        curConditions.sanityCheck_key = cellstr(conditions.sanityCheck_key(curRunIdx,:));
        curConditions.duration = conditions.duration(curRunIdx);
        numCurTrials = length(unique(curConditions.objRec_trialNum));
      
        for itrial = 1:numCurTrials
            
            % create an index for the current trial
            % remember, rows in onset_file.csv are NOT ordered
            % by trial number which is why need an index
            curTrialIdx = curConditions.objRec_trialNum==itrial;
            
            % check to make sure the current trial isn't supposed to be
            % excluded
            if curConditions.sanityCheck_model(curTrialIdx) <= excludereg
               
                % create trial-level directory
                trialDir = [runDir,filesep,'trial_',num2str(itrial,'%03d')];
                if ~isdir(trialDir)
                    mkdir(trialDir)
                end %if isdir(

                filename = [trialDir,filesep,'regs'];

                % figure out regressors that occur w/in the current run
                onsets = {}; names = {}; durations = {};

                % figure out durations (if using)
                if include_dur 
                    add_dur_cur_trial = curConditions.duration(curTrialIdx);
                    add_dur_all_other_trials = curConditions.duration(~curTrialIdx)';
                else
                    add_dur_cur_trial = 0;
                    add_dur_all_other_trials = zeros(size(curConditions.duration(~curTrialIdx)))';
                end
                
                % figure out onsets
                % set the first value to be a cell and then will group
                % values for each regressor together (this is a hacky solution)
                onsets = [num2cell(curConditions.Onset(curTrialIdx))];
                onsets = [onsets curConditions.Onset(~curTrialIdx)'];

                % figure out durations
                % use the same trick as w/ onsets to get information about
                % regressors to stay together
                durations = [num2cell(add_dur_cur_trial)];
                durations = [durations add_dur_all_other_trials]; 

                % append name of current regressor 
                cur_trial_name = ['Run_',num2str(cur_run_num,'%02d'),...
                                  '_Trial_',num2str(itrial,'%03d'),'_',curConditions.trialType{curTrialIdx},...
                                  '_EncVideo_',curConditions.VideoID{curTrialIdx}];
                names = [cellstr(cur_trial_name) cellstr(['NOT_' cur_trial_name])]; 

                % save out the regressor file
                if write_reg 
                    if exist(filename,'file')
                        regsOverwrite = input('Regs file already exists. Overwrite? (Y=1,N=0):');
                        if regsOverwrite
                            save(filename,'names','onsets','durations');
                        else
                            save([filename '+'],'names','onsets','durations'); 
                        end %if regsOverwrite
                    else
                        save(filename,'names','onsets','durations');
                    end %if exist(filename
                end %if write_reg
                
            end %if curConditions.sanityCheck_model(itrial)
            
        end %itrial
        
    end %irun
    
end %iSub