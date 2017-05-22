%This simple SPM regressor generation script assumes that there is a single
%column that codes conditions, and that all trials of interest are coded
%with a value less than conmax and all rows of no interest are coded with
%a value greater than conmax. Note that this script has no error detection
%and extremely limited flexibility for complex designs.
%
%M. Ritchey, Nov 2011
%
%Updated Jan 2012 by MR
%   - includes smart handling of cases in which not all trial types exist
%   for every run
%   - generates contrast vectors for each trial type and outputs to file
%Updated April 2012 by MR
%   - includes the option to specify a regressor file for padding contrast
%   vectors
%   - includes the option to weight contrasts by the proportional number of
%   trials per run, rather than weighting runs equally
%Updated June 2012 by MR
%   - includes the option to have duration set to zero or a behavioral
%   variable
%   - writes out a .mat file to the model directory with all of the model
%   specification values: model_info.mat
%   - includes sections for specifying run exceptions (e.g., which runs to
%   include, run naming)
%Updated August 2014 by HDZ

initialize_ABCDCon

%%%MODEL DEFINITIONS 
allModels(1).name = 'sanityCheck_model';
allModels(1).type = 'univariate';
allModels(1).contrasts = {[1 2 3 4 5 6] [1 2] [1 2;5 NaN] [1 2;6 NaN] [1 2;3 4]...
    [1 3] [1] [1;5] [1;6] [1;3]...
    [2 4] [2] [2;5] [2;6] [2;4]};
allModels(1).contrastNames = {'all' 'allRHits' 'allRHitsxCR' 'allRHitsxFA' 'allRHitsxFHits_Miss' ...
    'allBrown' 'brownRHits' 'brownRHitsxCR' 'brownRHitsxFA' 'brownRHitsxFHits_Miss'...
    'allGray' 'grayRHits' 'grayRHitsxCR' 'grayRHitsxFA' 'grayRHitsxFHits_Miss'};

allModels(2).name = 'objrecANDlocrec_model'; 

allModels(3).name = 'byRoombyHouse_model';
allModels(3).type = 'univariate';
allModels(3).contrasts = {[1] [2] [3] [4]...
    [5] [6] [7] [8]...
    [9]};
allModels(3).contrastNames = {'Brown_Rm1_RHit' 'Brown_Rm2_RHit' 'Gray_Rm1_RHit' 'Gray_Rm2_RHit' ...
    'Brown_Rm1_FHitsANDMisses' 'Brown_Rm2_FHitsANDMisses' 'Gray_Rm1_FHitsANDMisses' 'Gray_Rm2_FHitsANDMisses' ...
    'CR'};

% NB: this model is essentially redundant w/ 'stillsVSvideoBYhouse_model'
% except that the houses do not get their own regressors in this model (so
% we might think this model is doing a slightly worse job of capturing the
% data...?)
allModels(4).name = 'stillsVSvideo_model'; 
allModels(4).type = 'roomLocalizer';
allModels(4).contrasts = {[1] [2] [1 2] [1;2]};
allModels(4).contrastNames = {'StillsxBaseline', 'VideosxBaseline', 'StillsANDVideosxBaseline','StillxVideo'};

allModels(5).name = 'stillsVSvideoBYhouse_model';
allModels(5).type = 'roomLocalizer';
allModels(5).contrasts = {[1 2] [3 4] [1 2 3 4] [1 2; 3 4] [1;2] [3; 4] [1 3; 2 4]};
allModels(5).contrastNames = {'StillsxBaseline', 'VideosxBaseline', 'StillsANDVideosxBaseline', 'StillxVideo', 'BrownStillxGrayStill', 'BrownVideoxGrayVideo', 'BrownStillorVideoxGrayStillorVideo'};

allModels(6).name = 'stillsBYhouseNOVideoRegressors_model';
allModels(6).type = 'roomLocalizer';
allModels(6).contrasts = {[1 2] [1;2]};
allModels(6).contrastNames = {'StillsxBaseline', 'BrownStillxGrayStill'};

modelSelect = 0;
for imodel=1:size(allModels,2)
    curModel = allModels(imodel).name;
    curModelBase = strtok(curModel,'_');
    curModelType = allModels(imodel).type;
    contrasts = allModels(imodel).contrasts;
    contrastNames = allModels(imodel).contrastNames;
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

%%%MODEL OPTIONS
conLabel = '_cons';
write_reg = 1; %do you want to write out regressor mats?
write_con = 1; %do you want to write out contrasts?
include_dur = 1; %do you want to include durations? typically this is NO, but for roomLocalizer models we DO want to

%include motion regressors in the model? 1 for standard 6-param motion regs, 2 for motion+spikes
inc_motion = 2; 
spike_regs = 'spike_regs_rp.txt'; %used only if inc_motion==2 

%weight contrasts by trial numbers per run?
weight_runs = 0; %0 for no weighting (standard), 1 for weighting

%%%END OF MODIFICATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(outputDir)
    mkdir(outputDir);
end
save(strcat(outputDir,'model_info'));

for isub = 1:length(subjects)
    b.curSubj = ['s',num2str(subjects(isub),'%03d')]; 

    fprintf('\n----Working on subject %s----\n',b.curSubj)
    
    if ~isdir([outputDir,b.curSubj])
        mkdir(strcat(outputDir,b.curSubj));
    end
    
    if strcmp(curModelType, 'roomLocalizer')
        onsets_fname = 'room_localizer_onsets_file.csv';
    else
        onsets_fname = 'onsets_file.csv';
    end
    
    datafile = [rawBehavDir,b.curSubj,filesep,onsets_fname]; 
    
    conditions = tdfread(datafile,','); 
    % eliminate any rows where the model has a value of NaN
    % this is specific to `stillsBYhouseNOVideoRegressors_model` where
    % don't want to include video trials (which have been coded in R as
    % `NaN` in any of the regressors)
    % based on: http://www.mathworks.com/matlabcentral/newsreader/view_thread/236371
    % and also http://stackoverflow.com/questions/23932062/deleting-multiple-rows-from-all-fields-of-a-given-structure-array-using-matlab
    if strcmp(curModel,'stillsBYhouseNOVideoRegressors_model')
        non_na_rows = find(~isnan(conditions.stillsBYhouseNOVideoRegressors_model));
        % structfun(function, performed on structure, parameter (in this case, return
        % a structure w/ same field types as original), logical setting of
        % parameter)
        conditions = structfun(@(v) v(non_na_rows), conditions, 'Uniform', 0);
    end

    if ~isfield(conditions,curModel) 
        error('There is no column for the current model (%s) in onsets file (%s).\n',curModel,datafile);
    end %if 
    regs = conditions.(curModel);     
    regids = unique(regs); 
    tmp1 = cellstr(conditions.([curModelBase '_key']));
    for ireg=1:length(regids)
        % make sure regids are in the same order as regnames
        regnames{ireg} = cell2mat(unique(tmp1(regs==regids(ireg))));
    end %ireg=
    excludereg = 99;
    
    % replace values that may change between subjects
    % if there are any
    b.runs = {'run1','run2','run3','run4' };
    
    b = run_exceptions_ABCDCon(b);
    
    % set up things that are specific to the room_localizer run
    if strcmp(curModelType, 'roomLocalizer')
        b.runs = {'room_localizer'};
        numRecogTrials = 50; % 42 stills and 8 videos
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %contrast manager will store which trial types are in which column
    contrast_man = [];
    contrast_trinum = [];
    
    %%%Generate regressor files
    for iRun = 1:length(b.runs)
        
        cur_run = b.runs{iRun};
        
        if strcmp(curModelType, 'roomLocalizer')
            % since there is just one run, set up a column of 1s (ie, pull
            % all trials in this run). make sure this is a logical index so
            % that when use to pull out `regs` later it's getting values
            % in those rows
            curRunIdx = true(size(regs,1),1);
        else
            cur_run_num = str2double(strtok(cur_run, 'run'));        
            curRunIdx = conditions.BlockID==cur_run_num;
        end
        
        curfuncrun = [analMRIDir,b.curSubj,filesep,cur_run,filesep];
        filename = [outputDir,b.curSubj,filesep,curModel,'_',b.curSubj,'_',cur_run,'_regs'];
        
        % figure out regressors that occur w/in the current run
        % NB: this will throw a warning for subjects who have less than 4
        % runs (eg, s015)
        curRegs = regs(curRunIdx);
        if length(curRegs)~= (numRecogTrials/length(b.runs)) 
            warning('Check number of trials in %s.\n',cur_run)
        end %if
        curRegsIDs = unique(curRegs);
        numCurRegs = length(curRegsIDs);
        
        onsets = {}; names = {}; durations = {};
        for iReg = 1:numCurRegs
            %find trial types of interest
            if curRegs(iReg) <= excludereg
                rowids = find((regs==curRegsIDs(iReg)) & curRunIdx); 
            end

            %if they exist in this run, add the regressor information
            if rowids
                if include_dur %USE conditions.duration 
                    add_dur = conditions.duration(rowids)'; 
                else
                    add_dur = 0;% zeros(1,size(conditions.duration(curRunIdx),1));
                end
                onsets = [onsets conditions.Onset(rowids)']; 
                durations = [durations add_dur]; %append current value for durations onto prior value (b/c cellstr will put into next column)
               
                % append run number to regressor if running
                % byRoombyHouse_model, otherwise don't worry about that
                if(strcmp(curModel,'byRoombyHouse_model')==1)
                    names = [names strcat('run', cur_run_num, '_', regnames(regids==curRegsIDs(iReg)))];
                else
                    names = [names regnames(regids==curRegsIDs(iReg))]; %append name of current regressor (NB since embedded earlier in if rowids only adds when regressor exists for current run)
                end
                
                contrast_man = [contrast_man regids(regids==curRegsIDs(iReg))];
                contrast_trinum = [contrast_trinum length(onsets)]; %number of trials for current regressor
            else
                sprintf('Warning on %s: No trials for %s in %s\n',b.curSubj,regnames{regids==curRegs(iReg)},cur_run) 
            end
        end %iReg
        
        %save a set of regressors for each run
        if write_reg 
            if exist(filename,'file')
                regsOverwrite = input('Regs file already exists. Overwrite? (Y=1,N=0):');
                if regsOverwrite
                    save(filename,'names','onsets','durations');
                else
                    save([filename '+'],'names','onsets','durations'); %DEAL W WHEN ALREADY HAVE + FILE
                end %if regsOverwrite
            else
                save(filename,'names','onsets','durations');
            end %if exist(filename
        end
        
        %add motion regressor columns to contrast manager if desired
        if inc_motion==1 %pad run with 6 motion parameters
            contrast_man = [contrast_man repmat(excludereg+1,1,6)];
            contrast_trinum = [contrast_trinum repmat(0,1,6)];
        elseif inc_motion==2 %pad run with # columns in spike_regs
            curspikes = strcat(curfuncrun,spike_regs);
            spk = importdata(curspikes);
            spkcols = size(spk,2);
            contrast_man = [contrast_man repmat(excludereg+1,1,spkcols)];
            contrast_trinum = [contrast_trinum repmat(0,1,spkcols)];
        end
        
    end %iRun
    
    %%%Generate contrast vectors
    contrast_vectors = {};
    for iCon = 1:length(contrasts)
        curVector = zeros(size(contrast_man));
        curCon = contrasts{iCon};
        %find positively- and negatively-weighted columns
        for j=1:size(curCon,2)
            curVector(contrast_man==curCon(1,j)) = 1;
            if(size(curCon,1)>1)
                curVector(contrast_man==curCon(2,j)) = -1;
            end
        end
        %scale beta weights to sum to one
        if weight_runs==1 %weight runs by #trials per run
            trisums = abs(curVector.*contrast_trinum);
            curVector(curVector==1) = trisums(curVector==1)./(sum(contrast_trinum(curVector==1)));
            curVector(curVector==-1) = -1.*(trisums(curVector==-1)./(sum(contrast_trinum(curVector==-1))));
        else %equally weight all runs
            posweight = 1./(size(find(curVector==1),2));
            negweight = -1./(size(find(curVector==-1),2));
            curVector(curVector==1)=posweight;
            curVector(curVector==-1)=negweight;
        end
        contrast_vectors = [contrast_vectors curVector];
    end %iCon
    
    %save a single contrast file per subject
    filename = strcat(outputDir,b.curSubj,filesep,curModel,'_',b.curSubj,conLabel);
    if write_con
        save(filename,'contrast_vectors','contrastNames');
    end
    
end %iSub