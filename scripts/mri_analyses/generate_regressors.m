%This simple SPM regressor generation script assumes that there is a single
%column that codes conditions, and that all trials of interest are coded
%with a value less than conmax and all rows of no interest are coded with
%a value greater than conmax. Note that this script has no error detection
%and extremely limited flexibility for complex designs.
%
%M. Ritchey, Nov 2011
%
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

% enable selecting to analyze different models
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
include_dur = 0; %do you want to include durations? typically this is no

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
    
    onsets_fname = 'onsets_file.csv';
    
    datafile = [rawBehavDir,b.curSubj,filesep,onsets_fname]; 
    
    conditions = tdfread(datafile,','); 

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %contrast manager will store which trial types are in which column
    contrast_man = [];
    contrast_trinum = [];
    
    %%%Generate regressor files
    for iRun = 1:length(b.runs)
        
        cur_run = b.runs{iRun};

        cur_run_num = str2double(strtok(cur_run, 'run'));        
        curRunIdx = conditions.BlockID==cur_run_num;
        
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
               

                names = [names regnames(regids==curRegsIDs(iReg))]; %append name of current regressor (NB since embedded earlier in if rowids only adds when regressor exists for current run)
                
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