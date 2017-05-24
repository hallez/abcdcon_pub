% Loop across betas, and write out trial id information
%
% Halle R. Dimsdale-Zucker

%% Setup some basics
initialize_ABCDCon

%% Loop across subjects 
for isub=1:length(subjects)

    fprintf('\n----Working on s%s----\n',num2str(subjects(isub),'%03d'))

    b.curSubj = ['s' num2str(subjects(isub),'%03d')];
    b.analMRIDir = analMRIDir;
    b.modelDir = [analMRIDir,'multivariate_sanityCheck'];
    % if last subject overwrote b.runs, revert back 
    % //TODO figure out a better way of dealing w/ this 
    runbase = 'run';
    runnums = 1:4;
    nruns = length(runnums);
    b.runs = cell(1,nruns);
    for irun=1:nruns
        b.runs(irun) = cellstr([runbase,num2str(runnums(irun))]);
    end %for irun=

    % Run exceptions for subject-specific naming conventions
    b = run_exceptions_ABCDCon(b);

    %% Loop across runs
    for irun=1:length(b.runs)

        b.curRun = b.runs{irun};
        b.cur_run_num = str2double(strtok(b.curRun, 'run'));
        fprintf('Working on: %s\n',b.curRun);

        %% Loop across trials
        curRunDir = [b.modelDir,filesep,b.curSubj,filesep,b.curRun];

        curDirs = dir(curRunDir);
        trialCounter = 0;

        % figure out how many trials were in the current run (should always
        % be 63, but this is coded for generalizability) 
        for j = 1:size(curDirs,1)
            s=curDirs(j).name;
            if strfind(s,'trial_')
                trialCounter = trialCounter + 1;
            end %if strfind
        end %j=

        % check to ensure this is a reasonable amount of trials
        if trialCounter <= (numRecogTrials/nruns)
            numCurTrials = 1:trialCounter;
        else 
            error('Implausible number of trials for current run: %s.\n',b.curRun)
        end %if trialCounter

        % initialize a variable the size of the current number of trials
        pattern_trial_id = cell(1,trialCounter);

        for itrial=1:trialCounter

            b.curTrial = ['trial_',num2str(numCurTrials(itrial),'%03d')];

            % read in the beta_0001.img for the current trial and
            % run
            b.curTrial_dir = [b.modelDir,filesep,b.curSubj,filesep,b.curRun,filesep,b.curTrial];
            b.cur_beta=spm_vol([b.curTrial_dir,filesep,'beta_0001.img']);
            b.cur_beta.img = spm_read_vols(b.cur_beta);

            % figure out trial-specific label
            [TOKEN, REMAIN] = strtok(cellstr(b.cur_beta.descrip),'R');
            [TOKEN, REMAIN] = strtok(REMAIN,'*');
            pattern_trial_id(:,itrial) = TOKEN; 

        end %itrial=  

        % save out ids
        [path roi_name ext] = fileparts(b.cur_roi.fname);
        fname_ids_out = cellstr([b.modelDir,filesep,b.curSubj,filesep,'pattern_mtx_ids','_run_',num2str(b.cur_run_num,'%02d')]);
        save(fname_ids_out{:},'pattern_trial_id')

    end %irun=
                
end %isub

