%ContextABCD study
%Recognition script
%Halle Dimsdale-Zucker 

%% Intro stuff (like cleaning up, etc.)

initialize_ABCDCon
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock))); %sets random number generators to a random seed
cd(scriptsTask);

%% Flags, etc.

DEBUGGING_FLAG=0;

%% Script-specific variables

subject=input('Enter subject number: ');    
practice=input('Select 1 for practice, 0 for task: ');

% Font size and text formatting 
if strcmp(name, 'stim')
    %[upperLx upperLy lowerRx lowerRy] Where X coord increase going L-->R 
    % and Y coord increase Top-->Bottom (from: http://wiki.biac.duke.edu/biac:courses:ptb_2)  
    image_resize = [364 188 1530 836];
else
    image_resize = [64 188 1216 836];
end %if strcmp

% Instructions text -- NB some of these overwrite the initialize_ABCDCon
% values b/c in the MRI scanner
advancescreen = 'Please SQUEEZE the ball when you are ready to continue.\n';
behavStartScreen = 'If you have any questions, please ask the experimenter now.\n Otherwise, please hit ENTER when you are ready to begin.\n';
behavAdvanceScreen = 'Hit ENTER to advance to the next trial.\n';
btwnruns = 'Take a short break before we continue.\n Please lie still and wait while the scanner finishes running.\n';
endscreen = 'You are finished with this part.\n\n Please lie still and wait while the scanner finishes running.\n';
experimenterContinueScreen = 'Experimenter -- please hit ENTER to advance.\n';
getReadyScreen = 'The task is starting soon. Get ready!\n';
instruc = 'You will now be presented with items and asked\n to recall whether or not you saw the item earlier\n using the following scale:\n\n';
instruc2 = 'The object you see here may be a different size than\n when you studied it originally or presented\n from a slightly different viewpoint.\n\n This is NOT meant to trick you.\n\n It simply reflects that different people will be seeing\n the objects from different perspectives,\n and we have picked the most representative viewpoint.\n';
mem_q = 'Did you see this object presented in one of the houses earlier?\n';
MRIstartscreen = 'If you have any questions, please ask the experimenter now.\n Otherwise, please SQUEEZE the ball when you are ready to begin.\n';
scannerWaitScreen = 'The scanner is warming up. \n\n Please lie still and wait for the task to begin.\n';

% Timing
quesduration = 3.000*fast;
extraLen = 14.000*fast; 

%% Counterbalance info

if mod(subject,2)==1 %returns 1 when subject number is odd
    memscale = '(1)Remember details    (2)Only feels familiar    (3)New\n';
elseif mod(subject,2)==0
    memscale = '(1)New    (2)Only feels familiar    (3)Remember details\n';
end %if

%% Initialize file where will write data

fprintf('Setting up data file.\n');

%setup .dat file to write to later
datafilename = strcat(rawBehavDir, 's',num2str(subject,'%03d'),filesep,'ConABCD_objectRecog_s',num2str(subject,'%03d'),'.dat');  % name of data file to write to but not creating a file

if DEBUGGING_FLAG 
    datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectRecog_s',num2str(subject,'%03d'),'.dat'); 
    datafilepointer = fopen(datafilename,'wt'); 
elseif ~practice
    if subject <99 && fopen(datafilename, 'rt')~=-1 %r = open a file for reading. t = open in text mode. -1 is the value returned when fopen cannot open a file (e.g., if this file doens't exist yet--which is normally what you want bc it means you won't be overwriting anything!)
        tmp = input('Data file already exists! Override? (Y=1) ');
        tmp2 = input('Write over existing file? (Y=3,N=2) ');
        if tmp==1 && tmp2==2
            while fopen(datafilename,'rt')~=-1
                [a b c] = fileparts(datafilename);
                datafilename = strcat(a,filesep,b,'+',c); 
            end
            datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing..w = open for writing, discard existing contents. t = text mode
        elseif tmp==1 && tmp2==3
            % overwrite existing file
            datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectRecog_s',num2str(subject,'%03d'),'.dat'); 
            datafilepointer = fopen(datafilename,'wt'); 
        else
            fclose('all');
            error('Choose a different subject number.');
        end %if tmp==1 && tmp2==2
    else
        datafilepointer = fopen(datafilename,'wt'); 
    end % if subject<99 && fopen(....)
end %if

%% Set up sequence
if practice
    objrec.stim = {'AREDROSE8';'Ablackma3';'tricycle'};
    videos.Lure = [0;1;1];
    nobjects = 3;
    numblocks = 1;
    startblock = 1;
    numtrials = 3;
    objrec.jitter = ones(nobjects,1);
    save([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objrec_practiceSeq_',num2str(subject,'%03d')],'objrec','numblocks','numtrials','startblock','nobjects');
    extraLen = 1.000*fast;
elseif ~practice && subject < 99 && exist([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objrec_seq_',num2str(subject,'%03d'),'.mat'],'file')
    fprintf('Loading sequence for subject %03d.\n',subject)
    load([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objrec_seq_',num2str(subject,'%03d')]);
    startblock = 1; %RESET THIS VALUE IF NEED TO RESTART IN THE MIDDLE
else
    %% Read in stimulus info
    fprintf('Determining sequence for subject %03d.\n',subject)
    fprintf('Reading in stimnames using tdfread.\n');
    videos = tdfread([scriptsTask,'StimList_ByObject_040714.txt']);

    % reformat some fields w/in videos to be cells
    videos.MovieID = cellstr(videos.MovieID);
    videos.ObjectIDinMovie = cellstr(videos.ObjectIDinMovie); 
    videos.ObjectNumber = cellstr(videos.ObjectNumber);

    % load in objectEnc data
    fprintf('Loading in objenc data.\n')
    encFile = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'.dat');
    [aa bb cc] = fileparts(encFile);
    encFile_Alt = [aa,bb,'+',cc];
    if exist(encFile_Alt,'file')
        altfile = input('ObjEnc+ file exists. Load this file?: (Y=1,N=0)');
        if altfile
            encdata = tdfread(encFile_Alt,',');
        else
            encdata = tdfread(encFile,',');
        end %if 
    else 
        encdata = tdfread(encFile,',');
    end %if exist
    
    encdata.ObjectID = cellstr(encdata.ObjectID);
    
    %%  Number of blocks, trials, etc.
    fprintf('Figuring out numblocks, numtrials.\n')
    videoID = unique(videos.MovieID);
    numOldObjects = length(encdata.ObjectID);
    numNewObjects = round((numOldObjects/4))+2;
    nobjects = numOldObjects + numNewObjects;
    
    % check that have 20% new objects
    idealNew = (52/252);
    if numNewObjects/nobjects ~= idealNew
        warning('Incorrect percentage of new items, %.04f%% instead of %.04f%%\n',(numNewObjects/nobjects),idealNew)
    end % if
    
    numblocks = 4;
    startblock = 1;
    numtrials = nobjects/numblocks; 
    % check that numtrials is an integer
    if mod(numtrials,1)
        error('Number of trials per block is not a whole number.\n')
    end % if
    
    %% Set up jitter
    fprintf('Figuring out jitter.\n')
    % set jitter variables
    minjitter = 2.000;
    maxjitter = 8.000;
    idealmean = 4.000;

    redraw = 1;
    while redraw == 1
        % create an exponential distribution of potential jitters
        jitter = exprnd(minjitter,[numtrials,1]);

        % keep re-pulling this distribution until this distribution does not exceed
        % the max value you've set (NB: be careful not to make this value too low or it will fuck with the distribution)
        while max(jitter) >= maxjitter; 
            jitter = exprnd(minjitter,[numtrials,1]);
        end

        % reset the mean to be what you want
        jitdiff = idealmean - mean(jitter);
        jitter = jitter + jitdiff;
        jitter = jitter.*fast;

        % make sure jitter has optimal range
        idealJitRange = 6.500; % 6-7sec is ideal (think re HDR lag)
        jitRange = (max(jitter)) - (min(jitter));
        if ~(((idealJitRange - 1.000) <= jitRange) && (jitRange <= (idealJitRange + 1.000)))
            warning('Jitter range is %d which is not within ideal jitter range of %.02d +/- 1.\n',jitRange,idealJitRange)
            redraw = input('Redraw jitter? (Y=1, N=0): ');
        end % if
    end %while

    totalJitTime = sum(jitter);
    if round(totalJitTime)~=252
        warning('Total jitter (in seconds) is not equal to expected value of 252, instead, is %d\n',totalJitTime)
    end %if 
    h = figure;
    set(gcf,'Visible', 'off')
    hist(jitter,0.5:0.5:max(jitter))
    ylim([0 20])
    title(['Histogram of jitter lengths for s', num2str(subject,'%03d')])
    xlabel('Jitter length (in seconds)')
    ylabel('Frequency of each jitter length')
    print(h,'-dpdf',[rawBehavDir,'s',num2str(subject,'%03d'),filesep,'Jitter']);
    
    
    %% Randomization
    fprintf('Figuring out lures.\n')
    randObjectOrder = randperm(nobjects);

    % Figure out lures
    practiceObjects = videos.ObjectNumber(videos.Lure==2);
    allUnseenObjects = videos.ObjectNumber(~ismember(videos.ObjectNumber, encdata.ObjectID)); 
    possibleLures = allUnseenObjects(~ismember(allUnseenObjects,practiceObjects));
    randLureSelection_part1 = randperm(length(possibleLures));
    randLureSelection = randLureSelection_part1(1,1:numNewObjects);
    usedLures = possibleLures(randLureSelection);

    % Figure out all objects to be seen at recog
    allRecogObjects = vertcat(encdata.ObjectID,usedLures);
    if length(allRecogObjects)~=nobjects
        warning('Incorrect number of recog objects. Currently have %d, should have %d\n',length(allRecogObjects),nobjects)
    end %if

    % Randomize order of objects
    fprintf('Randomizing stimuli and jitter order.\n')
    allRecogObjects_Rand = cell(nobjects,1);
    VideoIDs_Rand = cell(nobjects,1);
    ObjectIDinMovie_Rand = cell(nobjects,1);
    ObjectNumber_Rand = cell(nobjects,1);
    CorrectResp_Rand = zeros(nobjects,1);
    for iobject=1:nobjects
        allRecogObjects_Rand{iobject} = allRecogObjects{randObjectOrder(iobject)};
        % THERE IS PROBABLY A BETTER WAY TO DEAL W/ THESE
        VideoIDs_Rand{iobject} = videos.MovieID{strcmp(videos.ObjectNumber,allRecogObjects{randObjectOrder(iobject)})};
        ObjectIDinMovie_Rand{iobject} = videos.ObjectIDinMovie{strcmp(videos.ObjectNumber,allRecogObjects{randObjectOrder(iobject)})};
        ObjectNumber_Rand{iobject} = videos.ObjectNumber{strcmp(videos.ObjectNumber,allRecogObjects{randObjectOrder(iobject)})};
        CorrectResp_Rand(iobject) = ~sum(strcmp(videos.ObjectNumber{strcmp(videos.ObjectNumber,allRecogObjects{randObjectOrder(iobject)})},possibleLures));
    end %for iobject=

    %split into runs and randomize jitters
    %trial orders for study and test
    for ii=1:numblocks
        objrec(ii).stim = allRecogObjects_Rand((ii.*numtrials)-(numtrials-1):ii.*numtrials);
        objrec(ii).VideoID = VideoIDs_Rand((ii.*numtrials)-(numtrials-1):ii.*numtrials);
        objrec(ii).ObjectIDinMovie = ObjectIDinMovie_Rand((ii.*numtrials)-(numtrials-1):ii.*numtrials);
        objrec(ii).ObjectNumber = ObjectNumber_Rand((ii.*numtrials)-(numtrials-1):ii.*numtrials);
        objrec(ii).CorrectResp = CorrectResp_Rand((ii.*numtrials)-(numtrials-1):ii.*numtrials);
        objrec(ii).jitter = jitter(randperm(size(jitter,1)));
    end
    
    %% Save randomization order before running
    fprintf('Saving out sequence.\n')
    save([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objrec_seq_',num2str(subject,'%03d')],'objrec','numblocks','numtrials','startblock','nobjects');
end % if practice 

%% Psychtoolbox stuff

fprintf('Setting up PTB stuff.\n');
%Check to see if has OpenGL compatibility; will abort if not compatible.
AssertOpenGL;  

% Dummy calls (prevents delays)
KbCheck; %Check the status of the keyboard. If any key is down (i.e., on/pressed), this will return the value 1. If not, 0.
WaitSecs(waitbtwninstruc); %Have script wait for (n) seconds (basically while matlab loads)
GetSecs; %report time back in a very accurate way

%swap from PTB's internal naming system to current keybord naming system
KbName('UnifyKeyNames'); 
%name the keys that subjects will be using
%name the keys that subjects will be using
enter = KbName('return');
escapeKey = KbName('ESCAPE');
oneResp = KbName('1!');
twoResp = KbName('2@');
threeResp = KbName('3#');
fourResp = KbName('4$');
trigger = KbName('5%');
enablekeys = [oneResp twoResp threeResp fourResp enter escapeKey];
keylist=zeros(1,256);%%create a list of 256 zeros
keylist(enablekeys)=1;%%set keys you interested in to 1

% Initialize PsychHID (MEX file that communicates with HID-compliant devices such as USB keyboards, etc.) and get list of devices
clear PsychHID;
if strcmp(name,'red');
    resp_device = [];
else
    devices = PsychHID('devices'); 
    %If Current Designs Keyboard is detected, make it the response device. Otherwise, use any keyboard.
    if ~isempty(find(strcmp({devices.manufacturer}, 'Current Designs, Inc.') & strcmp({devices.usageName}, 'Keyboard'),1));
        resp_device = find(strcmp({devices.manufacturer}, 'Current Designs, Inc.') & strcmp({devices.usageName}, 'Keyboard'));
    else
        resp_device = find(strcmp({devices.usageName}, 'Keyboard'));
    end
end %if 

%% Initialize the experiment 

%try/catch statement
try
    %% Set up PTB stuff (like the screen) 
    
    fprintf('Setting up PTB screen.\n');
    % Get screenNumber of stimulation display, and choose the maximum index, which is usually the right one.
    screens=Screen('Screens');
    screenNumber=max(screens);
    HideCursor;
    
    % Open a double buffered fullscreen window on the stimulation screen
    % 'screenNumber' and use 'gray' for color (NB: to truly use 'gray' this
    % would be [128 128 128]. [50 50 50] is a darker version of gray.
    % 'w' is the handle used to direct all drawing commands to that window
    % 'wRect' is a rectangle defining the size of the window. See "help PsychRects" for help on such rectangles
    [w, wRect]=Screen('OpenWindow',screenNumber, screenColor);
    [mx, my] = RectCenter(wRect);
    scene_x = 350;
    scene_y = 350;
    scene_mtx = [mx-scene_x/2,my-scene_y/2,mx+scene_x/2,my+scene_y/2]';
    if strcmp(name,'stim')
        ytext = (my*2)-250;  
    else
        ytext = (my*2)-100;
    end % if strcmp

    
    % Set text size
    Screen('TextSize', w, fontSize);
    
    Priority(MaxPriority(w));   % Set priority for script execution to realtime priority
    
    if strcmp(computer,'PCWIN') == 1
        ShowHideWinTaskbarMex(0);
    end
    
    % Initialize KbCheck and return to zero in case a button is pressed
    [KeyIsDown, endrt, KeyCode]=KbCheck;
    KeyIsDown = zeros(size(KeyIsDown));
    endrt = zeros(size(endrt));
    KeyCode = zeros(size(KeyCode));
    WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1 
    
    %% Give instructions
    
    fprintf('Giving instructions.\n');
    %pull the first screen of instructions
    if ~practice
        DrawFormattedText(w, [instruc linebreak memscale linebreak linebreak instruc2 linebreak MRIstartscreen linebreak experimenterContinueScreen], 'center', 'center', [255 255 255]);  % Write instructions
    elseif practice 
        DrawFormattedText(w, [instruc linebreak memscale linebreak linebreak instruc2 linebreak behavStartScreen], 'center', 'center', [255 255 255]);  % Write instructions
    end %if ~practice
    Screen('Flip', w);   % Update display to show instructions
    
    % Display instructions 'til KbCheck detects enter pressed
    %keyboard
    while (KeyCode(enter)==0)
        [KeyIsDown, RT, KeyCode]=KbCheck; 
        WaitSecs(waitbtwninstruc);    % Prevents overload  
    end
    
    %% Recog task
    
    fprintf('Starting the experiment.\n');
    
    %Clear out KBCheck, set things to 0, etc. 
    [KeyIsDown, endrt, KeyCode]=KbCheck;
    KeyIsDown = zeros(size(KeyIsDown));
    endrt = zeros(size(endrt));
    KeyCode = zeros(size(KeyCode));
    WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1
    
    for iblock=startblock:numblocks
%      for iblock=1 % use for testing of script
        
        if strcmp(name,'stim')
            DrawFormattedText(w, scannerWaitScreen, 'center', 'center', [255 255 255]);  % Write instructions
            Screen('Flip', w);   % Update display to show instructions

            %%% Wait for run to start
            %%% Change for FMRI to receive scanner trigger
            while (KeyCode(trigger)==0)
                [KeyIsDown, firstpulse, KeyCode]=KbCheck;
                WaitSecs(waitbtwninstruc);    % Prevents overload
            end
        
            %Clear out KBCheck, set things to 0, etc. 
            [KeyIsDown, endrt, KeyCode]=KbCheck;
            KeyIsDown = zeros(size(KeyIsDown));
            endrt = zeros(size(endrt));
            KeyCode = zeros(size(KeyCode));
            WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1

            % Wait til after the trigger to disable the trigger check
            RestrictKeysForKbCheck(enablekeys);

            %Added for KbQueue
            KbQueueCreate(resp_device,keylist);

            DrawFormattedText(w, getReadyScreen, 'center', 'center', [255 255 255]);  % Write fixation
        else
            %Clear out KBCheck, set things to 0, etc. 
            [KeyIsDown, firstpulse, KeyCode]=KbCheck;
            KeyIsDown = zeros(size(KeyIsDown));
            endrt = zeros(size(endrt));
            KeyCode = zeros(size(KeyCode));
            WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1
            
            % Wait til after the trigger to disable the trigger check
            RestrictKeysForKbCheck(enablekeys);

            %Added for KbQueue
            KbQueueCreate(resp_device,keylist);
            
            DrawFormattedText(w,fixation,'center','center',[255 255 255]);
        end %strcmp(...
        
        [VBLTimestamp expstart]=Screen('Flip', w);  % Update display to show +
        WaitSecs(4.000*fast);    % Display + for 4 seconds
        
        %Clear out KBCheck, set things to 0, etc. 
        [KeyIsDown, endrt, KeyCode]=KbCheck;
        KeyIsDown = zeros(size(KeyIsDown));
        endrt = zeros(size(endrt));
        KeyCode = zeros(size(KeyCode));
        WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1
        
        %initiate timecounter to keep trials on track
        timecounter = GetSecs;

       trialcounter = 0;
       for iobject=1:size(objrec(iblock).stim,1)
%         for iobject=1:2; %%USE FOR DEBUGGING
           trialcounter = trialcounter + 1;
           
                Screen('Close');    % Close to not overload

                %Initialize KbCheck and return to zero in case a button is pressed
                [KeyIsDown, endrt, KeyCode]=KbCheck;
                KeyIsDown = zeros(size(KeyIsDown));
                endrt = zeros(size(endrt));
                KeyCode = zeros(size(KeyCode));
                WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1

                %load image
                if ~practice
                    probe = [stimAloneDir,objrec(iblock).stim{iobject},file_ext]; 
                else
                    probe = [practicestimdir,objrec(iblock).stim{iobject},file_ext]; 
                end %if
                probe = imread(probe); 
                
                %present image
                showprobe = Screen('MakeTexture',w,probe);
                Screen('DrawTextures',w,showprobe,[],image_resize);
                DrawFormattedText(w,[mem_q, linebreak memscale],'center', ytext, [255 255 255]);
                [VBLtimestamp memqstart]=Screen('Flip',w);

                % initialize KbCheck and variables to avoid time dealys
                KeyIsDown = zeros(size(KeyIsDown));
                endrt = zeros(size(endrt));
                KeyCode = zeros(size(KeyCode));
                %Added for KbQueue
                KbQueueFlush();
                KbQueueStart();
                mem_resp = '';
                firstpress = [];

                while (GetSecs - memqstart) < quesduration
                    if ~KeyIsDown%(enablekeys)  % if key is pressed, stop recording response
                         [KeyIsDown firstpress]=KbQueueCheck(resp_device);
                         WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
                    end
                    mem_resp = KbName(find(firstpress,1,'last'));
                    if isempty(mem_resp) %added to deal w when subj doesn't make mem resp
                        mem_resp = 'xx';
                    end %
                    endrt = min(firstpress(find(firstpress)));
                end %while (GetSecs - memqstart) < quesduration
                
                DrawFormattedText(w, fixation, 'center', 'center', [255 255 255]);  % Write fixation to back buffer
                [VBLTimestamp jitterstart]=Screen('Flip', w);  % show fixation
                
                % continue recording during jittered ISI
                % increment timecounter and cut into jitter if necessary to
                % stay on track
                % timecounter is absolute time. computer time + time that
                % should have been spent on question + ideal amt of time on
                % jitter. gives you the clock time you want when jitter
                % ends. deals w/ if have lag in loading trials, etc. if
                % there;s a lag, cut this out of jitter. this holds
                % fixation cross on the screen until the next screen flips.
                % (no need to use waitsecs). 
                timecounter = timecounter + quesduration + objrec(iblock).jitter(iobject);
                while (GetSecs) < timecounter
                    if ~KeyIsDown%(enablekeys)  % if key is pressed, stop recording response
                        [KeyIsDown firstpress]=KbQueueCheck(resp_device);
                        WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
                    end
                    mem_resp = KbName(find(firstpress,1,'last'));
                    if isempty(mem_resp) %added to deal w when subj doesn't make mem resp
                        mem_resp = 'xx';
                    end %
                    endrt = min(firstpress(find(firstpress)));
                end
                
                % compute reaction time
                mem_rt = round(1000*(endrt-memqstart)); % computes RT for each trial. if wanted to know RT in absolute time, would use min(firstpress(find(firstpress))
                
                % Added for KbQueue
                if firstpress(escapeKey)
                    error('Escape key pressed');
                end

                if ~practice
                        %Write trial result to file: 
                        %first, create a header line for the output file
                        if iblock==1 && trialcounter==1
                            fprintf(datafilepointer,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
                                'SubjectID','VideoID','ObjectIDinMovie','ObjectID','ExptStart',...
                                'MemQuesStart','MemResp','CorrectMemResp','MemRT','JitterLength','Onset','BlockID','FirstPulseOnset');
                        end %if
                        %now write data to .dat file
                        fprintf('Saving data.\n');
                        fprintf(datafilepointer,'%03d,%s,%s,%s,%d,%d,%s,%d,%d,%d,%d,%d,%d\n',...
                            subject,...
                            objrec(iblock).VideoID{iobject},...
                            objrec(iblock).ObjectIDinMovie{iobject},...
                            objrec(iblock).ObjectNumber{iobject},...
                            expstart,...
                            memqstart,...
                            char(mem_resp),...
                            objrec(iblock).CorrectResp(iobject),...
                            mem_rt,...
                            objrec(iblock).jitter(iobject),...
                            memqstart - firstpulse,...
                            iblock,...
                            firstpulse);
                end %if ~practice

                if practice
                    DrawFormattedText(w,['Please explain your response to the experimenter.\n\n' behavAdvanceScreen],'center', 'center', [255 255 255]); 
                    Screen('Flip', w);   % Update display to show instructions
                    while (KeyCode(enter)==0)
                        [KeyIsDown, foo, KeyCode]=KbCheck;
                        WaitSecs(waitbtwninstruc);    % Prevents overload11
                    end
                    %Initialize KbCheck and return to zero in case a button is pressed
                    [KeyIsDown, endrt, KeyCode]=KbCheck;
                    KeyIsDown = zeros(size(KeyIsDown));
                    endrt = zeros(size(endrt));
                    KeyCode = zeros(size(KeyCode));
                    WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1
                    KbQueueFlush(); %Added for KbQueue
                end %if practice

                if iobject==nobjects && practice
                    DrawFormattedText(w, 'You have completed the practice.\n','center', 'center', [255 255 255]);
                    Screen('Flip', w);   % Update display to show instructions
                    while (KeyCode(enter)==0)
                        [KeyIsDown, foo, KeyCode]=KbCheck;
                        WaitSecs(waitbtwninstruc);    % Prevents overload
                    end
                elseif iobject==nobjects && ~practice
                    DrawFormattedText(w, endscreen, 'center', 'center', [255 255 255]);  % Write instructions
                    Screen('Flip', w);   % Update display to show instructions
                    while (KeyCode(enter)==0)
                        [KeyIsDown, foo, KeyCode]=KbCheck;
                        WaitSecs(waitbtwninstruc);    % Prevents overload
                    end
                end %if

                %Initialize KbCheck and return to zero in case a button is pressed
                [KeyIsDown, endrt, KeyCode]=KbCheck;
                KeyIsDown = zeros(size(KeyIsDown));
                endrt = zeros(size(endrt));
                KeyCode = zeros(size(KeyCode));
                WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1
                KbQueueFlush(); %Added for KbQueue
        end %iobject=
        
        %%% Add in extra timepoints for FMRI
        DrawFormattedText(w, fixation, 'center', 'center', [255 255 255]);  % Write fixation to back buffer
        [VBLTimestamp jitterstart]=Screen('Flip', w);  % show fixation
        
        % continue recording during jittered ISI
        while (GetSecs - jitterstart) < extraLen
            if ~KeyIsDown  % if key is pressed, stop recording response
                [KeyIsDown, endrt, KeyCode]=KbCheck;
                WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
            end
        end
        
        % Break until hit enter
        if iblock<numblocks
            DrawFormattedText(w, [btwnruns linebreak experimenterContinueScreen], 'center', 'center', [255 255 255]);  % Write instructions
        else
            DrawFormattedText(w, endscreen, 'center', 'center', [255 255 255]);  % Write instructions
        end
        
        Screen('Flip', w);   % Update display to show instructions
        
        while (KeyCode(enter)==0)
            [KeyIsDown, foo, KeyCode]=KbCheck;
            WaitSecs(waitbtwninstruc);    % Prevents overload
        end
        
        %reset key restriction to get trigger again
        RestrictKeysForKbCheck([]);
    end %iblock=
    %% Finish up %%
    
    % Cleanup at end of experiment
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    ShowHideWinTaskbarMex(1)
    
    % End of experiment:
    return;

catch %try
     % catch error in case something goes wrong in the 'try' part
    % Do same cleanup as at the end of a regular session
    
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    ShowHideWinTaskbarMex(1)
    
    psychrethrow(psychlasterror);   % Output the error message that describes the error
end