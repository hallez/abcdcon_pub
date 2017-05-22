%ContextABCD study
%Recognition script
%Halle Dimsdale-Zucker 


%% Intro stuff (like cleaning up, etc.)

initialize_ABCDCon
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock))); %sets random number generators to a random seed
cd(scriptsTask);

%% Flags
DEBUGGING_FLAG=0;

%% Script-specific variables

subject=input('Enter subject number: ');    % Ask for subject number input first
practice=input('Select 1 for practice, 0 for task: ');

%Instructions text
instruc = 'You will now be presented with items and asked\n to recall where you saw the item earlier\n using the following scale:\n\n';
instruc2 = 'If you cannot remember where you saw the object,\n you should make your best guess.\n';
location_q = 'Where did you see this object earlier?\n';
startscreen = 'If you have any questions, please ask the experimenter now.\n\n Otherwise, please hit ENTER to begin.\n';

%Timing
quesduration = 3.000*fast;

%% Counterbalance info

fprintf('Loading context enc CB info.\n');
load(strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_contextencCB_s',num2str(subject,'%03d')));

fprintf('Determining other CB.\n');
if practice==1 %&& mod(subject,2)==1
    name1 = 'Halle';
    name2 = 'Rob';
    locationscale = ['(1)' name1 'Rm1   '     '(2)' name1 'Rm2  '      '(3)' name2 'Rm1     '      '(4)' name2 'Rm2\n'];
elseif ~practice && strcmp(name1,'Alex') && strcmp(house1,'C:\Users\dynamic\Desktop\Halle\abcdcon_code\stimuli\Movies_Mp4Versions\Model_25_Brown_HiddenLayers_BirdsEyeAtEnd.mp4')
    locationscale = ['(1)' name1 'Rm1   '     '(2)' name1 'Rm2  '      '(3)' name2 'Rm1     '      '(4)' name2 'Rm2\n'];
elseif ~practice && strcmp(name1,'Jamie') && strcmp(house1,'C:\Users\dynamic\Desktop\Halle\abcdcon_code\stimuli\Movies_Mp4Versions\Model_25_Brown_HiddenLayers_BirdsEyeAtEnd.mp4')
    locationscale = ['(1)' name1 'Rm1   '     '(2)' name1 'Rm2  '      '(3)' name2 'Rm1     '      '(4)' name2 'Rm2\n'];
elseif ~practice && strcmp(name1,'Alex') && strcmp(house1,'C:\Users\dynamic\Desktop\Halle\abcdcon_code\stimuli\Movies_Mp4Versions\Model_25_Gray_HiddenLayers_BirdsEyeAtEnd.mp4')
    locationscale = ['(1)' name2 'Rm1   '     '(2)' name2 'Rm2  '      '(3)' name1 'Rm1     '      '(4)' name1 'Rm2\n'];
elseif ~practice && strcmp(name1,'Jamie') && strcmp(house1,'C:\Users\dynamic\Desktop\Halle\abcdcon_code\stimuli\Movies_Mp4Versions\Model_25_Gray_HiddenLayers_BirdsEyeAtEnd.mp4')
    locationscale = ['(1)' name2 'Rm1   '     '(2)' name2 'Rm2  '      '(3)' name1 'Rm1     '      '(4)' name1 'Rm2\n'];
end %if

%% Initialize file where will write data

fprintf('Setting up data file.\n');

%write to .dat files
datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_locationRecog_s',num2str(subject,'%03d'),'.dat');  % name of data file to write to but not creating a file

% read existing result file and prevent overwriting files from a previous subject (except for subjects > 99)
if DEBUGGING_FLAG 
    datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_locationRecog_s',num2str(subject,'%03d'),'.dat'); %overwrite existing file
    datafilepointer = fopen(datafilename,'wt');
elseif ~practice
    if subject <99 && fopen(datafilename, 'rt')~=-1; %r = open a file for reading. t = open in text mode. -1 is the value returned when fopen cannot open a file (e.g., if this file doens't exist yet--which is normally what you want bc it means you won't be overwriting anything!)
        tmp = input('Data file already exists! Override? (Y=1) ');
        tmp2 = input('Write over existing file? (Y=3,N=2) ');
        if tmp==1 && tmp2==2
            while fopen(datafilename,'rt')~=-1;
                [a b c] = fileparts(datafilename);
                datafilename = strcat(a,filesep,b,'+',c); %don't ever overwrite orig file--if a file for that subject number already exists and you override (ie, say it's okay to reuse that subject number), create the file but append '+' to the end. 
            end
            datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing..w = open for writing, discard existing contents. t = text mode
        elseif tmp==1 && tmp2==3
            %overwrite existing file
            datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_locationRecog_s',num2str(subject,'%03d'),'.dat'); %overwrite existing file
            datafilepointer = fopen(datafilename,'wt');
        else
            fclose('all');
            error('Choose a different subject number.');
        end
    else
        datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
    end
end %if

%% Save out sequence

if practice
    % "Read in" stimulus info
    load([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objenc_practiceSeq_',num2str(subject,'%03d')]);
    nobjects = length(stimlist.objInScene);
    videos.Lure = ones(1,nobjects);
    randObjectOrder = randperm(nobjects);

    % Save out sequence
    save([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_locrec_practiceSeq_',num2str(subject,'%03d')],'stimlist','videos','nobjects');
elseif ~practice && subject < 99 && exist([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_locrec_seq_',num2str(subject,'%03d'),'.mat'],'file')
    fprintf('Loading sequence for subject %03d.\n',subject)
    load([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_locrec_seq_',num2str(subject,'%03d')]);
else
    % Read in stimulus info
    fprintf('Reading in stimnames using tdfread.\n');
    videos = tdfread('StimList_ByObject_040714.txt');
    videos.MovieID = cellstr(videos.MovieID);
    videos.ObjectIDinMovie = cellstr(videos.ObjectIDinMovie); 
    videos.ObjectNumber = cellstr(videos.ObjectNumber);
    videoID = unique(videos.MovieID);
    if exist((strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'+.dat')),'file')
        altfile = input('ObjEnc+ file exists. Load this file?: (Y=1,N=0)');
        if altfile
            encdata = tdfread(strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'+.dat'),',');
        else
            encdata = tdfread(strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'.dat'),',');
        end %if
    else
        encdata = tdfread(strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'.dat'),',');
    end %if exist
    encdata.ObjectID = cellstr(encdata.ObjectID);
    nobjects=length(encdata.ObjectID);
    
    % Determining randomized object order
    randObjectOrder = randperm(nobjects);

    for iobject=1:nobjects
        stimlist{iobject} = encdata.ObjectID{randObjectOrder(iobject)};
    end %iobject=
    
    % Save out sequence 
    save([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_locrec_',num2str(subject,'%03d')],'encdata','videos','nobjects','stimlist','randObjectOrder','locationscale');
end %if

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
enter = KbName('return');
escapeKey = KbName('ESCAPE');
oneResp = KbName('1!');
twoResp = KbName('2@');
threeResp = KbName('3#');
fourResp = KbName('4$');
fiveResp = KbName('5%');
enablekeys = [oneResp twoResp threeResp fourResp fiveResp enter escapeKey];
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
    
    % Open a double buffered fullscreen window on the stimulation screen 'screenNumber' and use 'gray' for color
    % 'w' is the handle used to direct all drawing commands to that window
    % 'wRect' is a rectangle defining the size of the window. See "help PsychRects" for help on such rectangles
    [w, wRect]=Screen('OpenWindow',screenNumber, screenColor);
    [mx, my] = RectCenter(wRect);
    scene_x = 350;
    scene_y = 350;
    scene_mtx = [mx-scene_x/2,my-scene_y/2,mx+scene_x/2,my+scene_y/2]';
    ytext = (my*2)-100;
    
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
    DrawFormattedText(w, [instruc linebreak locationscale linebreak instruc2 linebreak startscreen], 'center', 'center', [255 255 255]);  % Write instructions
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
    
    trialcounter = 0;
    for iobject=1:nobjects;
%     for iobject=1:4; %use for debugging
        trialcounter = trialcounter + 1;
        
            %Wait til after the trigger to disable the trigger check
            RestrictKeysForKbCheck(enablekeys);

            %Added for KbQueue
            KbQueueCreate(resp_device,keylist);

            %Draw fixation cross
            DrawFormattedText(w, fixation, 'center', 'center', [255 255 255]);  % Write fixation
            [VBLTimestamp expstart]=Screen('Flip', w);  % Update display to show +
            WaitSecs(1.000);    % Display + 

            Screen('Close');    % Close to not overload

            %Initialize KbCheck and return to zero in case a button is pressed
            [KeyIsDown, endrt, KeyCode]=KbCheck;
            KeyIsDown = zeros(size(KeyIsDown));
            endrt = zeros(size(endrt));
            KeyCode = zeros(size(KeyCode));
            WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1

            %load image
            if ~practice
                probe = [stimAloneDir,stimlist{iobject},file_ext]; 
            else
                probe = [practicestimdir,stimlist.objAlone{randObjectOrder(iobject)},file_ext];
            end %if
            probe = imread(probe); 
            showprobe = Screen('MakeTexture',w,probe);
            %present image
            Screen('DrawTexture',w,showprobe,[],image_resize);
            DrawFormattedText(w,[location_q, linebreak locationscale],'center', ytext, [255 255 255]);
            [VBLtimestamp locationQstart]=Screen('Flip',w);

            % initialize KbCheck and variables to avoid time dealys
            KeyIsDown = zeros(size(KeyIsDown));
            endrt = zeros(size(endrt));
            KeyCode = zeros(size(KeyCode));
            %Added for KbQueue
            KbQueueFlush();
            KbQueueStart();
            location_resp = '';
            firstpress = [];

            while (GetSecs - locationQstart) < quesduration 
                if ~KeyIsDown%(enablekeys)  % if key is pressed, stop recording response
                    [KeyIsDown firstpress KeyCode]=KbQueueCheck(resp_device);
                    if KeyIsDown==1
                    end %KeyIsDown
                    WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
                end
                
                location_resp = KbName(find(firstpress,1,'last'));

                if ~isempty(location_resp) % advances once one of the allowable keys is pressed
                    if ismember(location_resp,{'1!','2@','3#','4$'})
                        break;
                    end %if ismember(
                end %~isempty
                
            end %while (GetSecs - memqstart) < quesduration
            
            if isempty(location_resp) %added to deal w when subj doesn't make mem resp
                location_resp = 'xx';
            end %

            endrt = min(firstpress(find(firstpress)));            
            location_rt=round(1000*(endrt-locationQstart)); % compute reaction time

        if practice && iobject==nobjects
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
        
        % Write trial result to file
            if ~practice
                % write header
                if trialcounter==1
                    fprintf(datafilepointer,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
                        'SubjectID','VideoID','ObjectIDinMovie','ObjectID','ExptStart',...
                        'LocationQuesStart','LocationResp','CorrectLocationResp','LocationRT',...
                        'Name1','House1','Name2','House2');
                end %if
                
                % write data
                fprintf('Saving data.\n');
                [path1 fname1 ext1] = fileparts(house1);
                [path2 fname2 ext2] = fileparts(house2);
                
                currObjMask = strcmp(videos.ObjectNumber,stimlist{iobject});
                
                fprintf(datafilepointer,'%03d,%s,%s,%s,%d,%d,%s,%d,%d,%s,%s,%s,%s\n',...
                subject,...
                videos.MovieID{currObjMask},...
                videos.ObjectIDinMovie{currObjMask},...
                videos.ObjectNumber{currObjMask},...
                expstart,...
                locationQstart,...
                char(location_resp),...
                videos.LocationID(currObjMask),...
                location_rt,...
                name1,...
                fname1,...
                name2,...
                fname2);
            end %if ~practice 
%         end %if videos.Lure(randObjectOrder(iobject))==1
    end %iobject=
    
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