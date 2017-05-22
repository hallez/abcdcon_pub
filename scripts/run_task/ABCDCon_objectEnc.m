% ContextABCD study
% encoding phase
% Halle Dimsdale-Zucker

%% Intro stuff (like cleaning up, etc.)

initialize_ABCDCon
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock))); %sets random number generators to a random seed
cd(scriptsTask);

%% Flags
DEBUGGING_FLAG=0;

%% Script-specific variables

fprintf('Setting up variables.\n');
subject=input('Enter subject number: ');    
practice=input('Select 1 for practice, 0 for task: '); 

% Instructions and on-screen text 
btwnvideos = 'Press ENTER to see the next 10 objects.\n';
instruc_location = 'Where did you see this object?\n';
instruc_seman = 'Is this object worth more than $50?\n';
location_scale = '1   2   3   4   5   6   7   8\n';

% Timing
sem_quesduration = 4.000*fast;
location_quesduration = 4.000*fast;

% Variables needed for movie
preloadsecs = 1; % for Screen('OpenMovie')
%abortit = 0; %flag to stop movie playback (just use when debugging)
rate = 1; %set movie playback rate (ie, play at normal speed)


%% Counterbalance info

fprintf('Loading context enc CB info.\n');
load(strcat(rawBehavDir,filesep,'s',num2str(subject,'%03d'),filesep,'ConABCD_contextencCB_s',num2str(subject,'%03d')));

%remember that house1 always goes w name1 and house2 w name2. house1 vs 2
%is random. name1 vs name2 is CB based on subject number (but all of this
%is done in the contextEnc script). 
if mod(subject,2)
    sem_scale = '(1)Yes             (2)No\n';
else
    sem_scale = '(1)No             (2)Yes\n';
end %if

% Instructions text
if ~practice
    instruc = 'You will now see objects in Alex and Jamie''s homes.\n Each video has 10 objects for you to learn.\n\n After the video, you will be asked\n two questions about each of the 10 objects:\n\n (1) Is the object worth more than $50?\n (2) Where in Alex or Jamie''s house was the object?\n';
else
    instruc = ['You will now see objects in Halle''s house\n so you can practice using the response keys.\n\n In the actual task, you will see a video of 10 objects\n in either Alex or Jamie''s house.\n\n For the practice, you will just see still pictures\n of objects in Halle''s house.\n\n You will be asked two questions\n about each of the 10 objects:\n\n (1) Is the object worth more than $50?\n',sem_scale, '\n(2) Where in Halle''s house was the object?\n\n'];
end %if

%% Initialize file where will write data

fprintf('Setting up data file.\n');

%write to .dat files
datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'.dat');  % name of data file to write to but not creating a file

% read existing result file and prevent overwriting files from a previous subject (except for subjects > 99)
if DEBUGGING_FLAG 
    datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'.dat'); %overwrite existing file
    datafilepointer = fopen(datafilename,'wt');
elseif ~practice 
    if subject <99 && fopen(datafilename, 'rt')~=-1 %r = open a file for reading. t = open in text mode. -1 is the value returned when fopen cannot open a file (e.g., if this file doens't exist yet--which is normally what you want bc it means you won't be overwriting anything!)
        tmp = input('Data file already exists! Override? (Y=1) ');
        tmp2 = input('Write over existing file? (Y=3,N=2) ');
        if tmp==1 && tmp2==2
            while fopen(datafilename,'rt')~=-1
                [a b c] = fileparts(datafilename);
                datafilename = strcat(a,filesep,b,'+',c); %don't ever overwrite orig file--if a file for that subject number already exists and you override (ie, say it's okay to reuse that subject number), create the file but append '+' to the end. 
            end
            datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing..w = open for writing, discard existing contents. t = text mode
        elseif tmp==1 && tmp2==3
            %overwrite existing file
            datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objectEnc_s',num2str(subject,'%03d'),'.dat'); %overwrite existing file
            datafilepointer = fopen(datafilename,'wt');
        else
            fclose('all');
            error('Choose a different subject number.');
        end %if tmp==1 && tmp2==2
    else 
        datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
    end %if subject <99 && fopen(datafilename, 'rt')~=-1
end %if

%% Save out sequence

% Load video and object stimuli names
videos = tdfread([scriptsTask 'StimList_ByObject_040714.txt']); 
videos.MovieID = cellstr(videos.MovieID);
videos.ObjectIDinMovie = cellstr(videos.ObjectIDinMovie);
videos.ObjectNumber = cellstr(videos.ObjectNumber);
videoID = unique(videos.MovieID);

if practice
    currMovieMask = strcmp(videos.MovieID,'Practice'); 
    currObjects = videos.ObjectIDinMovie(currMovieMask);
    numObj = length(currObjects);
    randObjectOrder = randperm(6);
    stimlist.objInScene = cell(length(randObjectOrder),1);
    stimlist.objAlone = cell(length(randObjectOrder),1);
    for iobject=1:length(randObjectOrder)
        stimlist.objInScene{iobject} = currObjects{randObjectOrder(iobject)};
        stimlist.objAlone{iobject} = videos.ObjectNumber{strcmp(stimlist.objInScene{iobject},videos.ObjectIDinMovie)};
    end %for iobject=
    save([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objenc_practiceSeq_',num2str(subject,'%03d')],'currMovieMask','currObjects','numObj','randObjectOrder','stimlist');
elseif ~practice && subject < 99 && exist([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objrec_seq_',num2str(subject,'%03d'),'.mat'],'file')
    fprintf('Loading sequence for subject %03d.\n',subject)
    load([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objenc_seq_',num2str(subject,'%03d')]);
else
    % Figure out number of stimuli, etc. 
    numEncObjects = 200;
    numvideos = numEncObjects/10;
    brownVideos = videoID(strncmp('Brown',videoID,5));
    grayVideos = videoID(strncmp('Gray',videoID,4));
    numBrownVideos = length(unique(brownVideos));
    numGrayVideos = length(unique(grayVideos));
    
    % Figure out which videos participant will see 
    brownPermutation = randperm(numBrownVideos,(numvideos/2));
    brownVideos_forEnc = brownVideos(brownPermutation);
    grayPermutation = randperm(numGrayVideos,(numvideos/2));
    grayVideos_forEnc = grayVideos(grayPermutation);
    
    [dir1 f1 ext1] = fileparts(house1);
    [dir2 f2 ext2] = fileparts(house2);

    randVideoOrder = randperm(numvideos/2);

    if mod(randi(100),2)
        firsthouse = f1;
    else
        firsthouse = f2;
    end %if mod(randi(100),2)

    brownCounter = 0;
    grayCounter = 0;
    
    for ivideo=1:numvideos
        % Figure out movie presentation order 
        if mod(ivideo,2) &&  strcmp(firsthouse,'Model_25_Brown_HiddenLayers_BirdsEyeAtEnd') %video is odd (ie, identity should match first house)
            %randomly select from Brown* videos
            brownCounter = brownCounter + 1;
            moviestring = brownVideos_forEnc{randVideoOrder(brownCounter)};
            housecolor = 'Brown';
        elseif mod(ivideo,2) &&  strcmp(firsthouse,'Model_25_Gray_HiddenLayers_BirdsEyeAtEnd') %video is odd (ie, identity should match first house)
            %randomly select from Gray* videos
            grayCounter = grayCounter + 1;
            moviestring = grayVideos_forEnc{randVideoOrder(grayCounter)};
            housecolor = 'Gray';
        elseif ~(mod(ivideo,2)) &&  strcmp(firsthouse,'Model_25_Brown_HiddenLayers_BirdsEyeAtEnd') %video is even (ie, should NOT match the identity of firsthouse)
            %randomly select from Gray* videos
            grayCounter = grayCounter + 1;
            moviestring = grayVideos_forEnc{randVideoOrder(grayCounter)};
            housecolor = 'Gray';
        elseif ~(mod(ivideo,2)) && strcmp(firsthouse,'Model_25_Gray_HiddenLayers_BirdsEyeAtEnd') %video is even (ie, should NOT match the identity of firsthouse)
            %randomly select from Brown* videos
            brownCounter = brownCounter + 1;
            moviestring = brownVideos_forEnc{randVideoOrder(brownCounter)};
            housecolor = 'Brown';
        end %if mod(ivideo,2)
        movienumber = moviestring(regexp(moviestring,'\d'));
        moviename = strcat('Model_25_', housecolor, 'Movie',movienumber,'_HiddenLayers.mp4');

        if strcmp(f1,'Model_25_Brown_HiddenLayers_BirdsEyeAtEnd') && strcmp(housecolor,'Brown')
            ownername = name1; 
        elseif strcmp(f2,'Model_25_Brown_HiddenLayers_BirdsEyeAtEnd') && strcmp(housecolor,'Brown')
            ownername = name2;
        elseif strcmp(f1,'Model_25_Gray_HiddenLayers_BirdsEyeAtEnd') && strcmp(housecolor,'Gray')
            ownername = name1; 
        elseif strcmp(f2,'Model_25_Gray_HiddenLayers_BirdsEyeAtEnd') && strcmp(housecolor,'Gray')
            ownername = name2; 
        end % if strcmp
        
        stimlist(ivideo).moviename = moviename;
        stimlist(ivideo).name = ownername;
        stimlist(ivideo).moviestring = moviestring;
        
    end %ivideo=
    
    % Save out sequence
    save([rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_objenc_seq_',num2str(subject,'%03d')],'stimlist','videos','numvideos','numEncObjects');

end %if practice

%% Psychtoolbox stuff

fprintf('Setting up PTB stuff.\n');

AssertOpenGL; %Check to see if has OpenGL compatibility; will abort if not compatible.

% Dummy calls (prevents delays)
KbCheck; %Check the status of the keyboard. If any key is down (i.e., on/pressed), this will return the value 1. If not, 0.
WaitSecs(waitbtwninstruc); %Have script wait for (n) seconds (basically while matlab loads)
GetSecs; %report time back in a very accurate way

% Keyboard stuff
KbName('UnifyKeyNames'); %swap from PTB's internal naming system to current keybord naming system
%name the keys that subjects will be using
%name the keys that subjects will be using
enter = KbName('return');
escapeKey = KbName('ESCAPE');
oneResp = KbName('1!');
twoResp = KbName('2@');
threeResp = KbName('3#');
fourResp = KbName('4$');
fiveResp = KbName('5%');
sixResp = KbName('6^');
sevenResp = KbName('7&');
eightResp = KbName('8*');
enablekeys = [oneResp twoResp threeResp fourResp fiveResp sixResp sevenResp eightResp enter escapeKey];
keylist=zeros(1,256); %create a list of 256 zeros
keylist(enablekeys)=1; %set keys you interested in to 1

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

%% Initialize the experiment :D

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
    
    Screen('TextSize', w, fontSize);    % Set text size
    
    Priority(MaxPriority(w));   % Set priority for script execution to realtime priority
    
    if strcmp(computer,'PCWIN') == 1
        ShowHideWinTaskbarMex(0);
    end
    
    % Initialize KbCheck and return to zero in case a button is pressed
    [KeyIsDown, endrt, KeyCode]=KbCheck;
    KeyIsDown = zeros(size(KeyIsDown));
    endrt = zeros(size(endrt));
    KeyCode = zeros(size(KeyCode));
    WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1 CHANGE BACK TO 1 IF HAVE PROBLEMS
        
    %% Give instructions
    
    fprintf('Giving instructions.\n');
    %pull the first screen of instructions
    DrawFormattedText(w, [instruc linebreak advancescreen], 'center', 'center', [255 255 255]);  % Write instructions in white (255 255 255) font
    Screen('Flip', w);   % Update display to show instructions
    
   
    % Display instructions until KbCheck detects enter pressed
    while (KeyCode(enter)==0)
        [KeyIsDown, RT, KeyCode]=KbCheck; 
        WaitSecs(waitbtwninstruc);    % Prevents overload  
    end 
   
    %% Loop through blocks of trials
    
    fprintf('Starting the experiment.\n');
    
    %Clear out KBCheck, set things to 0, etc. 
    [KeyIsDown, endrt, KeyCode]=KbCheck;
    KeyIsDown = zeros(size(KeyIsDown));
    endrt = zeros(size(endrt));
    KeyCode = zeros(size(KeyCode));
    WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1 
    
    if practice
        for iobject=1:6
%         for iobject=1 %use for debugging
            if ~strcmp(currObjects{randObjectOrder(iobject)},'PracticeLure01')
                
                RestrictKeysForKbCheck(enablekeys); 

                KbQueueCreate(resp_device,keylist); % Added for KbQueue

                %Fixation cross
                DrawFormattedText(w, fixation, 'center', 'center', [255 255 255]);  % Write fixation
                [VBLTimestamp expstart]=Screen('Flip', w);  % Update display to show +
                WaitSecs(1.000);    % Display + for x seconds

                Screen('Close');    % Close to not overload

                if ~strcmp(videos.ObjectNumber{iobject},'AREDROSE8')
                %SEMANTIC JUDGMENT
                %load image
                probe = [practicestimdir, stimlist.objInScene{iobject}, '.png'];
                probe = imread(probe); 
                showprobe = Screen('MakeTexture',w,probe);
                %present image
                Screen('DrawTexture',w,showprobe,[],image_resize); 
                DrawFormattedText(w,[instruc_seman sem_scale],'center',ytext,[255 255 255]); 
                [VBLtimestamp semQstart]=Screen('Flip',w); 

                % initialize KbCheck and variables to avoid time dealys
                KeyIsDown = zeros(size(KeyIsDown));
                endrt = zeros(size(endrt));
                KeyCode = zeros(size(KeyCode));
                %Added for KbQueue
                KbQueueFlush();
                KbQueueStart();
                sem_resp = '';
                firstpress = [];

                while (GetSecs - semQstart) < sem_quesduration
                    % FOR AN EXAMPLE OF WAITING FOR A SPECIFIC KEY TO BE PRESSED, SEE: https://raw.githubusercontent.com/Psychtoolbox-3/Psychtoolbox-3/beta/Psychtoolbox/PsychDemos/KbQueueDemo.m
                    if ~KeyIsDown%(enablekeys)  % if key is pressed, stop recording response
                        [KeyIsDown firstpress]=KbQueueCheck(resp_device);
                        WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
                    end

                    sem_resp = KbName(find(firstpress,1,'last'));

                    if isempty(sem_resp) %added to deal w when subj doesn't make mem resp
                        sem_resp = 'xx';
                    end %

                    sem_endrt = min(firstpress(find(firstpress)));            
                    sem_rt=round(1000*(sem_endrt-semQstart)); % compute reaction time 
                end %while

                % initialize KbCheck and variables to avoid time dealys
                KeyIsDown = zeros(size(KeyIsDown));
                endrt = zeros(size(endrt));
                KeyCode = zeros(size(KeyCode));
                %Added for KbQueue
                KbQueueFlush();
                KbQueueStart();
                firstpress = [];

                %LOCATION JUDGMENT
                probe = [practicestimdir, 'HRZApartment_practice.png'];
                probe = imread(probe); 
                showprobe = Screen('MakeTexture',w,probe);
                %present image
                Screen('DrawTexture',w,showprobe,[],image_resize); 
                DrawFormattedText(w,instruc_location,'center',ytext,[255 255 255]); 
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

                while (GetSecs - locationQstart) < location_quesduration
                    if ~KeyIsDown%(enablekeys)  % if key is pressed, stop recording response
                        [KeyIsDown firstpress]=KbQueueCheck(resp_device);
                        WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
                    end

                    location_resp = KbName(find(firstpress,1,'last'));

                    if isempty(location_resp) %added to deal w when subj doesn't make mem resp
                        location_resp = 'xx';
                    end %

                    location_endrt = min(firstpress(find(firstpress)));            
                    location_rt=round(1000*(location_endrt-locationQstart)); % compute reaction time 
                end %while
                end %if
            end %if
            
            if iobject==6
                DrawFormattedText(w, 'You have completed the practice.\n','center', 'center', [255 255 255]);
                Screen('Flip',w);
                % Display instructions until KbCheck detects enter pressed
                while (KeyCode(enter)==0)
                    [KeyIsDown, RT, KeyCode]=KbCheck; 
                    WaitSecs(waitbtwninstruc);    % Prevents overload  
                end
            end %if iobject==numObj
            
        end %iobject=
    else
        for ivideo=1:numvideos; 
%          for ivideo=1; %use for debugging
            RestrictKeysForKbCheck(enablekeys);  

            KbQueueCreate(resp_device,keylist); % Added for KbQueue

            %Fixation cross
            DrawFormattedText(w, fixation, 'center', 'center', [255 255 255]);  % Write fixation
            [VBLTimestamp expstart]=Screen('Flip', w);  % Update display to show +
            WaitSecs(1.000);    % Display + for x seconds

            timecounter = GetSecs; %initiate timecounter to keep trials on track

            Screen('Close');    % Close to not overload

            %Loading screen
            DrawFormattedText(w,loadscreen,'center', 'center', [255 255 255]);
            Screen('Flip',w);
            WaitSecs(preloadsecs+.005);

            %Initialize KbCheck and return to zero in case a button is pressed
            [KeyIsDown, endrt, KeyCode]=KbCheck;
            KeyIsDown = zeros(size(KeyIsDown));
            endrt = zeros(size(endrt));
            KeyCode = zeros(size(KeyCode));
            WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1

            %Load movie into memory and onto Screen
            %IF HAVE ISSUES SETTING THIS UP, TYPE help GStreamer (ERROR
            %MESSAGE WILL TELL YOU TO DO THIS AS WELL). FOLLOW INSTRUCTIONS
            %TO DOWNLOAD A VERSION OF GSTREAMER (As of 5/21/14 old versions
            %are located here:
            %http://gstreamer.freedesktop.org/src/gstreamer/). Install and
            %unpack the tarball (tar xzvf <.tar.gz file> as per
            %instructions here:
            %http://www.haskell.org/haskellwiki/How_to_unpack_a_tar_file_in_Windows).
            %If get an error re needing a C compiler, follow instructions
            %here:
            %http://stackoverflow.com/questions/6394755/how-to-install-gcc-on-windows-7-machine
            %including going here
            %http://sourceforge.net/projects/mingw-w64/files/ and
            %downloading MinGw. May need to set mingw/bin/gcc.exe onto
            %path. Can do by typing env in the control panel search bar and
            %going to the "Edit environment variables" section. Paste the
            %path to the gcc.exe file into the PATH variable and
            %restart/logoff that user for gcc to be seen on the path. 
            [movie movieduration fps imgw imgh] = Screen('OpenMovie', w, [videosDir stimlist(ivideo).moviename], [], preloadsecs);
            Screen('PlayMovie',movie,rate,0); % plays movie at rate 1 (ie, play at normal speed and play just once)

            % Show basic info about movie:
            fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', stimlist(ivideo).moviename, movieduration, fps, imgw, imgh);

            % Playback loop:
            tex = 0;
            while (KeyCode(enter)==0)

                Screen('TextSize', w, fontSizeBig);    % Set text size to be larger for house labels
    %                 while (abortit<2) %nice for debugging b/c can hit esc to get out of movie, but elim for actual expt
                    while (1) 
                        if (abs(rate)>0) 
                            % Return next frame in movie, in sync with current playback
                            % time and sound.
                            % tex either the texture handle or zero if no new frame is
                            % ready yet. pts = Presentation timestamp in seconds.
                            [tex pts] = Screen('GetMovieImage', w, movie, 1);

                            if tex<0
                            % No. This means that the end of this movie is reached.
                            % We exit the loop and prepare next movie:
                            break;
                            end; %if tex<0

                            if (tex>0)
                                % Yes. Draw the new texture immediately to screen:
                                %dstRect = CenterRect(ScaleRect(Screen('Rect',
                                %myMovieImageTexture), 3,3 ) , Screen('Rect',
                                %myWindow)); %original version
                                %dstRect = -2240       -1108        3520
                                %2132 when use 3,3
                                %dstRect =  160         242        1120         782 use 0.5,0.5
                                dstRect = CenterRect(ScaleRect(Screen('Rect', tex), .6,.6 ) , Screen('Rect',w));
                                Screen('DrawTexture', w, tex, [], dstRect);
                                DrawFormattedText(w,[stimlist(ivideo).name possession],'center',100,[255 255 255]);

                                % Update display:
                                Screen('Flip', w);

                                % Release texture:
                                Screen('Close', tex);
                            end; %if (tex>0)
                        end %if (abs(rate)>0)
                    end %while (1)


                    Screen('CloseMovie',movie);
    %                 Check for abortion by user:
    %                         abortit=0;
    %                         [keyIsDown,secs,keyCode]=KbCheck; %#ok<ASGLU>
    %                         if (keyIsDown==1 && keyCode(esc)) 
    %                             % Set the abort-demo flag.
    %                             abortit=2;
    %                             break;
    %                         end; %if
    %                 end %while (abortit<2)

                Screen('TextSize', w, fontSize);    % Reset text size back to normal

                %Pull objects from current movie, randomize their order
                currMovieMask = strcmp(videos.MovieID,stimlist(ivideo).moviestring);
                currObjects = videos.ObjectIDinMovie(currMovieMask);
                numObj = length(currObjects);
                randObjectOrder = randperm(10); 

                for iobject=1:numObj
                    %SEMANTIC JUDGMENT
                    %load image
                    probe = [stimSceneDir, currObjects{randObjectOrder(iobject)}, '.png'];
                    probe = imread(probe); 
                    showprobe = Screen('MakeTexture',w,probe);
                    %present image
                    Screen('DrawTexture',w,showprobe,[],dstRect); 
                    DrawFormattedText(w,[instruc_seman sem_scale],'center',ytext,[255 255 255]); 
                    [VBLtimestamp semQstart]=Screen('Flip',w); 

                    % initialize KbCheck and variables to avoid time dealys
                    KeyIsDown = zeros(size(KeyIsDown));
                    endrt = zeros(size(endrt));
                    KeyCode = zeros(size(KeyCode));
                    %Added for KbQueue
                    KbQueueFlush();
                    KbQueueStart();
                    sem_resp = '';
    %                 sem_endrt = semQstart; 
                    firstpress = [];

                    while (GetSecs - semQstart) < sem_quesduration
                        if ~KeyIsDown  % if key is pressed, stop recording response
                            [KeyIsDown firstpress]=KbQueueCheck(resp_device);
                            WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
                        end

                        sem_resp = KbName(find(firstpress,1,'last'));

                        if isempty(sem_resp) %added to deal w when subj doesn't make mem resp
                            sem_resp = 'xx';
                        end %

                        sem_endrt = min(firstpress(find(firstpress)));            
                        sem_rt=round(1000*(sem_endrt-semQstart)); % compute reaction time 
                    end %while

                    % initialize KbCheck and variables to avoid time dealys
                    KeyIsDown = zeros(size(KeyIsDown));
                    endrt = zeros(size(endrt));
                    KeyCode = zeros(size(KeyCode));
                    %Added for KbQueue
                    KbQueueFlush();
                    KbQueueStart();
                    firstpress = [];

                    %LOCATION JUDGMENT
                    %load image
                    if strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie1_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie2_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie3_HiddenLayers.mp4')||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie4_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie5_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie6_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie7_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie8_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie9_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie10_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie11_HiddenLayers.mp4') ||...
                            strcmp(stimlist(ivideo).moviename,'Model_25_BrownMovie12_HiddenLayers.mp4')
                        probe = [stimSceneDir, 'Model_25_Brown_HiddenLayers_Quadrants0001.png'];
                    else
                        probe = [stimSceneDir, 'Model_25_Gray_HiddenLayers_Quadrants0001.png'];
                    end %if

                    probe = imread(probe); 
                    showprobe = Screen('MakeTexture',w,probe);
                    %present image
                    Screen('DrawTexture',w,showprobe,[],dstRect); 
                    DrawFormattedText(w,instruc_location,'center',ytext,[255 255 255]); 
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

                    while (GetSecs - locationQstart) < location_quesduration
                        if ~KeyIsDown  % if key is pressed, stop recording response
                            [KeyIsDown firstpress]=KbQueueCheck(resp_device);
                            WaitSecs(waitbtwninstruc);    % wait 1 ms before checking the keyboard again to prevent overload
                        end

                        location_resp = KbName(find(firstpress,1,'last'));

                        if isempty(location_resp) %added to deal w when subj doesn't make mem resp
                            location_resp = 'xx';
                        end %

                        location_endrt = min(firstpress(find(firstpress)));            
                        location_rt=round(1000*(location_endrt-locationQstart)); % compute reaction time 
                    end %while


                     if ~practice 
                        %Write trial result to file: 
                        %first, create a header line for the output file
                        if ivideo==1 && iobject==1 
                            fprintf(datafilepointer,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
                                'SubjectID','VideoID','ObjectIDinMovie','ObjectID','ExptStart','SemQuesStart','SemResp','SemRT',...
                                'LocationQuesStart','LocationResp','CorrectLocation','LocationRT');
                        end %if
                        %now write data to .dat file
                        fprintf('Saving data.\n');
                        tmp_objID2 = videos.ObjectNumber(strcmp(videos.ObjectIDinMovie,currObjects{randObjectOrder(iobject)}));
                        fprintf(datafilepointer,'%03d,%s,%s,%s,%d,%d,%s,%d,%d,%s,%d,%d\n',...
                            subject,...
                            stimlist(ivideo).moviestring,...
                            currObjects{randObjectOrder(iobject)},...
                            tmp_objID2{1},...
                            expstart,...
                            semQstart,...
                            char(sem_resp),...
                            sem_rt,...
                            locationQstart,...
                            char(location_resp),...
                            videos.LocationID(strcmp(videos.ObjectIDinMovie,currObjects{randObjectOrder(iobject)})),...
                            location_rt)
                     end %if ~practice
                end %iobject=

                if ivideo==numvideos
                    DrawFormattedText(w,endscreen,'center','center',[255 255 255]);
                else
                    DrawFormattedText(w,btwnvideos,'center','center',[255 255 255]);
                end %if

                Screen('Flip',w);
                % Display instructions until KbCheck detects enter pressed
                while (KeyCode(enter)==0)
                    [KeyIsDown, RT, KeyCode]=KbCheck; 
                    WaitSecs(waitbtwninstruc);    % Prevents overload  
                end
            end %while (KeyCode(enter)==0)

            %Initialize KbCheck and return to zero in case a button is pressed
            [KeyIsDown, endrt, KeyCode]=KbCheck;
            KeyIsDown = zeros(size(KeyIsDown));
            endrt = zeros(size(endrt));
            KeyCode = zeros(size(KeyCode));
            WaitSecs(waitbtwninstruc); %add this to clear out responses so KbCheck will be 0 rather than 1
            KbQueueFlush(); %Added for KbQueue


        end %ivideo=
    end %if practice
        
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

