% ContextABCD study
% Context encoding phase
% Halle Dimsadale-Zucker

%% Intro stuff 

initialize_ABCDCon
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock))); %sets random number generators to a random seed
cd(scriptsTask);

%% Flags
DEBUGGING_FLAG=0;

%% Script-specific variables 

fprintf('Setting up variables.\n');
subject=input('Enter subject number: '); 

% Experiment-specific variables
numhouses = 2;

% Instructions and on-screen text 
instruc = 'You will now see Alex and Jamie''s homes.\n You will be asked to draw a map of each house\n so pay close attention to the houses''\n layouts and landmarks.\n';
instruc_view2_p1 = 'You will now see';
instruc_view2_p2 ='again.\n';

%% Counterbalance info

fprintf('Figuring out counterbalancing.\n');

if mod(subject,2) %returns 1 when subject number is odd 
    name1 = 'Alex';
    name2 = 'Jamie';
else
    name1 = 'Jamie';
    name2 = 'Alex';
end

%randomize house1 vs house2 
if mod(randi(100),2)
    house1 = [videosDir 'Model_25_Brown_HiddenLayers_BirdsEyeAtEnd.mp4'];
    house2 = [videosDir 'Model_25_Gray_HiddenLayers_BirdsEyeAtEnd.mp4'];
else
    house1 = [videosDir 'Model_25_Gray_HiddenLayers_BirdsEyeAtEnd.mp4'];
    house2 = [videosDir 'Model_25_Brown_HiddenLayers_BirdsEyeAtEnd.mp4'];
end %if 

[dir1 file1 ext1] = fileparts(house1);
[dir2 file2 ext2] = fileparts(house2);

%% Initialize file where will write data

fprintf('Setting up data file.\n');

if ~exist(strcat(rawBehavDir,'s',num2str(subject,'%03d')),'dir')
    mkdir(strcat(rawBehavDir,'s',num2str(subject,'%03d')));
end

%write to .dat files
datafilename = strcat(rawBehavDir, 's',num2str(subject,'%03d'),filesep,'ConABCD_contextencCB_s',num2str(subject,'%03d'),'.dat');  % name of data file to write to but not creating a file

% read existing result file and prevent overwriting files from a previous subject (except for subjects > 99)
if DEBUGGING_FLAG 
    datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_contextencCB_s',num2str(subject,'%03d'),'.dat'); %overwrite existing file
    datafilepointer = fopen(datafilename,'wt');
else
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
            datafilename = strcat(rawBehavDir,'s',num2str(subject,'%03d'),filesep,'ConABCD_contextencCB_s',num2str(subject,'%03d'),'.dat'); %overwrite existing file
            datafilepointer = fopen(datafilename,'wt');
        else
            fclose('all');
            error('Choose a different subject number.');
        end
    else
        datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
    end %if subjects < 99 & fopen....
end %if

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
enter = KbName('return');
escapeKey = KbName('ESCAPE');
oneResp = KbName('1!');
twoResp = KbName('2@');
threeResp = KbName('3#');
fourResp = KbName('4$');
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
    WaitSecs(1); %add this to clear out responses so KbCheck will be 0 rather than 1 CHANGE BACK TO 1 IF HAVE PROBLEMS
        
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
    
    %start with a fixation cross
    [KeyIsDown, endrt, KeyCode]=KbCheck;
    KeyIsDown = zeros(size(KeyIsDown));
    endrt = zeros(size(endrt));
    KeyCode = zeros(size(KeyCode));
    WaitSecs(1); %add this to clear out responses so KbCheck will be 0 rather than 1 CHANGE BACK TO 1 IF HAVE PROBLEMS
    
    DrawFormattedText(w, fixation, 'center', 'center', [255 255 255]);  % Make fixation screen 
    Screen('Flip', w);   % Update display to show fix cross
    
    for irep=1:2;
       for ihouse=1:numhouses; 
    %     for ihouse=1; %use for piloting

            [KeyIsDown, endrt, KeyCode] = KbCheck; % basically zeroing out KeyIsDown, etc. so that will accurately collect and time subj's responses 
            WaitSecs(1); % Prevents overload

            RestrictKeysForKbCheck(enablekeys); % Wait til after the trigger to disable the trigger check

            KbQueueCreate(resp_device,keylist); % Added for KbQueue

            %Fixation cross
            DrawFormattedText(w, fixation, 'center', 'center', [255 255 255]);  % Write fixation
            [VBLTimestamp expstart]=Screen('Flip', w);  % Update display to show +
            WaitSecs(1);    % Display + for 4 seconds CHANGE BACK TO 1 IF HAVE PROBLEMS

            preloadsecs = 1; % for Screen('OpenMovie')

            Screen('Close');    % Close to not overload

            %Loading screen
            DrawFormattedText(w,loadscreen,'center', 'center', [255 255 255]);
            Screen('Flip',w);
            WaitSecs(2.000);

    %         abortit = 0; %flag to stop movie playback (just use when debugging)
            rate = 1; %set movie playback rate (ie, play at normal speed)

            if ihouse==1
                moviename = house1;
            else
                moviename = house2;
            end %if


            [movie movieduration fps imgw imgh] = Screen('OpenMovie', w, moviename, [], preloadsecs);

            Screen('PlayMovie',movie,rate,0); % plays movie at rate 1 (ie, play at normal speed and play just once)

            % Show basic info about movie:
            fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);

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
                                if ihouse==1
                                    DrawFormattedText(w,[name1 possession],'center',100,[255 255 255]);
                                elseif ihouse==2
                                    DrawFormattedText(w,[name2 possession],'center',100,[255 255 255]);
                                end %if

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

                %Advance to second viewing
                if ihouse==1 && irep==1
                    DrawFormattedText(w,['Press ENTER to see ' name2 possession],'center','center',[255 255 255]);
                elseif ihouse==2 && irep==1
                    DrawFormattedText(w, [instruc_view2_p1 [' '] name1 possession [' '] instruc_view2_p2 advancescreen], 'center', 'center', [255 255 255]);
                elseif ihouse==1 && irep==2
                    DrawFormattedText(w, [instruc_view2_p1 [' '] name2 possession [' '] instruc_view2_p2 advancescreen], 'center', 'center', [255 255 255]);
                elseif ihouse==2 && irep==2
                    DrawFormattedText(w,endscreen,'center','center',[255 255 255]);
                end %if
                Screen('Flip',w);
                % Display instructions until KbCheck detects enter pressed
                while (KeyCode(enter)==0)
                    [KeyIsDown, RT, KeyCode]=KbCheck; 
                    WaitSecs(waitbtwninstruc);    % Prevents overload  
                end
    %             
            end %while (KeyCode(enter)==0)

            %Initialize KbCheck and return to zero in case a button is pressed
            [KeyIsDown, endrt, KeyCode]=KbCheck;
            KeyIsDown = zeros(size(KeyIsDown));
            endrt = zeros(size(endrt));
            KeyCode = zeros(size(KeyCode));
            WaitSecs(1); %add this to clear out responses so KbCheck will be 0 rather than 1

            % Write trial result to file:
            % Create headerline
            [aa bb cc] = fileparts(datafilename);
            if ihouse==1 && irep==1
                fprintf(datafilepointer, '%s,%s,%s,%s,%s,%s,%s\n', ...
                    'SubjectID','name1','name2','house1','house2','house1fpath','house2fpath');
            end %if
    %             
            fprintf('Saving data');
            fprintf(datafilepointer,'%03d,%s,%s,%s,%s,%s,%s\n',...
                subject,...
                name1,...
                name2,...
                file1,...
                file2,...
                house1,...
                house2);
            save(strcat(aa,'/',bb),'name1','name2','house1','house2'); %saves to .mat file

            KbQueueFlush(); %Added for KbQueue   

       end %ihouse=
    end %irep=
    
    %% Finish up %%
    
    % Cleanup at end of experiment
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    ShowHideWinTaskbarMex(1)
    fprintf('Name1: %s, House1: %s\n',name1,house1)
    
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
    
    fprintf('Name1: %s, House1: %s\n',name1,house1)
    psychrethrow(psychlasterror);   % Output the error message that describes the error
end
