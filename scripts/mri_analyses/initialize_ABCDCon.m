% Initialize script for ABCDCon directories
% Author: Halle R. Dimsadle-Zucker

fprintf('\nRunning initialize_ABCDCon.\n')

% clean up
fprintf('Cleaning up before we begin.\n')
clear all;
close all;
fclose('all');
clc;

% initialize directories - assumes all paths are relative to a project_dir
% which is ../../ relative to this script
fprintf('Initializing directories.\n')
project_dir = fullfile(pwd,'..','..');
name = getComputerName(); %//TODO use getComputerName that's in vendor 
addpath(genpath(fullfile(project_dir,'vendor','yamlmatlab','0.4.3')));
config = ReadYaml(fullfile(project_dir,'config.yml'));
addpath(genpath(fullfile(project_dir,'vendor')));

rawMRIDir = [fullfile(project_dir,config.directories.raw_mri) filesep];
analMRIDir = [fullfile(project_dir,config.directories.analyzed_mri) filesep];

% stimulus directories
practicestimdir = [fullfile(project_dir,config.directories.stimuli.practice) filesep];
stimSceneDir = [fullfile(project_dir,config.directories.stimuli.object_in_scene) filesep];
stimAloneDir = [fullfile(project_dir,config.directories.stimuli.object_alone) filesep];
videosDir = [fullfile(project_dir,config.directories.stimuli.house_video) filesep];
houseStillsDir = [fullfile(project_dir,config.directories.stimuli.house_still) filesep];

% script directories
scriptsMRIAnal = [fullfile(project_dir,config.directories.mri_analysis_scripts) filesep];
scriptsBehavAnal = [fullfile(project_dir,config.directories.behavioral_analysis_scripts) filesep];
scriptsTask = [fullfile(project_dir,config.directories.run_experiment_scripts) filesep];

% data directories
rawBehavDir = [fullfile(project_dir,config.directories.raw_behavioral) filesep];
analBehavDir = [fullfile(project_dir,config.directories.analyzed_behavioral) filesep];

%% Variables for running task scripts
% instructions and text screens
advancescreen = 'Please hit ENTER to continue.\n';
endscreen = 'You are finished with this part.\n\n Please notify the experimenter.\n';
file_ext = '.png';
fixation = '+';
loadscreen = 'Please wait. The next video is loading.';
possession = '''s house';

% font and visual settings
fontSize = 25;
fontSizeBig = 40;
image_resize = [64 188 1216 836];
linebreak = '\n';
screenColor = [50 50 50];
ytext = 750;

% timing
fast = 1.000; %set fast to a proportion to speed up testing, otherwise, set to 1.000
waitbtwninstruc = 0.001*fast; %CHANGE BACK TO 1 IF HAVE PROBLEMS

%% Variables for analysis scripts

subjects = [1:2,6:14,18:21,23:30];  

runbase = 'run';
runnums = 1:4;
nruns = length(runnums);
b.runs = cell(1,nruns);
for irun=1:nruns
    b.runs(irun) = cellstr([runbase,num2str(runnums(irun))]);
end %for irun=

numRecogTrials = 252;

rawLocalizerString = '001_localizer_32ch';
convertedLocalizerString = 'localizer_32ch_0001';

nii = '.nii';
mat = '.mat';

fprintf('Done initializing. Returning to current script.\n\n')
