function [Subjects, StudyFolder] = ABL_OpenSim_Setup_Batch(path)
% Batch Processing Script for restructuring ABL .anc and .trc files for OpenSim Compatability

% INPUT: .anc and .trc files
% OUTPUT: .trc and .mot

% Created by: Billy Clark & Ricky Pimentel
% Applied Biomechanics Laboratory - UNC Chapel Hill
% 2020

%% Setup Directories
clear; close all; clc; warning off; dbstop if error;
CurrPath = dir; % get current folder and add to path
addpath(genpath(CurrPath(1).folder));
addpath(genpath('C:\Users\richa\Documents\Packages\BOS'));
addpath(genpath('C:\Users\richa\Documents\Packages\MoCapTools'));
% path = 'D:\UNC_ABL\FpMetabolics_2020\SubjData';
path = 'D:\UNC_ABL\AgeFatigue';
% path = 'C:\Users\richa\Documents\Bioengineering\Projects\GT_Modeling\BatchProcessing';

%% Modular Settings
Settings.ConvertForces = 'No'; % Convert .anb forces to .mot files?
Settings.PlotVoltage = 'No'; % GRF/Forces plotting settings
Settings.PlotGRFs = 'No';
Settings.PlotMoments = 'No';

% convert TRC marker trajectories
Settings.ConvertTrajectories = 'No'; % convert marker trajectory files (.TRCs)

Settings.CalculateHJC = 'Yes'; % calculate hip joint center from functional trials
Settings.HipJointCenter = 'Yes'; % use static files with hip joint centers
Settings.PrintResults = 'No'; % print results in command window

Settings.MassFromStatic = 'InfoSheet'; %'InfoSheet';
Settings.IdentifyCrossover = 'Yes'; % identify crossover steps when walking on treadmill

%% Set location of directories
% Set subject location
if exist('path','var')
    StudyFolder = path;
else
    StudyFolder = uigetdir(CurrPath(1).folder, 'Select Folder Containing Subject Data');
end
SubjPath = dir(StudyFolder);
SubjPath(1:2) = []; 
addpath(genpath(StudyFolder));

% identify subject and trials within folder
NumSubjFolders = sum([SubjPath.isdir]==1); % identify number of subjects
SubjFolders = [SubjPath.isdir];

% Loop Through Subjects% initialize subject structure
Subjects(NumSubjFolders).name = [];
Subjects(NumSubjFolders).Folders = [];
Subjects(NumSubjFolders).Demo = [];
Subjects(NumSubjFolders).Trials = [];
SC = 1;

% locate cal file
FPcal_file = 'C:\Users\richa\Documents\Packages\BOS\Functions\FP_Calibrations\forcepla.cal';
if ~ exist(FPcal_file, 'file') % if default file doesnt exist, select a new file
    [FPcal_file, calfilepath] = uigetfile('.cal', 'Select a force plate calibration file');
    addpath(genpath(calfilepath));
    clearvars calfilepath
end

%% Loop through each subject
for SubjLoop = find(SubjFolders==1) % first two items are garbage
    
    tic; % start processing timer for each subject
    
    % Define subject filename, folder, path, and directory
    Subjects(SC).name = SubjPath(SubjLoop).name;
    Subjects(SC).Folders.folder = SubjPath(SubjLoop).folder;
    Subjects(SC).Folders.path = strcat(Subjects(SC).Folders.folder, '\', Subjects(SC).name);
    Subjects(SC).Folders.dir = dir(Subjects(SC).Folders.path);
    
    % display which subject is being processed
    disp(' ');
    disp(strcat('PROCESSING__', Subjects(SC).name));
    disp(' ');
    
    % if non-existent, create OpenSim Folder to store data
    Subjects(SC).Folders.OpenSimFolder = strcat(Subjects(SC).Folders.path, '\OpenSim');
    if ~exist(Subjects(SC).Folders.OpenSimFolder, 'dir')
        mkdir(Subjects(SC).Folders.OpenSimFolder);
    end
    
    %% Set up Loop through all Trials of Subject
    % identify of number of trials
    NumTrials = sum(contains({Subjects(SC).Folders.dir.name}, '.trc'));
    Trials(NumTrials).name = [];   % intialize Trial sub-structure
    Trials(NumTrials).subject = [];
    Trials(NumTrials).folder = [];
    
    TC = 1;  % trial loop counter
    for j = find(contains({Subjects(SC).Folders.dir.name}, '.trc')) % identify TRC files to loop through
        
        Trials(TC).name = Subjects(SC).Folders.dir(j).name(1:end-4); % save trial name
        TrialName = Subjects(SC).Folders.dir(j).name(1:end-4); % save trial name as variable
        Trials(TC).subject = Subjects(SC).name; % save subject name
        Trials(TC).folder = Subjects(SC).Folders.path; % save folder name
        
        % if static trial - label accordingly
        if contains(Trials(TC).name,'Static') || ...
                contains(Trials(TC).name,'static')
            StaticFile = Trials(TC).name;
            Trials(TC).type = 'static';
        elseif contains(Trials(TC).name,'hjc') || ...
                contains(Trials(TC).name,'HJC')
            Trials(TC).type = 'hjc';
        else
            Trials(TC).type = 'walking';
        end
        
        % define specific filenames 
        Trials(TC).files.OpenSimTRC = strcat( TrialName, '_OpenSim.trc');
        Trials(TC).files.OpenSimGRF = strcat( TrialName, '_OpenSimGRF.mot');
        Trials(TC).files.TRC = strcat( TrialName, '.trc');
        Trials(TC).files.ANC = strcat( TrialName, '.anc');
        % if not copying over to opensim files, rename as original files
        if exist(Trials(TC).files.OpenSimTRC, 'file') == 0
            if strcmp(Settings.ConvertTrajectories, 'No')
                Trials(TC).files.OpenSimTRC = Trials(TC).files.TRC;
            end
        end
        if exist(Trials(TC).files.OpenSimGRF, 'file') == 0
            if strcmp(Settings.ConvertForces, 'No')
                if exist(Trials(TC).files.ANC, 'file')
                    Trials(TC).files.OpenSimGRF = Trials(TC).files.ANC;
                else
                    Trials(TC).files.OpenSimGRF = [Trials(TC).files.TRC(1:end-4) '.mot'];
                end
            end
        end
        
        TC = TC + 1; % on to next trial (next TRC file)
    end
    
    % delete empty trials
    ToDel = arrayfun(@(Trials) isempty(Trials.name),Trials);
    Trials(ToDel) = [];
    
    % identify static & walking trials
    WalkingTrials = contains({Trials.type}, 'walking');
    StaticTrial = contains({Trials.type}, 'static');
    
    if sum(StaticTrial) > 1 % if too many static trials - user selects one
        ps = 'Multiple Static Trials Detected, Select the one to use to scale the opensim model (no exo)'; 
        StaticTrials = {Trials(StaticTrial).name};
        s = listdlg('ListString', StaticTrials, 'PromptString', ps);  
        staticname = StaticTrials{s};
        StaticTrial = contains({Trials.name}, staticname);
    end
    clearvars TrialCounter
    
    %% Calculate Hip Joint Centers if desired
    if strcmp(Settings.CalculateHJC, 'Yes')
        %         for i = 1:length(Trials)
        %             filesAddHJC{i} = [Trials(i).name '.trc'];
        %         end
        
        % add HJCs to static trial
        Osim.findHJC(strcat(StaticFile, '.trc'), Trials(StaticTrial).folder)
    end
    close all;
    
    %% Convert Force Files (MOT)
    if strcmp(Settings.ConvertForces, 'Yes')
        disp(' '); 
        disp('Converting Force Files'); 
        
        % Pre-load FPcal_file to speed up process
        [FPcal_data.S, FPcal_data.pos, FPcal_data.origin, FPcal_data.R] = load_fpcal(FPcal_file);
        
        for i = 1:length(Trials)
            cd(Trials(i).folder); % change to folder with trials
            convertFPdata_OpenSim(Trials(i).files.ANC, FPcal_data, Trials(i).folder, Settings);
        end
        
        disp('All .anc files converted to .mot');
    end
    clearvars AbovePath ToCut FileDir Files2move i j newFile
    cd(CurrPath(1).folder);
    
    
    %% Convert Marker Trajectories
    if strcmp(Settings.ConvertTrajectories, 'Yes')
        disp(' '); 
        disp('Converting Marker Files'); 
        
        %% Make any changes to TRC File
        for i = 1:length(Trials) % loop through files to convert
            
            disp(['Converting   ' Trials(i).name '   TRC data']);
            
            % Copy the TRC files to a OpenSim folder to prevent overwriting
            inFile = strcat(Trials(i).folder,'\', Trials(i).files.TRC);
            cd(Subjects(SC).Folders.path); 
            addpath(genpath(Subjects(SC).Folders.OpenSimFolder));  % add output to path
            trc = Osim.readTRC(inFile);
            TRC = table2array(trc);
            
            %% check for gaps or undefined data in original marker file
            [~, NumColumn] = size(TRC); % get matrix dimensions
            
            % save times of trial
            Trials(i).Times.TRC = TRC(:,1);
            TRCnew = zeros(size(TRC));
            
            %% Rotate Coordinate System
            ROTATE = 1;
            if ROTATE == 1 % same as method above, flip Z and Y
                TRCnew(:,1) = TRC(:,1);
                for n = 2:3:NumColumn
                    TRCnew(:, n) = TRC(:, n); %X
                    TRCnew(:, n+1) = TRC(:, n+2); % old Z to new Y
                    TRCnew(:, n+2) = -TRC(:, n+1); % negative old Y to new Z
                end
            else
                TRCnew = TRC;
            end
            clearvars n ii ij
            
            %% Write data to new TRC file
            Headers = trc.Properties.VariableNames;
            Headers{1} = 'Header';
            trcNew = array2table(TRCnew, 'VariableNames', Headers);
            NewFile = strcat(Subjects(SC).Folders.OpenSimFolder, '\', Trials(i).files.OpenSimTRC);
            Osim.writeTRC(trcNew, 'FilePath', NewFile);
            
        end
    end
    clearvars e eF eS eW file i ans Angle
    
    
    %% Get subject mass from static trial
    if strcmp(Settings.MassFromStatic, 'Yes')
        i = StaticTrial;
        disp(' ');
        disp('Getting subject mass from static trial');
        Trials(i).GRF = Osim.readMOT(strcat(Trials(i).folder, '\OpenSim\', Trials(i).files.OpenSimGRF));
%         Trials(i).GRF = LoadGRF(strcat(Trials(i).folder, '\OpenSim\', Trials(i).files.OpenSimGRF), 0); % load grfs
        Cols = {'ground_force_vz','x1_ground_force_vz'}; % define vGRF columns
%         Cols = {'ground_force_vy','1_ground_force_vy'}; % define vGRF columns
        vGRFs(:,1) = Trials(i).GRF.(Cols{1});
        vGRFs(:,2) = Trials(i).GRF.(Cols{2});
        AvgvGRF = mean(sum(vGRFs,2)); % sum and average
        Subjects(SC).Demo.mass = AvgvGRF ./ 9.81; % calculate mass from newtons
        clearvars Cols vGRFs AvgvGRF
        disp(' ');
    elseif strcmp(Settings.MassFromStatic, 'InfoSheet')
        fn = [Subjects(SC).name '_Info.csv']; 
        T = readtable(fn);
        r = contains(T{:,1}, 'Mass');
        Subjects(SC).Demo.mass = str2double(char(T{r, 2}));
    end
    
    %% Identify Crossover steps
    % Crossover analysis
    
    if strcmp(Settings.IdentifyCrossover, 'Yes')
        disp(' ');
        disp('Identifying Crossover Steps');
        
        % run crossover analysis
        [Trials] = CrossoverAnalysis(Trials, WalkingTrials, Subjects(SC).Demo.mass);
        close; % close figure
        
    end
    
    %% Finalize iteration for subject loop
    % remove some variables to save storage space
%     Fields2Del = {'GRF'};
%     if isfield(Trials, Fields2Del)
%         Trials = rmfield(Trials, Fields2Del);
%     end
    
    % add Trials struct to Subjects struct
    Subjects(SC).Trials = Trials;
    
    % end of loop cleaning and message
    T = toc;
    disp(' ');
    disp(strcat([Subjects(SC).name, ' processed in:', ' ', num2str(round(T)), ' seconds']));
    disp(' ');
    disp(' ');
    clearvars TrialCounter Trials FData TRCnew TRC col ans d e eF eS eW filename2 i ii ij ind j k m n Str ...
        Ind Lind Linds Lval Margin Rind Rinds Rval ROTATE Src Dst val StartMetrics Gaps CrossoverFig ...
        AllLefts AllRights AllSteps ColumnHeaderStart LStrikes RStrikes TrialTimes Window
    
    % iterate to next subject
    SC = SC + 1;
    
end

%% save Subjects data after processing loop
save(strcat(StudyFolder, '\Subjects.mat'), 'Subjects');

disp(' ');
disp('ALL SUBJECTS PROCESSED');

end
