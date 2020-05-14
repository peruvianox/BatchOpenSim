function [Subjects, StudyFolder] = ABL_OpenSim_Setup_Batch(path)
% Batch Processing Script for restructuring ABL .anc and .trc files for OpenSim Compatability
% INPUT: .anc and .trc
% OUTPUT: .trc and .mot
% Notes: must be accompanied by FPcal_file, covertFPdata_OpenSim function

% Created by: Billy Clark & Ricky Pimentel 
% Applied Biomechanics Laboratory - UNC Chapel Hill
% 2020

%% Setup Directories
clear; close all; clc;
CurrPath = dir; % get current folder and add to path
addpath(genpath(CurrPath(1).folder));
warning off
dbstop if error

%% Modular Settings
Settings.ConvertForces = 'No'; % Convert .anb forces to .mot or .forces files?
Settings.ConvertForceType = 'OpenSim'; % other options: 'Forces', 'Both';
Settings.PlotVoltage = 'No'; % GRF/Forces plotting settings
Settings.PlotGRFs = 'No'; 
Settings.PlotMoments = 'No'; 
Settings.CompareGRFandForces = 'No'; % compare converted force files (forces vs opensim)

% convert TRC marker trajectories?
Settings.ConvertTrajectories = 'Yes'; % convert marker trajectory files (.TRCs)
Settings.ConvertTrajectoryType = 'OpenSim'; % no other type currently

% Crossover analysis
Settings.IdentifyCrossover = 'Yes'; % treadmill walking - identify crossover steps
Settings.ExcludeCrossover = 'Yes'; % exclude crossover steps from analysis, or make sure they are not in the time range for analysis
Settings.SaveCrossoverPlot = 'Yes'; % plot crossover step analysis figures

% Time Parsing setup
% number of windows to extract from trial (ie best biofeedback, best ultrasound tracking)
% if 0, do not apply stride times and just pull first strides from whole trial
Settings.TimeParse.NumWindows = 1;
Settings.TimeParse.NumStrides = 2; % number of strides to extract from window (1,2,3,...)

Settings.AddTorso = 'Yes'; % add torso points to only lowerbody walking data
Settings.PlotTorso = 'No'; % plot generalized torso curve

Settings.HipJointCenter = 'Yes'; % use static files with hip joint centers
Settings.PrintResults = 'No'; % print results in command window


Settings.Billy = 0; % Billy set to 1
% % matches speeds when adding torso
% % sets stride parsing to right side only
% % uses these strides for identifying trial start/stop times
% Ind = contains({CurrPath.name}, 'Stride_Times_BILLY.xlsx');
% StrideFile = strcat(CurrPath(Ind).folder, '\', CurrPath(Ind).name);

%% Set location of directories
% Set subject location
if exist('path','var')
    StudyFolder = path; 
else
    StudyFolder = uigetdir(CurrPath(1).folder, 'Select Folder Containing Subject Data');
end
SubjPath = dir(StudyFolder);
addpath(genpath(StudyFolder));

% identify subject and trials within folder
NumSubjFolders = sum([SubjPath.isdir]==1)-2; % identify number of subjects
SubjFolders = [SubjPath.isdir];

% Loop Through Subjects% initialize subject structure
Subjects(NumSubjFolders).name = [];
Subjects(NumSubjFolders).Folders = [];
Subjects(NumSubjFolders).Demo = [];
Subjects(NumSubjFolders).Trials = [];
SubjCounter = 1;


% locate cal file
FPcal_file = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\2_anc_to_forces\forcepla_Nov2019_updated.cal'; % Ricky's path

% Billy change for your path to cal file
if Settings.Billy == 1
    FPcal_file = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CodeLibrary\2_anc_to_forces\forcepla_April2018.cal';
end

if ~ exist(FPcal_file, 'file') % if default file doesnt exist, select a new file
    [FPcal_file, calfilepath] = uigetfile('.cal', 'Select a force plate calibration file');
    addpath(genpath(calfilepath));
    clearvars calfilepath
end


%% Loop through each subject
for SubjLoop = 3:length(SubjFolders) % first two items are garbage
    
    if SubjFolders(SubjLoop) == 0
        continue % skip iteration bc it isnt a subject folder
    end
    
    tic; % start time counter
    
    % Define subject filename, folder, path, and directory
    Subjects(SubjCounter).name = SubjPath(SubjLoop).name;
    Subjects(SubjCounter).Folders.folder = SubjPath(SubjLoop).folder;
    Subjects(SubjCounter).Folders.path = strcat(Subjects(SubjCounter).Folders.folder, '\', Subjects(SubjCounter).name);
    Subjects(SubjCounter).Folders.dir = dir(Subjects(SubjCounter).Folders.path);
    
    % display which subject is being processed
    disp(' ');
    disp(strcat('PROCESSING__', Subjects(SubjCounter).name));
    disp(' ');
    
    % if non-existent, create OpenSim Folder to store data
    Subjects(SubjCounter).Folders.OpenSimFolder = strcat(Subjects(SubjCounter).Folders.path, '\OpenSim');
    if ~exist(Subjects(SubjCounter).Folders.OpenSimFolder, 'dir')
        mkdir(Subjects(SubjCounter).Folders.OpenSimFolder);
    end
    
    % identify of number of trials
    NumTrials = sum(contains({Subjects(SubjCounter).Folders.dir.name}, '.trc'));
    Trials(NumTrials).name = [];   % intialize Trial sub-structure
    Trials(NumTrials).subject = [];
    Trials(NumTrials).folder = [];
    
    
    %% Set up Individual Trial Loops
    TrialCounter = 1;  % trial loop counter
    for j = find(contains({Subjects(SubjCounter).Folders.dir.name}, '.trc')) % identify TRC files to loop through
        
        Str = strsplit(Subjects(SubjCounter).Folders.dir(j).name(1:end-4), '_'); 
        Trials(TrialCounter).name = Str{1}; % save trial name
        TrialName = Subjects(SubjCounter).Folders.dir(j).name(1:end-4); % save trial name as variable
        Trials(TrialCounter).subject = Subjects(SubjCounter).name; % save subject name
        Trials(TrialCounter).folder = Subjects(SubjCounter).Folders.path; % save folder name
        
        % if static trial - label accordingly
        if contains(Trials(TrialCounter).name,'Static') || ...
                contains(Trials(TrialCounter).name,'static')
            Static_name = Trials(TrialCounter).name;
            Trials(TrialCounter).type = 'static';
        else
            Trials(TrialCounter).type = 'walking';
        end
        
        % define specific filenames
        Trials(TrialCounter).files.OpenSimTRC = strcat( TrialName, '_OpenSim.trc');
        Trials(TrialCounter).files.OpenSimGRF = strcat( TrialName, '_OpenSimGRF.mot');
        Trials(TrialCounter).files.Forces = strcat( TrialName, '.forces');
        Trials(TrialCounter).files.TRC = strcat( TrialName, '.trc');
        Trials(TrialCounter).files.ANC = strcat( TrialName, '.anc');
        
        TrialCounter = TrialCounter + 1; % on to next trial (next TRC file)
    end
    WalkingTrials = contains({Trials.type}, 'walking');
    StaticTrials = ~WalkingTrials;
    clearvars TrialCounter
    
    %% Convert Force Files (MOT)
    if strcmp(Settings.ConvertForces, 'Yes')
        
        % Pre-load FPcal_file to speed up process
        [FPcal_data.S, FPcal_data.pos, FPcal_data.origin, FPcal_data.R] = load_fpcal(FPcal_file);
        
        for i = 1:length(Trials)
            % two different MOT conversion type depending on application
            if strcmp(Settings.ConvertForceType, 'Forces') || strcmp(Settings.ConvertForceType, 'Both')
                %                 [mass, FPused, files, directoryf] =
                convertFPdata(Trials(i).files.ANC, FPcal_file, Trials(i).folder);
                
            elseif strcmp(Settings.ConvertForceType, 'OpenSim') || strcmp(Settings.Settings.ConvetForceType, 'Both')
                
                cd(Trials(i).folder); % change to folder with trials
                
                %                 [trials(i).mass, trials(i).FPused, trials(i).files, trials(i).directoryf, trials(i).Analog] = ...
                convertFPdata_OpenSim(Trials(i).files.ANC, FPcal_data, Trials(i).folder, 0, Settings);
                
            end
        end
        
        disp(' ');
        disp('All .FORCE files converted to .MOT');
    end
    clearvars AbovePath ToCut FileDir Files2move i j newFile
    cd(CurrPath(1).folder);
    
    
    %% Convert Marker Trajectories
    if strcmp(Settings.ConvertTrajectories, 'Yes')
        
        %% Make any changes to the TRC File
        ROTATE = 1; % setting to rotate coordinate system
        Angle=deg2rad(-90); % angle by which to rotate
        
        Gaps(NumTrials).Frames = []; % initialize gaps structure
        Gaps(NumTrials).Times = [];
        Gaps(NumTrials).Log = [];
        
        for i = 1:length(Trials) % loop through files to convert
            disp(['...Converting ' Trials(i).name]);
            
            % Copy the TRC files to a OpenSim folder to prevent overwriting
            copyfile(strcat(Trials(i).folder,'\', Trials(i).files.TRC), ...
                strcat(Subjects(SubjCounter).Folders.OpenSimFolder, '\',  Trials(i).files.OpenSimTRC));
            addpath(Subjects(SubjCounter).Folders.OpenSimFolder);  % add output to path
            file = strcat(Subjects(SubjCounter).Folders.OpenSimFolder, '\',  Trials(i).files.OpenSimTRC);
            
            % load TRC data
            data = importdata(file,'\t',5);
            TRC = data.data;
            TRCnew = TRC;
            
            if isnan(TRCnew(1,2))  %delete first row if necessary
                TRC(1,:) = [];
                TRCnew(1,:) = [];
            end
            if isnan(TRCnew(2,1))  %delete first column if necessary
                TRC(:,1) = [];
                TRCnew(:,1) = [];
            end
            
            
            %% check for gaps or undefined data in original marker file
            if sum(sum(isnan(TRCnew))) > 0
                
                % log as trial with gaps in structure and logical
                Gaps(i).Log = 1; 
                Trials(i).type = 'GAPS';
                Lines = find(sum(isnan(TRCnew), 2));
                Gaps(i).Trial = Trials(i).name;
                Gaps(i).Times = TRCnew(Lines, 2);
                Gaps(i).Frames = TRCnew(Lines, 1);
                
            else
                Gaps(i).Log = 0; 
                
                [NumRow, NumColumn] = size(TRC); % get matrix dimensions
                % save times of trial
                Trials(i).Times.TRC = TRC(:,2);
                
                %% Rotate Coordinate System
                % in new ABL lab (MEJ 10306) always need to rotate coordinate system.
                if ROTATE == 1 % same as method above, flip Z and Y
                    n = 3;
                    for ii = 1:(size(TRC, 2) - 2) / 3
                        TRCnew(:, n) = TRC(:, n); %X
                        TRCnew(:, n+1) = TRC(:, n+2); % old Z to new Y
                        TRCnew(:, n+2) = -TRC(:, n+1); % negative old Y to new Z
                        n = n+3;
                    end
                end
                clearvars n ii ij
                
                %% Write data to new TRC file
                % Correct Headers
                fid = fopen(file,'r'); %Read TRC files to determine correct headers
                FData = textscan(fid,'%q'); % Looks at headers
                FData = FData{1,1}; % Reorganizes
                FData{4,1} = file; % Renames file name
                fclose(fid);
                
                % find header locations in file
                StartMetrics = find(strcmp(FData,'DataRate'), 1);
                LastCategory = find(strcmp(FData,'OrigNumFrames'), 1);
                ColumnHeaderStart = find(strcmp(FData,'Frame#'), 1);
                
                fid = fopen(file,'w'); % Open TRC files to write
                % Copies over header info
                fprintf(fid,'%s\t', FData{1:StartMetrics-1,1});
                fprintf(fid,'\n'); % new line
                fprintf(fid,'%s\t', FData{StartMetrics:LastCategory,1});
                fprintf(fid,'\n'); % new line
                fprintf(fid,'%s\t', FData{LastCategory+1:ColumnHeaderStart-1,1});
                fprintf(fid,'\n'); % new line
                fprintf(fid,'%s\t', FData{ColumnHeaderStart:ColumnHeaderStart+1,1});
                
                % Copies over column headers with two tabs in between each marker
                NCol = (NumColumn-2) / 3;
                for m=1:NCol
                    fprintf(fid,'%s\t', FData{m+ColumnHeaderStart+1,1});
                    fprintf(fid,'\t'); % tab
                    fprintf(fid,'\t'); % tab
                end
                
                fprintf(fid,'\n');  % new line and
                fprintf(fid,'\t');  % tab over twice to pass the Frame and Time columns
                fprintf(fid,'\t');
                
                % Labels x, y, and z columns
                d = 0;
                for m = 1:NCol
                    d = d+1;
                    fprintf(fid,'%s\t', ['X',num2str(d)]);
                    fprintf(fid,'%s\t', ['Y',num2str(d)]);
                    fprintf(fid,'%s\t', ['Z',num2str(d)]);
                end
                
                % Adds a space between xyz and data.
                fprintf(fid,'\n');
                fprintf(fid,'\n');
                
                % Inputs new/old marker x, y, z info
                for ii = 1:NumRow
                    fprintf(fid,'%3.4f\t',TRCnew(ii,:)); fprintf(fid,'\n');
                end
                
                clearvars data FData fid TRC TRCnew ii d m NCol NumRow NumColumn
                fclose('all');
                
                % This seems pointless, but OpenSim won't recognize the file unless it is opened and saved again
                e = actxserver('excel.application');
                eW = e.Workbooks;
                eF = eW.Open(strcat(Subjects(SubjCounter).Folders.OpenSimFolder, '\', Trials(i).files.OpenSimTRC)); % open file
                eS = eF.ActiveSheet;
                eF.Save;
                eF.Close; % close the file
                e.Quit; % close Excel entirely
                
            end
        end
        
    end
    clearvars e eF eS eW file i ans Angle
    
    
    %% If Gaps present in any data, log and display error to correct
    if sum([Gaps.Log]) > 0
        %                 Msg = [Trials(i).subject ,'    ', Trials(i).name, '    ', 'has gaps or invalid data'];
        i = 1;
        Msg = ['Gaps or invalid points present in ' Trials(i).subject, '. See "Gaps" structure for details'];
        error(Msg);
    else
        disp('All New TRC Files Saved');
    end
    
    
    %% Compare GRFs and Forces
    % check to make sure .Forces and OpenSimGRF.mot have the same results
    if strcmp(Settings.CompareGRFandForces, 'Yes')
        disp(' ');
        disp('Comparing GRFs and COPs from OpenSim conversion');
        
        for i = find(WalkingTrials)
            % load GRFs
            Trials(i).GRF = LoadGRF(Trials(i).files.OpenSimGRF, 0);
            
            % Compare forces
            [Trials(i).Forces, Trials(i).GRF] = CompareGRFandCOPs(Trials(i).files.Forces, Trials(i).GRF, 'Yes');
        end
    end
    
    %% Get subject mass from static trial
    i = find(StaticTrials);
    disp(' ');
    disp('Getting subject mass from static trial');
    Trials(i).GRF = LoadGRF(strcat(Trials(i).folder, '\OpenSim\', Trials(i).files.OpenSimGRF), 0); % load grfs
    Cols = {'ground_force_vy','1_ground_force_vy'}; % define vGRF columns
    vGRFs = Trials(i).GRF.AllData(:,contains(Trials(i).GRF.ColHeaders,Cols)); % identify columns from full dataset
    AvgvGRF = mean(sum(vGRFs,2)); % sum and average
    Subjects(SubjCounter).Demo.mass = AvgvGRF ./ 9.81; % calculate mass from newtons
    clearvars Cols vGRFs AvgvGRF
    disp(' ');
    
    %% Identify Crossover steps
    if strcmp(Settings.IdentifyCrossover, 'Yes')
        disp(' ');
        disp('Identifying Crossover Steps');
        if Settings.Billy == 0
            
            % initialize figure for speed
            CrossoverFig = figure('Position',[100 100 1200 800]);
            
            % crossover analysis
            [Trials] = CrossoverAnalysis(Trials, WalkingTrials, Subjects(SubjCounter).Demo.mass, Settings, CrossoverFig);
            
            close; % close figure
            
        else
            
            for i = find(WalkingTrials)
                % load GRF data if the field doesnt exist in Trials structure
                % or if the structure is empty for Trial of interest
                if ~isfield(Trials(i), 'GRF')
                    Trials(i).GRF = LoadGRF(Trials(i).files.OpenSimGRF, 0);
                elseif isempty(Trials(i).GRF)
                    Trials(i).GRF = LoadGRF(Trials(i).files.OpenSimGRF, 0);
                end
                
                % get temporal spatial data
                Trials(i).TSData = TreadmillTempSpatData(Trials(i).GRF);
            end
            
        end
    end
    
    
    %% Ensure no crossover steps are in window
    % ensure GRF and TSData fields exist in the Trials structure
    for i = find(WalkingTrials)
        if isempty(Trials(i).GRF) % load GRF data
            Trials(i).GRF = LoadGRF(Trials(i).files.OpenSimGRF);
        end
        
        if isfield(Trials, 'TSData') == 0 % compute temporal spatial gait events
            Trials(i).TSData = TreadmillTempSpatData(Trials(i).GRF);
        elseif isempty(Trials(i).TSData)
            Trials(i).TSData = TreadmillTempSpatData(Trials(i).GRF);
        end
    end
    
    if strcmp(Settings.ExcludeCrossover, 'Yes')
        disp(' ');
        disp('Excluding crossover steps from stride analysis windows');
        
        if Settings.TimeParse.NumWindows == 1
            
            if Settings.Billy == 1
                
                % selecting strides for for billy's study
                %                 StrideFile = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\BatchOpenSim\Stride_Times_BILLY.xlsx';
                if exist(StrideFile, 'file') == 0 % locate stride times file
                    [Stridefile, StrideFilePath] = uigetfile('Select file containing stride times of interest');
                    StrideFile = strcat(StrideFilePath, '\', Stridefile);
                end
                [~, ~, BillyStrides] = xlsread(StrideFile, 'Strike1');
                BillyTimes = BillyStrides(1:9, 2:end);
                col = find(contains(BillyTimes(1,:), Subjects(SubjCounter).name));
                TrialTimes = horzcat(BillyTimes(:,1), BillyTimes(:,col));
                Window = 1;
                Margin = 0.1;
                
                for i = find(WalkingTrials)
                    
                    Str = strsplit(Trials(i).name, '_');
                    
                    for j = 1:length(TrialTimes)
                        if strcmp(TrialTimes{j,1}, Str{1})
                            
                            Trials(i).Times.StartWindow = TrialTimes{j,2};
                            
                            % double check timing by looking for foot strike at that instance
                            [Rval, Rind] = min(abs(Trials(i).TSData.R_Strike(:,2) - Trials(i).Times.StartWindow));
                            
                            % get next window of foot strikes
                            Rinds = Rind:Rind + Window;
                            
                            % apply number of strides and sidedness
                            % subtract error and margin for start
                            Trials(i).Times.Start = Trials(i).TSData.R_Strike(Rinds,2) - Rval - Margin;
                            % add error and margin for end
                            Trials(i).Times.End = Trials(i).TSData.R_Strike(Rinds+1,2) + Rval + Margin;
                            
                            break
                        end
                        
                    end
                    
                end
                
                
            else
                % SPECIFIC TO METABOLICS BIOFEEDBACK DATA PROCESSING
                % determine if crossover steps occur in 10 step window for trials
                [~, ~, BestStrides] = xlsread('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\BestStrides.xlsx');
                BestTimes = BestStrides(3:end, 16:end);
                col = contains(BestTimes(1,:), Subjects(SubjCounter).name);
                TrialTimes = horzcat(BestTimes(:,2), BestTimes(:,col));
                Window = 10;
                Margin = 0.1; % margin on either side of the gait event for start/finish for each window
                NumStrides = Settings.TimeParse.NumStrides; % set number of strides to extract for RRA and later on
                
                % find first NumStrides consecutive strides in the window that dont have a crossover step
                
                % link times to trials
                for i = find(WalkingTrials)
                    
                    Str = strsplit(Trials(i).name, '_');
                    for j = 1:length(TrialTimes)
                        if strcmp(TrialTimes{j,1}, Str{1})
                            Trials(i).Times.StartWindow = TrialTimes{j,2};
                            % double check timing by looking for foot strike at that instance
                            [Lval, Lind] = min(abs(Trials(i).TSData.L_Strike(:,2) - Trials(i).Times.StartWindow));
                            [Rval, Rind] = min(abs(Trials(i).TSData.R_Strike(:,2) - Trials(i).Times.StartWindow));
                            [~, ind] = min([Lval, Rval]); % find which side steps first
                            
                            % get next 10 foot strikes
                            LStrikes = Lind:Lind + Window;
                            RStrikes = Rind:Rind + Window;
                            
                            % determine if any are crossover steps
                            AllLefts = ~ismember(LStrikes, Trials(i).Cross.L_Stride);
                            AllRights = ~ismember(RStrikes, Trials(i).Cross.R_Stride);
                            AllSteps = AllLefts + AllRights == 2;
                            
                            % get first 3 non-crossover steps
                            Linds = LStrikes(find(AllSteps, NumStrides));
                            Rinds = RStrikes(find(AllSteps, NumStrides));
                            
                            % apply number of strides and sidedness
                            if ind == 1 % if left strike is closer to start window
                                % subtract error and margin for start
                                Trials(i).Times.Start = Trials(i).TSData.L_Strike(Linds,2) - Lval - Margin;
                                % add error and margin for end
                                Trials(i).Times.End = Trials(i).TSData.R_Strike(Rinds+2,2) + Lval + Margin;
                            else % if right strike is closer to start window
                                % subtract error and margin for start
                                Trials(i).Times.Start = Trials(i).TSData.R_Strike(Rinds,2) - Rval - Margin;
                                % add error and margin for end
                                Trials(i).Times.End = Trials(i).TSData.L_Strike(Linds+2,2) + Rval + Margin;
                            end
                            
                            % make sure stride times arent before trial start or after trial ends
                            if  Trials(i).Times.Start < Trials(i).Times.TRC(1)
                                Trials(i).Times.Start = Trials(i).Times.TRC(1);
                            end
                            
                            if  Trials(i).Times.End > Trials(i).Times.TRC(end)
                                Trials(i).Times.End= Trials(i).Times.TRC(end);
                            end
                            
                            if strcmp(Settings.PrintResults, 'Yes') % print results of stride times
                                disp(strcat('Trial : ' , Trials(i).name));
                                for j = 1:NumStrides
                                    disp(strcat('Stride : ', num2str(j), '      Start : ', num2str(Trials(i).Times.Start(j)),...
                                        '      End : ', num2str(Trials(i).Times.End(j))));
                                end
                                disp(' ');
                            end
                            
                            break
                        end
                        
                    end
                    
                end % walking trial loop
                
            end
        end
    end 
    
    
    %% Add Torso
    if strcmp(Settings.AddTorso, 'Yes')
        disp(' ');
        disp('Adding Torso Motion to files');
        
        Trials = AddTorso(Trials, WalkingTrials, Settings);
        
        % move non-torso trials to non-torso folder
        mkdir(strcat(Subjects(SubjCounter).Folders.OpenSimFolder, '\NoTorso'));
        
        for i = 1:length(Trials)
            Src = strcat(Subjects(SubjCounter).Folders.OpenSimFolder, '\', Trials(i).files.OpenSimTRC);
            Dst = strcat(Subjects(SubjCounter).Folders.OpenSimFolder, '\NoTorso\', Trials(i).files.OpenSimTRC);
            movefile(Src, Dst);
%             Trials(i).files.OpenSimAddTorso = strcat(Trials(i).files.OpenSimTRC(1:end-4), 'AddTorso.trc');
        end
    end
    
    %% Finalize iteration for subject loop
    
    % save TRC marker data
    for i = 1:length(Trials)
        Trials(i).TRC.data = Trials(i).data;
        Trials(i).TRC.textdata = Trials(i).textdata;
        Trials(i).TRC.colheaders = Trials(i).colheaders;
    end
    
    % subject structure update
    Fields2Del = {'GRF', 'FiltData', 'AllCycles','Cycles2','Export', 'NewData','MedFiltData',...
        'Cycles','data','textdata','colheaders', 'AddTorsoFullFileName'};
    Trials = rmfield(Trials, Fields2Del); % remove GRFs to save space
    Subjects(SubjCounter).Trials = Trials;
    
    
    T = toc;
    disp(' ');
    disp(strcat([Subjects(SubjCounter).name, ' processed in:', ' ', num2str(round(T)), ' seconds']));
    disp(' ');
    disp(' ');
    
    clearvars TrialCounter Trials FData TRCnew TRC col ans d e eF eS eW filename2 i ii ij ind j k m n Str ...
        Ind Lind Linds Lval Margin Rind Rinds Rval ROTATE Src Dst val StartMetrics Gaps CrossoverFig ...
        AllLefts AllRights AllSteps ColumnHeaderStart LStrikes RStrikes TrialTimes Window
    
    % iterate to next subject
    SubjCounter = SubjCounter + 1;
    
end

%% save Subjects
save(strcat(StudyFolder, '\Subjects.mat'), 'Subjects');

disp(' ');
disp('ALL SUBJECTS PROCESSED');

end






