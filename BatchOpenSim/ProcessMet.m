% Process CMC Metabolics Results

clear;
clc;
dbstop if error
cd('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts');
addpath(genpath('BatchOpenSim'));
addpath(genpath('CodeLibrary'));

% load subjext files
subjectPath = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\exp';
if ~exist(subjectPath, 'dir')
    subjectPath = uigetdir(CurrFolder, 'Select Folder Containing Subject Data');
end
addpath(genpath(subjectPath));
load('Subjects.mat');

% define metabolics folder
MetabolicsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\MetabolicReports';
addpath(genpath(MetabolicsFolder));
MetabolicsDir = dir(MetabolicsFolder);

% add walking speed to demographics
if isfield(Subjects, 'PrefWalkSpeed') == 0
    Demo.File = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Metabolics_Feedback_Demographics';
    [~,~,Demo.Raw] = xlsread(Demo.File);
    Demo.WalkSpeedCol = strcmp(Demo.Raw(1,:), 'Pref Speed (m/s)');
    
    % loop through subjects
    for subj = 1:length(Subjects)
        SubjRow = find(strcmp(Demo.Raw(:,1), Subjects(subj).name),1);
        Subjects(subj).Demo.WalkSpeed = Demo.Raw{SubjRow, Demo.WalkSpeedCol}; % walk speed in m/s
        Subjects(subj).PrefWalkSpeed =  Demo.Raw{SubjRow, Demo.WalkSpeedCol};
    end
end

%% Extract Metabolic data
clc;
IsStoFile = ~[MetabolicsDir.isdir];
dbstop if error
j = 1;
Nfiles =  sum([MetabolicsDir.isdir]==0);
Met(Nfiles).Subject = [];
for i = 1:length(MetabolicsDir)
    if IsStoFile(i) == 1
        
        % get stride times from GRF for each trial
        Str = strsplit(MetabolicsDir(i).name, '_');
        Met(j).Subject = Str{1};
        
        TrialName = Str{2};
        Subj = contains({Subjects.name}, Str{1});
        Trial = contains({Subjects(Subj).Trials.name}, Str{2});
        TSData = Subjects(Subj).Trials(Trial).TSData;
        
        % create structure
        Met(j).Trial = TrialName;
        % label as left or right side
        if contains(MetabolicsDir(i).name, 'Left')
            Met(j).Side = 'Left';
        else
            Met(j).Side = 'Right';
        end
        
        % get muscle and joint metabolics
        [Met(j).UMB, Met(j).BHAR, Met(j).Data] = GetMuscleMetabolics(MetabolicsDir(i).name, TSData, MetabolicsFolder);
        
        % get cycle times for computing integrals later on
        Met(j).Time = Met(j).Data.Interp(:,1);
        
        % extract gait events for each side
        if strcmp(Met(j).Side, 'Left')
            [~, Ind] = min(abs(Met(j).Time(1) - TSData.L_Strike(:,2)));
            Met(j).GC_Start = TSData.L_Strike(Ind,2);
            [~, Met(j).GC_StartInd] = min(abs(Met(j).Time - TSData.L_Strike(Ind,2)));
            Met(j).GC_End = TSData.L_Strike(Ind+1,2);
            [~, Met(j).GC_EndInd] = min(abs(Met(j).Time - TSData.L_Strike(Ind+1,2)));
            Met(j).GC_Off = TSData.L_Off(Ind,2);
            [~, Met(j).GC_OffInd] = min(abs(Met(j).Time - TSData.L_Off(Ind,2)));
            
            OppInd = find(TSData.R_Strike(:, 2) > TSData.L_Strike(Ind, 2), 1);
            Met(j).GC_2DS = TSData.R_Strike(OppInd, 2);
            [~,Met(j).GC_2DSInd] =  min(abs(Met(j).Time - TSData.R_Strike(OppInd,2)));
            Met(j).GC_1DS = TSData.R_Off(OppInd-1, 2);
            [~, Met(j).GC_1DSInd] =  min(abs(Met(j).Time - TSData.R_Off(OppInd-1,2)));
            
        else
            % for right side
            [~, Ind] = min(abs(Met(j).Time(1) - TSData.R_Strike(:,2)));
            Met(j).GC_Start = TSData.R_Strike(Ind,2);
            [~, Met(j).GC_StartInd] = min(abs(Met(j).Time - TSData.R_Strike(Ind,2)));
            Met(j).GC_End = TSData.R_Strike(Ind+1,2);
            [~, Met(j).GC_EndInd] = min(abs(Met(j).Time - TSData.R_Strike(Ind+1,2)));
            Met(j).GC_Off = TSData.R_Off(Ind,2);
            [~, Met(j).GC_OffInd] = min(abs(Met(j).Time - TSData.R_Off(Ind,2)));
            
            OppInd = find(TSData.L_Strike(:, 2) > TSData.R_Strike(Ind, 2), 1);
            Met(j).GC_2DS = TSData.L_Strike(OppInd, 2);
            [~,Met(j).GC_2DSInd] =  min(abs(Met(j).Time - TSData.L_Strike(OppInd,2)));
            Met(j).GC_1DS = TSData.L_Off(OppInd-1, 2);
            [~, Met(j).GC_1DSInd] =  min(abs(Met(j).Time - TSData.L_Off(OppInd-1,2)));
        end
        % create logical of TS times for each file
        Z = zeros(length(Met(j).Time),1);                       % stance
        Z(Met(j).GC_StartInd: Met(j).GC_OffInd-1) = 1;
        Met(j).LogTimes.Stance = logical(Z);
        Z = zeros(length(Met(j).Time),1);                       % 1DS
        Z(Met(j).GC_StartInd: Met(j).GC_1DSInd-1) = 1;
        Met(j).LogTimes.DS1 = logical(Z);
        Z = zeros(length(Met(j).Time),1);                       % single support
        Z(Met(j).GC_1DSInd: Met(j).GC_2DSInd-1) = 1;
        Met(j).LogTimes.SingSup = logical(Z);
        Z = zeros(length(Met(j).Time),1);                       % 2DS
        Z(Met(j).GC_2DSInd: Met(j).GC_OffInd-1) = 1;
        Met(j).LogTimes.DS2 = logical(Z);
        Z = zeros(length(Met(j).Time),1);                       % swing
        Z(Met(j).GC_OffInd: Met(j).GC_EndInd) = 1;
        Met(j).LogTimes.Swing = logical(Z);
        Z = zeros(length(Met(j).Time),1);                       % stride
        Z(Met(j).GC_StartInd: Met(j).GC_EndInd) = 1;
        Met(j).LogTimes.Stride = logical(Z);
        
        clearvars Zeros Ind OppInd
        j = j + 1;
    end
end

clearvars i j k Str Subj Trial TrialName IsStoFile Z

%% Check other Outputs

% Muscle Forces (actuators)
ResultsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\ActuatorForces';
% [Actuators] = CheckActuators(ResultsFolder, Subjects);
[Actuators, Subjects] = CheckActuators(ResultsFolder, Subjects, 'Yes', 'Yes');

% muscle controls (activation)
ResultsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\MuscleControls';
[Activations, Subjects] = CheckActivations(ResultsFolder);

% kinematics
ResultsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\Kinematics';
[Kinematics, Subjects] = CheckKinematics(ResultsFolder);



%% Combine left and right Met trials for the same subject into structure
clc;
for i = 1:length(Subjects) % loop through subjects
    
    for j = 1:length(Subjects(i).Trials) % loop through trials
        
        %     get first right/left stride
        subj = contains({Met.Subject}, Subjects(i).name);
        Str = strsplit(Subjects(i).Trials(j).name, '_');
        trial = contains({Met.Trial}, Str{1});
        left =  contains({Met.Side}, 'Left');
        right = contains({Met.Side}, 'Right');
        
        SubjTrialLeft = subj + trial + left == 3;
        SubjTrialRight = subj + trial + right == 3;
        
        if sum(SubjTrialLeft) > 0
            Subjects(i).Trials(j).Met.Left = Met(SubjTrialLeft);
        end
        if sum(SubjTrialRight) > 0
            Subjects(i).Trials(j).Met.Right = Met(SubjTrialRight);
        end
        
        clearvars SubjTrialLeft SubjTrialRight left right Str subj trial
    end
end

clearvars MetabolicsDir MetabolicsFolder subjectPath TSData i j Nfiles Z ans

%% process Met files
clc;
for i = 1:length(Subjects) % loop through subjects
    
    for j = 1:length(Subjects(i).Trials)
        
        % skip trial/subject if no field or if an empty structure
        if isfield(Subjects(i).Trials(j), 'Met') == 0
            continue
        end
        if isempty(Subjects(i).Trials(j).Met)
            continue
        end
        
        % average left and right strides and integrate over periods of
        % interest
        clearvars Array Mat
        % Hip
        Subjects(i).Trials(j).Met.Hip.Columns = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Hip.Columns;
        % Umberger
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Hip.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Hip.right_parsed;
        Subjects(i).Trials(j).Met.Hip.UMB_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Hip.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Hip.right_total_parsed;
        Subjects(i).Trials(j).Met.Hip.UMB_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Hip.UMB_Avg_Sum = sum(Subjects(i).Trials(j).Met.Hip.UMB_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Hip.UMB_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Hip.UMB_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Hip.UMB_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Hip.left_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Hip.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Hip.left_total_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Hip.right_total_parsed, ...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        % Bhargava
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Hip.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Hip.right_parsed;
        Subjects(i).Trials(j).Met.Hip.BHAR_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Hip.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Hip.right_total_parsed;
        Subjects(i).Trials(j).Met.Hip.BHAR_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Hip.BHAR_Avg_Sum = sum(Subjects(i).Trials(j).Met.Hip.BHAR_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Hip.BHAR_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Hip.BHAR_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Hip.BHAR_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Hip.left_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Hip.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Hip.left_total_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Hip.right_total_parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        % Knee
        Subjects(i).Trials(j).Met.Knee.Columns = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Knee.Columns;
        % Umberger
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Knee.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Knee.right_parsed;
        Subjects(i).Trials(j).Met.Knee.UMB_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Knee.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Knee.right_total_parsed;
        Subjects(i).Trials(j).Met.Knee.UMB_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Knee.UMB_Avg_Sum = sum(Subjects(i).Trials(j).Met.Knee.UMB_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Knee.UMB_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Knee.UMB_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Knee.UMB_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Knee.left_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Knee.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Knee.left_total_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Knee.right_total_parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        % Bhargava
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Knee.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Knee.right_parsed;
        Subjects(i).Trials(j).Met.Knee.BHAR_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Knee.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Knee.right_total_parsed;
        Subjects(i).Trials(j).Met.Knee.BHAR_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Knee.BHAR_Avg_Sum = sum(Subjects(i).Trials(j).Met.Knee.BHAR_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Knee.BHAR_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Knee.BHAR_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Knee.BHAR_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Knee.left_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Knee.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Knee.left_total_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Knee.right_total_parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        
        % Ankle
        Subjects(i).Trials(j).Met.Ankle.Columns = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Ankle.Columns;
        % Umberger
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Ankle.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Ankle.right_parsed;
        Subjects(i).Trials(j).Met.Ankle.UMB_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Ankle.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Ankle.right_total_parsed;
        Subjects(i).Trials(j).Met.Ankle.UMB_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Ankle.UMB_Avg_Sum = sum(Subjects(i).Trials(j).Met.Ankle.UMB_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Ankle.UMB_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Ankle.UMB_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Ankle.UMB_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Ankle.left_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Ankle.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Ankle.left_total_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Ankle.right_total_parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        % Bhargava
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Ankle.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Ankle.right_parsed;
        Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Ankle.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Ankle.right_total_parsed;
        Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_Sum = sum(Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Ankle.BHAR_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Ankle.left_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Ankle.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Ankle.left_total_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Ankle.right_total_parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        % Trunk
        Subjects(i).Trials(j).Met.Trunk.Columns = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Trunk.Columns;
        % Umberger
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Trunk.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Trunk.right_parsed;
        Subjects(i).Trials(j).Met.Trunk.UMB_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Trunk.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Trunk.right_total_parsed;
        Subjects(i).Trials(j).Met.Trunk.UMB_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Trunk.UMB_Avg_Sum = sum(Subjects(i).Trials(j).Met.Trunk.UMB_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Trunk.UMB_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Trunk.UMB_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Trunk.UMB_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Trunk.left_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Trunk.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Trunk.left_total_parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Trunk.right_total_parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        % Bhargava
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Trunk.left_parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Trunk.right_parsed;
        Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Trunk.left_total_parsed;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Trunk.right_total_parsed;
        Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_Sum = sum(Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_total_parsed);
        clearvars Array Mat
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Trunk.BHAR_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Trunk.left_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Trunk.right_parsed, ...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Trunk.left_total_parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Trunk.right_total_parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        % Total
        % Umberger
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Total.parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Total.parsed;
        Subjects(i).Trials(j).Met.Total.UMB_parsed = mean(Array, 3);
        clearvars Array
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Total.UMB_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Total.parsed, ...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Total.parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        % Bhargava
        Array(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Total.parsed;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Total.parsed;
        Subjects(i).Trials(j).Met.Total.BHAR_parsed = mean(Array, 3);
        clearvars Array
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Total.BHAR_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Total.parsed, ...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Total.parsed,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
    end
end


%% Save data into structures by each joint and normalize to body mass
k = 0;
m = 0;
n = 0;
o = 0;
p = 0;
clc;
Hip(5).Condition = [];
Knee(5).Condition = [];
Ankle(5).Condition = [];
Total(5).Condition = [];
Muscles(5).Condition = [];

for i = 1:length(Subjects) % loop through subjects
    
    for j = 1:length(Subjects(i).Trials)
        % skip if static trial
        if contains(Subjects(i).Trials(j).name, 'static')
            %         if j > 5
            continue
        end
        
        % skip trial/subject if no field or if an empty structure
        if isfield(Subjects(i).Trials(j), 'Met') == 0
            continue
        end
        if isempty(Subjects(i).Trials(j).Met)
            continue
        end
        
        % save in arrays for averaging across subjects
        Str = Subjects(i).Trials(j).name(1:4);
        if strcmp(Str, 'Fm20')
            trial = 2;
            k = k + 1;
            ind = k;
            Hip(trial).Condition = 'Fm20';
            Knee(trial).Condition = 'Fm20';
            Ankle(trial).Condition = 'Fm20';
            Total(trial).Condition = 'Fm20';
            Muscles(trial).Condition = 'Fm20';
        elseif strcmp(Str, 'Fm40')
            trial = 1;
            m = m + 1;
            ind = m;
            Hip(trial).Condition = 'Fm40';
            Knee(trial).Condition = 'Fm40';
            Ankle(trial).Condition = 'Fm40';
            Total(trial).Condition = 'Fm40';
            Muscles(trial).Condition = 'Fm40';
        elseif strcmp(Str, 'Fp20')
            trial = 4;
            n = n + 1;
            ind = n;
            Hip(trial).Condition = 'Fp20';
            Knee(trial).Condition = 'Fp20';
            Ankle(trial).Condition ='Fp20';
            Total(trial).Condition = 'Fp20';
            Muscles(trial).Condition = 'Fp20';
        elseif strcmp(Str, 'Fp40')
            trial = 5;
            o = o + 1;
            ind = o;
            Hip(trial).Condition = 'Fp40';
            Knee(trial).Condition =  'Fp40';
            Ankle(trial).Condition = 'Fp40';
            Total(trial).Condition = 'Fp40';
            Muscles(trial).Condition = 'Fp40';
        elseif strcmp(Str, 'Norm')
            trial = 3;
            p = p + 1;
            ind = p;
            Hip(trial).Condition = 'Norm';
            Knee(trial).Condition = 'Norm';
            Ankle(trial).Condition = 'Norm';
            Total(trial).Condition = 'Norm';
            Muscles(trial).Condition = 'Norm';
        end
        
        % group muscles by joint
        % hip
        Hip(trial).UMB_SubjectsIncluded{ind} = Subjects(i).name;
        Hip(trial).UMB_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Avg_parsed / Subjects(i).Demo.mass;
        Hip(trial).UMB_Total(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Avg_total_parsed / Subjects(i).Demo.mass;
        Hip(trial).UMB_Columns = Subjects(i).Trials(j).Met.Hip.Columns;
        Hip(trial).UMB_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Avg_total_Sum / Subjects(i).Demo.mass;
        Hip(trial).UMB_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Avg_Sum / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).UMB_Work.Stride_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.UMB_Energy.Stride.Avg_Muscles_Int / Subjects(i).Demo.mass;
        
        Hip(trial).BHAR_SubjectsIncluded{ind} = Subjects(i).name;
        Hip(trial).BHAR_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Avg_parsed / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Total(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Avg_total_parsed / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Columns = Subjects(i).Trials(j).Met.Hip.Columns;
        Hip(trial).BHAR_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Avg_total_Sum / Subjects(i).Demo.mass;
        Hip(trial).BHAR_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Avg_Sum / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Hip(trial).BHAR_Work.Stride_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Hip.BHAR_Energy.Stride.Avg_Muscles_Int / Subjects(i).Demo.mass;
        
        % knee
        Knee(trial).UMB_SubjectsIncluded{ind} = Subjects(i).name;
        Knee(trial).UMB_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Avg_parsed / Subjects(i).Demo.mass;
        Knee(trial).UMB_Total(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Avg_total_parsed / Subjects(i).Demo.mass;
        Knee(trial).UMB_Columns = Subjects(i).Trials(j).Met.Knee.Columns;
        Knee(trial).UMB_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Avg_total_Sum / Subjects(i).Demo.mass;
        Knee(trial).UMB_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Avg_Sum / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).UMB_Work.Stride_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.UMB_Energy.Stride.Avg_Muscles_Int / Subjects(i).Demo.mass;
        
        Knee(trial).BHAR_SubjectsIncluded{ind} = Subjects(i).name;
        Knee(trial).BHAR_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Avg_parsed / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Total(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Avg_total_parsed / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Columns = Subjects(i).Trials(j).Met.Knee.Columns;
        Knee(trial).BHAR_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Avg_total_Sum / Subjects(i).Demo.mass;
        Knee(trial).BHAR_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Avg_Sum / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Knee(trial).BHAR_Work.Stride_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Knee.BHAR_Energy.Stride.Avg_Muscles_Int / Subjects(i).Demo.mass;
        
        % ankle
        Ankle(trial).UMB_SubjectsIncluded{ind} = Subjects(i).name;
        Ankle(trial).UMB_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Avg_parsed / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Total(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Avg_total_parsed / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Columns = Subjects(i).Trials(j).Met.Ankle.Columns;
        Ankle(trial).UMB_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Avg_total_Sum / Subjects(i).Demo.mass;
        Ankle(trial).UMB_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Avg_Sum / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).UMB_Work.Stride_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.UMB_Energy.Stride.Avg_Muscles_Int / Subjects(i).Demo.mass;
        
        Ankle(trial).BHAR_SubjectsIncluded{ind} = Subjects(i).name;
        Ankle(trial).BHAR_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_parsed / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Total(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_total_parsed / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Columns = Subjects(i).Trials(j).Met.Ankle.Columns;
        Ankle(trial).BHAR_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_total_Sum / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Avg_Sum / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Ankle(trial).BHAR_Work.Stride_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Ankle.BHAR_Energy.Stride.Avg_Muscles_Int / Subjects(i).Demo.mass;
        
        % trunk
        %         Trunk(trial).UMB_SubjectsIncluded{ind} = Subjects(i).name;
        %         Trunk(trial).UMB_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Avg_parsed / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Total(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Avg_total_parsed / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Columns = Subjects(i).Trials(j).Met.Trunk.Columns;
        %         Trunk(trial).UMB_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Avg_total_Sum / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Avg_Sum / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %        Trunk(trial).UMB_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).UMB_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.UMB_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %
        %         Trunk(trial).BHAR_SubjectsIncluded{ind} = Subjects(i).name;
        %         Trunk(trial).BHAR_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_parsed / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Total(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_total_parsed / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Columns = Subjects(i).Trials(j).Met.Trunk.Columns;
        %         Trunk(trial).BHAR_SumTot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_total_Sum / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_SumMusc(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Avg_Sum / Subjects(i).Demo.mass;
        %          Trunk(trial).BHAR_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.Stance_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.Stance.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.DS1_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.DS1.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.SingSup_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.SingSup.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.DS2_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.DS2.Avg_Muscles_Int / Subjects(i).Demo.mass;
        %          Trunk(trial).BHAR_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        %         Trunk(trial).BHAR_Work.Swing_Muscles(:,:,ind) = Subjects(i).Trials(j).Met.Trunk.BHAR_Energy.Swing.Avg_Muscles_Int / Subjects(i).Demo.mass;
        
        
        % total
        Total(trial).UMB_SubjectsIncluded{ind} = Subjects(i).name;
        Total(trial).UMB_Total(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_parsed / Subjects(i).Demo.mass;
        Total(trial).UMB_Total_Avg(:,:,ind) =   mean(Total(trial).UMB_Total(:,:,ind));
        Total(trial).UMB_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        
        Total(trial).BHAR_SubjectsIncluded{ind} = Subjects(i).name / Subjects(i).Demo.mass;
        Total(trial).BHAR_Total(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_parsed / Subjects(i).Demo.mass;
        Total(trial).BHAR_Total_Avg(:,:,ind) =   mean(Total(trial).BHAR_Total(:,:,ind));
        Total(trial).BHAR_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).BHAR_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).BHAR_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).BHAR_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).BHAR_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).BHAR_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).BHAR_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        
        % extract muscle metabolics and average left and right sides
        for musc = 3:48
            Muscles(trial).UMB_Cols{musc-2} = Subjects(i).Trials(j).Met.Left.UMB.ParsedMuscles(musc).name;
            Muscles(trial).UMB_Muscles(:,musc-2, ind) = mean([...
                Subjects(i).Trials(j).Met.Left.UMB.ParsedMuscles(musc).data.left_parsed,...
                Subjects(i).Trials(j).Met.Right.UMB.ParsedMuscles(musc).data.right_parsed],2) / Subjects(i).Demo.mass;
        end
        
        % sum similar muscles over stance phase for a smaller muscle matrix
        
        % 4 muscles with multiple lines of action
        % Glute med
        Gmed_Ind = contains([Muscles(trial).UMB_Cols], 'glut_med');
        Gmed = sum(Muscles(trial).UMB_Muscles(:,Gmed_Ind, ind), 2);
        % Glute min
        Gmin_Ind = contains([Muscles(trial).UMB_Cols], 'glut_min');
        Gmin = sum(Muscles(trial).UMB_Muscles(:,Gmin_Ind, ind), 2);
        % Glute max
        Gmax_Ind = contains([Muscles(trial).UMB_Cols], 'glut_max');
        Gmax = sum(Muscles(trial).UMB_Muscles(:,Gmax_Ind, ind), 2);
        % Add Mag
        Amag_Ind = contains([Muscles(trial).UMB_Cols], 'add_mag');
        Amag = sum(Muscles(trial).UMB_Muscles(:,Amag_Ind, ind), 2);
        
        % copy over all muscles
        Muscles(trial).UMB_CombMuscles(:,1:46,ind) = Muscles(trial).UMB_Muscles(:,:,ind);
        
        % remove 4x3=12 muscles with muscles with multiple lines of action
        LOAM = sum([Gmed_Ind', Gmin_Ind', Gmax_Ind', Amag_Ind'], 2) > 0;
        Muscles(trial).UMB_CombCols = Muscles(trial).UMB_Cols;
        Muscles(trial).UMB_CombCols{end+1} = 'glut_med';
        Muscles(trial).UMB_CombCols{end+1} = 'glut_min';
        Muscles(trial).UMB_CombCols{end+1} = 'glut_max';
        Muscles(trial).UMB_CombCols{end+1} = 'add_mag';
        
        % add in summed 4 muscles
        Muscles(trial).UMB_CombMuscles(:,47:50,ind) = [Gmed, Gmin, Gmax, Amag];
        
    end
end

% delete summed muscles (1,2,3s)
LOAM_ = logical([LOAM; 0; 0; 0; 0]);
for i = 1:length(Muscles)
    Muscles(i).UMB_CombCols(LOAM_) = [];
    Muscles(i).UMB_CombMuscles(:,LOAM_, :) = [];
end

clearvars n o p Str i j k m ans musc trial ind LOAM LOAM_...
    Amag Amag_Ind Gmax Gmax_Ind Gmin Gmin_Ind Gmed Gmed_Ind

%% Generate overall averages and standard devs per joint
clc;
for trial = [1 2 3 4 5]
    Hip(trial).UMB_MusclesAvg = mean(Hip(trial).UMB_Muscles, 3);
    Hip(trial).UMB_MusclesStd = std(Hip(trial).UMB_Muscles, 0, 3);
    Hip(trial).UMB_TotalAvg = mean(Hip(trial).UMB_Total, 3);
    Hip(trial).UMB_TotalStd = std(Hip(trial).UMB_Total, 0, 3);
    Hip(trial).UMB_Sum_Avg = mean(Hip(trial).UMB_SumTot, 3);
    Hip(trial).UMB_Sum_Std =std(Hip(trial).UMB_SumTot, 0, 3);
    Hip(trial).BHAR_MusclesAvg = mean(Hip(trial).BHAR_Muscles, 3);
    Hip(trial).BHAR_MusclesStd = std(Hip(trial).BHAR_Muscles, 0, 3);
    Hip(trial).BHAR_TotalAvg = mean(Hip(trial).BHAR_Total, 3);
    Hip(trial).BHAR_TotalStd = std(Hip(trial).BHAR_Total, 0, 3);
    Hip(trial).BHAR_Sum_Avg = mean(Hip(trial).BHAR_SumTot, 3);
    Hip(trial).BHAR_Sum_Std =std(Hip(trial).BHAR_SumTot, 0, 3);
    
    Knee(trial).UMB_MusclesAvg = mean(Knee(trial).UMB_Muscles, 3);
    Knee(trial).UMB_MusclesStd = std(Knee(trial).UMB_Muscles, 0, 3);
    Knee(trial).UMB_TotalAvg = mean(Knee(trial).UMB_Total, 3);
    Knee(trial).UMB_TotalStd = std(Knee(trial).UMB_Total, 0, 3);
    Knee(trial).UMB_Sum_Avg = mean(Knee(trial).UMB_SumTot, 3);
    Knee(trial).UMB_Sum_Std =std(Knee(trial).UMB_SumTot, 0, 3);
    Knee(trial).BHAR_MusclesAvg = mean(Knee(trial).BHAR_Muscles, 3);
    Knee(trial).BHAR_MusclesStd = std(Knee(trial).BHAR_Muscles, 0, 3);
    Knee(trial).BHAR_TotalAvg = mean(Knee(trial).BHAR_Total, 3);
    Knee(trial).BHAR_TotalStd = std(Knee(trial).BHAR_Total, 0, 3);
    Knee(trial).BHAR_Sum_Avg = mean(Knee(trial).BHAR_SumTot, 3);
    Knee(trial).BHAR_Sum_Std =std(Knee(trial).BHAR_SumTot, 0, 3);
    
    Ankle(trial).UMB_MusclesAvg = mean(Ankle(trial).UMB_Muscles, 3);
    Ankle(trial).UMB_MusclesStd = std(Ankle(trial).UMB_Muscles, 0, 3);
    Ankle(trial).UMB_TotalAvg = mean(Ankle(trial).UMB_Total, 3);
    Ankle(trial).UMB_TotalStd = std(Ankle(trial).UMB_Total, 0, 3);
    Ankle(trial).UMB_Sum_Avg = mean(Ankle(trial).UMB_SumTot, 3);
    Ankle(trial).UMB_Sum_Std =std(Ankle(trial).UMB_SumTot, 0, 3);
    Ankle(trial).BHAR_MusclesAvg = mean(Ankle(trial).BHAR_Muscles, 3);
    Ankle(trial).BHAR_MusclesStd = std(Ankle(trial).BHAR_Muscles, 0, 3);
    Ankle(trial).BHAR_TotalAvg = mean(Ankle(trial).BHAR_Total, 3);
    Ankle(trial).BHAR_TotalStd = std(Ankle(trial).BHAR_Total, 0, 3);
    Ankle(trial).BHAR_Sum_Avg = mean(Ankle(trial).BHAR_SumTot, 3);
    Ankle(trial).BHAR_Sum_Std =std(Ankle(trial).BHAR_SumTot, 0, 3);
    
    %     Trunk(trial).UMB_MusclesAvg = mean(Trunk(trial).UMB_Muscles, 3);
    %     Trunk(trial).UMB_MusclesStd = std(Trunk(trial).UMB_Muscles, 0, 3);
    %     Trunk(trial).UMB_TotalAvg = mean(Trunk(trial).UMB_Total, 3);
    %     Trunk(trial).UMB_TotalStd = std(Trunk(trial).UMB_Total, 0, 3);
    %     Trunk(trial).UMB_Sum_Avg = mean(Trunk(trial).UMB_SumTot, 3);
    %     Trunk(trial).UMB_Sum_Std =std(Trunk(trial).UMB_SumTot, 0, 3);
    %     Trunk(trial).BHAR_MusclesAvg = mean(Trunk(trial).BHAR_Muscles, 3);
    %     Trunk(trial).BHAR_MusclesStd = std(Trunk(trial).BHAR_Muscles, 0, 3);
    %     Trunk(trial).BHAR_TotalAvg = mean(Trunk(trial).BHAR_Total, 3);
    %     Trunk(trial).BHAR_TotalStd = std(Trunk(trial).BHAR_Total, 0, 3);
    %     Trunk(trial).BHAR_Sum_Avg = mean(Trunk(trial).BHAR_SumTot, 3);
    %     Trunk(trial).BHAR_Sum_Std =std(Trunk(trial).BHAR_SumTot, 0, 3);
    
    Total(trial).UMB_TotalAvg = mean(Total(trial).UMB_Total, 3);
    Total(trial).UMB_TotalStd = std(Total(trial).UMB_Total, 0, 3);
    Total(trial).BHAR_TotalAvg = mean(Total(trial).BHAR_Total, 3);
    Total(trial).BHAR_TotalStd = std(Total(trial).BHAR_Total, 0, 3);
    
end


%% Correlate modelled and measured metabolic cost
clearvars R n o p
% load measured met cost
[~, ~, Raw] = xlsread('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Subject_NetMetabolics.xlsx');
c = 0;
TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray');  rgb('DarkSlateGray'); rgb('Black')];
SubjColors = colormap(jet(length(Subjects)));
TrialSymbols = {'-', '-', 'o', '+','+'};
TxtFnt = 15;
AxFnt = 12;
MkrSz = 12;
TitleFnt = 16;
close all; clc;
CorrVal = figure('Position', [100 100 1000 800]);

for subj = 1:length(Subjects) % loop through subjects
    row = contains(Raw(:,1), Subjects(subj).name);
    for trial = 1:length(Subjects(subj).Trials)  % loop through trials
        Conditions = Raw(2,2:7);
        for i = 1:length(Conditions)  % loop through conditions
            
            if contains(Subjects(subj).Trials(trial).name, Conditions{i})
                c = c + 1;
                Subjects(subj).Trials(trial).AvgMetCost = Raw{row, i+1};
                MetCost(c).Subj = Subjects(subj).name;
                MetCost(c).Trial = Subjects(subj).Trials(trial).name;
%                 MetCost(c).Measured = Subjects(subj).Trials(trial).AvgMetCost / Subjects(subj).PrefWalkSpeed;
                MetCost(c).Measured = Subjects(subj).Trials(trial).AvgMetCost;
                MetCost(c).UmbModelled = Subjects(subj).Trials(trial).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
                MetCost(c).BharModelled = Subjects(subj).Trials(trial).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
                
                % plot correlations
                %                 subplot(235);
                %                 hold on;
                %                 if contains(Subjects(subj).Trials(trial).name, 'Fm40')
                %                     %                     text(MetCost(c).Measured, MetCost(c).UmbModelled, '\bf -', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).UmbModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Fm20')
                %                     %                     text(MetCost(c).Measured, MetCost(c).UmbModelled, '-', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).UmbModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Norm')
                %                     %                     text(MetCost(c).Measured, MetCost(c).UmbModelled, 'o', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).UmbModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Fp20')
                %                     %                     text(MetCost(c).Measured, MetCost(c).UmbModelled, '+', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).UmbModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Fp40')
                %                     %                     text(MetCost(c).Measured, MetCost(c).UmbModelled, '\bf +', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).UmbModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 end
                %
                %                 subplot(236); hold on;
                %                 if contains(Subjects(subj).Trials(trial).name, 'Fm40')
                %                     %                     text(MetCost(c).Measured, MetCost(c).BharModelled, '\bf -', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).BharModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Fm20')
                %                     %                     text(MetCost(c).Measured, MetCost(c).BharModelled, '-', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).BharModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Norm')
                %                     %                     text(MetCost(c).Measured, MetCost(c).BharModelled, 'o', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).BharModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Fp20')
                %                     %                     text(MetCost(c).Measured, MetCost(c).BharModelled, '+', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).BharModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 elseif contains(Subjects(subj).Trials(trial).name, 'Fp40')
                %                     %                     text(MetCost(c).Measured, MetCost(c).BharModelled, '\bf +', 'Color', SubjColors(subj, :), 'FontSize', TxtFnt);
                %                     plot(MetCost(c).Measured, MetCost(c).BharModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %                 end
                %
                %                 % correlate models
                %                 subplot(234);          hold on;
                %                 plot(MetCost(c).BharModelled, MetCost(c).UmbModelled, '.', 'Color', SubjColors(subj, :), 'MarkerSize', MkrSz);
                %
            end
        end
    end
end

MetCostTable = struct2table(MetCost);
[R.Umb,P.Umb] = corr([MetCost.Measured]', [MetCost.UmbModelled]');
[R.Bhar,P.Bhar] = corr([MetCost.Measured]', [MetCost.BharModelled]');
[R.BharUmb,P.BharUmb] = corr([MetCost.UmbModelled]', [MetCost.BharModelled]');

% Identify trials
Cond_Ind(:,1) = contains([MetCostTable.Trial], 'Fm40');
Cond_Ind(:,2) = contains([MetCostTable.Trial], 'Fm20');
Cond_Ind(:,3) = contains([MetCostTable.Trial], 'Norm');
Cond_Ind(:,4) = contains([MetCostTable.Trial], 'Fp20');
Cond_Ind(:,5) = contains([MetCostTable.Trial], 'Fp40');

for cond = 1:5
    subplot(235); hold on;
    plot(MetCostTable.Measured(Cond_Ind(:,cond)), MetCostTable.UmbModelled(Cond_Ind(:,cond)),...
        '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
    subplot(236); hold on;
    plot(MetCostTable.Measured(Cond_Ind(:,cond)), MetCostTable.BharModelled(Cond_Ind(:,cond)),...
        '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
    subplot(234); hold on;
    plot(MetCostTable.UmbModelled(Cond_Ind(:,cond)), MetCostTable.BharModelled(Cond_Ind(:,cond)),...
        '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
end

subplot(234);
x = linspace(0, 15); 
plot(x,x, '-k'); 
ylim([7.5 15]);
xlim([7.5 15]);
title('Inter-Model', 'FontSize', TitleFnt);
ylabel('Bhargava');
xlabel('Umberger');
ax = gca;
ax.FontSize = AxFnt;
text(11,9, ['R = ' num2str(round(R.BharUmb, 4))]);
% text(11,8.5, ['P = ' num2str(P.BharUmb, 4)]);
text(11, 8.5, ['P < 0.001']);

subplot(235);
% x = linspace(0, 15); 
% plot(x,x, '-k'); 
ylim([7.5 15]);
xlim([2 10]);
title('Umberger Model', 'FontSize', TitleFnt);
xlabel('Measured');
ylabel('Simulated');
ax = gca;
ax.FontSize = AxFnt;
text(7, 9, ['R = ' num2str(round(R.Umb, 4))]);
% text(7, 8.5, ['P = ' num2str(P.Umb)]);
text(7, 8.5, ['P < 0.001']);

subplot(236);
% x = linspace(0, 15); 
% plot(x,x, '-k'); 
ylim([7.5 15]);
xlim([2 10]);
title('Bhargava Model', 'FontSize', TitleFnt);
xlabel('Measured');
ylabel('Simulated');
ax = gca;
ax.FontSize = AxFnt;
text(7, 9, ['R = ' num2str(round(R.Bhar, 4))]);
% text(7, 8.5, ['P = ' num2str(P.Bhar)]);
text(7, 8.5, ['P < 0.001']);


%% Total metabolic cost by condition
% figure; hold on;
% MkrSz = 10;
DotMargin = 0.2;
LW = 1.25;
X = [1 2 3 4 5];

% loop through subjects to plot by color
% for subj = 1:length(Subjects)
%     % define values
%     ind = contains({Subjects(subj).Trials.name}, 'Fm40');
%     %     cond(1) = mean(Subjects(subj).Trials(ind).Met.Total.UMB_parsed) / Subjects(subj).Demo.mass;
%     Mescond(1) = Subjects(subj).Trials(ind).AvgMetCost;
%     UMBcond(1) = Subjects(subj).Trials(ind).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     BHARcond(1) = Subjects(subj).Trials(ind).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     ind = contains({Subjects(subj).Trials.name}, 'Fm20');
%     %     cond(2) = mean(Subjects(subj).Trials(ind).Met.Total.UMB_parsed) / Subjects(subj).Demo.mass;
%     Mescond(2) = Subjects(subj).Trials(ind).AvgMetCost;
%     UMBcond(2) = Subjects(subj).Trials(ind).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     BHARcond(2) = Subjects(subj).Trials(ind).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     ind = contains({Subjects(subj).Trials.name}, 'Norm');
%     %     cond(3) = mean(Subjects(subj).Trials(ind).Met.Total.UMB_parsed) / Subjects(subj).Demo.mass;
%     Mescond(3) = Subjects(subj).Trials(ind).AvgMetCost;
%     UMBcond(3) = Subjects(subj).Trials(ind).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     BHARcond(3) = Subjects(subj).Trials(ind).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     ind = contains({Subjects(subj).Trials.name}, 'Fp20');
%     %     cond(4) = mean(Subjects(subj).Trials(ind).Met.Total.UMB_parsed) / Subjects(subj).Demo.mass;
%     Mescond(4) = Subjects(subj).Trials(ind).AvgMetCost;
%     UMBcond(4) = Subjects(subj).Trials(ind).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     BHARcond(4) = Subjects(subj).Trials(ind).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     ind = contains({Subjects(subj).Trials.name}, 'Fp40');
%     %     cond(5) = mean(Subjects(subj).Trials(ind).Met.Total.UMB_parsed) / Subjects(subj).Demo.mass;
%     Mescond(5) = Subjects(subj).Trials(ind).AvgMetCost;
%     UMBcond(5) = Subjects(subj).Trials(ind).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     BHARcond(5) = Subjects(subj).Trials(ind).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
%     % plot
%     subplot(231); hold on;
%     plot(X, Mescond, '-s', 'MarkerSize', MkrSz, 'Color', SubjColors(subj, :), 'MarkerFaceColor', SubjColors(subj, :));
%     subplot(232); hold on;
%     plot(X, UMBcond, '-s', 'MarkerSize', MkrSz, 'Color', SubjColors(subj, :), 'MarkerFaceColor', SubjColors(subj, :));
%     subplot(233); hold on;
%     plot(X, BHARcond, '-s', 'MarkerSize', MkrSz, 'Color', SubjColors(subj, :), 'MarkerFaceColor', SubjColors(subj, :));
% end

% Identify trials
Cond_Ind(:,1) = contains([MetCostTable.Trial], 'Fm40');
Cond_Ind(:,2) = contains([MetCostTable.Trial], 'Fm20');
Cond_Ind(:,3) = contains([MetCostTable.Trial], 'Norm');
Cond_Ind(:,4) = contains([MetCostTable.Trial], 'Fp20');
Cond_Ind(:,5) = contains([MetCostTable.Trial], 'Fp40');

for cond = 1:5
    subplot(231); hold on;
    plot(cond - DotMargin, MetCostTable.Measured(Cond_Ind(:,cond)), '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
    
    h = boxplot2(MetCostTable.Measured(Cond_Ind(:,cond)),  cond);
    i = cond;
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    
    subplot(232); hold on;
    plot(cond - DotMargin, MetCostTable.UmbModelled(Cond_Ind(:,cond)), '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
    
    h = boxplot2(MetCostTable.UmbModelled(Cond_Ind(:,cond)),  cond);
    i = cond;
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    
    subplot(233); hold on;
    plot(cond - DotMargin, MetCostTable.BharModelled(Cond_Ind(:,cond)), '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
    
    h = boxplot2(MetCostTable.BharModelled(Cond_Ind(:,cond)),  cond);
    i = cond;
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    
end

subplot(231);
xlim([0.5 5.5]);
title('Empirical', 'FontSize', TitleFnt);
ax = gca;
ax.XTick = [ 1 2 3 4 5];
ax.XTickLabel = {'-40%', '-20%','Norm','+20%','+40%'};
ax.FontSize = AxFnt;
ylabel('Metabolic Cost (W / kg)');

subplot(232);
xlim([0.5 5.5]);
ylim([7 17]);
title('Umberger Model', 'FontSize', TitleFnt);
ax = gca;
ax.XTick = [ 1 2 3 4 5];
ax.XTickLabel = {'-40%', '-20%','Norm','+20%','+40%'};
ax.FontSize = AxFnt;
ylabel('Metabolic Cost (W / kg)');

subplot(233);
xlim([0.5 5.5]);
ylim([7 17]);
title('Bhargava Model', 'FontSize', TitleFnt);
ax = gca;
ax.XTick = [ 1 2 3 4 5];
ax.XTickLabel = {'-40%', '-20%','Norm','+20%','+40%'};
ax.FontSize = AxFnt;
ylabel('Metabolic Cost (W / kg)');

subplotsqueeze(CorrVal, 1.05);
saveas(CorrVal, 'CorrVal.png');
saveas(CorrVal, 'CorrVal.pdf');

%% plot metabolic profile across the gait cycle
clc; close all;
yPk = 400; % set y peak
LW = 2; % set line width
TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray'); rgb('DarkSlateGray'); rgb('Black')];
HipColors = [rgb('PowderBlue'); rgb('DeepSkyBlue'); rgb('RoyalBlue');  rgb('Navy'); rgb('MidnightBlue')];
% KneeColors = [rgb('Orange'); rgb('DarkOrange'); rgb('Coral'); rgb('Chocolate'); rgb('Sienna')];
KneeColors = [rgb('SandyBrown'); rgb('Orange'); rgb('DarkOrange'); rgb('Chocolate'); rgb('Sienna')];
AnkleColors = [rgb('PaleGreen'); rgb('LimeGreen'); rgb('ForestGreen'); rgb('Green'); rgb('DarkGreen')];
FntSz = 12; 
JointEnergyCostsBothModels = figure('Position',[50 50 1200 800]);
TxtHt = 11.25; 
Stride = 1:100;
c = 0;
L = 0;
yPk = 12;
TotalyPk = 30;
for i = 1:5
    c = c + 1;
    
    % total
    subplot(421); hold on; % grid on;
    plot(Stride, Total(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', TotalColors(c,:));
    %     plot(Stride, Total(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
    ylabel('W / kg'); ylim([0 TotalyPk]);
    
    
    subplot(422); hold on; %grid on;
    %     plot(Stride, Total(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
    plot(Stride, Total(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', TotalColors(c,:));
    %     ylabel('W / kg');
    ylim([0 TotalyPk]);
    
    
    % Hip
    subplot(423); hold on;% grid on;
    plot(Stride, Hip(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', HipColors(c,:));
    %     plot(Stride, Hip(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
    ylabel('W / kg'); ylim([0 yPk]);
   
    subplot(424); hold on;% grid on;
    %     plot(Stride, Hip(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
    plot(Stride, Hip(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', HipColors(c,:));
    %     ylabel('W / kg');
    ylim([0 yPk]);
    
    % Knee
    subplot(425); hold on; %grid on;
    plot(Stride, Knee(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', KneeColors(c,:));
    %      plot(Stride, Knee(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
    ylabel('W / kg'); ylim([0 yPk]);
    
    subplot(426); hold on; %grid on;
    %     plot(Stride, Knee(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
    plot(Stride, Knee(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', KneeColors(c,:));
    %     ylabel('W / kg');
    ylim([0 yPk]);
    
    % Ankle
    subplot(427); hold on; %grid on;
    plot(Stride, Ankle(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', AnkleColors(c,:));
    %     plot(Stride, Ankle(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
    ylabel('W / kg'); ylim([0 yPk]);
    xlabel('% Gait Cycle');
    
    subplot(428); hold on; %grid on;
    %     plot(Stride, Ankle(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
    plot(Stride, Ankle(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', AnkleColors(c,:));
    %     ylabel('W / kg');
    ylim([0 yPk]);
    xlabel('% Gait Cycle');
    
%     L = L + 1;
%     %     CurveNameUMB{L} = strcat(Hip(i).Condition, ' UMB');
%     %     CurveNameBHAR{L} = strcat(Hip(i).Condition, ' BHAR');
%     CurveName{L} = Hip(i).Condition;
end

CurveName = {'-40','-20','Norm','+20','+40'}; 
subplot(421);
    text(50, 28, '\bf Total UMB', 'Color',TotalColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
ax = gca;
ax.XTick = [0 25 50 75 100];
ax.XTickLabel = [];
ax.FontSize = FntSz; 

subplot(422);
     text(50, 28, '\bf Total BHAR', 'Color',TotalColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% legend(CurveName, 'Location','Best');
ax = gca;
ax.XTick = [0 25 50 75 100];
ax.XTickLabel = [];
ax.YTick = [0 10 20 30];
ax.YTickLabel = [];
ax.FontSize = FntSz; 

subplot(423);
    text(50, TxtHt, '\bf Hip UMB', 'Color',HipColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
ax = gca;
ax.XTick = [0 25 50 75 100];
ax.XTickLabel = [];
ax.YTick = [0 3 6 9 12];
ax.FontSize = FntSz; 

subplot(424);
     text(50, TxtHt, '\bf Hip BHAR', 'Color',HipColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% legend(CurveName, 'Location','Best');
ax = gca;
ax.XTick = [0 25 50 75 100];
ax.XTickLabel = [];
ax.YTick = [0 3 6 9 12];
ax.YTickLabel = [];
ax.FontSize = FntSz; 

subplot(425);
   text(50, TxtHt,'\bf Knee UMB', 'Color',KneeColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
ax = gca;
ax.XTick = [0 25 50 75 100];
ax.YTick = [0 3 6 9 12];
ax.XTickLabel = [];
ax.FontSize = FntSz; 

subplot(426);
 text(50, TxtHt,'\bf Knee BHAR', 'Color',KneeColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% legend(CurveName, 'Location','Best');
ax = gca;
ax.XTick = [0 25 50 75 100];
ax.XTickLabel = [];
ax.YTick = [0 3 6 9 12];
ax.YTickLabel = [];
ax.FontSize = FntSz; 

subplot(427);
text(50, TxtHt, '\bf Ankle UMB', 'Color',AnkleColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
ax = gca;
ax.XTick = [0 25 50 75 100];
% % ax.XTickLabel = [];
ax.YTick = [0 3 6 9 12];
ax.FontSize = FntSz; 

subplot(428);
text(50, TxtHt,'\bf Ankle BHAR', 'Color',AnkleColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% legend(CurveName, 'Location','Best');
ax = gca;
ax.XTick = [0 25 50 75 100];
% ax.XTickLabel = [];
ax.YTick = [0 3 6 9 12];
ax.YTickLabel = [];
ax.FontSize = FntSz; 

subplotsqueeze(JointEnergyCostsBothModels, 1.2);
saveas(JointEnergyCostsBothModels, 'JointEnergyCostsBothModels.png');
saveas(JointEnergyCostsBothModels, 'JointEnergyCostsBothModels.pdf');

clearvars ax AxFnt c i ind int L LW Mescond MkrSz subj trial TxtFnt TotalyPk UMBcond BHARcond X yPk


%% Stats for metabolic cost
clc; clearvars Table table Meas MeasTable
Conditions = {'Hip1','Hip2','Hip3','Hip4','Hip5', 'Ankle1','Ankle2','Ankle3','Ankle4','Ankle5'};
k = 1;
for j = 1:length(Subjects)
    
    % Total
    DS1TotalTable(k).y1 = Total(1).UMB_Work.DS1_Tot(1,1,j);
    DS1TotalTable(k).y2 = Total(2).UMB_Work.DS1_Tot(1,1,j);
    DS1TotalTable(k).y3 = Total(3).UMB_Work.DS1_Tot(1,1,j);
    DS1TotalTable(k).y4 = Total(4).UMB_Work.DS1_Tot(1,1,j);
    DS1TotalTable(k).y5 = Total(5).UMB_Work.DS1_Tot(1,1,j);
    
    SingSupTotalTable(k).Subject = Subjects(j).name;
    SingSupTotalTable(k).Joint = 'Total';
    SingSupTotalTable(k).y1 = Total(1).UMB_Work.SingSup_Tot(1,1,j);
    SingSupTotalTable(k).y2 = Total(2).UMB_Work.SingSup_Tot(1,1,j);
    SingSupTotalTable(k).y3 = Total(3).UMB_Work.SingSup_Tot(1,1,j);
    SingSupTotalTable(k).y4 = Total(4).UMB_Work.SingSup_Tot(1,1,j);
    SingSupTotalTable(k).y5 = Total(5).UMB_Work.SingSup_Tot(1,1,j);
    
    DS2TotalTable(k).Subject = Subjects(j).name;
    DS2TotalTable(k).Joint = 'Total';
    DS2TotalTable(k).y1 = Total(1).UMB_Work.DS2_Tot(1,1,j);
    DS2TotalTable(k).y2 = Total(2).UMB_Work.DS2_Tot(1,1,j);
    DS2TotalTable(k).y3 = Total(3).UMB_Work.DS2_Tot(1,1,j);
    DS2TotalTable(k).y4 = Total(4).UMB_Work.DS2_Tot(1,1,j);
    DS2TotalTable(k).y5 = Total(5).UMB_Work.DS2_Tot(1,1,j);
    
    StanceTotalTable(k).Subject = Subjects(j).name;
    StanceTotalTable(k).Joint = 'Total';
    StanceTotalTable(k).y1 = Total(1).UMB_Work.Stance_Tot(1,1,j);
    StanceTotalTable(k).y2 = Total(2).UMB_Work.Stance_Tot(1,1,j);
    StanceTotalTable(k).y3 = Total(3).UMB_Work.Stance_Tot(1,1,j);
    StanceTotalTable(k).y4 = Total(4).UMB_Work.Stance_Tot(1,1,j);
    StanceTotalTable(k).y5 = Total(5).UMB_Work.Stance_Tot(1,1,j);
    
    SwingTotalTable(k).Subject = Subjects(j).name;
    SwingTotalTable(k).Joint = 'Total';
    SwingTotalTable(k).y1 = Total(1).UMB_Work.Swing_Tot(1,1,j);
    SwingTotalTable(k).y2 = Total(2).UMB_Work.Swing_Tot(1,1,j);
    SwingTotalTable(k).y3 = Total(3).UMB_Work.Swing_Tot(1,1,j);
    SwingTotalTable(k).y4 = Total(4).UMB_Work.Swing_Tot(1,1,j);
    SwingTotalTable(k).y5 = Total(5).UMB_Work.Swing_Tot(1,1,j);
    
    StrideTotalTable(k).Subject = Subjects(j).name;
    StrideTotalTable(k).Joint = 'Total';
    StrideTotalTable(k).y1 = Total(1).UMB_Work.Stride_Tot(1,1,j);
    StrideTotalTable(k).y2 = Total(2).UMB_Work.Stride_Tot(1,1,j);
    StrideTotalTable(k).y3 = Total(3).UMB_Work.Stride_Tot(1,1,j);
    StrideTotalTable(k).y4 = Total(4).UMB_Work.Stride_Tot(1,1,j);
    StrideTotalTable(k).y5 = Total(5).UMB_Work.Stride_Tot(1,1,j);
    
    % Hip
    DS1HipTable(k).y1 = Hip(1).UMB_Work.DS1_Tot(1,1,j);
    DS1HipTable(k).y2 = Hip(2).UMB_Work.DS1_Tot(1,1,j);
    DS1HipTable(k).y3 = Hip(3).UMB_Work.DS1_Tot(1,1,j);
    DS1HipTable(k).y4 = Hip(4).UMB_Work.DS1_Tot(1,1,j);
    DS1HipTable(k).y5 = Hip(5).UMB_Work.DS1_Tot(1,1,j);
    
    SingSupHipTable(k).Subject = Subjects(j).name;
    SingSupHipTable(k).Joint = 'Hip';
    SingSupHipTable(k).y1 = Hip(1).UMB_Work.SingSup_Tot(1,1,j);
    SingSupHipTable(k).y2 = Hip(2).UMB_Work.SingSup_Tot(1,1,j);
    SingSupHipTable(k).y3 = Hip(3).UMB_Work.SingSup_Tot(1,1,j);
    SingSupHipTable(k).y4 = Hip(4).UMB_Work.SingSup_Tot(1,1,j);
    SingSupHipTable(k).y5 = Hip(5).UMB_Work.SingSup_Tot(1,1,j);
    
    DS2HipTable(k).Subject = Subjects(j).name;
    DS2HipTable(k).Joint = 'Hip';
    DS2HipTable(k).y1 = Hip(1).UMB_Work.DS2_Tot(1,1,j);
    DS2HipTable(k).y2 = Hip(2).UMB_Work.DS2_Tot(1,1,j);
    DS2HipTable(k).y3 = Hip(3).UMB_Work.DS2_Tot(1,1,j);
    DS2HipTable(k).y4 = Hip(4).UMB_Work.DS2_Tot(1,1,j);
    DS2HipTable(k).y5 = Hip(5).UMB_Work.DS2_Tot(1,1,j);
    
    StanceHipTable(k).Subject = Subjects(j).name;
    StanceHipTable(k).Joint = 'Hip';
    StanceHipTable(k).y1 = Hip(1).UMB_Work.Stance_Tot(1,1,j);
    StanceHipTable(k).y2 = Hip(2).UMB_Work.Stance_Tot(1,1,j);
    StanceHipTable(k).y3 = Hip(3).UMB_Work.Stance_Tot(1,1,j);
    StanceHipTable(k).y4 = Hip(4).UMB_Work.Stance_Tot(1,1,j);
    StanceHipTable(k).y5 = Hip(5).UMB_Work.Stance_Tot(1,1,j);
    
    SwingHipTable(k).Subject = Subjects(j).name;
    SwingHipTable(k).Joint = 'Hip';
    SwingHipTable(k).y1 = Hip(1).UMB_Work.Swing_Tot(1,1,j);
    SwingHipTable(k).y2 = Hip(2).UMB_Work.Swing_Tot(1,1,j);
    SwingHipTable(k).y3 = Hip(3).UMB_Work.Swing_Tot(1,1,j);
    SwingHipTable(k).y4 = Hip(4).UMB_Work.Swing_Tot(1,1,j);
    SwingHipTable(k).y5 = Hip(5).UMB_Work.Swing_Tot(1,1,j);
    
    StrideHipTable(k).Subject = Subjects(j).name;
    StrideHipTable(k).Joint = 'Hip';
    StrideHipTable(k).y1 = Hip(1).UMB_Work.Stride_Tot(1,1,j);
    StrideHipTable(k).y2 = Hip(2).UMB_Work.Stride_Tot(1,1,j);
    StrideHipTable(k).y3 = Hip(3).UMB_Work.Stride_Tot(1,1,j);
    StrideHipTable(k).y4 = Hip(4).UMB_Work.Stride_Tot(1,1,j);
    StrideHipTable(k).y5 = Hip(5).UMB_Work.Stride_Tot(1,1,j);
    
    % knee
    DS1KneeTable(k).y1 = Knee(1).UMB_Work.DS1_Tot(1,1,j);
    DS1KneeTable(k).y2 = Knee(2).UMB_Work.DS1_Tot(1,1,j);
    DS1KneeTable(k).y3 = Knee(3).UMB_Work.DS1_Tot(1,1,j);
    DS1KneeTable(k).y4 = Knee(4).UMB_Work.DS1_Tot(1,1,j);
    DS1KneeTable(k).y5 = Knee(5).UMB_Work.DS1_Tot(1,1,j);
    
    SingSupKneeTable(k).Subject = Subjects(j).name;
    SingSupKneeTable(k).Joint = 'Knee';
    SingSupKneeTable(k).y1 = Knee(1).UMB_Work.SingSup_Tot(1,1,j);
    SingSupKneeTable(k).y2 = Knee(2).UMB_Work.SingSup_Tot(1,1,j);
    SingSupKneeTable(k).y3 = Knee(3).UMB_Work.SingSup_Tot(1,1,j);
    SingSupKneeTable(k).y4 = Knee(4).UMB_Work.SingSup_Tot(1,1,j);
    SingSupKneeTable(k).y5 = Knee(5).UMB_Work.SingSup_Tot(1,1,j);
    
    DS2KneeTable(k).Subject = Subjects(j).name;
    DS2KneeTable(k).Joint = 'Knee';
    DS2KneeTable(k).y1 = Knee(1).UMB_Work.DS2_Tot(1,1,j);
    DS2KneeTable(k).y2 = Knee(2).UMB_Work.DS2_Tot(1,1,j);
    DS2KneeTable(k).y3 = Knee(3).UMB_Work.DS2_Tot(1,1,j);
    DS2KneeTable(k).y4 = Knee(4).UMB_Work.DS2_Tot(1,1,j);
    DS2KneeTable(k).y5 = Knee(5).UMB_Work.DS2_Tot(1,1,j);
    
    StanceKneeTable(k).Subject = Subjects(j).name;
    StanceKneeTable(k).Joint = 'Knee';
    StanceKneeTable(k).y1 = Knee(1).UMB_Work.Stance_Tot(1,1,j);
    StanceKneeTable(k).y2 = Knee(2).UMB_Work.Stance_Tot(1,1,j);
    StanceKneeTable(k).y3 = Knee(3).UMB_Work.Stance_Tot(1,1,j);
    StanceKneeTable(k).y4 = Knee(4).UMB_Work.Stance_Tot(1,1,j);
    StanceKneeTable(k).y5 = Knee(5).UMB_Work.Stance_Tot(1,1,j);
    
    SwingKneeTable(k).Subject = Subjects(j).name;
    SwingKneeTable(k).Joint = 'Knee';
    SwingKneeTable(k).y1 = Knee(1).UMB_Work.Swing_Tot(1,1,j);
    SwingKneeTable(k).y2 = Knee(2).UMB_Work.Swing_Tot(1,1,j);
    SwingKneeTable(k).y3 = Knee(3).UMB_Work.Swing_Tot(1,1,j);
    SwingKneeTable(k).y4 = Knee(4).UMB_Work.Swing_Tot(1,1,j);
    SwingKneeTable(k).y5 = Knee(5).UMB_Work.Swing_Tot(1,1,j);
    
    StrideKneeTable(k).Subject = Subjects(j).name;
    StrideKneeTable(k).Joint = 'Knee';
    StrideKneeTable(k).y1 = Knee(1).UMB_Work.Stride_Tot(1,1,j);
    StrideKneeTable(k).y2 = Knee(2).UMB_Work.Stride_Tot(1,1,j);
    StrideKneeTable(k).y3 = Knee(3).UMB_Work.Stride_Tot(1,1,j);
    StrideKneeTable(k).y4 = Knee(4).UMB_Work.Stride_Tot(1,1,j);
    StrideKneeTable(k).y5 = Knee(5).UMB_Work.Stride_Tot(1,1,j);
    
    % ankle
    DS1AnkleTable(k).Subject = Subjects(j).name;
    DS1AnkleTable(k).Joint = 'Ankle';
    DS1AnkleTable(k).y1 = Ankle(1).UMB_Work.DS1_Tot(1,1,j);
    DS1AnkleTable(k).y2 = Ankle(2).UMB_Work.DS1_Tot(1,1,j);
    DS1AnkleTable(k).y3 = Ankle(3).UMB_Work.DS1_Tot(1,1,j);
    DS1AnkleTable(k).y4 = Ankle(4).UMB_Work.DS1_Tot(1,1,j);
    DS1AnkleTable(k).y5 = Ankle(5).UMB_Work.DS1_Tot(1,1,j);
    
    SingSupAnkleTable(k).Subject = Subjects(j).name;
    SingSupAnkleTable(k).Joint = 'Ankle';
    SingSupAnkleTable(k).y1 = Ankle(1).UMB_Work.SingSup_Tot(1,1,j);
    SingSupAnkleTable(k).y2 = Ankle(2).UMB_Work.SingSup_Tot(1,1,j);
    SingSupAnkleTable(k).y3 = Ankle(3).UMB_Work.SingSup_Tot(1,1,j);
    SingSupAnkleTable(k).y4 = Ankle(4).UMB_Work.SingSup_Tot(1,1,j);
    SingSupAnkleTable(k).y5 = Ankle(5).UMB_Work.SingSup_Tot(1,1,j);
    
    DS2AnkleTable(k).Subject = Subjects(j).name;
    DS2AnkleTable(k).Joint = 'Ankle';
    DS2AnkleTable(k).y1 = Ankle(1).UMB_Work.DS2_Tot(1,1,j);
    DS2AnkleTable(k).y2 = Ankle(2).UMB_Work.DS2_Tot(1,1,j);
    DS2AnkleTable(k).y3 = Ankle(3).UMB_Work.DS2_Tot(1,1,j);
    DS2AnkleTable(k).y4 = Ankle(4).UMB_Work.DS2_Tot(1,1,j);
    DS2AnkleTable(k).y5 = Ankle(5).UMB_Work.DS2_Tot(1,1,j);
    
    StanceAnkleTable(k).Subject = Subjects(j).name;
    StanceAnkleTable(k).Joint = 'Ankle';
    StanceAnkleTable(k).y1 = Ankle(1).UMB_Work.Stance_Tot(1,1,j);
    StanceAnkleTable(k).y2 = Ankle(2).UMB_Work.Stance_Tot(1,1,j);
    StanceAnkleTable(k).y3 = Ankle(3).UMB_Work.Stance_Tot(1,1,j);
    StanceAnkleTable(k).y4 = Ankle(4).UMB_Work.Stance_Tot(1,1,j);
    StanceAnkleTable(k).y5 = Ankle(5).UMB_Work.Stance_Tot(1,1,j);
    
    SwingAnkleTable(k).Subject = Subjects(j).name;
    SwingAnkleTable(k).Joint = 'Ankle';
    SwingAnkleTable(k).y1 = Ankle(1).UMB_Work.Swing_Tot(1,1,j);
    SwingAnkleTable(k).y2 = Ankle(2).UMB_Work.Swing_Tot(1,1,j);
    SwingAnkleTable(k).y3 = Ankle(3).UMB_Work.Swing_Tot(1,1,j);
    SwingAnkleTable(k).y4 = Ankle(4).UMB_Work.Swing_Tot(1,1,j);
    SwingAnkleTable(k).y5 = Ankle(5).UMB_Work.Swing_Tot(1,1,j);
    
    StrideAnkleTable(k).Subject = Subjects(j).name;
    StrideAnkleTable(k).Joint = 'Ankle';
    StrideAnkleTable(k).y1 = Ankle(1).UMB_Work.Stride_Tot(1,1,j);
    StrideAnkleTable(k).y2 = Ankle(2).UMB_Work.Stride_Tot(1,1,j);
    StrideAnkleTable(k).y3 = Ankle(3).UMB_Work.Stride_Tot(1,1,j);
    StrideAnkleTable(k).y4 = Ankle(4).UMB_Work.Stride_Tot(1,1,j);
    StrideAnkleTable(k).y5 = Ankle(5).UMB_Work.Stride_Tot(1,1,j);
    
    k = k+1;
end

% Total
DS1Totaltable = struct2table(DS1TotalTable);
[rm.DS1Total] = fitrm(DS1Totaltable, 'y1-y5 ~ 1');
ranovatbl.DS1Total = ranova(rm.DS1Total);
MC.DS1Total = multcompare(rm.DS1Total, 'Time', 'ComparisonType', 'lsd');

SingSupTotaltable = struct2table(SingSupTotalTable);
[rm.SingSupTotal] = fitrm(SingSupTotaltable, 'y1-y5~1');
[ranovatbl.SingSupTotal] = ranova(rm.SingSupTotal);
MC.SingSupTotal = multcompare(rm.SingSupTotal, 'Time', 'ComparisonType', 'lsd');

DS2Totaltable = struct2table(DS2TotalTable);
[rm.DS2Total] = fitrm(DS2Totaltable, 'y1-y5~1');
[ranovatbl.DS2Total] = ranova(rm.DS2Total);
MC.DS2Total = multcompare(rm.DS2Total, 'Time', 'ComparisonType', 'lsd');

SwingTotaltable = struct2table(SwingTotalTable);
[rm.SwingTotal] = fitrm(SwingTotaltable, 'y1-y5~1');
[ranovatbl.SwingTotal] = ranova(rm.SwingTotal);
MC.SwingTotal = multcompare(rm.SwingTotal, 'Time', 'ComparisonType', 'lsd');

StanceTotaltable = struct2table(StanceTotalTable);
[rm.StanceTotal] = fitrm(StanceTotaltable, 'y1-y5~1');
[ranovatbl.StanceTotal] = ranova(rm.StanceTotal);
MC.StanceTotal = multcompare(rm.StanceTotal, 'Time', 'ComparisonType', 'lsd');

StrideTotaltable = struct2table(StrideTotalTable);
[rm.StrideTotal] = fitrm(StrideTotaltable, 'y1-y5~1');
[ranovatbl.StrideTotal] = ranova(rm.StrideTotal);
MC.StrideTotal = multcompare(rm.StrideTotal, 'Time', 'ComparisonType', 'lsd');

% Hip
DS1Hiptable = struct2table(DS1HipTable);
[rm.DS1Hip] = fitrm(DS1Hiptable, 'y1-y5 ~ 1');
ranovatbl.DS1Hip = ranova(rm.DS1Hip);
MC.DS1Hip = multcompare(rm.DS1Hip, 'Time', 'ComparisonType', 'lsd');

SingSupHiptable = struct2table(SingSupHipTable);
[rm.SingSupHip] = fitrm(SingSupHiptable, 'y1-y5~1');
[ranovatbl.SingSupHip] = ranova(rm.SingSupHip);
MC.SingSupHip = multcompare(rm.SingSupHip, 'Time', 'ComparisonType', 'lsd');

DS2Hiptable = struct2table(DS2HipTable);
[rm.DS2Hip] = fitrm(DS2Hiptable, 'y1-y5~1');
[ranovatbl.DS2Hip] = ranova(rm.DS2Hip);
MC.DS2Hip = multcompare(rm.DS2Hip, 'Time', 'ComparisonType', 'lsd');

SwingHiptable = struct2table(SwingHipTable);
[rm.SwingHip] = fitrm(SwingHiptable, 'y1-y5~1');
[ranovatbl.SwingHip] = ranova(rm.SwingHip);
MC.SwingHip = multcompare(rm.SwingHip, 'Time', 'ComparisonType', 'lsd');

StanceHiptable = struct2table(StanceHipTable);
[rm.StanceHip] = fitrm(StanceHiptable, 'y1-y5~1');
[ranovatbl.StanceHip] = ranova(rm.StanceHip);
MC.StanceHip = multcompare(rm.StanceHip, 'Time', 'ComparisonType', 'lsd');

StrideHiptable = struct2table(StrideHipTable);
[rm.StrideHip] = fitrm(StrideHiptable, 'y1-y5~1');
[ranovatbl.StrideHip] = ranova(rm.StrideHip);
MC.StrideHip = multcompare(rm.StrideHip, 'Time', 'ComparisonType', 'lsd');

% Knee
DS1Kneetable = struct2table(DS1KneeTable);
[rm.DS1Knee] = fitrm(DS1Kneetable, 'y1-y5 ~ 1');
ranovatbl.DS1Knee = ranova(rm.DS1Knee);
MC.DS1Knee = multcompare(rm.DS1Knee, 'Time', 'ComparisonType', 'lsd');

SingSupKneetable = struct2table(SingSupKneeTable);
[rm.SingSupKnee] = fitrm(SingSupKneetable, 'y1-y5~1');
[ranovatbl.SingSupKnee] = ranova(rm.SingSupKnee);
MC.SingSupKnee = multcompare(rm.SingSupKnee, 'Time', 'ComparisonType', 'lsd');

DS2Kneetable = struct2table(DS2KneeTable);
[rm.DS2Knee] = fitrm(DS2Kneetable, 'y1-y5~1');
[ranovatbl.DS2Knee] = ranova(rm.DS2Knee);
MC.DS2Knee = multcompare(rm.DS2Knee, 'Time', 'ComparisonType', 'lsd');

SwingKneetable = struct2table(SwingKneeTable);
[rm.SwingKnee] = fitrm(SwingKneetable, 'y1-y5~1');
[ranovatbl.SwingKnee] = ranova(rm.SwingKnee);
MC.SwingKnee = multcompare(rm.SwingKnee, 'Time', 'ComparisonType', 'lsd');

StanceKneetable = struct2table(StanceKneeTable);
[rm.StanceKnee] = fitrm(StanceKneetable, 'y1-y5~1');
[ranovatbl.StanceKnee] = ranova(rm.StanceKnee);
MC.StanceKnee = multcompare(rm.StanceKnee, 'Time', 'ComparisonType', 'lsd');

StrideKneetable = struct2table(StrideKneeTable);
[rm.StrideKnee] = fitrm(StrideKneetable, 'y1-y5~1');
[ranovatbl.StrideKnee] = ranova(rm.StrideKnee);
MC.StrideKnee = multcompare(rm.StrideKnee, 'Time', 'ComparisonType', 'lsd');

% Ankle
DS1Ankletable = struct2table(DS1AnkleTable);
[rm.DS1Ankle] = fitrm(DS1Ankletable, 'y1-y5~1');
[ranovatbl.DS1Ankle] = ranova(rm.DS1Ankle);
MC.DS1Ankle = multcompare(rm.DS1Ankle, 'Time', 'ComparisonType', 'lsd');

SingSupAnkletable = struct2table(SingSupAnkleTable);
[rm.SingSupAnkle] = fitrm(SingSupAnkletable, 'y1-y5~1');
[ranovatbl.SingSupAnkle] = ranova(rm.SingSupAnkle);
MC.SingSupAnkle = multcompare(rm.SingSupAnkle, 'Time', 'ComparisonType', 'lsd');

DS2Ankletable = struct2table(DS2AnkleTable);
[rm.DS2Ankle] = fitrm(DS2Ankletable, 'y1-y5~1');
[ranovatbl.DS2Ankle] = ranova(rm.DS2Ankle);
MC.DS2Ankle = multcompare(rm.DS2Ankle, 'Time', 'ComparisonType', 'lsd');

SwingAnkletable = struct2table(SwingAnkleTable);
[rm.SwingAnkle] = fitrm(SwingAnkletable, 'y1-y5~1');
[ranovatbl.SwingAnkle] = ranova(rm.SwingAnkle);
MC.SwingAnkle = multcompare(rm.SwingAnkle, 'Time', 'ComparisonType', 'lsd');

StanceAnkletable = struct2table(StanceAnkleTable);
[rm.StanceAnkle] = fitrm(StanceAnkletable, 'y1-y5~1');
[ranovatbl.StanceAnkle] = ranova(rm.StanceAnkle);
MC.StanceAnkle = multcompare(rm.StanceAnkle, 'Time', 'ComparisonType', 'lsd');

StrideAnkletable = struct2table(StrideAnkleTable);
[rm.StrideAnkle] = fitrm(StrideAnkletable, 'y1-y5~1');
[ranovatbl.StrideAnkle] = ranova(rm.StrideAnkle);
MC.StrideAnkle = multcompare(rm.StrideAnkle, 'Time', 'ComparisonType', 'lsd');



%% Effect Size Statistics
% eta squared
cols = {'y1','y2','y3','y4','y5'}; 
% total
MES.DS1Total = mes1way(DS1Totaltable{:,cols} , 'eta2' ); 
MES.DS2Total = mes1way(DS2Totaltable{:,cols} , 'eta2' ); 
MES.SingSupTotal = mes1way(SingSupTotaltable{:,cols} , 'eta2' ); 
MES.SwingTotal = mes1way(SwingTotaltable{:,cols} , 'eta2' ); 
MES.StrideTotal = mes1way(StrideTotaltable{:,cols} , 'eta2' ); 

% hip
MES.DS1Hip = mes1way(DS1Hiptable{:,cols} , 'eta2' ); 
MES.DS2Hip = mes1way(DS2Hiptable{:,cols} , 'eta2' ); 
MES.SingSupHip = mes1way(SingSupHiptable{:,cols} , 'eta2' ); 
MES.SwingHip = mes1way(SwingHiptable{:,cols} , 'eta2' ); 
MES.StrideHip = mes1way(StrideHiptable{:,cols} , 'eta2' ); 

% knee
MES.DS1Knee = mes1way(DS1Kneetable{:,cols} , 'eta2' ); 
MES.DS2Knee = mes1way(DS2Kneetable{:,cols} , 'eta2' ); 
MES.SingSupKnee = mes1way(SingSupKneetable{:,cols} , 'eta2' ); 
MES.SwingKnee = mes1way(SwingKneetable{:,cols} , 'eta2' ); 
MES.StrideKnee = mes1way(StrideKneetable{:,cols} , 'eta2' ); 

% ankle
MES.DS1Ankle = mes1way(DS1Ankletable{:,cols} , 'eta2' ); 
MES.DS2Ankle = mes1way(DS2Ankletable{:,cols} , 'eta2' ); 
MES.SingSupAnkle = mes1way(SingSupAnkletable{:,cols} , 'eta2' ); 
MES.SwingAnkle = mes1way(SwingAnkletable{:,cols} , 'eta2' ); 
MES.StrideAnkle = mes1way(StrideAnkletable{:,cols} , 'eta2' ); 

% cohens D
cols = {'y1','y2','y4','y5'}; 
Norm = 'y3';
% total
MC.DS1Total_Cd = (nanmean(DS1Totaltable{:,cols}) - nanmean(DS1Totaltable{:,Norm}) )...
    / nanstd(DS1Totaltable{:,Norm}); 
MC.SingSupTotal_Cd = (nanmean(SingSupTotaltable{:,cols}) - nanmean(SingSupTotaltable{:,Norm}) )...
    / nanstd(SingSupTotaltable{:,Norm}); 
MC.DS2Total_Cd = (nanmean(DS2Totaltable{:,cols}) - nanmean(DS2Totaltable{:,Norm}) )...
    / nanstd(DS2Totaltable{:,Norm}); 
MC.SwingTotal_Cd = (nanmean(SwingTotaltable{:,cols}) - nanmean(SwingTotaltable{:,Norm}) )...
    / nanstd(SwingTotaltable{:,Norm}); 
MC.StrideTotal_Cd = (nanmean(StrideTotaltable{:,cols}) - nanmean(StrideTotaltable{:,Norm}) )...
    / nanstd(StrideTotaltable{:,Norm}); 

% hip
MC.DS1Hip_Cd = (nanmean(DS1Hiptable{:,cols}) - nanmean(DS1Hiptable{:,Norm}) )...
    / nanstd(DS1Hiptable{:,Norm}); 
MC.SingSupHip_Cd = (nanmean(SingSupHiptable{:,cols}) - nanmean(SingSupHiptable{:,Norm}) )...
    / nanstd(SingSupHiptable{:,Norm}); 
MC.DS2Hip_Cd = (nanmean(DS2Hiptable{:,cols}) - nanmean(DS2Hiptable{:,Norm}) )...
    / nanstd(DS2Hiptable{:,Norm}); 
MC.SwingHip_Cd = (nanmean(SwingHiptable{:,cols}) - nanmean(SwingHiptable{:,Norm}) )...
    / nanstd(SwingHiptable{:,Norm}); 
MC.StrideHip_Cd = (nanmean(StrideHiptable{:,cols}) - nanmean(StrideHiptable{:,Norm}) )...
    / nanstd(StrideHiptable{:,Norm}); 

% knee
MC.DS1Knee_Cd = (nanmean(DS1Kneetable{:,cols}) - nanmean(DS1Kneetable{:,Norm}) )...
    / nanstd(DS1Kneetable{:,Norm}); 
MC.SingSupKnee_Cd = (nanmean(SingSupKneetable{:,cols}) - nanmean(SingSupKneetable{:,Norm}) )...
    / nanstd(SingSupKneetable{:,Norm}); 
MC.DS2Knee_Cd = (nanmean(DS2Kneetable{:,cols}) - nanmean(DS2Kneetable{:,Norm}) )...
    / nanstd(DS2Kneetable{:,Norm}); 
MC.SwingKnee_Cd = (nanmean(SwingKneetable{:,cols}) - nanmean(SwingKneetable{:,Norm}) )...
    / nanstd(SwingKneetable{:,Norm}); 
MC.StrideKnee_Cd = (nanmean(StrideKneetable{:,cols}) - nanmean(StrideKneetable{:,Norm}) )...
    / nanstd(StrideKneetable{:,Norm}); 

% ankle
MC.DS1Ankle_Cd = (nanmean(DS1Ankletable{:,cols}) - nanmean(DS1Ankletable{:,Norm}) )...
    / nanstd(DS1Ankletable{:,Norm}); 
MC.SingSupAnkle_Cd = (nanmean(SingSupAnkletable{:,cols}) - nanmean(SingSupAnkletable{:,Norm}) )...
    / nanstd(SingSupAnkletable{:,Norm}); 
MC.DS2Ankle_Cd = (nanmean(DS2Ankletable{:,cols}) - nanmean(DS2Ankletable{:,Norm}) )...
    / nanstd(DS2Ankletable{:,Norm}); 
MC.SwingAnkle_Cd = (nanmean(SwingAnkletable{:,cols}) - nanmean(SwingAnkletable{:,Norm}) )...
    / nanstd(SwingAnkletable{:,Norm}); 
MC.StrideAnkle_Cd = (nanmean(StrideAnkletable{:,cols}) - nanmean(StrideAnkletable{:,Norm}) )...
    / nanstd(StrideAnkletable{:,Norm}); 

%% Create figure of joint metabolic costs ABSTRACT FIGURE
% clc;
% close all;
% HipColors = [rgb('PowderBlue'); rgb('DeepSkyBlue'); rgb('RoyalBlue');  rgb('MediumBlue'); rgb('MidnightBlue')];
% AnkleColors = [rgb('PaleGreen'); rgb('LimeGreen'); rgb('ForestGreen');  rgb('Green'); rgb('DarkGreen')];
% MkrSz = 15;
% SmlMkr = 8;
%
% JointEnergyCostsUMB = figure('Position',[50 50 1100 300]);
% hold on;
% x = [1 2 3 4 5];
% LW = 1.25;
% DotMargin = 0.32;
% MeanMargin = 0.12;
% SingSupOffset = 6;
% DS2Offset = 12;
% StanceOffset = 18;
% SwingOffset = 24;
% TxtFnt = 12;
% CapEndSz = 0.1;
%
% % hip
% % subplot(121)
% for i = 1:5
%
%     % initial dub sup
%     xRel = x(i)-MeanMargin;
%     h = boxplot2(Hip(i).UMB_Work.DS1_Tot,  xRel);
%     h.box.Color = HipColors(i,:);
%     h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);
%     h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);
%     h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % sing sup
%     xRel = x(i)-MeanMargin+SingSupOffset;
%     h = boxplot2(Hip(i).UMB_Work.SingSup_Tot,  xRel);
%     h.box.Color = HipColors(i,:);
%     h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);
%     h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);
%     h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % secondary dub sup
%     xRel = x(i)-MeanMargin+DS2Offset;
%     h = boxplot2(Hip(i).UMB_Work.DS2_Tot,  xRel);
%     h.box.Color = HipColors(i,:);
%     h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);
%     h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);
%     h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % stance
%     xRel = x(i)-MeanMargin+StanceOffset;
%     h = boxplot2(Hip(i).UMB_Work.Stance_Tot,  xRel);
%     h.box.Color = HipColors(i,:);
%     h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);
%     h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);
%     h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % stance
%     xRel = x(i)-MeanMargin+SwingOffset;
%     h = boxplot2(Hip(i).UMB_Work.Swing_Tot,  xRel);
%     h.box.Color = HipColors(i,:);
%     h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);
%     h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);
%     h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     for k = 1:length(Hip(i).UMB_Work.Stance_Tot)
%
%         % Umberger
%         plot(x(i)-DotMargin, Hip(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%         plot(x(i)-DotMargin+SingSupOffset , Hip(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%         plot(x(i)-DotMargin+DS2Offset , Hip(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%         plot(x(i)-DotMargin+StanceOffset, Hip(i).UMB_Work.Stance_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%         plot(x(i)-DotMargin+SwingOffset, Hip(i).UMB_Work.Swing_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%     end
% end
%
% % Ankle
% % subplot(122);
% for i = 1:5
%     % initial dub sup
%     xRel = x(i)+MeanMargin;
%     h = boxplot2(Ankle(i).UMB_Work.DS1_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);
%     h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);
%     h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);
%     h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     % single sup
%     xRel = x(i)+MeanMargin+SingSupOffset;
%     h = boxplot2(Ankle(i).UMB_Work.SingSup_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);
%     h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);
%     h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);
%     h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     % initial dub sup
%     xRel = x(i)+MeanMargin+DS2Offset;
%     h = boxplot2(Ankle(i).UMB_Work.DS2_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);
%     h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);
%     h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);
%     h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     % stance
%     xRel = x(i)+MeanMargin+StanceOffset;
%     h = boxplot2(Ankle(i).UMB_Work.Stance_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);
%     h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);
%     h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);
%     h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     % swing
%     xRel = x(i)+MeanMargin+SwingOffset;
%     h = boxplot2(Ankle(i).UMB_Work.Swing_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);
%     h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);
%     h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);
%     h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;
%     h.med.LineWidth = LW;
%     h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';
%     h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;
%     h.ladj.LineWidth = LW;
%     h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     for k = 1:length(Ankle(i).UMB_Work.Stance_Tot)
%
%         % Umberger
%         plot(x(i)+DotMargin, Ankle(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%         plot(x(i)+DotMargin+SingSupOffset , Ankle(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%         plot(x(i)+DotMargin+DS2Offset , Ankle(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%         plot(x(i)+DotMargin+StanceOffset, Ankle(i).UMB_Work.Stance_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%         plot(x(i)+DotMargin+SwingOffset, Ankle(i).UMB_Work.Swing_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%
%     end
%
% end
% ypos = 3.3;
% %     title('Metabolic Cost Per Stride Phase');
% text(3,ypos, {'Double Support', 'Leading Leg'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
% text(9,ypos+0.1, {'Single Support'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
% text(15,ypos, {'Double Support', 'Trailing Leg'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
% text(21,ypos+0.1, {'Stance Phase'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
% text(27,ypos+0.1, {'Swing Phase'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
%
% TickLabel = {'-40%', '-20%','Norm','+20%','+40%'};
% ylabel('Matabolic Power (W / kg)');
% xlim([0 30]);
% ylim([0 3.5]);
% vline(18, 'k');
% ax = gca;
% ax.XTick = [1 2 3 4 5 7 8 9 10 11 13 14 15 16 17 19 20 21 22 23 25 26 27 28 29];
% ax.XTickLabels = TickLabel;
% ax.YTick = [0 0.5 1 1.5 2 2.5 3 3.5];
% ax.YTickLabels = {'0', '0.5', '1.0' , '1.5', '2.0', '2.5', '3.0', '3.5'};
% ax.FontSize = 9;
%
% text(1, 2.75, '\bf HIP MUSCULATURE', 'Color', rgb('RoyalBlue'), 'FontSize',12);
% text(1, 2.5, '\bf ANKLE MUSCULATURE', 'Color', rgb('ForestGreen'), 'FontSize',12);
%
% % add asterisks for sig diff conditions and % differences for
% % HipHt = 2.5;
% HipPos = 3;
% AnklePos = 3;
% Asterisk = 20;
% SigFnt = 8;
% Rounding = 2;
% AsteriskMargin = mean([MeanMargin, DotMargin]);
% % hip
% if ranovatbl.DS1Hip.pValue(1) < 0.05
%     % -40 condition
%     text(1-AsteriskMargin, 2.0, '*', 'Color', HipColors(1,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(1).UMB_Work.DS1_Tot]) - mean([Hip(3).UMB_Work.DS1_Tot])), Rounding);
%     %        text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     % -20 condition
%     text(2-AsteriskMargin, 2.0, '*', 'Color', HipColors(2,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(2).UMB_Work.DS1_Tot]) - mean([Hip(3).UMB_Work.DS1_Tot])), Rounding);
%     %        text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(2-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     % + 20 condition
%     text(4-AsteriskMargin, 2.0, '*', 'Color', HipColors(4,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(4).UMB_Work.DS1_Tot]) - mean([Hip(3).UMB_Work.DS1_Tot])), Rounding);
%     %        text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(4-AsteriskMargin, 1.95, [num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
% if ranovatbl.SingSupHip.pValue(1) < 0.05
%     % + 40 condition
%     text(11-AsteriskMargin, 2.5, '*', 'Color', HipColors(5,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(5).UMB_Work.DS1_Tot]) - mean([Hip(3).UMB_Work.DS1_Tot])), Rounding);
%     %        text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(11-AsteriskMargin, 2.45, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
% if ranovatbl.DS2Hip.pValue(1) < 0.05
%     % + 20 condition
%     text(16-AsteriskMargin, 2.0, '*', 'Color', HipColors(4,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(4).UMB_Work.DS1_Tot]) - mean([Hip(3).UMB_Work.DS1_Tot])), Rounding);
%     %        text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(16-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     % + 40 condition
%     text(17-AsteriskMargin, 2.5, '*', 'Color', HipColors(5,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(5).UMB_Work.DS1_Tot]) - mean([Hip(3).UMB_Work.DS1_Tot])), Rounding);
%     %        text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(17-AsteriskMargin, 2.45, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
% if ranovatbl.StanceHip.pValue(1) < 0.05
%     % -40 condition
%     text(19-AsteriskMargin, 3, '*', 'Color', HipColors(1,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(1).UMB_Work.Stance_Tot]) - mean([Hip(3).UMB_Work.Stance_Tot])), Rounding);
%     %        text(19-AsteriskMargin, 2.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(19-AsteriskMargin, 2.95, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%
%     % +20 condition
%     text(22-AsteriskMargin, 3, '*', 'Color', HipColors(4,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(4).UMB_Work.Stance_Tot]) - mean([Hip(3).UMB_Work.Stance_Tot])), Rounding);
%     %        text(19-AsteriskMargin, 2.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(22-AsteriskMargin, 2.95, [num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%
%     %          text(22-AsteriskMargin, 3, '*', 'Color', HipColors(4,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     %          diff = round((mean([Hip(4).UMB_Work.Stance_Tot]) - mean([Hip(3).UMB_Work.Stance_Tot])), Rounding);
%     %        text(22-AsteriskMargin, 2.95, [num2str(diff)], 'Color', HipColors(4,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
% if ranovatbl.SwingHip.pValue(1) < 0.05
%     % + 40 condition
%     text(29-AsteriskMargin, 2.5, '*', 'Color', HipColors(5,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Hip(5).UMB_Work.DS1_Tot]) - mean([Hip(3).UMB_Work.DS1_Tot])), Rounding);
%     %        text(1-AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', HipColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(29-AsteriskMargin, 2.45, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
%
% % Ankle
% if ranovatbl.DS1Ankle.pValue(1) < 0.05
%     %     text(AnklePos, 1.35, '*', 'Color', rgb('ForestGreen'), 'FontSize',20);
% end
% if ranovatbl.SingSupAnkle.pValue(1) < 0.05
%     % + 20 condition
%     text(10+AsteriskMargin, 2, '*', 'Color', AnkleColors(4,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Ankle(4).UMB_Work.SingSup_Tot]) - mean([Ankle(3).UMB_Work.SingSup_Tot])), Rounding);
%     %        text(10+AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', AnkleColors(4,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %    text(10+AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     % + 40 condition
%     text(11+AsteriskMargin, 2, '*', 'Color', AnkleColors(5,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Ankle(5).UMB_Work.SingSup_Tot]) - mean([Ankle(3).UMB_Work.SingSup_Tot])), Rounding);
%     %        text(11+AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', AnkleColors(5,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %   text(11+AsteriskMargin, 1.95, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
% if ranovatbl.DS2Ankle.pValue(1) < 0.05
%     %     -40 condition
%     text(13+AsteriskMargin, 2, '*', 'Color', AnkleColors(1,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Ankle(1).UMB_Work.DS2_Tot]) - mean([Ankle(3).UMB_Work.DS2_Tot])), Rounding);
%     %        text(13+AsteriskMargin, 1.95, [num2str(diff)], 'Color', AnkleColors(1,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %     text(13+AsteriskMargin, 1.95, [num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
% if ranovatbl.StanceAnkle.pValue(1) < 0.05
%     %  +20 condition
%     text(22+AsteriskMargin, 0.5, '*', 'Color', AnkleColors(4,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Ankle(4).UMB_Work.DS2_Tot]) - mean([Ankle(3).UMB_Work.DS2_Tot])), Rounding);
%     %        text(22+AsteriskMargin, 0.45, ['+' num2str(diff)], 'Color', AnkleColors(4,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %   text(22+AsteriskMargin, 0.45, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     % + 40 condition
%     text(23+AsteriskMargin, 0.5, '*', 'Color', AnkleColors(5,:), 'FontSize',Asterisk, 'HorizontalAlignment', 'center');
%     diff = round((mean([Ankle(5).UMB_Work.DS2_Tot]) - mean([Ankle(3).UMB_Work.DS2_Tot])), Rounding);
%     %        text(23+AsteriskMargin, 0.45, ['+' num2str(diff)], 'Color', AnkleColors(5,:), 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
%     %  text(23+AsteriskMargin, 0.45, ['+' num2str(diff)], 'Color', 'k', 'FontSize',SigFnt, 'HorizontalAlignment', 'center');
% end
% if ranovatbl.SwingAnkle.pValue(1) < 0.05
%     %     text(AnklePos+24, 1.8, '*', 'Color', rgb('ForestGreen'), 'FontSize',Asterisk);
% end
%
% subplotsqueeze(JointEnergyCostsUMB, 1.1);
% % saveas(JointEnergyCostsUMB, 'JointEnergyCostsUMB.png');
% % saveas(JointEnergyCostsUMB, 'JointEnergyCostsUMB.pdf');
%
% clearvars i j k L m n o p Raw row h ax c diff x X

%% Stride-Stance-Swing Metabolic Costs
clc;
close all;

TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray');  rgb('DarkSlateGray'); rgb('Black')];
HipColors = [rgb('PowderBlue'); rgb('DeepSkyBlue'); rgb('RoyalBlue');  rgb('Navy'); rgb('MidnightBlue')];
% KneeColors = [rgb('Orange'); rgb('DarkOrange'); rgb('Coral');  rgb('Chocolate'); rgb('Sienna')];
KneeColors = [rgb('SandyBrown'); rgb('Orange'); rgb('DarkOrange'); rgb('Chocolate'); rgb('Sienna')];
AnkleColors = [rgb('PaleGreen'); rgb('LimeGreen'); rgb('ForestGreen');  rgb('Green'); rgb('DarkGreen')];

CycleEnergyCostsUMB = figure('Position',[50 50 1000 800]);
hold on;
x = [1 2 3 4 5];
LW = 1.25;
DotMargin = 0.2;
MeanMargin = 0;
StrideOffset = 0;
% StanceOffset = 6;
DS1Offset = 6;
SingSupOffset = 12;
DS2Offset = 18;
SwingOffset = 24;

TxtFnt = 12;
CapEndSz = 0.1;
SmlMkr = 7;

% Total
subplot(411); hold on;
for i = 1:5
    % stride
    xRel = x(i)-MeanMargin+StrideOffset;
    h = boxplot2(Total(i).UMB_Work.Stride_Tot,  xRel);
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % initial dub sup
    xRel = x(i)-MeanMargin + DS1Offset;
    h = boxplot2(Total(i).UMB_Work.DS1_Tot,  xRel);
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % sing sup
    xRel = x(i)-MeanMargin + SingSupOffset;
    h = boxplot2(Total(i).UMB_Work.SingSup_Tot,  xRel);
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % secondary dub sup
    xRel = x(i)-MeanMargin+DS2Offset;
    h = boxplot2(Total(i).UMB_Work.DS2_Tot,  xRel);
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % swing
    xRel = x(i)-MeanMargin+SwingOffset;
    h = boxplot2(Total(i).UMB_Work.Swing_Tot,  xRel);
    h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
    h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    
    
    for k = 1:length(Total(i).UMB_Work.Stance_Tot)
        % Umberger
        plot(x(i)-DotMargin+StrideOffset, Total(i).UMB_Work.Stride_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
        plot(x(i)-DotMargin+DS1Offset, Total(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
        plot(x(i)-DotMargin+SingSupOffset , Total(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
        plot(x(i)-DotMargin+DS2Offset , Total(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
        %             plot(x(i)-DotMargin+StanceOffset, Total(i).UMB_Work.Stance_Tot(k), '.', 'MarkerSize',SmlMkr,...
        %             'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
        plot(x(i)-DotMargin+SwingOffset, Total(i).UMB_Work.Swing_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
    end
end

% hip
subplot(412); hold on;
for i = 1:5
    
    % stride
    xRel = x(i)-MeanMargin+StrideOffset;
    h = boxplot2(Hip(i).UMB_Work.Stride_Tot,  xRel);
    h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
    h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % initial dub sup
    xRel = x(i)-MeanMargin+DS1Offset;
    h = boxplot2(Hip(i).UMB_Work.DS1_Tot,  xRel);
    h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
    h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % sing sup
    xRel = x(i)-MeanMargin+SingSupOffset;
    h = boxplot2(Hip(i).UMB_Work.SingSup_Tot,  xRel);
    h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
    h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % secondary dub sup
    xRel = x(i)-MeanMargin+DS2Offset;
    h = boxplot2(Hip(i).UMB_Work.DS2_Tot,  xRel);
    h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
    h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % swing
    xRel = x(i)-MeanMargin+SwingOffset;
    h = boxplot2(Hip(i).UMB_Work.Swing_Tot,  xRel);
    h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
    h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    for k = 1:length(Hip(i).UMB_Work.Stance_Tot)
        % Umberger
        plot(x(i)-DotMargin+StrideOffset, Hip(i).UMB_Work.Stride_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
        plot(x(i)-DotMargin+DS1Offset, Hip(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
        plot(x(i)-DotMargin+SingSupOffset , Hip(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
        plot(x(i)-DotMargin+DS2Offset , Hip(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
        %        plot(x(i)-DotMargin+StanceOffset, Hip(i).UMB_Work.Stance_Tot(k), '.', 'MarkerSize',SmlMkr,...
        %             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
        plot(x(i)-DotMargin+SwingOffset, Hip(i).UMB_Work.Swing_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
    end
end

% Knee
subplot(413); hold on;
for i = 1:5
    % stride
    xRel = x(i)-MeanMargin+StrideOffset;
    h = boxplot2(Knee(i).UMB_Work.Stride_Tot,  xRel);
    h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
    h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % initial dub sup
    xRel = x(i)-MeanMargin+DS1Offset;
    h = boxplot2(Knee(i).UMB_Work.DS1_Tot,  xRel);
    h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
    h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % sing sup
    xRel = x(i)-MeanMargin+SingSupOffset;
    h = boxplot2(Knee(i).UMB_Work.SingSup_Tot,  xRel);
    h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
    h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % secondary dub sup
    xRel = x(i)-MeanMargin+DS2Offset;
    h = boxplot2(Knee(i).UMB_Work.DS2_Tot,  xRel);
    h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
    h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    % swing
    xRel = x(i)-MeanMargin+SwingOffset;
    h = boxplot2(Knee(i).UMB_Work.Swing_Tot,  xRel);
    h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
    h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';
    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;
    h.ladj.LineWidth = LW;
    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
    
    for k = 1:length(Knee(i).UMB_Work.Stance_Tot) % plot individual points
        % Umberger
        plot(x(i)-DotMargin+StrideOffset, Knee(i).UMB_Work.Stride_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
        plot(x(i)-DotMargin+DS1Offset, Knee(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
        plot(x(i)-DotMargin+SingSupOffset , Knee(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
        plot(x(i)-DotMargin+DS2Offset , Knee(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
        %            plot(x(i)-DotMargin+StanceOffset, Knee(i).UMB_Work.Stance_Tot(k), '.', 'MarkerSize',SmlMkr,...
        %             'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
        plot(x(i)-DotMargin+SwingOffset, Knee(i).UMB_Work.Swing_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
    end
end

% Ankle
subplot(414); hold on;
for i = 1:5
    % BOXPLOT
    % stride
    xRel = x(i)+MeanMargin+StrideOffset;
    h = boxplot2(Ankle(i).UMB_Work.Stride_Tot, xRel);
    h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
    h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
    
    % initial dub sup
    xRel = x(i)+MeanMargin+DS1Offset;
    h = boxplot2(Ankle(i).UMB_Work.DS1_Tot, xRel);
    h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
    h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
    
    % single sup
    xRel = x(i)+MeanMargin+SingSupOffset;
    h = boxplot2(Ankle(i).UMB_Work.SingSup_Tot, xRel);
    h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
    h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
    
    % secondary dub sup
    xRel = x(i)+MeanMargin+DS2Offset;
    h = boxplot2(Ankle(i).UMB_Work.DS2_Tot, xRel);
    h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
    h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
    
    % swing
    xRel = x(i)+MeanMargin+SwingOffset;
    h = boxplot2(Ankle(i).UMB_Work.Swing_Tot, xRel);
    h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
    h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
    h.out.Marker = 'none';
    h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
    h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
    h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
    h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;
    h.uadj.LineWidth = LW;
    h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
    h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
    
    for k = 1:length(Ankle(i).UMB_Work.Stance_Tot) % plot individual points
        % Umberger
        plot(x(i)-DotMargin+StrideOffset, Ankle(i).UMB_Work.Stride_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
        plot(x(i)-DotMargin+DS1Offset, Ankle(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
        plot(x(i)-DotMargin+SingSupOffset , Ankle(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
        plot(x(i)-DotMargin+DS2Offset , Ankle(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
        %           plot(x(i)-DotMargin+StanceOffset, Ankle(i).UMB_Work.Stance_Tot(k), '.', 'MarkerSize',SmlMkr,...
        %             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
        plot(x(i)-DotMargin+SwingOffset, Ankle(i).UMB_Work.Swing_Tot(k), '.', 'MarkerSize',SmlMkr,...
            'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
    end
    
end

% AXES and LABELS
% phase titles
subplot(4,1,1);
ypos = 21;
% text(StanceOffset + 3,ypos+0.1, {'\bf Stance Phase'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(StrideOffset + 3,ypos+0.1, {'\bf Entire Stride'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(DS1Offset + 3,ypos+0.1, {'\bf Double Support'; '\bf Leading Leg'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(SingSupOffset + 3,ypos+0.1, {'\bf Single Limb Support'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(DS2Offset + 3,ypos+0.1, {'\bf Double Support'; '\bf Trailing Leg'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(SwingOffset + 3,ypos+0.1, {'\bf Swing Phase'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');

% condition labels
for j = 1:4
    subplot(4,1,j);
    %     TickLabel = {'-40%', '-20%','Norm','+20%','+40%'};
    ylabel('W / kg');
    xlim([0.5 29.5]);
    ax = gca;
    ax.XTick = [1 2 3 4 5 7 8 9 10 11 13 14 15 16 17 19 20 21 22 23 25 26 27 28 29];
    ax.XTickLabels = [];
    ax.FontSize = 10;
end

subplot(4,1,4);
TickLabel = {'-40', '-20','Norm','+20','+40'};
ax = gca;
ax.XTick = [1 2 3 4 5 7 8 9 10 11 13 14 15 16 17 19 20 21 22 23 25 26 27 28 29];
ax.XTickLabels = TickLabel;
ax.YTick = [0 1 2 3];
ax.FontSize = 10;

% other axis
MuscLabelCent = 15;
subplot(4,1,1);
ylim([0 18]);
vline(6, 'k');
text(MuscLabelCent, 17, '\bf ALL MUSCLES', 'Color', TotalColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');
subplot(4,1,2);
ylim([0 5]);
vline(6, 'k');
text(MuscLabelCent, 4.2, '\bf HIP MUSCLES', 'Color',HipColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');
subplot(4,1,3);
ylim([0 4]);
vline(6, 'k');
text(MuscLabelCent, 2.7, '\bf KNEE MUSCLES', 'Color', KneeColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');
subplot(4,1,4);
ylim([0 3]);
vline(6, 'k');
text(MuscLabelCent, 2.2, '\bf ANKLE MUSCLES', 'Color', AnkleColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');

% label right axis as kcal \ kg
% yyaxis right;
% ylim([0 4 / 4184]);
% ylabel('kcal / kg');
% set(gca,'ycolor','k');
% set(gcf,'color','w');

% ASTERISKS
alpha = 0.05;
Asterisk = 20; 
% total
subplot(411);
TopMargin = 0.5;
% if ranovatbl.StanceTotal.pValue(1) < alpha
%     multcomp = MC.StanceTotal(MC.StanceTotal.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Total(Log).UMB_Work.Stance_Tot]) + TopMargin;
%             text(Log+StanceOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
if ranovatbl.SwingTotal.pValue(1) < alpha
    multcomp = MC.SwingTotal(MC.SwingTotal.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Total(Log).UMB_Work.Swing_Tot]) + TopMargin;
            text(Log+SwingOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.StrideTotal.pValue(1) < alpha
    multcomp = MC.StrideTotal(MC.StrideTotal.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Total(Log).UMB_Work.Stride_Tot]) + TopMargin;
            text(Log+StrideOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS1Total.pValue(1) < alpha
    multcomp = MC.DS1Total(MC.DS1Total.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Total(Log).UMB_Work.DS1_Tot]) + TopMargin;
            text(Log+DS1Offset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.SingSupTotal.pValue(1) < alpha
    multcomp = MC.SingSupTotal(MC.SingSupTotal.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Total(Log).UMB_Work.SingSup_Tot]) + TopMargin;
            text(Log+SingSupOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS2Total.pValue(1) < alpha
    multcomp = MC.DS2Total(MC.DS2Total.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Total(Log).UMB_Work.DS2_Tot]) + TopMargin;
            text(Log+DS2Offset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end

% hip
subplot(412);
TopMargin = 0.25;
% if ranovatbl.StanceHip.pValue(1) < alpha
%     multcomp = MC.StanceHip(MC.StanceHip.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Hip(Log).UMB_Work.Stance_Tot]) + TopMargin;
%             text(Log+StanceOffset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
if ranovatbl.SwingHip.pValue(1) < alpha
    multcomp = MC.SwingHip(MC.SwingHip.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Hip(Log).UMB_Work.Swing_Tot]) + TopMargin;
            text(Log+SwingOffset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.StrideHip.pValue(1) < alpha
    multcomp = MC.StrideHip(MC.StrideHip.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Hip(Log).UMB_Work.Stride_Tot]) + TopMargin;
            text(Log+StrideOffset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS1Hip.pValue(1) < alpha
    multcomp = MC.DS1Hip(MC.DS1Hip.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Hip(Log).UMB_Work.DS1_Tot]) + TopMargin;
            text(Log+DS1Offset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.SingSupHip.pValue(1) < alpha
    multcomp = MC.SingSupHip(MC.SingSupHip.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Hip(Log).UMB_Work.SingSup_Tot]) + TopMargin;
            text(Log+SingSupOffset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS2Hip.pValue(1) < alpha
    multcomp = MC.DS2Hip(MC.DS2Hip.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Hip(Log).UMB_Work.DS2_Tot]) + TopMargin;
            text(Log+DS2Offset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end

% Knee
subplot(413);
TopMargin = 0.25;
% if ranovatbl.StanceKnee.pValue(1) < alpha
%     multcomp = MC.StanceKnee(MC.StanceKnee.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Knee(Log).UMB_Work.Stance_Tot]) + TopMargin;
%             text(Log+StanceOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
if ranovatbl.SwingKnee.pValue(1) < alpha
    multcomp = MC.SwingKnee(MC.SwingKnee.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Knee(Log).UMB_Work.Swing_Tot]) + TopMargin;
            text(Log+SwingOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.StrideKnee.pValue(1) < alpha
    multcomp = MC.StrideKnee(MC.StrideKnee.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Knee(Log).UMB_Work.Stride_Tot]) + TopMargin;
            text(Log+StrideOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS1Knee.pValue(1) < alpha
    multcomp = MC.DS1Knee(MC.DS1Knee.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Knee(Log).UMB_Work.DS1_Tot]) + TopMargin;
            text(Log+DS1Offset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.SingSupKnee.pValue(1) < alpha
    multcomp = MC.SingSupKnee(MC.SingSupKnee.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Knee(Log).UMB_Work.SingSup_Tot]) + TopMargin;
            text(Log+SingSupOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS2Knee.pValue(1) < alpha
    multcomp = MC.DS2Knee(MC.DS2Knee.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Knee(Log).UMB_Work.DS2_Tot]) + TopMargin;
            text(Log+DS2Offset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end

% Ankle
subplot(414);
TopMargin = 0.25;
% if ranovatbl.StanceAnkle.pValue(1) < alpha
%     multcomp = MC.StanceAnkle(MC.StanceAnkle.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Ankle(Log).UMB_Work.Stance_Tot]) + TopMargin;
%             text(Log+StanceOffset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
if ranovatbl.SwingAnkle.pValue(1) < alpha
    multcomp = MC.SwingAnkle(MC.SwingAnkle.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Ankle(Log).UMB_Work.Swing_Tot]) + TopMargin;
            text(Log+SwingOffset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.StrideAnkle.pValue(1) < alpha
    multcomp = MC.StrideAnkle(MC.StrideAnkle.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Ankle(Log).UMB_Work.Stride_Tot]) + TopMargin;
            text(Log+StrideOffset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS1Ankle.pValue(1) < alpha
    multcomp = MC.DS1Ankle(MC.DS1Ankle.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Ankle(Log).UMB_Work.DS1_Tot]) + TopMargin;
            text(Log+DS1Offset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.SingSupAnkle.pValue(1) < alpha
    multcomp = MC.SingSupAnkle(MC.SingSupAnkle.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Ankle(Log).UMB_Work.SingSup_Tot]) + TopMargin;
            text(Log+SingSupOffset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end
if ranovatbl.DS2Ankle.pValue(1) < alpha
    multcomp = MC.DS2Ankle(MC.DS2Ankle.Time_1==3,:);
    for i = 1:4
        if multcomp.pValue(i) < alpha
            Log = multcomp.Time_2(i);
            y = max([Ankle(Log).UMB_Work.DS2_Tot]) + TopMargin;
            text(Log+DS2Offset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
        end
    end
end

% 
% subplot(4,1,1);
% ylim([0 2.5]);
subplot(4,1,2);
ylim([0 4.5]);
subplot(4,1,3);
ylim([0 3]);
subplot(4,1,4);
ylim([0 2.5]);

subplotsqueeze(CycleEnergyCostsUMB, 1.12);
saveas(CycleEnergyCostsUMB, 'CycleEnergyCostsUMB.png');
saveas(CycleEnergyCostsUMB, 'CycleEnergyCostsUMB.pdf');

%% Stance Phase Metabolic Costs
%
% SmlMkr = 8;
% StanceEnergyCostsUMB = figure('Position',[10 10 800 800]);
% x = [1 2 3 4 5];
% LW = 1.25;
% DotMargin = 0.15;
% MeanMargin = 0;
% SingSupOffset = 6;
% DS1Offset = 0;
% DS2Offset = 12;
% TxtFnt = 12;
% CapEndSz = 0.1;
%
% % Total
% subplot(411); hold on;
% for i = 1:5
%     % initial dub sup
%     xRel = x(i)-MeanMargin;
%     h = boxplot2(Total(i).UMB_Work.DS1_Tot,  xRel);
%     h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
%     h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % sing sup
%     xRel = x(i)-MeanMargin+SingSupOffset;
%     h = boxplot2(Total(i).UMB_Work.SingSup_Tot,  xRel);
%     h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
%     h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % secondary dub sup
%     xRel = x(i)-MeanMargin+DS2Offset;
%     h = boxplot2(Total(i).UMB_Work.DS2_Tot,  xRel);
%     h.box.Color = TotalColors(i,:);    h.ladj.Color = TotalColors(i,:);
%     h.lwhis.Color = TotalColors(i,:);    h.med.Color = TotalColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = TotalColors(i,:);    h.uwhis.Color = TotalColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     for k = 1:length(Total(i).UMB_Work.Stance_Tot)
%         % Umberger
%         plot(x(i)-DotMargin, Total(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
%         plot(x(i)-DotMargin+SingSupOffset , Total(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
%         plot(x(i)-DotMargin+DS2Offset , Total(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', TotalColors(i,:), 'MarkerEdgeColor', TotalColors(i,:));
%     end
% end
%
% % hip
% subplot(412); hold on;
% for i = 1:5
%
%     % initial dub sup
%     xRel = x(i)-MeanMargin;
%     h = boxplot2(Hip(i).UMB_Work.DS1_Tot,  xRel);
%     h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % sing sup
%     xRel = x(i)-MeanMargin+SingSupOffset;
%     h = boxplot2(Hip(i).UMB_Work.SingSup_Tot,  xRel);
%     h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % secondary dub sup
%     xRel = x(i)-MeanMargin+DS2Offset;
%     h = boxplot2(Hip(i).UMB_Work.DS2_Tot,  xRel);
%     h.box.Color = HipColors(i,:);    h.ladj.Color = HipColors(i,:);
%     h.lwhis.Color = HipColors(i,:);    h.med.Color = HipColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = HipColors(i,:);    h.uwhis.Color = HipColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     for k = 1:length(Hip(i).UMB_Work.Stance_Tot)
%         % Umberger
%         plot(x(i)-DotMargin, Hip(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%         plot(x(i)-DotMargin+SingSupOffset , Hip(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%         plot(x(i)-DotMargin+DS2Offset , Hip(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', HipColors(i,:), 'MarkerEdgeColor', HipColors(i,:));
%     end
% end
%
% % Knee
% subplot(413); hold on;
% for i = 1:5
%     % BOXPLOT
%     % initial dub sup
%     xRel = x(i)-MeanMargin;
%     h = boxplot2(Knee(i).UMB_Work.DS1_Tot,  xRel);
%     h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
%     h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % sing sup
%     xRel = x(i)-MeanMargin+SingSupOffset;
%     h = boxplot2(Knee(i).UMB_Work.SingSup_Tot,  xRel);
%     h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
%     h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     % secondary dub sup
%     xRel = x(i)-MeanMargin+DS2Offset;
%     h = boxplot2(Knee(i).UMB_Work.DS2_Tot,  xRel);
%     h.box.Color = KneeColors(i,:);    h.ladj.Color = KneeColors(i,:);
%     h.lwhis.Color = KneeColors(i,:);    h.med.Color = KneeColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = KneeColors(i,:);    h.uwhis.Color = KneeColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz  xRel+CapEndSz];
%
%     for k = 1:length(Knee(i).UMB_Work.Stance_Tot) % plot individual points
%         % Umberger
%         plot(x(i)-DotMargin, Knee(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
%         plot(x(i)-DotMargin+SingSupOffset , Knee(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
%         plot(x(i)-DotMargin+DS2Offset , Knee(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', KneeColors(i,:), 'MarkerEdgeColor', KneeColors(i,:));
%     end
% end
%
% % Ankle
% subplot(414); hold on;
% for i = 1:5
%     % initial dub sup
%     xRel = x(i)+MeanMargin;
%     h = boxplot2(Ankle(i).UMB_Work.DS1_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     % single sup
%     xRel = x(i)+MeanMargin+SingSupOffset;
%     h = boxplot2(Ankle(i).UMB_Work.SingSup_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     % initial dub sup
%     xRel = x(i)+MeanMargin+DS2Offset;
%     h = boxplot2(Ankle(i).UMB_Work.DS2_Tot, xRel);
%     h.box.Color = AnkleColors(i,:);    h.ladj.Color = AnkleColors(i,:);
%     h.lwhis.Color = AnkleColors(i,:);    h.med.Color = AnkleColors(i,:);
%     h.out.Marker = 'none';
%     h.uadj.Color = AnkleColors(i,:);    h.uwhis.Color = AnkleColors(i,:);
%     h.box.LineWidth = LW;    h.med.LineWidth = LW;    h.lwhis.LineWidth = LW;
%     h.uwhis.LineStyle = '-';    h.lwhis.LineStyle = '-';
%     h.uwhis.LineWidth = LW;    h.ladj.LineWidth = LW;    h.uadj.LineWidth = LW;
%     h.ladj.XData = [xRel-CapEndSz xRel+CapEndSz];
%     h.uadj.XData = [xRel-CapEndSz xRel+CapEndSz];
%
%     for k = 1:length(Ankle(i).UMB_Work.Stance_Tot) % plot individual points
%         % Umberger
%         plot(x(i)-DotMargin, Ankle(i).UMB_Work.DS1_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%         plot(x(i)-DotMargin+SingSupOffset , Ankle(i).UMB_Work.SingSup_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%         plot(x(i)-DotMargin+DS2Offset , Ankle(i).UMB_Work.DS2_Tot(k), '.', 'MarkerSize',SmlMkr,...
%             'MarkerFaceColor', AnkleColors(i,:), 'MarkerEdgeColor', AnkleColors(i,:));
%     end
% end
%
% % AXES AND LABELS
% subplot(4,1,1);
% ypos = 7;
% text(3,ypos, {'\bf Double Support', '\bf Leading Leg'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
% text(9,ypos+0.1, {'\bf Single Support'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
% text(15,ypos, {'\bf Double Support', '\bf Trailing Leg'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
%
% for j = 1:4
%     subplot(4,1,j);
%     ylabel('W / kg');
%     xlim([0.5 17.5]);
%     ax = gca;
%     ax.XTick = [1 2 3 4 5 7 8 9 10 11 13 14 15 16 17];
%     ax.XTickLabels = [];
%     ax.FontSize = 10;
% end
%
% subplot(4,1,4);
% TickLabel = {'-40%', '-20%','Norm','+20%','+40%'};
% ax.XTickLabels = TickLabel;
%
% subplot(4,1,1);
% ylim([0 6]);
% text(9, 5.6, '\bf ALL MUSCLES', 'Color', TotalColors(3,:), 'FontSize',12,'HorizontalAlignment', 'center');
% subplot(4,1,2);
% ylim([0 2]);
% text(9, 1.9, '\bf HIP MUSCLES', 'Color',HipColors(3,:), 'FontSize',12,'HorizontalAlignment', 'center');
% subplot(4,1,3);
% ylim([0 1.5]);
% text(9, 1.42, '\bf KNEE MUSCLES', 'Color', KneeColors(3,:), 'FontSize',12,'HorizontalAlignment', 'center');
% subplot(4,1,4);
% ylim([0 1.5]);
% text(9, 1.42, '\bf ANKLE MUSCLES', 'Color', AnkleColors(3,:), 'FontSize',12,'HorizontalAlignment', 'center');
%
% % ASTERISKS
% alpha = 0.05;
% Asterisk = 20;
%
% % total
% subplot(411);
% TopMargin = 0.4;
% if ranovatbl.DS1Total.pValue(1) < alpha
%     multcomp = MC.DS1Total(MC.DS1Total.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Total(Log).UMB_Work.DS1_Tot]) + TopMargin;
%             text(Log+DS1Offset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.SingSupTotal.pValue(1) < alpha
%     multcomp = MC.SingSupTotal(MC.SingSupTotal.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Total(Log).UMB_Work.SingSup_Tot]) + TopMargin;
%             text(Log+SingSupOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.DS2Total.pValue(1) < alpha
%     multcomp = MC.DS2Total(MC.DS2Total.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Total(Log).UMB_Work.DS2_Tot]) + TopMargin;
%             text(Log+DS2Offset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% % hip
% subplot(412);
% TopMargin = 0.25;
% if ranovatbl.DS1Hip.pValue(1) < alpha
%     multcomp = MC.DS1Hip(MC.DS1Hip.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Hip(Log).UMB_Work.DS1_Tot]) + TopMargin;
%             text(Log+DS1Offset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.SingSupHip.pValue(1) < alpha
%     multcomp = MC.SingSupHip(MC.SingSupHip.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Hip(Log).UMB_Work.SingSup_Tot]) + TopMargin;
%             text(Log+SingSupOffset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.DS2Hip.pValue(1) < alpha
%     multcomp = MC.DS2Hip(MC.DS2Hip.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Hip(Log).UMB_Work.DS2_Tot]) + TopMargin;
%             text(Log+DS2Offset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% % Knee
% subplot(413);
% TopMargin = 0.25;
% if ranovatbl.DS1Knee.pValue(1) < alpha
%     multcomp = MC.DS1Knee(MC.DS1Knee.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Knee(Log).UMB_Work.DS1_Tot]) + TopMargin;
%             text(Log+DS1Offset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.SingSupKnee.pValue(1) < alpha
%     multcomp = MC.SingSupKnee(MC.SingSupKnee.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Knee(Log).UMB_Work.SingSup_Tot]) + TopMargin;
%             text(Log+SingSupOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.DS2Knee.pValue(1) < alpha
%     multcomp = MC.DS2Knee(MC.DS2Knee.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Knee(Log).UMB_Work.DS2_Tot]) + TopMargin;
%             text(Log+DS2Offset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% % Ankle
% subplot(414);
% TopMargin = 0.25;
% if ranovatbl.DS1Ankle.pValue(1) < alpha
%     multcomp = MC.DS1Ankle(MC.DS1Ankle.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Ankle(Log).UMB_Work.DS1_Tot]) + TopMargin;
%             text(Log+DS1Offset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.SingSupAnkle.pValue(1) < alpha
%     multcomp = MC.SingSupAnkle(MC.SingSupAnkle.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Ankle(Log).UMB_Work.SingSup_Tot]) + TopMargin;
%             text(Log+SingSupOffset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
% if ranovatbl.DS2Ankle.pValue(1) < alpha
%     multcomp = MC.DS2Ankle(MC.DS2Ankle.Time_1==3,:);
%     for i = 1:4
%         if multcomp.pValue(i) < alpha
%             Log = multcomp.Time_2(i);
%             y = max([Ankle(Log).UMB_Work.DS2_Tot]) + TopMargin;
%             text(Log+DS2Offset, y, '*', 'Color', AnkleColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
%         end
%     end
% end
%
% subplotsqueeze(StanceEnergyCostsUMB, 1.12);
% saveas(StanceEnergyCostsUMB, 'StanceEnergyCostsUMB.png');
% saveas(StanceEnergyCostsUMB, 'StanceEnergyCostsUMB.pdf');
%
% clearvars i j k L m n o p Raw row h ax c diff x X


%% Muscle Plots
% MusclePlots = figure('Position',[20 20 800 800]);
LW = 1.5;
[~, NumMusc,NumSub] = size(Muscles(1).UMB_CombMuscles);
% MuscColors = colormap(jet(NumMusc));
j = 1;
for cond = 1:length(Muscles)
    %     subplot(5,1,cond); hold on;
    %   get average and std of muscles over gait cycle
    for musc = 1:NumMusc
        Muscles(cond).UMB_CombMuscleAvg(:,musc) = mean(Muscles(cond).UMB_CombMuscles(:,musc, :), 3);
        Muscles(cond).UMB_CombMuscleStd(:,musc) = std(Muscles(cond).UMB_CombMuscles(:,musc, :), 0, 3);
        
        
        %     find major contributors from the 38 muscles
        if mean( Muscles(cond).UMB_CombMuscleAvg(:,musc)) > 0.1
            IncludedMuscles(musc) = 1;
            %                  plot(Muscles(cond).UMB_CombMuscleAvg(:,musc), 'Color', MuscColors(musc,:),'LineWidth',LW);
        end
        
        
        for subj = 1:length(Subjects)
            MuscleMet(j).Cond = Muscles(cond).Condition;
            MuscleMet(j).CondNum = cond;
            MuscleMet(j).Subject = Subjects(subj).name;
            Str = char(Muscles(cond).UMB_Cols(musc));
            MuscleMet(j).(Str) = Muscles(cond).UMB_CombMuscles(:,musc, subj);
            
            j = j+1; 
        end
    end
    
    %   ylim([0 6]);
    %   ylabel('W / kg');
    %   title(Muscles(cond).Condition);
end
% subplot(5,1,5); hold on;
% for musc = 1:NumMusc
%     if IncludedMuscles(musc) == 1
%     annotation('textbox', [0.025, 0.1+musc/45, 0, 0], 'string', ...
%         strcat('\bf ',Muscles(cond).UMB_CombCols{musc}) , 'Color',MuscColors(musc,:));
%     end
% end
% xlabel('% Gait Cycle');

%% Create SPM structure
clearvars MuscleMet
j = 1; 
for cond = 1:length(Muscles)
    for subj = 1:length(Subjects)
        MuscleMet.Cond{j} = Muscles(cond).Condition;
        MuscleMet.CondNum(j) = cond;
        MuscleMet.Subject{j} = Subjects(subj).name;
        for musc = 1:length(Muscles(cond).UMB_CombCols)
            Str = char(Muscles(cond).UMB_CombCols(musc));
            MuscleMet.(Str)(j, :) = Muscles(cond).UMB_CombMuscles(:,musc, subj)';
        end
        j = j+1;
    end
end

% export Muscle Table
save('MuscleMetCosts.mat','MuscleMet'); 

% run SPM in Python

% {'semimem','semiten','bifemlh','bifemsh','sar','add_long','add_brev','tfl','pect','grac','iliacus','psoas','quad_fem','gem','peri','rect_fem','vas_med','vas_int','vas_lat','med_gas','lat_gas','soleus','tib_post','flex_dig','flex_hal','tib_ant','per_brev','per_long','per_tert','ext_dig','ext_hal','ercspn','intobl','extobl','glut_med','glut_min','glut_max','add_mag'}

%% 
MetSPM = importdata('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\MuscleMetCostsSPM.csv');

% determine top 18 muscles to plot
[SUM, I] = sort(sum(Muscles(3).UMB_CombMuscleAvg, 1))
SortNames = Muscles(3).UMB_CombCols(I)
%%
clc; close all;
LW = 1.5; % set line width
% TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray');  rgb('DarkSlateGray'); rgb('Black')];
TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
FntSz = 9;

MuscleMetCosts = figure('Position',[50 50 1000 800]);
Stride = 1:100;
yPk = 2;
for i = 1:5
%     % erector spinae
%     subplot(6,3,1); hold on;
%     Ind = contains([Muscles(i).UMB_CombCols], 'ercspn');
%     plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
%     title('Erector Spinae');
%     ylabel('W / kg'); ylim([ 0 yPk]);
%     ax = gca;
%     ax.XTick = [0 25 50 75 100];
%     ax.XTickLabel = [];
%     ax.FontSize = FntSz;
%     % SPM shading
%     if i == 5
%         SPMshade('ercspn', MetSPM, 2)
%     end
    
    % psoas
    subplot(6,3,1); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'psoas');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Psoas');
        ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('psoas', MetSPM, 2)
    end
    
    % iliacus
    subplot(6,3,2); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'iliacus');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Iliacus');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('iliacus', MetSPM, 2)
    end
    
    % glut max
    subplot(6,3,6); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'glut_max');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Gluteus Maximus');
%     ylabel('W / kg'); 
    ylim([ 0 3]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('glut_max', MetSPM, 3)
    end
    
    % glut med
    subplot(6,3,5); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'glut_med');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Gluteus Medius');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca; 
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('glut_med', MetSPM, 2)
    end
    
    %         % add_mag
    %     subplot(6,3,6); hold on;
    %     Ind = contains([Muscles(i).UMB_CombCols], 'add_mag');
    %     plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    %     title('Adductor Mag');
    %     ylabel('W / kg'); ylim([ 0 yPk]);
    
    % rec_fem
    subplot(6,3,3); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'rect_fem');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Rectus Femoris');
    %     ylabel('W / kg'); 
    ylim([ 0 3]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.YTick = [0 1 2 3];
    ax.XTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('rect_fem', MetSPM, 3)
    end
    
      % sar
    subplot(6,3,4); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'sar');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Sartorius');
        ylabel('W / kg'); 
    ylim([ 0 2]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.YTick = [0 1 2 3];
    ax.XTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('sar', MetSPM, 3)
    end
    
    
    % vas_med
    subplot(6,3,7); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'vas_med');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Vastus Medialis');
    ylabel('W / kg'); ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('vas_med', MetSPM, 2)
    end
    
    % vas_int
    subplot(6,3,8); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'vas_int');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Vastus Intermedialis');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('vas_int', MetSPM, 2)
    end
    
    % vas_lat
    subplot(6,3,9); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'vas_lat');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Vastus Lateralis');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('vas_lat', MetSPM, 2)
    end
    
    % bifemlh
    subplot(6,3,10); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'bifemlh');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Biceps Femoris - Long Head');
    ylabel('W / kg'); ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('bifemlh', MetSPM, 2)
    end
    
    % bifemsh
    subplot(6,3,11); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'bifemsh');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Biceps Femoris - Short Head');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('bifemsh', MetSPM, 2)
    end
    
    % semimem
    subplot(6,3,13); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'semimem');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Semimembranosus');
    ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('semimem', MetSPM, 2)
    end
    
%     % semiten
%     subplot(6,3,13); hold on;
%     Ind = contains([Muscles(i).UMB_CombCols], 'semiten');
%     plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
%     title('Semitendinosus');
%     ylabel('W / kg'); ylim([ 0 yPk]);
%     ax = gca;
%     ax.XTick = [0 25 50 75 100];
%     ax.XTickLabel = [];
%     ax.FontSize = FntSz;
%     % SPM shading
%     if i == 5
%         SPMshade('semiten', MetSPM, 2)
%     end
    
     % add_long
    subplot(6,3,12); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'add_long');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Adductor Longus');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('add_long', MetSPM, 2)
    end
    
    % soleus
    subplot(6,3,18); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'soleus');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Soleus');
    %     ylabel('W / kg'); 
    ylim([ 0 6]);
%     ax = gca;
%     ax.XTick = [0 25 50 75 100];
%     ax.XTickLabel = [];
%     ax.FontSize = FntSz;
        ax = gca; 
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = {'0' '25' '50' '75' '100'};
%     ax.YTickLabel = [];
    xlabel('% Gait Cycle');
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('soleus', MetSPM, 6)
    end
    
    % med_gas
    subplot(6,3,17); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'med_gas');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Medial Gastrocnemius');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
%     ax = gca;
%     ax.XTick = [0 25 50 75 100];
%     ax.XTickLabel = [];
%     ax.YTickLabel = [];
%     ax.FontSize = FntSz;
    ax = gca; 
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = {'0' '25' '50' '75' '100'};
    ax.YTickLabel = [];
    xlabel('% Gait Cycle');
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('med_gas', MetSPM, 2)
    end
    
    % lat_gas
    subplot(6,3,16); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'lat_gas');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Lateral Gastrocnemius');
    ylabel('W / kg'); ylim([ 0 yPk]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = {'0' '25' '50' '75' '100'};
    xlabel('% Gait Cycle');
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('lat_gas', MetSPM, 2)
    end
    
    % tib_ant
    subplot(6,3,14); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'tib_ant');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Tibialis Anterior');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
%     ax = gca;
%     ax.XTick = [0 25 50 75 100];
%     ax.XTickLabel = {'0' '25' '50' '75' '100'};
        ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
%     xlabel('% Gait Cycle');
%     ax.YTickLabel = [];
%     ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('tib_ant', MetSPM, 2)
    end
    
    % ext_dig
    subplot(6,3,15); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'ext_dig');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Extensor Digitorum');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
        ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
%     ax = gca; 
%     ax.XTick = [0 25 50 75 100];
%     ax.XTickLabel = {'0' '25' '50' '75' '100'};
%     ax.YTickLabel = [];
%     xlabel('% Gait Cycle');
%     ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('ext_dig', MetSPM, 2)
    end
end

subplot(631);
Conditions = {'-40', '-20','Norm','+20','+40'}; 
legend(Conditions, 'Location','Best');
% 
% subplot(6,3,10);
% Conditions = {'-40', '-20','Norm','+20','+40'}; 
% legend(Conditions, 'Location','Best');

subplotsqueeze(MuscleMetCosts, 1.12);

saveas(MuscleMetCosts, 'MuscleMetCosts.png'); 
saveas(MuscleMetCosts, 'MuscleMetCosts.pdf'); 


