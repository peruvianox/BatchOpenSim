% Process CMC Metabolics Results after batch processing
% create figures and run stats

clear;
clc;
% dbstop if error
% cd('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts');
cd('C:\Users\richa\Documents\Packages\OpenSim\Scripts');
addpath(genpath('BatchOpenSim'));
addpath(genpath('CodeLibrary'));

% load subjext files
% subjectPath = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\exp';
subjectPath = 'D:\UNC_ABL\FpMetabolics_2019';
if ~exist(subjectPath, 'dir')
    subjectPath = uigetdir(CurrFolder, 'Select Folder Containing Subject Data');
end
addpath(genpath(subjectPath));
load('Subjects.mat');

% define metabolics folder
% MetabolicsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\MetabolicReports';
MetabolicsFolder = 'C:\Users\richa\Documents\Packages\OpenSim\Scripts\CMC_Results\MetabolicReports';
addpath(genpath(MetabolicsFolder));
MetabolicsDir = dir(MetabolicsFolder);

% add walking speed to demographics
if isfield(Subjects, 'PrefWalkSpeed') == 0
%     Demo.File = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Metabolics_Feedback_Demographics';
    Demo.File = 'C:\Users\richa\Documents\Packages\OpenSim\Metabolics_Feedback_Demographics';
    
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
[Struct, Subjects] = ExtractData(MetabolicsFolder, Subjects, 'met');

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
         Array(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Total.left_total_parsedSum;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Total.right_total_parsedSum;
        Subjects(i).Trials(j).Met.Total.UMB_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.UMB.ParsedJoints.Total.left_total_parsedSum;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.UMB.ParsedJoints.Total.right_total_parsedSum;
        Subjects(i).Trials(j).Met.Total.UMB_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Total.UMB_Avg_Sum = sum(Subjects(i).Trials(j).Met.Total.UMB_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Total.UMB_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Total.UMB_Avg_total_parsed);
        clearvars Array
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Total.UMB_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Total.parsed, ...
            Subjects(i).Trials(j).Met.Left.UMB.InterpJoints.Total.left_total_parsedSum, Subjects(i).Trials(j).Met.Right.UMB.InterpJoints.Total.right_total_parsedSum,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
        
        % Bhargava
            Array(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Total.left_total_parsedSum;
        Array(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Total.right_total_parsedSum;
        Subjects(i).Trials(j).Met.Total.BHAR_Avg_parsed = mean(Array, 3);
        Mat(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Total.left_total_parsedSum;
        Mat(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Total.right_total_parsedSum;
        Subjects(i).Trials(j).Met.Total.BHAR_Avg_total_parsed = mean(Mat, 3);
        Subjects(i).Trials(j).Met.Total.BHAR_Avg_Sum = sum(Subjects(i).Trials(j).Met.Total.BHAR_Avg_parsed); % sums
        Subjects(i).Trials(j).Met.Total.BHAR_Avg_total_Sum = sum(Subjects(i).Trials(j).Met.Total.BHAR_Avg_total_parsed);
        clearvars Array
        % Integrate for Energy Expenditure
        Subjects(i).Trials(j).Met.Total.BHAR_Energy = Integrate4EnergyCost(...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Total.parsed, ...
            Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Total.left_total_parsedSum, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Total.right_total_parsedSum,...
            Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
       
        
%         Array(:,:,1) = Subjects(i).Trials(j).Met.Left.BHAR.ParsedJoints.Total.parsed;
%         Array(:,:,2) = Subjects(i).Trials(j).Met.Right.BHAR.ParsedJoints.Total.parsed;
%         Subjects(i).Trials(j).Met.Total.BHAR_parsed = mean(Array, 3);
%         clearvars Array
%         % Integrate for Energy Expenditure
%         Subjects(i).Trials(j).Met.Total.BHAR_Energy = Integrate4EnergyCost(...
%             Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Total.parsed, ...
%             Subjects(i).Trials(j).Met.Left.BHAR.InterpJoints.Total.parsed, Subjects(i).Trials(j).Met.Right.BHAR.InterpJoints.Total.parsed,...
%             Subjects(i).Trials(j).Met.Left.LogTimes, Subjects(i).Trials(j).Met.Right.LogTimes, Subjects(i).PrefWalkSpeed);
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
        Total(trial).UMB_Total(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Avg_parsed / Subjects(i).Demo.mass;
        Total(trial).UMB_Total_Avg(:,:,ind) =   mean(Total(trial).UMB_Total(:,:,ind));
        Total(trial).UMB_Work.Stride_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.Stance_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.Stance.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.DS1_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.DS1.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.SingSup_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.SingSup.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.DS2_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.DS2.Avg_Total_Int / Subjects(i).Demo.mass;
        Total(trial).UMB_Work.Swing_Tot(:,:,ind) = Subjects(i).Trials(j).Met.Total.UMB_Energy.Swing.Avg_Total_Int / Subjects(i).Demo.mass;
        
        Total(trial).BHAR_SubjectsIncluded{ind} = Subjects(i).name / Subjects(i).Demo.mass;
        Total(trial).BHAR_Total(:,:,ind) = Subjects(i).Trials(j).Met.Total.BHAR_Avg_parsed / Subjects(i).Demo.mass;
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

    %% get gait event averages parsed as % of cycle
    
    for subj = 1:length(Subjects)
        for trial = 1:length(Subjects(subj).Trials)
            if strcmp(Subjects(subj).Trials(trial).type, 'static')
                continue
            end
            
            NumFutureStrides = 5;
            
            L_TrlStart = Subjects(subj).Trials(trial).Met.Left.Time(1);
            R_TrlStart = Subjects(subj).Trials(trial).Met.Right.Time(1);
            
            TrlStart = min([L_TrlStart, R_TrlStart]); 
            
            [~, L_StartInd] = min(abs(Subjects(subj).Trials(trial).TSData.L_Strike(:,2) - L_TrlStart));
            [~, L_EndInd] = min(abs(Subjects(subj).Trials(trial).TSData.L_Off(:,2) - L_TrlStart));
            if L_EndInd<L_StartInd
                L_EndInd = L_EndInd + 1;
            end
            L_End = Subjects(subj).Trials(trial).TSData.L_Strike(L_StartInd+NumFutureStrides,2);
            
            [~, R_StartInd] = min(abs(Subjects(subj).Trials(trial).TSData.R_Strike(:,2) - R_TrlStart));
            [~, R_EndInd] = min(abs(Subjects(subj).Trials(trial).TSData.R_Off(:,2) - R_TrlStart));
            if R_EndInd<R_StartInd
                R_EndInd = R_EndInd + 1;
            end
              R_End = Subjects(subj).Trials(trial).TSData.R_Strike(R_StartInd+NumFutureStrides,2);
            
            TrlEnd = max([L_End R_End]); 
            
            ind1 = Subjects(subj).Trials(trial).TSData.L_Strike(:,2) > TrlStart;
            ind2 = Subjects(subj).Trials(trial).TSData.L_Strike(:,2) < TrlEnd;
            L_Strikes = Subjects(subj).Trials(trial).TSData.L_Strike(ind1+ind2 ==2, 2);
            
            ind1 = Subjects(subj).Trials(trial).TSData.L_Off(:,2) > TrlStart;
            ind2 = Subjects(subj).Trials(trial).TSData.L_Off(:,2) < TrlEnd;
            L_Offs = Subjects(subj).Trials(trial).TSData.L_Off(ind1+ind2 ==2, 2);
            
            ind1 = Subjects(subj).Trials(trial).TSData.R_Strike(:,2) > TrlStart;
            ind2 = Subjects(subj).Trials(trial).TSData.R_Strike(:,2) < TrlEnd;
            R_Strikes = Subjects(subj).Trials(trial).TSData.R_Strike(ind1+ind2 ==2, 2);
            
            ind1 = Subjects(subj).Trials(trial).TSData.R_Off(:,2) > TrlStart;
            ind2 = Subjects(subj).Trials(trial).TSData.R_Off(:,2) < TrlEnd;
            R_Offs = Subjects(subj).Trials(trial).TSData.R_Off(ind1+ind2 ==2, 2);
            
%             L_Strikes = Subjects(subj).Trials(trial).TSData.L_Strike(L_StartInd:L_StartInd+NumFutureStrides-1,2);
%             L_Offs = Subjects(subj).Trials(trial).TSData.L_Off(L_EndInd:L_EndInd+NumFutureStrides-1,2);
%             
%             R_Strikes = Subjects(subj).Trials(trial).TSData.R_Strike(R_StartInd:R_StartInd+NumFutureStrides-1,2);
%             R_Offs = Subjects(subj).Trials(trial).TSData.R_Off(R_EndInd:R_EndInd+NumFutureStrides-1,2);
            
            % make sure no gait events are missing
%             CheckGaitEvents
            
           [Subjects(subj).Trials(trial).TSData.L_Times, Subjects(subj).Trials(trial).TSData.R_Times] ...
               = GaitEventParse(L_Strikes, L_Offs, R_Strikes, R_Offs);
           
           % calculate average gait event times bilaterally
           % and save in TSData structure
           Subjects(subj).Trials(trial).TSData.Off_Avg = nanmean(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).FootOff_p, ...
               Subjects(subj).Trials(trial).TSData.R_Times(:).FootOff_p]);
           Subjects(subj).Trials(trial).TSData.Off_std = nanstd(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).FootOff_p,...
               Subjects(subj).Trials(trial).TSData.R_Times(:).FootOff_p]);
           
             Subjects(subj).Trials(trial).TSData.OppOff_Avg = nanmean(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).OppFootOff_p, ...
               Subjects(subj).Trials(trial).TSData.R_Times(:).OppFootOff_p]);
           Subjects(subj).Trials(trial).TSData.OppOff_std = nanstd(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).OppFootOff_p,...
               Subjects(subj).Trials(trial).TSData.R_Times(:).OppFootOff_p]);
           
             Subjects(subj).Trials(trial).TSData.OppOn_Avg = nanmean(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).OppFootOn_p, ...
               Subjects(subj).Trials(trial).TSData.R_Times(:).OppFootOn_p]);
           Subjects(subj).Trials(trial).TSData.OppOn_std = nanstd(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).OppFootOn_p, ...
               Subjects(subj).Trials(trial).TSData.R_Times(:).OppFootOn_p]);
           
     Subjects(subj).Trials(trial).TSData.StrideTime_Avg = nanmean(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).StrideTime, ...
               Subjects(subj).Trials(trial).TSData.R_Times(:).StrideTime]);
           Subjects(subj).Trials(trial).TSData.StrideTime_std = nanstd(...
               [Subjects(subj).Trials(trial).TSData.L_Times(:).StrideTime...
               Subjects(subj).Trials(trial).TSData.R_Times(:).StrideTime]);
        end
    end
    
    clearvars ind1 ind2 L_Offs R_Offs L_Strikes R_Strikes L_TrlStart L_End R_TrlStart R_End TrlStart TrlEnd ans...
        L_StartInd R_StartInd L_EndInd R_EndInd n offs opp_offs strikes opp_strikes sf ToDel GCs i j Start End Mat ...
        r_offs r_strikes l_offs l_strikes Iffs Strujes Opp_Offs Offs Opp_Off Opp_Strikes R_Times L_Times Off_Ind
    
    %% Create table of stride timing data
    StrideFreq = struct.empty(0); 
    
    for subj = 1:length(Subjects)
        for trial = 1:length(Subjects(subj).Trials)
            
            if strcmp(Subjects(subj).Trials(trial).type, 'static')
                continue
            end
            
             StrideFreq(trial).Cond = Subjects(subj).Trials(trial).name;
            
            StrideFreq(trial).StrideTimeAvg(subj) = Subjects(subj).Trials(trial).TSData.StrideTime_Avg;
            StrideFreq(trial).StrideTimeStd(subj) = Subjects(subj).Trials(trial).TSData.StrideTime_std;
            sf = StrideFreq(trial).StrideTimeAvg(subj);
            
            % leading limb dub support
            StrideFreq(trial).DS1Avg(subj) = Subjects(subj).Trials(trial).TSData.OppOff_Avg;
            StrideFreq(trial).DS1Avg_T(subj) = StrideFreq(trial).DS1Avg(subj) * sf;
            
            % single support
            StrideFreq(trial).SingSupAvg(subj) = Subjects(subj).Trials(trial).TSData.OppOn_Avg - Subjects(subj).Trials(trial).TSData.OppOff_Avg;
            StrideFreq(trial).SingSupAvg_T(subj) = StrideFreq(trial).SingSupAvg(subj) * sf;
     
            % trailing limb dub support
            StrideFreq(trial).DS2Avg(subj) = Subjects(subj).Trials(trial).TSData.Off_Avg - Subjects(subj).Trials(trial).TSData.OppOn_Avg;
            StrideFreq(trial).DS2Avg_T(subj) = StrideFreq(trial).DS2Avg(subj) * sf;
            
            % swing
            StrideFreq(trial).SwingAvg(subj) = 1 - Subjects(subj).Trials(trial).TSData.Off_Avg;
            StrideFreq(trial).SwingAvg_T(subj) = StrideFreq(trial).SwingAvg(subj) * sf;
            
        end
    end
    
    % calculate subject - level group averages and standard deviations
    SF = struct.empty(0);
    t = 0; 
    for trial = [2 1 5 3 4]
        t = t + 1; 
        
         SF(t).Cond =  StrideFreq(trial).Cond;
        SF(t).StrideTimeAvg = mean(StrideFreq(trial).StrideTimeAvg);
        SF(t).StrideTimeStd = std(StrideFreq(trial).StrideTimeAvg);
        
        % Initial dub sup
        SF(t).DS1Avg = mean(StrideFreq(trial).DS1Avg);
        SF(t).DS1Std = std(StrideFreq(trial).DS1Avg);
        SF(t).DS1Avg_T = mean(StrideFreq(trial).DS1Avg_T);
        SF(t).DS1Std_T = std(StrideFreq(trial).DS1Avg_T);
        
        % sing sup
        SF(t).SingSupAvg = mean(StrideFreq(trial).SingSupAvg);
        SF(t).SingSupStd = std(StrideFreq(trial).SingSupAvg);
        SF(t).SingSupAvg_T = mean(StrideFreq(trial).SingSupAvg_T);
        SF(t).SingSupStd_T = std(StrideFreq(trial).SingSupAvg_T);
        
        % double support trailling limb
        SF(t).DS2Avg = mean(StrideFreq(trial).DS2Avg);
        SF(t).DS2Std = std(StrideFreq(trial).DS2Avg);
        SF(t).DS2Avg_T = mean(StrideFreq(trial).DS2Avg_T);
        SF(t).DS2Std_T = std(StrideFreq(trial).DS2Avg_T);
        
        % swing
        SF(t).SwingAvg = mean(StrideFreq(trial).SwingAvg);
        SF(t).SwingStd = std(StrideFreq(trial).SwingAvg);
        SF(t).SwingAvg_T = mean(StrideFreq(trial).SwingAvg_T);
        SF(t).SwingStd_T = std(StrideFreq(trial).SwingAvg_T);
    end
    
    % export data to excel table
    Data = zeros(10);
    T = [1 3 5 7 9];
    P = [2 4 6 8 10];
    
    % Stride time
    Data(1,T) = [SF.StrideTimeAvg];
    Data(1,P) = 0;
    Data(2,T) = [SF.StrideTimeStd];
    Data(2,P) = 0;
    
    % lead limb double support
    Data(3,T) = [SF.DS1Avg_T];
    Data(3,P) = 100 * [SF.DS1Avg];
    Data(4,T) = [SF.DS1Std_T];
    Data(4,P) = 100 * [SF.DS1Std];
    
    % lead limb double support
    Data(5,T) = [SF.SingSupAvg_T];
    Data(5,P) = 100 * [SF.SingSupAvg];
    Data(6,T) = [SF.SingSupStd_T];
    Data(6,P) = 100 * [SF.SingSupStd];
    
    % lead limb double support
    Data(7,T) = [SF.DS2Avg_T];
    Data(7,P) = 100 * [SF.DS2Avg];
    Data(8,T) = [SF.DS2Std_T];
    Data(8,P) = 100 * [SF.DS2Std];
    
    % lead limb double support
    Data(9,T) = [SF.SwingAvg_T];
    Data(9,P) = 100 *  [SF.SwingAvg];
    Data(10,T) = [SF.SwingStd_T];
    Data(10,P) = 100 * [SF.SwingStd];

    % write to table
    % copy and past it to table 2
    

    
    
%% Correlate modelled and measured metabolic cost
clearvars R n o p
% load measured met cost
% [~, ~, Raw] = xlsread('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Subject_NetMetabolics.xlsx');
[~, ~, Raw] = xlsread('C:\Users\richa\Documents\Packages\OpenSim\Subject_NetMetabolics.xlsx');
c = 0;
% TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray');  rgb('DarkSlateGray'); rgb('Black')];
TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];


SubjColors = colormap(jet(length(Subjects)));
TrialSymbols = {'-', '-', 'o', '+','+'};
TxtFnt = 12;
AxFnt = 9;
MkrSz = 12;
TitleFnt = 16;
close all; clc;
CorrVal = figure('Position', [100 100 1000 800]);
Settings.AddStandingMetRate = 'Yes'; 

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
                MetCost(c).Measured = Subjects(subj).Trials(trial).AvgMetCost;
                MetCost(c).Standing =  Raw{row, 8};
                %                 MetCost(c).Measured = Subjects(subj).Trials(trial).AvgMetCost + MetCost(c).Standing;
                
                MetCost(c).UmbModelled = Subjects(subj).Trials(trial).Met.Total.UMB_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
                MetCost(c).BharModelled = Subjects(subj).Trials(trial).Met.Total.BHAR_Energy.Stride.Avg_Total_Int / Subjects(subj).Demo.mass;
                
                GaitEvents(c).Subj = Subjects(subj).name;
                GaitEvents(c).Trial = Subjects(subj).Trials(trial).name;
                GaitEvents(c).Off_Avg = Subjects(subj).Trials(trial).TSData.Off_Avg;
                GaitEvents(c).Off_std = Subjects(subj).Trials(trial).TSData.Off_std;
                GaitEvents(c).OppOff_Avg = Subjects(subj).Trials(trial).TSData.OppOff_Avg;
                GaitEvents(c).OppOff_std = Subjects(subj).Trials(trial).TSData.OppOff_std;
                GaitEvents(c).OppOn_Avg = Subjects(subj).Trials(trial).TSData.OppOn_Avg;
                GaitEvents(c).OppOn_std = Subjects(subj).Trials(trial).TSData.OppOn_std;
                
%                 R_Start = Subjects(subj).Trials(trial).Met.Right.GC_Start;
                R_Off = (Subjects(subj).Trials(trial).Met.Right.GC_Off - Subjects(subj).Trials(trial).Met.Right.GC_Start) / ...
                    (Subjects(subj).Trials(trial).Met.Right.GC_End - Subjects(subj).Trials(trial).Met.Right.GC_Start);
                R_OppOff = (Subjects(subj).Trials(trial).Met.Right.GC_1DS - Subjects(subj).Trials(trial).Met.Right.GC_Start) / ...
                    (Subjects(subj).Trials(trial).Met.Right.GC_End - Subjects(subj).Trials(trial).Met.Right.GC_Start);
                R_OppOn = (Subjects(subj).Trials(trial).Met.Right.GC_2DS - Subjects(subj).Trials(trial).Met.Right.GC_Start) / ...
                    (Subjects(subj).Trials(trial).Met.Right.GC_End - Subjects(subj).Trials(trial).Met.Right.GC_Start);
                
                      L_Off = (Subjects(subj).Trials(trial).Met.Left.GC_Off - Subjects(subj).Trials(trial).Met.Left.GC_Start) / ...
                    (Subjects(subj).Trials(trial).Met.Left.GC_End - Subjects(subj).Trials(trial).Met.Left.GC_Start);
                L_OppOff = (Subjects(subj).Trials(trial).Met.Left.GC_1DS - Subjects(subj).Trials(trial).Met.Left.GC_Start) / ...
                    (Subjects(subj).Trials(trial).Met.Left.GC_End - Subjects(subj).Trials(trial).Met.Left.GC_Start);
                L_OppOn = (Subjects(subj).Trials(trial).Met.Left.GC_2DS - Subjects(subj).Trials(trial).Met.Left.GC_Start) / ...
                    (Subjects(subj).Trials(trial).Met.Left.GC_End - Subjects(subj).Trials(trial).Met.Left.GC_Start);
                 
                 
                
                 GaitEvents(c).Off_Avg1 = mean([L_Off R_Off]); 
%                 GaitEvents(c).Off_std1 = Subjects(subj).Trials(trial).TSData.Off_std;
                GaitEvents(c).OppOff_Avg1 = mean([L_OppOff R_OppOff]); 
%                 GaitEvents(c).OppOff_std1 = Subjects(subj).Trials(trial).TSData.OppOff_std;
                GaitEvents(c).OppOn_Avg1 = mean([L_OppOn R_OppOn]); 
%                 GaitEvents(c).OppOn_std1 = Subjects(subj).Trials(trial).TSData.OppOn_std;
                
                
            end
        end
    end
end

clearvars R P
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
    subplot(335); hold on;
    plot(MetCostTable.Measured(Cond_Ind(:,cond)), MetCostTable.UmbModelled(Cond_Ind(:,cond)),...
        '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
    subplot(336); hold on;
    plot(MetCostTable.Measured(Cond_Ind(:,cond)), MetCostTable.BharModelled(Cond_Ind(:,cond)),...
        '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
    subplot(334); hold on;
    plot(MetCostTable.UmbModelled(Cond_Ind(:,cond)), MetCostTable.BharModelled(Cond_Ind(:,cond)),...
        '.', 'MarkerSize', MkrSz, 'Color', TotalColors(cond, :), 'MarkerFaceColor',TotalColors(cond, :));
end

% generate linear trendlines
x = linspace(0, 15); 
Trend.Umb_fit = polyfit(MetCostTable.Measured, MetCostTable.UmbModelled, 1); 
Trend.Umb_val = polyval(Trend.Umb_fit, x); 

Trend.Bhar_fit = polyfit(MetCostTable.Measured, MetCostTable.BharModelled, 1); 
Trend.Bhar_val = polyval(Trend.Bhar_fit, x); 

Trend.Both_fit = polyfit(MetCostTable.UmbModelled, MetCostTable.BharModelled, 1); 
Trend.Both_val = polyval(Trend.Both_fit, x); 

LW = 1; 

% plot correlations
subplot(334);
x = linspace(0, 15); 
plot(x,x, '--k', 'LineWidth',LW, 'Color',rgb('LightGray')); 
plot(x,Trend.Both_val, '-k', 'LineWidth',LW, 'Color',rgb('Gray')); 
ylim([2 10]);
xlim([2 10]);
title('Inter-Model', 'FontSize', TitleFnt);
ylabel('Bhargava (W/kg)');
xlabel('Umberger (W/kg)');
ax = gca;
ax.FontSize = AxFnt;
text(10, 4, ['R = ' num2str(round(R.BharUmb, 4))], 'HorizontalAlignment', 'right');
% text(11,8.5, ['P = ' num2str(P.BharUmb, 4)]);
text(10, 3.25, ['P < 0.001'], 'HorizontalAlignment', 'right');
text(10, 2.5, ['y = ', num2str(round(Trend.Both_fit(1),2)), 'x + ', num2str(round(Trend.Both_fit(2),2))], 'HorizontalAlignment', 'right');
text(1, 11, '\bf D', 'HorizontalAlignment','center')

subplot(335);
plot(x,x, '--k', 'LineWidth',LW, 'Color',rgb('LightGray')); 
plot(x,Trend.Umb_val, '-k', 'LineWidth',LW, 'Color',rgb('Gray')); 
ylim([2 10]);
xlim([2 10]);
title('Umberger Model', 'FontSize', TitleFnt);
xlabel('Measured (W/kg)');
ylabel('Simulated (W/kg)');
ax = gca;
ax.FontSize = AxFnt;
text(10, 4, ['R = ' num2str(round(R.Umb, 4))], 'HorizontalAlignment', 'right');
% text(7, 8.5, ['P = ' num2str(P.Umb)]);
text(10, 3.25, ['P < 0.001'], 'HorizontalAlignment', 'right');
text(10, 2.5, ['y = ', num2str(round(Trend.Umb_fit(1),2)), 'x + ', num2str(round(Trend.Umb_fit(2),2))], 'HorizontalAlignment', 'right');
text(1, 11, '\bf E', 'HorizontalAlignment','center')

subplot(336);
plot(x,x, '--k', 'LineWidth',LW, 'Color',rgb('LightGray')); 
plot(x,Trend.Bhar_val, '-k',  'LineWidth',LW, 'Color',rgb('Gray')); 
ylim([2 10]);
xlim([2 10]);
title('Bhargava Model', 'FontSize', TitleFnt);
xlabel('Measured (W/kg)');
ylabel('Simulated (W/kg)');
ax = gca;
ax.FontSize = AxFnt;
text(10, 4, ['R = ' num2str(round(R.Bhar, 4))], 'HorizontalAlignment', 'right');
% text(7, 8.5, ['P = ' num2str(P.Bhar)]);
text(10, 3.25, ['P < 0.001'], 'HorizontalAlignment', 'right');
text(10, 2.5, ['y = ', num2str(round(Trend.Bhar_fit(1),2)), 'x + ', num2str(round(Trend.Bhar_fit(2),2))], 'HorizontalAlignment', 'right');
text(1, 11, '\bf F', 'HorizontalAlignment','center')


     %% Stats for metabolic cost
clc; clearvars Table table Meas MeasTable
Conditions = {'Hip1','Hip2','Hip3','Hip4','Hip5', 'Ankle1','Ankle2','Ankle3','Ankle4','Ankle5'};
k = 1;

for j = 1:length(Subjects)
    
    % stride time - in seconds
    ST(k).y1 = StrideFreq(2).StrideTimeAvg(j) ;
    ST(k).y2 = StrideFreq(1).StrideTimeAvg(j) ;
    ST(k).y3 = StrideFreq(5).StrideTimeAvg(j) ;
    ST(k).y4 = StrideFreq(3).StrideTimeAvg(j) ;
    ST(k).y5 = StrideFreq(4).StrideTimeAvg(j) ;
    
    % Leading Limb Dub Sup - in % GC
    DS1(k).y1 = StrideFreq(2).DS1Avg(j);
    DS1(k).y2 = StrideFreq(1).DS1Avg(j);
    DS1(k).y3 = StrideFreq(5).DS1Avg(j);
    DS1(k).y4 = StrideFreq(3).DS1Avg(j);
    DS1(k).y5 = StrideFreq(4).DS1Avg(j);
    
    %Single sup  - in % GC
    SingSup(k).y1 = StrideFreq(2).SingSupAvg(j);
    SingSup(k).y2 = StrideFreq(1).SingSupAvg(j);
    SingSup(k).y3 = StrideFreq(5).SingSupAvg(j);
    SingSup(k).y4 = StrideFreq(3).SingSupAvg(j);
    SingSup(k).y5 = StrideFreq(4).SingSupAvg(j);
    
    % trailing limb dub sup - in % GC
    DS2(k).y1 = StrideFreq(2).DS2Avg(j);
    DS2(k).y2 = StrideFreq(1).DS2Avg(j);
    DS2(k).y3 = StrideFreq(5).DS2Avg(j);
    DS2(k).y4 = StrideFreq(3).DS2Avg(j);
    DS2(k).y5 = StrideFreq(4).DS2Avg(j);
    
    % swing - in % GC
    Swing(k).y1 = StrideFreq(2).SwingAvg(j);
    Swing(k).y2 = StrideFreq(1).SwingAvg(j);
    Swing(k).y3 = StrideFreq(5).SwingAvg(j);
    Swing(k).y4 = StrideFreq(3).SwingAvg(j);
    Swing(k).y5 = StrideFreq(4).SwingAvg(j);
    
    % ---------------------------------------------------------------------
    % overall metabolic cost
    MetabolicCost(k).y1 = Subjects(j).Trials(2).AvgMetCost; 
    MetabolicCost(k).y2 = Subjects(j).Trials(1).AvgMetCost; 
    MetabolicCost(k).y3 = Subjects(j).Trials(5).AvgMetCost; 
    MetabolicCost(k).y4 = Subjects(j).Trials(3).AvgMetCost; 
    MetabolicCost(k).y5 = Subjects(j).Trials(4).AvgMetCost; 
    
    % ---------------------------------------------------------------------
    % Total
    DS1TotalTable(k).Subject = Subjects(j).name;
    DS1TotalTable(k).Joint = 'Total';
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
    DS1HipTable(k).Subject = Subjects(j).name;
    DS1HipTable(k).Joint = 'Hip';
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
    DS1KneeTable(k).Subject = Subjects(j).name;
    DS1KneeTable(k).Joint = 'Knee';
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

% Normality Tests
% Ahmed BenSaïda (2020). Shapiro-Wilk and Shapiro-Francia normality tests. 
% (https://www.mathworks.com/matlabcentral/fileexchange/13964-shapiro-wilk-and-shapiro-francia-normality-tests), 
% MATLAB Central File Exchange. Retrieved December 15, 2020.
addpath(genpath('BatchOpenSim/PostProcessingFxns/swtest')); 
n = 5;
mkr = 12; 
StrideTimeFig = figure('Position',[100 100 1000 400]); 
txt = 18; 

% stride times
StrideTimetable = struct2table(ST); 
[rm.StrideTime] = fitrm(StrideTimetable, 'y1-y5 ~ 1');
ranovatbl.StrideTime = ranova(rm.StrideTime);
MC.StrideTime = multcompare(rm.StrideTime, 'Time', 'ComparisonType', 'lsd');
disp('Stride Time Avg = '); 
mean(table2array(StrideTimetable), 1)
disp('Stride Time Std = '); 
std(table2array(StrideTimetable), 1)
H = zeros(1, n); 
for i = 1:n
    H(i) = swtest(table2array(StrideTimetable(:,i)));
end
disp('Stride Time Normality = '); 
disp(H);
x = ones(1,12); 
subplot(151); hold on; 
plot(x, [ST.y1], '.', 'MarkerSize',mkr, 'Color',TotalColors(1,:)); 
plot(x*2, [ST.y2], '.', 'MarkerSize',mkr, 'Color',TotalColors(2,:)); 
plot(x*3, [ST.y3], '.', 'MarkerSize',mkr, 'Color',TotalColors(3,:)); 
plot(x*4, [ST.y4], '.', 'MarkerSize',mkr, 'Color',TotalColors(4,:)); 
plot(x*5, [ST.y5], '.', 'MarkerSize',mkr, 'Color',TotalColors(5,:)); 
title('Stride TIme'); 
xlim([0.5 5.5]); 
ylim([0.5 1.5]); 
ylabel('s'); 
ax = gca; 
ax.XTickLabel = {'-40','-20','Norm','+20','+40'};
ax.YTick = [0.5 0.75 1 1.25 1.5]; 
text(1,max([ST.y1]) + 0.1, '*', 'FontSize', txt, 'Color',TotalColors(1,:), 'HorizontalAlignment','center'); 
text(2,max([ST.y2]) + 0.1, '*', 'FontSize', txt, 'Color',TotalColors(2,:), 'HorizontalAlignment','center'); 
text(5,max([ST.y5]) + 0.1, '*', 'FontSize', txt, 'Color',TotalColors(5,:), 'HorizontalAlignment','center'); 

DS1table = struct2table(DS1); 
[rm.DS1] = fitrm(DS1table, 'y1-y5 ~ 1');
ranovatbl.DS1 = ranova(rm.DS1);
MC.DS1 = multcompare(rm.DS1, 'Time', 'ComparisonType', 'lsd');
disp('DS1 Time Avg = '); 
mean(table2array(DS1table), 1)
disp('DS1 Time Std = '); 
std(table2array(DS1table), 1)
H = zeros(1, n); 
for i = 1:n
    H(i) = swtest(table2array(DS1table(:,i)));
end
disp('DS1 Time Normality = '); 
disp(H);
subplot(152); hold on; 
plot(x, [DS1.y1], '.', 'MarkerSize',mkr, 'Color',TotalColors(1,:)); 
plot(x*2, [DS1.y2], '.', 'MarkerSize',mkr, 'Color',TotalColors(2,:)); 
plot(x*3, [DS1.y3], '.', 'MarkerSize',mkr, 'Color',TotalColors(3,:)); 
plot(x*4, [DS1.y4], '.', 'MarkerSize',mkr, 'Color',TotalColors(4,:)); 
plot(x*5, [DS1.y5], '.', 'MarkerSize',mkr, 'Color',TotalColors(5,:)); 
text(4,max([DS1.y4]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(4,:), 'HorizontalAlignment','center'); 
text(5,max([DS1.y5]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(5,:), 'HorizontalAlignment','center'); 
title({'Double Support'; 'Leading Limb'}); 
xlim([0.5 5.5]); 
ylim([0.1 0.25]); 
ylabel('% Gait Cycle'); 
ax = gca; 
ax.XTickLabel = {'-40','-20','Norm','+20','+40'};
ax.YTickLabel = [10 15 20 25];


SingSuptable = struct2table(SingSup); 
[rm.SingSup] = fitrm(SingSuptable, 'y1-y5 ~ 1');
ranovatbl.SingSup = ranova(rm.SingSup);
MC.SingSup = multcompare(rm.SingSup, 'Time', 'ComparisonType', 'lsd');
disp('SingSup Time Avg = '); 
mean(table2array(SingSuptable), 1)
disp('SingSupTime Std = '); 
std(table2array(SingSuptable), 1)
H = zeros(1, n); 
for i = 1:n
    H(i) = swtest(table2array(SingSuptable(:,i)));
end
disp('SingSup Time Normality = '); 
disp(H);
subplot(153); hold on; 
plot(x, [SingSup.y1], '.', 'MarkerSize',mkr, 'Color',TotalColors(1,:)); 
plot(x*2, [SingSup.y2], '.', 'MarkerSize',mkr, 'Color',TotalColors(2,:)); 
plot(x*3, [SingSup.y3], '.', 'MarkerSize',mkr, 'Color',TotalColors(3,:)); 
plot(x*4, [SingSup.y4], '.', 'MarkerSize',mkr, 'Color',TotalColors(4,:)); 
plot(x*5, [SingSup.y5], '.', 'MarkerSize',mkr, 'Color',TotalColors(5,:)); 
text(4,max([SingSup.y4]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(4,:), 'HorizontalAlignment','center'); 
text(5,max([SingSup.y5]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(5,:), 'HorizontalAlignment','center'); 
title('Single Limb Support'); 
xlim([0.5 5.5]); 
ylim([0.25 0.4]); 
ylabel('% Gait Cycle'); 
ax = gca; 
ax.XTickLabel = {'-40','-20','Norm','+20','+40'};
ax.YTickLabel = [25 30 35 40];

DS2table = struct2table(DS2); 
[rm.DS2] = fitrm(DS2table, 'y1-y5 ~ 1');
ranovatbl.DS2 = ranova(rm.DS2);
MC.DS2 = multcompare(rm.DS2, 'Time', 'ComparisonType', 'lsd');
disp('DS2 Time Avg = '); 
mean(table2array(DS2table), 1)
disp('DS2 Time Std = '); 
std(table2array(DS2table), 1)
H = zeros(1, n); 
for i = 1:n
    H(i) = swtest(table2array(DS2table(:,i)));
end
disp('DS2 Time Normality = '); 
disp(H);
subplot(154); hold on; 
plot(x, [DS2.y1], '.', 'MarkerSize',mkr, 'Color',TotalColors(1,:)); 
plot(x*2, [DS2.y2], '.', 'MarkerSize',mkr, 'Color',TotalColors(2,:)); 
plot(x*3, [DS2.y3], '.', 'MarkerSize',mkr, 'Color',TotalColors(3,:)); 
plot(x*4, [DS2.y4], '.', 'MarkerSize',mkr, 'Color',TotalColors(4,:)); 
plot(x*5, [DS2.y5], '.', 'MarkerSize',mkr, 'Color',TotalColors(5,:)); 
text(4,max([DS2.y4]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(4,:), 'HorizontalAlignment','center'); 
text(5,max([DS2.y5]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(5,:), 'HorizontalAlignment','center'); 
title({'Double Support'; 'Trailing Limb'}); 
xlim([0.5 5.5]); 
ylim([0.1 0.25]); 
ylabel('% Gait Cycle'); 
ax = gca; 
ax.XTickLabel = {'-40','-20','Norm','+20','+40'};
ax.YTickLabel = [10 15 20 25];

Swingtable = struct2table(Swing); 
[rm.Swing] = fitrm(Swingtable, 'y1-y5 ~ 1');
ranovatbl.Swing = ranova(rm.Swing);
MC.Swing = multcompare(rm.Swing, 'Time', 'ComparisonType', 'lsd');
disp('Swing Time Avg = '); 
mean(table2array(Swingtable), 1)
disp('Swing Time Std = '); 
std(table2array(Swingtable), 1)
H = zeros(1, n); 
for i = 1:n
    H(i) = swtest(table2array(Swingtable(:,i)));
end
disp('Swing Time Normality = '); 
disp(H);
subplot(155); hold on; 
plot(x, [Swing.y1], '.', 'MarkerSize',mkr, 'Color',TotalColors(1,:)); 
plot(x*2, [Swing.y2], '.', 'MarkerSize',mkr, 'Color',TotalColors(2,:)); 
plot(x*3, [Swing.y3], '.', 'MarkerSize',mkr, 'Color',TotalColors(3,:)); 
plot(x*4, [Swing.y4], '.', 'MarkerSize',mkr, 'Color',TotalColors(4,:)); 
plot(x*5, [Swing.y5], '.', 'MarkerSize',mkr, 'Color',TotalColors(5,:)); 
text(4,max([Swing.y4]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(4,:), 'HorizontalAlignment','center'); 
text(5,max([Swing.y5]) + 0.015, '*', 'FontSize', txt, 'Color',TotalColors(5,:), 'HorizontalAlignment','center'); 
title('Swing*'); 
xlim([0.5 5.5]); 
ylim([0.25 0.4]); 
ylabel('% Gait Cycle'); 
ax = gca; 
ax.XTickLabel = {'-40','-20','Norm','+20','+40'};
ax.YTickLabel = [25 30 35 40];

% ---------------------------------------------------------------
% met cost
MetabolicCostTable = struct2table(MetabolicCost); 
[rm.MetCost] = fitrm(MetabolicCostTable, 'y1-y5 ~ 1');
ranovatbl.MetCost = ranova(rm.MetCost);
MC.MetCost = multcompare(rm.MetCost, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 1:n
    H(i) = swtest(table2array(MetabolicCostTable(:,i)));
end
disp('MetCost Normality = '); 
disp(H); 

% Total
DS1Totaltable = struct2table(DS1TotalTable);
[rm.DS1Total] = fitrm(DS1Totaltable, 'y1-y5 ~ 1');
ranovatbl.DS1Total = ranova(rm.DS1Total);
MC.DS1Total = multcompare(rm.DS1Total, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS1Totaltable(:,i)));
end
disp('DS1 Total Normality = '); 
disp(H(3:7)); 

SingSupTotaltable = struct2table(SingSupTotalTable);
[rm.SingSupTotal] = fitrm(SingSupTotaltable, 'y1-y5~1');
[ranovatbl.SingSupTotal] = ranova(rm.SingSupTotal);
MC.SingSupTotal = multcompare(rm.SingSupTotal, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SingSupTotaltable(:,i)));
end
disp('SingSup Total Normality = '); 
disp(H(3:7)); 

DS2Totaltable = struct2table(DS2TotalTable);
[rm.DS2Total] = fitrm(DS2Totaltable, 'y1-y5~1');
[ranovatbl.DS2Total] = ranova(rm.DS2Total);
MC.DS2Total = multcompare(rm.DS2Total, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS2Totaltable(:,i)));
end
disp('DS2 Total Normality = '); 
disp(H(3:7)); 

SwingTotaltable = struct2table(SwingTotalTable);
[rm.SwingTotal] = fitrm(SwingTotaltable, 'y1-y5~1');
[ranovatbl.SwingTotal] = ranova(rm.SwingTotal);
MC.SwingTotal = multcompare(rm.SwingTotal, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SwingTotaltable(:,i)));
end
disp('Swing Total Normality = '); 
disp(H(3:7)); 

StanceTotaltable = struct2table(StanceTotalTable);
[rm.StanceTotal] = fitrm(StanceTotaltable, 'y1-y5~1');
[ranovatbl.StanceTotal] = ranova(rm.StanceTotal);
MC.StanceTotal = multcompare(rm.StanceTotal, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StanceTotaltable(:,i)));
end
disp('Stance Total Normality = '); 
disp(H(3:7)); 

StrideTotaltable = struct2table(StrideTotalTable);
[rm.StrideTotal] = fitrm(StrideTotaltable, 'y1-y5~1');
[ranovatbl.StrideTotal] = ranova(rm.StrideTotal);
MC.StrideTotal = multcompare(rm.StrideTotal, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StrideTotaltable(:,i)));
end
disp('Stride Total Normality = '); 
disp(H(3:7)); 

% Hip
DS1Hiptable = struct2table(DS1HipTable);
[rm.DS1Hip] = fitrm(DS1Hiptable, 'y1-y5 ~ 1');
ranovatbl.DS1Hip = ranova(rm.DS1Hip);
MC.DS1Hip = multcompare(rm.DS1Hip, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS1Hiptable(:,i)));
end
disp('DS1 Hip Normality = '); 
disp(H(3:7)); 

SingSupHiptable = struct2table(SingSupHipTable);
[rm.SingSupHip] = fitrm(SingSupHiptable, 'y1-y5~1');
[ranovatbl.SingSupHip] = ranova(rm.SingSupHip);
MC.SingSupHip = multcompare(rm.SingSupHip, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SingSupHiptable(:,i)));
end
disp('SingSup Hip Normality = '); 
disp(H(3:7)); 

DS2Hiptable = struct2table(DS2HipTable);
[rm.DS2Hip] = fitrm(DS2Hiptable, 'y1-y5~1');
[ranovatbl.DS2Hip] = ranova(rm.DS2Hip);
MC.DS2Hip = multcompare(rm.DS2Hip, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS2Hiptable(:,i)));
end
disp('DS2 Hip Normality = '); 
disp(H(3:7)); 

SwingHiptable = struct2table(SwingHipTable);
[rm.SwingHip] = fitrm(SwingHiptable, 'y1-y5~1');
[ranovatbl.SwingHip] = ranova(rm.SwingHip);
MC.SwingHip = multcompare(rm.SwingHip, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SwingHiptable(:,i)));
end
disp('Swing Hip Normality = '); 
disp(H(3:7)); 

StanceHiptable = struct2table(StanceHipTable);
[rm.StanceHip] = fitrm(StanceHiptable, 'y1-y5~1');
[ranovatbl.StanceHip] = ranova(rm.StanceHip);
MC.StanceHip = multcompare(rm.StanceHip, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StanceHiptable(:,i)));
end
disp('Stance Hip Normality = '); 
disp(H(3:7)); 

StrideHiptable = struct2table(StrideHipTable);
[rm.StrideHip] = fitrm(StrideHiptable, 'y1-y5~1');
[ranovatbl.StrideHip] = ranova(rm.StrideHip);
MC.StrideHip = multcompare(rm.StrideHip, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StrideHiptable(:,i)));
end
disp('Stride Hip Normality = '); 
disp(H(3:7)); 

% Knee
DS1Kneetable = struct2table(DS1KneeTable);
[rm.DS1Knee] = fitrm(DS1Kneetable, 'y1-y5 ~ 1');
ranovatbl.DS1Knee = ranova(rm.DS1Knee);
MC.DS1Knee = multcompare(rm.DS1Knee, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS1Kneetable(:,i)));
end
disp('DS1 Knee Normality = '); 
disp(H(3:7)); 

SingSupKneetable = struct2table(SingSupKneeTable);
[rm.SingSupKnee] = fitrm(SingSupKneetable, 'y1-y5~1');
[ranovatbl.SingSupKnee] = ranova(rm.SingSupKnee);
MC.SingSupKnee = multcompare(rm.SingSupKnee, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SingSupKneetable(:,i)));
end
disp('SingSup Knee Normality = '); 
disp(H(3:7)); 

DS2Kneetable = struct2table(DS2KneeTable);
[rm.DS2Knee] = fitrm(DS2Kneetable, 'y1-y5~1');
[ranovatbl.DS2Knee] = ranova(rm.DS2Knee);
MC.DS2Knee = multcompare(rm.DS2Knee, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS2Kneetable(:,i)));
end
disp('DS2 Knee Normality = '); 
disp(H(3:7)); 

SwingKneetable = struct2table(SwingKneeTable);
[rm.SwingKnee] = fitrm(SwingKneetable, 'y1-y5~1');
[ranovatbl.SwingKnee] = ranova(rm.SwingKnee);
MC.SwingKnee = multcompare(rm.SwingKnee, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SwingKneetable(:,i)));
end
disp('Swing Knee Normality = '); 
disp(H(3:7)); 

StanceKneetable = struct2table(StanceKneeTable);
[rm.StanceKnee] = fitrm(StanceKneetable, 'y1-y5~1');
[ranovatbl.StanceKnee] = ranova(rm.StanceKnee);
MC.StanceKnee = multcompare(rm.StanceKnee, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StanceKneetable(:,i)));
end
disp('Stance Knee Normality = '); 
disp(H(3:7)); 

StrideKneetable = struct2table(StrideKneeTable);
[rm.StrideKnee] = fitrm(StrideKneetable, 'y1-y5~1');
[ranovatbl.StrideKnee] = ranova(rm.StrideKnee);
MC.StrideKnee = multcompare(rm.StrideKnee, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StrideKneetable(:,i)));
end
disp('StrideKnee Normality = '); 
disp(H(3:7)); 

% Ankle
DS1Ankletable = struct2table(DS1AnkleTable);
[rm.DS1Ankle] = fitrm(DS1Ankletable, 'y1-y5~1');
[ranovatbl.DS1Ankle] = ranova(rm.DS1Ankle);
MC.DS1Ankle = multcompare(rm.DS1Ankle, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS1Ankletable(:,i)));
end
disp('DS1 Ankle Normality = '); 
disp(H(3:7)); 

SingSupAnkletable = struct2table(SingSupAnkleTable);
[rm.SingSupAnkle] = fitrm(SingSupAnkletable, 'y1-y5~1');
[ranovatbl.SingSupAnkle] = ranova(rm.SingSupAnkle);
MC.SingSupAnkle = multcompare(rm.SingSupAnkle, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SingSupAnkletable(:,i)));
end
disp('SingSup Ankle Normality = '); 
disp(H(3:7)); 

DS2Ankletable = struct2table(DS2AnkleTable);
[rm.DS2Ankle] = fitrm(DS2Ankletable, 'y1-y5~1');
[ranovatbl.DS2Ankle] = ranova(rm.DS2Ankle);
MC.DS2Ankle = multcompare(rm.DS2Ankle, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(DS2Ankletable(:,i)));
end
disp('DS2 Ankle Normality = '); 
disp(H(3:7)); 

SwingAnkletable = struct2table(SwingAnkleTable);
[rm.SwingAnkle] = fitrm(SwingAnkletable, 'y1-y5~1');
[ranovatbl.SwingAnkle] = ranova(rm.SwingAnkle);
MC.SwingAnkle = multcompare(rm.SwingAnkle, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(SwingAnkletable(:,i)));
end
disp('Swing Ankle Normality = '); 
disp(H(3:7)); 

StanceAnkletable = struct2table(StanceAnkleTable);
[rm.StanceAnkle] = fitrm(StanceAnkletable, 'y1-y5~1');
[ranovatbl.StanceAnkle] = ranova(rm.StanceAnkle);
MC.StanceAnkle = multcompare(rm.StanceAnkle, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StanceAnkletable(:,i)));
end
disp('Stance Ankle Normality = '); 
disp(H(3:7)); 

StrideAnkletable = struct2table(StrideAnkleTable);
[rm.StrideAnkle] = fitrm(StrideAnkletable, 'y1-y5~1');
[ranovatbl.StrideAnkle] = ranova(rm.StrideAnkle);
MC.StrideAnkle = multcompare(rm.StrideAnkle, 'Time', 'ComparisonType', 'lsd');
H = zeros(1, n); 
for i = 3:n+2
    H(i) = swtest(table2array(StrideAnkletable(:,i)));
end
disp('Stride Ankle Normality = '); 
disp(H(3:7)); 

% save stride time fig
print('StrideTimeFig','-djpeg','-r300');
% close(StrideTimeFig); 

%% Non parametric tests
% stride measures
[rmNP.StrideTime_p, ~, rmNP.StrideTime] = friedman(table2array(StrideTimetable),1, 'off'); 
MCNP.StrideTime = multcompare(rmNP.StrideTime, 'Display','off','CType', 'lsd');
[ rmNP.DS1_p, ~, rmNP.DS1] = friedman(table2array(DS1table),1, 'off'); 
MCNP.DS1 = multcompare(rmNP.DS1, 'Display','off','CType', 'lsd');
[ rmNP.SingSup_p, ~, rmNP.SingSup] = friedman(table2array(SingSuptable),1, 'off'); 
MCNP.SingSup = multcompare(rmNP.SingSup, 'Display','off','CType', 'lsd');
[ rmNP.DS2_p, ~, rmNP.DS2] = friedman(table2array(DS2table),1, 'off'); 
MCNP.DS2 = multcompare(rmNP.DS2, 'Display','off','CType', 'lsd');
[ rmNP.Swing_p, ~, rmNP.Swing] = friedman(table2array(Swingtable),1, 'off'); 
MCNP.Swing = multcompare(rmNP.Swing, 'Display','off','CType', 'lsd');

% met cost
[ rmNP.MetCost_p, ~, rmNP.MetCost] = friedman(table2array(MetabolicCostTable),1, 'off'); 
MCNP.MetCost = multcompare(rmNP.MetCost, 'Display','off','CType', 'lsd');

% joint level powers
% total
[ rmNP.DS1Total_p, ~, rmNP.DS1Total] = friedman(table2array(DS1Totaltable(:,3:7)),1, 'off'); 
MCNP.DS1Total = multcompare(rmNP.DS1Total, 'Display','off','CType', 'lsd');
[ rmNP.SingSupTotal_p, ~, rmNP.SingSupTotal] = friedman(table2array(SingSupTotaltable(:,3:7)),1, 'off'); 
MCNP.SingSupTotal = multcompare(rmNP.SingSupTotal, 'Display','off','CType', 'lsd');
[ rmNP.DS2Total_p, ~, rmNP.DS2Total] = friedman(table2array(DS2Totaltable(:,3:7)),1, 'off'); 
MCNP.DS2Total = multcompare(rmNP.DS2Total, 'Display','off','CType', 'lsd');
[ rmNP.SwingTotal_p, ~, rmNP.SwingTotal] = friedman(table2array(SwingTotaltable(:,3:7)),1, 'off'); 
MCNP.SwingTotal = multcompare(rmNP.SwingTotal, 'Display','off','CType', 'lsd');
[ rmNP.StanceTotal_p, ~, rmNP.StanceTotal] = friedman(table2array(StanceTotaltable(:,3:7)),1, 'off'); 
MCNP.StanceTotal = multcompare(rmNP.StanceTotal, 'Display','off','CType', 'lsd');
[ rmNP.StrideTotal_p, ~, rmNP.StrideTotal] = friedman(table2array(StrideTotaltable(:,3:7)),1, 'off'); 
MCNP.StrideTotal = multcompare(rmNP.StrideTotal, 'Display','off','CType', 'lsd');
% hip
[ rmNP.DS1Hip_p, ~, rmNP.DS1Hip] = friedman(table2array(DS1Hiptable(:,3:7)),1, 'off'); 
MCNP.DS1Hip = multcompare(rmNP.DS1Hip, 'Display','off','CType', 'lsd');
[ rmNP.SingSupHip_p, ~, rmNP.SingSupHip] = friedman(table2array(SingSupHiptable(:,3:7)),1, 'off'); 
MCNP.SingSupHip = multcompare(rmNP.SingSupHip, 'Display','off','CType', 'lsd');
[ rmNP.DS2Hip_p, ~, rmNP.DS2Hip] = friedman(table2array(DS2Hiptable(:,3:7)),1, 'off'); 
MCNP.DS2Hip = multcompare(rmNP.DS2Hip, 'Display','off','CType', 'lsd');
[ rmNP.SwingHip_p, ~, rmNP.SwingHip] = friedman(table2array(SwingHiptable(:,3:7)),1, 'off'); 
MCNP.SwingHip = multcompare(rmNP.SwingHip, 'Display','off','CType', 'lsd');
[ rmNP.StanceHip_p, ~, rmNP.StanceHip] = friedman(table2array(StanceHiptable(:,3:7)),1, 'off'); 
MCNP.StanceHip = multcompare(rmNP.StanceHip, 'Display','off','CType', 'lsd');
[ rmNP.StrideHip_p, ~, rmNP.StrideHip] = friedman(table2array(StrideHiptable(:,3:7)),1, 'off'); 
MCNP.StrideHip = multcompare(rmNP.StrideHip, 'Display','off','CType', 'lsd');

% knee
[ rmNP.DS1Knee_p, ~, rmNP.DS1Knee] = friedman(table2array(DS1Kneetable(:,3:7)),1, 'off'); 
MCNP.DS1Knee = multcompare(rmNP.DS1Knee, 'Display','off','CType', 'lsd');
[ rmNP.SingSupKnee_p, ~, rmNP.SingSupKnee] = friedman(table2array(SingSupKneetable(:,3:7)),1, 'off'); 
MCNP.SingSupKnee = multcompare(rmNP.SingSupKnee, 'Display','off','CType', 'lsd');
[ rmNP.DS2Knee_p, ~, rmNP.DS2Knee] = friedman(table2array(DS2Kneetable(:,3:7)),1, 'off'); 
MCNP.DS2Knee = multcompare(rmNP.DS2Knee, 'Display','off','CType', 'lsd');
[ rmNP.SwingKnee_p, ~, rmNP.SwingKnee] = friedman(table2array(SwingKneetable(:,3:7)),1, 'off'); 
MCNP.SwingKnee = multcompare(rmNP.SwingKnee, 'Display','off','CType', 'lsd');
[ rmNP.StanceKnee_p, ~, rmNP.StanceKnee] = friedman(table2array(StanceKneetable(:,3:7)),1, 'off'); 
MCNP.StanceKnee = multcompare(rmNP.StanceKnee, 'Display','off','CType', 'lsd');
[ rmNP.StrideKnee_p, ~, rmNP.StrideKnee] = friedman(table2array(StrideKneetable(:,3:7)),1, 'off'); 
MCNP.StrideKnee = multcompare(rmNP.StrideKnee, 'Display','off','CType', 'lsd');

% ankle
[ rmNP.DS1Ankle_p, ~, rmNP.DS1Ankle] = friedman(table2array(DS1Ankletable(:,3:7)),1, 'off'); 
MCNP.DS1Ankle = multcompare(rmNP.DS1Ankle, 'Display','off','CType', 'lsd');
[ rmNP.SingSupAnkle_p, ~, rmNP.SingSupAnkle] = friedman(table2array(SingSupAnkletable(:,3:7)),1, 'off'); 
MCNP.SingSupAnkle = multcompare(rmNP.SingSupAnkle, 'Display','off','CType', 'lsd');
[ rmNP.DS2Ankle_p, ~, rmNP.DS2Ankle] = friedman(table2array(DS2Ankletable(:,3:7)),1, 'off'); 
MCNP.DS2Ankle = multcompare(rmNP.DS2Ankle, 'Display','off','CType', 'lsd');
[ rmNP.SwingAnkle_p, ~, rmNP.SwingAnkle] = friedman(table2array(SwingAnkletable(:,3:7)),1, 'off'); 
MCNP.SwingAnkle = multcompare(rmNP.SwingAnkle, 'Display','off','CType', 'lsd');
[ rmNP.StanceAnkle_p, ~, rmNP.StanceAnkle] = friedman(table2array(StanceAnkletable(:,3:7)),1, 'off'); 
MCNP.StanceAnkle = multcompare(rmNP.StanceAnkle, 'Display','off','CType', 'lsd');
[ rmNP.StrideAnkle_p, ~, rmNP.StrideAnkle] = friedman(table2array(StrideAnkletable(:,3:7)),1, 'off'); 
MCNP.StrideAnkle = multcompare(rmNP.StrideAnkle, 'Display','off','CType', 'lsd');

%% Effect Size Statistics
% eta squared
cols = {'y1','y2','y3','y4','y5'}; 

% met cost
MES.MetCost = mes1way(MetabolicCostTable{:,cols} , 'eta2' ); 

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



%% Total metabolic cost by condition
% figure; hold on;
% MkrSz = 10;
DotMargin = 0.2;
LW = 1.25;
X = [1 2 3 4 5];

% Identify trials
Cond_Ind(:,1) = contains([MetCostTable.Trial], 'Fm40');
Cond_Ind(:,2) = contains([MetCostTable.Trial], 'Fm20');
Cond_Ind(:,3) = contains([MetCostTable.Trial], 'Norm');
Cond_Ind(:,4) = contains([MetCostTable.Trial], 'Fp20');
Cond_Ind(:,5) = contains([MetCostTable.Trial], 'Fp40');

for cond = 1:5
    subplot(331); hold on;
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
    
    subplot(332); hold on;
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
    
    subplot(333); hold on;
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

subplot(331);
xlim([0.5 5.5]);
ylim([2 10]); 
title('Empirical', 'FontSize', TitleFnt);
ax = gca;
ax.XTick = [ 1 2 3 4 5];
ax.XTickLabel = {'-40', '-20','Norm','+20','+40'};
ax.FontSize = AxFnt;
ylabel('W / kg');

subplot(332);
xlim([0.5 5.5]);
ylim([2 10]);
title('Umberger Model', 'FontSize', TitleFnt);
ax = gca;
ax.XTick = [ 1 2 3 4 5];
ax.XTickLabel = {'-40', '-20','Norm','+20','+40'};
ax.FontSize = AxFnt;
ylabel('W / kg');

subplot(333);
xlim([0.5 5.5]);
ylim([2 10]);
title('Bhargava Model', 'FontSize', TitleFnt);
ax = gca;
ax.XTick = [ 1 2 3 4 5];
ax.XTickLabel = {'-40', '-20','Norm','+20','+40'};
ax.FontSize = AxFnt;
ylabel('W / kg');

% measured met cost ANOVA and post hoc
MesCostTbl = table(MetCostTable.Measured(Cond_Ind(:,1)),... 
    MetCostTable.Measured(Cond_Ind(:,2)),... 
    MetCostTable.Measured(Cond_Ind(:,3)),... 
    MetCostTable.Measured(Cond_Ind(:,4)),... 
    MetCostTable.Measured(Cond_Ind(:,5)),... 
    'VariableNames', {'y1','y2','y3','y4','y5'});
[rm.MesCost] = fitrm(MesCostTbl, 'y1-y5 ~ 1');
ranovatbl.MesCost = ranova(rm.MesCost);
MC.MesCost = multcompare(rm.MesCost, 'Time', 'ComparisonType', 'lsd');

% umberger met cost ANOVA and post hoc
UmbCostTbl = table(MetCostTable.UmbModelled(Cond_Ind(:,1)),... 
    MetCostTable.UmbModelled(Cond_Ind(:,2)),... 
    MetCostTable.UmbModelled(Cond_Ind(:,3)),... 
    MetCostTable.UmbModelled(Cond_Ind(:,4)),... 
    MetCostTable.UmbModelled(Cond_Ind(:,5)),... 
    'VariableNames', {'y1','y2','y3','y4','y5'});
[rm.UmbCost] = fitrm(UmbCostTbl, 'y1-y5 ~ 1');
ranovatbl.UmbCost = ranova(rm.UmbCost);
MC.UmbCost = multcompare(rm.UmbCost, 'Time', 'ComparisonType', 'lsd');

% bhargava met cost ANOVA and post hoc
BharCostTbl = table(MetCostTable.BharModelled(Cond_Ind(:,1)),... 
    MetCostTable.BharModelled(Cond_Ind(:,2)),... 
    MetCostTable.BharModelled(Cond_Ind(:,3)),... 
    MetCostTable.BharModelled(Cond_Ind(:,4)),... 
    MetCostTable.BharModelled(Cond_Ind(:,5)),... 
    'VariableNames', {'y1','y2','y3','y4','y5'});
[rm.BharCost] = fitrm(BharCostTbl, 'y1-y5 ~ 1');
ranovatbl.BharCost = ranova(rm.BharCost);
MC.BharCost = multcompare(rm.BharCost, 'Time', 'ComparisonType', 'lsd');

% plot asterisks above 
alpha = 0.05;
Asterisk = 20; 
TopMargin = 0.5;
% measured
subplot(331); 
y = max([MetCostTable.Measured(Cond_Ind(:,1))]) + TopMargin;
text(1, y, '*', 'Color', TotalColors(1,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
y = max([MetCostTable.Measured(Cond_Ind(:,2))]) + TopMargin;
text(2, y, '*', 'Color', TotalColors(2,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
y = max([MetCostTable.Measured(Cond_Ind(:,4))]) + TopMargin;
text(4, y, '*', 'Color', TotalColors(4,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
y = max([MetCostTable.Measured(Cond_Ind(:,5))]) + TopMargin;
text(5, y, '*', 'Color', TotalColors(5,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
text(0, 11, '\bf A', 'HorizontalAlignment','center')

% umberger
subplot(332); 
y = max([MetCostTable.UmbModelled(Cond_Ind(:,1))]) + TopMargin;
text(1, y, '*', 'Color', TotalColors(1,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
y = max([MetCostTable.UmbModelled(Cond_Ind(:,5))]) + TopMargin;
text(5, y, '*', 'Color', TotalColors(5,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
text(0, 11, '\bf B', 'HorizontalAlignment','center')

% bhargava
subplot(333); 
y = max([MetCostTable.BharModelled(Cond_Ind(:,1))]) + TopMargin;
text(1, y, '*', 'Color', TotalColors(1,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
y = max([MetCostTable.BharModelled(Cond_Ind(:,5))]) + TopMargin;
text(5, y, '*', 'Color', TotalColors(5,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
text(0, 11, '\bf C', 'HorizontalAlignment','center')

%% show total metabolic cost across the gait cycle
LW = 2; 
subplot(3,3,[7 8 9]); 
TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
HipColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
KneeColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
AnkleColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];

FntSz = 12; 
% JointEnergyCostsBothModels = figure('Position',[50 50 1200 800]);
TxtHt = 11.25; 
Stride = 1:100;
c = 0;
L = 0;
yPk = 12;
TotalyPk = 30;

GaitEventTable = struct2table(GaitEvents); 

Cond_Ind(:,1) = contains([GaitEventTable.Trial], 'Fm40');
Cond_Ind(:,2) = contains([GaitEventTable.Trial], 'Fm20');
Cond_Ind(:,3) = contains([GaitEventTable.Trial], 'Norm');
Cond_Ind(:,4) = contains([GaitEventTable.Trial], 'Fp20');
Cond_Ind(:,5) = contains([GaitEventTable.Trial], 'Fp40');

for i = [2 4 1 5 3]
    
    % total
    hold on; % grid on;
    plot(Stride, Total(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', TotalColors(i,:));
    plot(Stride, Total(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', TotalColors(i,:));
    
end

for i = [2 4 1 5 3]
    H = vline(100*nanmean(GaitEventTable.Off_Avg1(Cond_Ind(:,i))), '-');
    H.LineWidth = LW/2;
    H.Color =  TotalColors(i,:);
    
        H = vline(100*nanmean(GaitEventTable.OppOff_Avg1(Cond_Ind(:,i))), '--');
    H.LineWidth = LW/2;
    H.Color =  TotalColors(i,:);
    
        H = vline(100*nanmean(GaitEventTable.OppOn_Avg1(Cond_Ind(:,i))), '--');
    H.LineWidth = LW/2;
    H.Color =  TotalColors(i,:);
    
end
    ylabel('W / kg'); ylim([0 15]);
    xlabel('% Gait Cycle'); 
    title('Total Modelled Metabolic Power'); 
    text(-5, 17, '\bf G', 'HorizontalAlignment','center')
    
    text(100, 14, ['\bf' char(8211) ' Umberger'], 'HorizontalAlignment','right','FontSize',TxtFnt);
    text(100, 12, '\bf - - Bhargava', 'HorizontalAlignment','right','FontSize',TxtFnt);
    %             text(100, 10, '\bf | Gait Events', 'HorizontalAlignment','right','FontSize',TxtFnt);

subplotsqueeze(CorrVal, 1.05);
% saveas(CorrVal, 'CorrVal.png');
% print('StrideTimeFig','-dpng','-r300');
print('CorrVal','-djpeg','-r300');

%% plot BOTH MODELS across the gait cycle 
% clc; close all;
% yPk = 400; % set y peak
% LW = 2; % set line width
% % TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray'); rgb('DarkSlateGray'); rgb('Black')];
% % HipColors = [rgb('PowderBlue'); rgb('DeepSkyBlue'); rgb('RoyalBlue');  rgb('Navy'); rgb('MidnightBlue')];
% % % KneeColors = [rgb('Orange'); rgb('DarkOrange'); rgb('Coral'); rgb('Chocolate'); rgb('Sienna')];
% % KneeColors = [rgb('SandyBrown'); rgb('Orange'); rgb('DarkOrange'); rgb('Chocolate'); rgb('Sienna')];
% % AnkleColors = [rgb('PaleGreen'); rgb('LimeGreen'); rgb('ForestGreen'); rgb('Green'); rgb('DarkGreen')];
% 
% TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
% HipColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
% KneeColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
% AnkleColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
% 
% FntSz = 12; 
% JointEnergyCostsBothModels = figure('Position',[50 50 1200 800]);
% TxtHt = 10; 
% Stride = 1:100;
% c = 0;
% L = 0;
% yPk = 10;
% TotalyPk = 15;
% for i = [2 4 1 5 3]
%     c = i;
%     
%     % total
%     subplot(421); hold on; % grid on;
%     plot(Stride, Total(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', TotalColors(c,:));
%     %     plot(Stride, Total(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
%     ylabel('W / kg'); ylim([0 TotalyPk]);
%     
%     
%     subplot(422); hold on; %grid on;
%     %     plot(Stride, Total(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
%     plot(Stride, Total(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', TotalColors(c,:));
%     %     ylabel('W / kg');
%     ylim([0 TotalyPk]);
%     
%     
%     % Hip
%     subplot(423); hold on;% grid on;
%     plot(Stride, Hip(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', HipColors(c,:));
%     %     plot(Stride, Hip(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
%     ylabel('W / kg'); ylim([0 yPk]);
%    
%     subplot(424); hold on;% grid on;
%     %     plot(Stride, Hip(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
%     plot(Stride, Hip(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', HipColors(c,:));
%     %     ylabel('W / kg');
%     ylim([0 yPk]);
%     
%     % Knee
%     subplot(425); hold on; %grid on;
%     plot(Stride, Knee(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', KneeColors(c,:));
%     %      plot(Stride, Knee(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
%     ylabel('W / kg'); ylim([0 yPk]);
%     
%     subplot(426); hold on; %grid on;
%     %     plot(Stride, Knee(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
%     plot(Stride, Knee(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', KneeColors(c,:));
%     %     ylabel('W / kg');
%     ylim([0 yPk]);
%     
%     % Ankle
%     subplot(427); hold on; %grid on;
%     plot(Stride, Ankle(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', AnkleColors(c,:));
%     %     plot(Stride, Ankle(i).BHAR_TotalAvg, '--r', 'LineWidth', LW, 'Color', Colors(c,:));
%     ylabel('W / kg'); ylim([0 yPk]);
%     xlabel('% Gait Cycle');
%     
%     subplot(428); hold on; %grid on;
%     %     plot(Stride, Ankle(i).UMB_TotalAvg, 'r', 'LineWidth', LW, 'Color', Colors(c,:));
%     plot(Stride, Ankle(i).BHAR_TotalAvg, '-r', 'LineWidth', LW, 'Color', AnkleColors(c,:));
%     %     ylabel('W / kg');
%     ylim([0 yPk]);
%     xlabel('% Gait Cycle');
%     
% %     L = L + 1;
% %     %     CurveNameUMB{L} = strcat(Hip(i).Condition, ' UMB');
% %     %     CurveNameBHAR{L} = strcat(Hip(i).Condition, ' BHAR');
%     CurveName{i} = Hip(i).Condition;
% end
% 
% % CurveName = {'-40','-20','Norm','+20','+40'}; 
% subplot(421);
%     text(50, 15, '\bf Total UMB', 'Color',TotalColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% ax.XTickLabel = [];
% ax.FontSize = FntSz; 
% 
% subplot(422);
%      text(50, 15, '\bf Total BHAR', 'Color',TotalColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best');
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% ax.XTickLabel = [];
% % ax.YTick = [0 10 20 30];
% ax.YTickLabel = [];
% ax.FontSize = FntSz; 
% 
% subplot(423);
%     text(50, TxtHt, '\bf Hip UMB', 'Color',HipColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% ax.XTickLabel = [];
% ax.YTick = [0 2 4 6 8 10];
% ax.FontSize = FntSz; 
% 
% subplot(424);
%      text(50, TxtHt, '\bf Hip BHAR', 'Color',HipColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best');
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% ax.XTickLabel = [];
% ax.YTick = [0 2 4 6 8 10];
% ax.YTickLabel = [];
% ax.FontSize = FntSz; 
% 
% subplot(425);
%    text(50, TxtHt,'\bf Knee UMB', 'Color',KneeColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% ax.YTick = [0 2 4 6 8 10];
% ax.XTickLabel = [];
% ax.FontSize = FntSz; 
% 
% subplot(426);
%  text(50, TxtHt,'\bf Knee BHAR', 'Color',KneeColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best');
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% ax.XTickLabel = [];
% ax.YTick = [0 2 4 6 8 10];
% ax.YTickLabel = [];
% ax.FontSize = FntSz; 
% 
% subplot(427);
% text(50, TxtHt, '\bf Ankle UMB', 'Color',AnkleColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best', 'FontSize',FntSz - 3);
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% % % ax.XTickLabel = [];
% ax.YTick = [0 2 4 6 8 10];
% ax.FontSize = FntSz; 
% 
% subplot(428);
% text(50, TxtHt,'\bf Ankle BHAR', 'Color',AnkleColors(3,:), 'HorizontalAlignment','center', 'FontSize',FntSz);
% % legend(CurveName, 'Location','Best');
% ax = gca;
% ax.XTick = [0 25 50 75 100];
% % ax.XTickLabel = [];
% ax.YTick = [0 2 4 6 8 10];
% ax.YTickLabel = [];
% ax.FontSize = FntSz; 
% 
% subplotsqueeze(JointEnergyCostsBothModels, 1.2);
% saveas(JointEnergyCostsBothModels, 'JointEnergyCostsBothModels.png');
% % saveas(JointEnergyCostsBothModels, 'JointEnergyCostsBothModels.pdf');
% 
% clearvars ax AxFnt c i ind int L LW Mescond MkrSz subj trial TxtFnt TotalyPk UMBcond BHARcond X yPk


%% Stride-Stance-Swing Metabolic Costs
clc;
close all;

TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
HipColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
KneeColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
AnkleColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];

CycleEnergyCostsUMB = figure('Position',[50 50 1000 800]);
hold on;
x = [1 2 3 4 5];
LW = 1.25;
DotMargin = 0.3;
MeanMargin = 0;
StrideOffset = 0;
% StanceOffset = 6;
DS1Offset = 6;
SingSupOffset = 12;
DS2Offset = 18;
SwingOffset = 24;

TxtFnt = 10;
CapEndSz = 0.1;
SmlMkr = 8;

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
ypos = 11;
% text(StanceOffset + 3,ypos+0.1, {'\bf Stance Phase'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(StrideOffset + 3,ypos+0.1, {'\bf Entire Stride'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(DS1Offset + 3,ypos+0.1, {'\bf Leading Limb'; '\bf Double Support'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(SingSupOffset + 3,ypos+0.1, {'\bf Single Support'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(DS2Offset + 3,ypos+0.1, {'\bf Trailing Limb'; '\bf Double Support'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(SwingOffset + 3,ypos+0.1, {'\bf Swing Phase'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
text(36,ypos+0.1, {'\bf Gait Cycle'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');

% ticks label
for j = 1:4
    subplot(4,1,j);
    %     TickLabel = {'-40%', '-20%','Norm','+20%','+40%'};
    ylabel('W / kg');
%     xlim([0.5 29.5]);
    ax = gca;
  ax.XTick = [1 2 3 4 5 7 8 9 10 11 13 14 15 16 17 19 20 21 22 23 25 26 27 28 29 ...
    30 33 36 39 42];
    ax.XTickLabels = [];
    ax.FontSize = 7.5;
end

% figure sub labels A,B,C,D and Tick label at bottom
x = -1.5; 
TxtFnt = 12; 
subplot(4,1,1);
text(x,11, {'\bf A'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
subplot(4,1,2); 
text(x,4.6, {'\bf B'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
subplot(4,1,3); 
text(x,2.7, {'\bf C'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
subplot(4,1,4);
text(x,2.3, {'\bf D'}, 'FontSize', TxtFnt, 'HorizontalAlignment', 'center');
TickLabel = {'-40', '-20','Norm','+20','+40', '-40', '-20','Norm','+20','+40',...
    '-40', '-20','Norm','+20','+40', '-40', '-20','Norm','+20','+40',...
    '-40', '-20','Norm','+20','+40', '0%', '25%','50%','75%','100%'};
ax = gca;
ax.XTick = [1 2 3 4 5 7 8 9 10 11 13 14 15 16 17 19 20 21 22 23 25 26 27 28 29 ...
    30 33 36 39 42];
ax.XTickLabels = TickLabel;
ax.YTick = [0 1 2 3];
ax.FontSize = 7.5;


% Set y axis boundaries
MuscLabelCent = 15;
subplot(4,1,1);
ylim([0 10]);
vline(6, 'k');
text(MuscLabelCent, 9.3, '\bf ALL MUSCLES', 'Color', TotalColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');
subplot(4,1,2);
ylim([0 4.5]);
vline(6, 'k');
text(MuscLabelCent, 4.2, '\bf HIP MUSCLES', 'Color',HipColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');
subplot(4,1,3);
ylim([0 2.5]);
vline(6, 'k');
text(MuscLabelCent, 2.3, '\bf KNEE MUSCLES', 'Color', KneeColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');
subplot(4,1,4);
ylim([0 2.2]);
vline(6, 'k');
text(MuscLabelCent, 2, '\bf ANKLE MUSCLES', 'Color', AnkleColors(3,:), 'FontSize',12, 'HorizontalAlignment', 'center');


%% Add ASTERISKS
alpha = 0.05;
Asterisk = 20; 
% total
subplot(411);
TopMargin = 0.5;
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
Log = 4;
y = max([Total(Log).UMB_Work.DS1_Tot]) + TopMargin;
text(Log+DS1Offset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
  Log = 5;
y = max([Total(Log).UMB_Work.DS1_Tot]) + TopMargin;
text(Log+DS1Offset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');

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
Log = 1;
y = max([Total(Log).UMB_Work.SingSup_Tot]) + TopMargin;
text(Log+SingSupOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
Log = 4;
y = max([Total(Log).UMB_Work.SingSup_Tot]) + TopMargin;
text(Log+SingSupOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
Log = 5;
y = max([Total(Log).UMB_Work.SingSup_Tot]) + TopMargin;
text(Log+SingSupOffset, y, '*', 'Color', TotalColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');

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
Log = 1;
y = max([Hip(Log).UMB_Work.DS1_Tot]) + TopMargin;
text(Log+DS1Offset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
Log = 2;
y = max([Hip(Log).UMB_Work.DS1_Tot]) + TopMargin;
text(Log+DS1Offset, y, '*', 'Color', HipColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');

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
TopMargin = 0.2;
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
Log = 1;
y = max([Knee(Log).UMB_Work.SingSup_Tot]) + TopMargin;
text(Log+SingSupOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
Log = 4;
y = max([Knee(Log).UMB_Work.SingSup_Tot]) + TopMargin;
text(Log+SingSupOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');
Log = 5;
y = max([Knee(Log).UMB_Work.SingSup_Tot]) + TopMargin;
text(Log+SingSupOffset, y, '*', 'Color', KneeColors(Log,:), 'FontSize', Asterisk, 'HorizontalAlignment', 'center');

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
TopMargin = 0.2;
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


%% plot joint level curves across the gait cycle to the right
Stride = linspace(30,42); 

% add curves for each condition
for i = [2 4 1 5 3]
    
    % total
    subplot(411); yyaxis right
    plot(Stride, Total(i).UMB_TotalAvg,  '-r', 'LineWidth', LW, 'Color', TotalColors(i,:));
    
    %     % Hip
    subplot(412); yyaxis right
    plot(Stride, Hip(i).UMB_TotalAvg,  '-r', 'LineWidth', LW, 'Color', TotalColors(i,:));
    
    %     % Knee
    subplot(413); yyaxis right
    plot(Stride, Knee(i).UMB_TotalAvg,  '-r', 'LineWidth', LW, 'Color', TotalColors(i,:));
    
    %     % Ankle
    subplot(414); yyaxis right
    plot(Stride, Ankle(i).UMB_TotalAvg,  '-r', 'LineWidth', LW, 'Color', TotalColors(i,:));

end
w = 1.5; 

% draw thick vertical likes and adjust axes
subplot(411); 
xlim([0 42]); 
v = vline(30, '-k'); 
v.LineWidth = w;
ax = gca; 
ax.YColor = 'k';

subplot(412); 
xlim([0 42]); 
ax = gca; 
ax.YColor = 'k';
ax.YTick = [0 1 2 3 4];
ax.YLim = [0 4]; 
v = vline(30, '-k'); 
v.LineWidth = w;

subplot(413); 
xlim([0 42]); 
v = vline(30, '-k'); 
v.LineWidth = w;
ax = gca; 
ax.YColor = 'k';
ax.YTick = [0 1 2 3];

subplot(414); 
xlim([0 42]); 
v = vline(30, '-k'); 
v.LineWidth = w;
ax = gca; 
ax.YColor = 'k';

subplotsqueeze(CycleEnergyCostsUMB, 1.12);
% saveas(CycleEnergyCostsUMB, 'CycleEnergyCostsUMB.png');
% saveas(CycleEnergyCostsUMB, 'CycleEnergyCostsUMB.pdf');
print('CycleEnergyCostsUMB','-djpeg','-r300');


%% plot joint level curves 
% Stride = 1:100;
% 
% % create subaxes
% ax_Total = axes('Position',[.30 .85 .6 .08]);
% ax_Hip = axes('Position',[.30 .62 .6 .08]);
% ax_Knee = axes('Position',[.30 .41 .6 .08]);
% ax_Ankle = axes('Position',[.30 .19 .6 .08]);
% 
% % add curves for each condition
% for i = [2 4 1 5 3]
%     
%     % total
%     hold(ax_Total, 'on');
%     plot(ax_Total, Stride, Total(i).UMB_TotalAvg,  'r', 'LineWidth', LW, 'Color', TotalColors(i,:));
% %     plot(ax_Total, Stride, Total(i).UMB_TotalAvg,  'r', 'LineWidth', LW, 'Color', TotalColors(i,:));
% %     h = vline(100*SF(i).DS1Avg,'r-');%, 'Color',  TotalColors(i,:));
% %     h.color =   TotalColors(i,:);
%     
%     % Hip
%     hold(ax_Hip, 'on');
%     plot(ax_Hip, Stride, Hip(i).UMB_TotalAvg,  'r', 'LineWidth', LW, 'Color', TotalColors(i,:));
%     
%     % Knee
%     hold(ax_Knee, 'on');
%     plot(ax_Knee, Stride, Knee(i).UMB_TotalAvg,  'r', 'LineWidth', LW, 'Color', TotalColors(i,:));
%     
%     % Ankle
%     hold(ax_Ankle, 'on');
%     plot(ax_Ankle, Stride, Ankle(i).UMB_TotalAvg,  'r', 'LineWidth', LW, 'Color', TotalColors(i,:));
%     
% end
% 
% % set X axis ticks
% ax_Total.XTick = [0 25 50 75 100];
% ax_Hip.XTick = [0 25 50 75 100];
% ax_Knee.XTick = [0 25 50 75 100];
% ax_Ankle.XTick = [0 25 50 75 100];
% 
% % additional y labels
% text(-35, 49, '\bf All Muscles', 'Rotation', 90, 'Fontsize',14); 
% text(-35, 30, '\bf Hip Muscles', 'Rotation', 90, 'Fontsize',14); 
% text(-35, 12, '\bf Knee Muscles', 'Rotation', 90, 'Fontsize',14); 
% text(-35, -7, '\bf Ankle Muscles', 'Rotation', 90, 'Fontsize',14); 
% 
% subplotsqueeze(CycleEnergyCostsUMB, 1.12);
% saveas(CycleEnergyCostsUMB, 'CycleEnergyCostsUMB.png');
% % saveas(CycleEnergyCostsUMB, 'CycleEnergyCostsUMB.pdf');

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

%% export Muscle Table
save('MuscleMetCosts.mat','MuscleMet'); 

% run SPM in Python

% {'semimem','semiten','bifemlh','bifemsh','sar','add_long','add_brev','tfl','pect','grac','iliacus','psoas','quad_fem','gem','peri','rect_fem','vas_med','vas_int','vas_lat','med_gas','lat_gas','soleus','tib_post','flex_dig','flex_hal','tib_ant','per_brev','per_long','per_tert','ext_dig','ext_hal','ercspn','intobl','extobl','glut_med','glut_min','glut_max','add_mag'}

%% Import muscle table after python
MetSPM = importdata('C:\Users\richa\Documents\Packages\OpenSim\Scripts\MuscleMetCostsSPM.csv');

% determine top 18 muscles to plot
[SUM, I] = sort(sum(Muscles(3).UMB_CombMuscleAvg, 1)); 
SortNames = Muscles(3).UMB_CombCols(I); 

%% PLOT MUSCLE METABOLICS
clc; close all;
LW = 2; % set line width
% TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray');  rgb('DarkSlateGray'); rgb('Black')];
TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
FntSz = 9;

MuscleMetCosts = figure('Position',[50 50 1200 1000]);
Stride = 1:100;
yPk = 2;
for i = [2 4 1 5 3]
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
    subplot(7,3,1); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'psoas');
    Pl(i).C = plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Psoas');
        ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
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
    subplot(7,3,2); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'iliacus');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Iliacus');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('iliacus', MetSPM, 2)
    end
    
     % rec_fem
    subplot(7,3,3); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'rect_fem');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Rectus Femoris');
    %     ylabel('W / kg'); 
    ylim([ 0 3]);
         xlim([0 100]); 
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
    
          % glut min
    subplot(7,3,4); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'glut_min');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Gluteus Minimus');
    ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca; 
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('glut_min', MetSPM, 2)
    end

    
    % glut med
    subplot(7,3,5); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'glut_med');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Gluteus Medius');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
    xlim([0 100]); 
    ax = gca; 
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('glut_med', MetSPM, 2)
    end
    
    
 % glut max
    subplot(7,3,6); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'glut_max');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Gluteus Maximus');
%     ylabel('W / kg'); 
    ylim([ 0 3]);
        xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%        ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('glut_max', MetSPM, 3)
    end
    
   
    
      % sar
    subplot(7,3,7); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'sar');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Sartorius');
        ylabel('W / kg'); 
    ylim([ 0 2]);
         xlim([0 100]); 
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
    
    % add_mag
    subplot(7,3,8); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'add_mag');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Adductor Magnus');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('add_mag', MetSPM, 2)
    end
    
     % add_long
    subplot(7,3,9); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'add_long');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Adductor Longus');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('add_long', MetSPM, 2)
    end
    
    
   
    
    % vas_lat
    subplot(7,3,10); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'vas_lat');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Vastus Lateralis');
        ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('vas_lat', MetSPM, 2)
    end
    
    % vas_med
    subplot(7,3,11); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'vas_med');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Vastus Medialis');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('vas_med', MetSPM, 2)
    end
    
    % vas_int
    subplot(7,3,12); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'vas_int');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Vastus Intermedius');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('vas_int', MetSPM, 2)
    end
    
    
    
     % bifemlh
    subplot(7,3,13); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'bifemlh');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Biceps Femoris - Long Head');
    ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('bifemlh', MetSPM, 2)
    end
    
   
    
    % semimem
    subplot(7,3,14); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'semimem');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Semimembranosus');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
    ax.YTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('semimem', MetSPM, 2)
    end
    
%     % semiten
    subplot(7,3,15); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'semiten');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Semitendinosus');
%     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('semiten', MetSPM, 2)
    end
    
     % bifemsh
    subplot(7,3,16); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'bifemsh');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Biceps Femoris - Short Head');
        ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTickLabel = [];
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('bifemsh', MetSPM, 2)
    end
    
      % tib_ant
    subplot(7,3,17); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'tib_ant');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Tibialis Anterior');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
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
    subplot(7,3,18); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'ext_dig');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Extensor Digitorum');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
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
    
    
    % lat_gas
    subplot(7,3,19); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'lat_gas');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Lateral Gastrocnemius');
    ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = {'0' '25' '50' '75' '100'};
    xlabel('% Gait Cycle');
    ax.FontSize = FntSz;
    % SPM shading
    if i == 5
        SPMshade('lat_gas', MetSPM, 2)
    end
    
    % med_gas
    subplot(7,3,20); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'med_gas');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Medial Gastrocnemius');
    %     ylabel('W / kg'); 
    ylim([ 0 yPk]);
         xlim([0 100]); 
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
    
       % soleus
    subplot(7,3,21); hold on;
    Ind = contains([Muscles(i).UMB_CombCols], 'soleus');
    plot(Stride, Muscles(i).UMB_CombMuscleAvg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Soleus');
    %     ylabel('W / kg'); 
    ylim([ 0 6]);
         xlim([0 100]); 
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
  
end

subplot(731);
Conditions = {'+40', '+20','Norm','-20','-40'}; 
legend([Pl(5).C, Pl(4).C, Pl(3).C, Pl(2).C, Pl(1).C], Conditions, 'Location','Best');
% Conditions = {'-40', '-20','Norm','+20','+40'}; 
% legend(Conditions, 'Location','Best');

% subplot(7,3,10);
% Conditions = {'-40', '-20','Norm','+20','+40'}; 
% legend(Conditions, 'Location','Best');

subplotsqueeze(MuscleMetCosts, 1.12);

% saveas(MuscleMetCosts, 'MuscleMetCosts.png'); 
% saveas(MuscleMetCosts, 'MuscleMetCosts.pdf'); 
print('MuscleMetCosts','-djpeg','-r300');


%% Check other Outputs

ResultsFolder = 'C:\Users\richa\Documents\Packages\OpenSim\Scripts\CMC_Results\ActuatorForces';
[Actuators, Subjects] = CheckActuators(ResultsFolder, Subjects, 'No', 'No', 'Yes');

ResultsFolder = 'C:\Users\richa\Documents\Packages\OpenSim\Scripts\CMC_Results\MuscleControls';
[Activations, Subjects] = CheckActivations(ResultsFolder, Subjects, 'No', 'Yes');

% kinematics
ResultsFolder = 'C:\Users\richa\Documents\Packages\OpenSim\Scripts\CMC_Results\Kinematics';
[Kinematics, Subjects] = CheckKinematics(ResultsFolder, Subjects, 'No','Yes');


