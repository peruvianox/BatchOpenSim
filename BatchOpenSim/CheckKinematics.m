function [Kinematics, Subjects] = CheckKinematics(ResultsFolder)

% Ensure Muscle Kinematics are within the threshold


%% Identify files
if exist('ResultsFolder','var') == 0
    ResultsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\Kinematics';
end

KinematicsFolder = ResultsFolder;
KinematicsDir = dir(KinematicsFolder);
addpath(genpath(KinematicsFolder));
addpath(genpath('CodeLibrary'));
addpath(genpath('BatchOpenSim'));

% add subject data
if exist('Subjects', 'var') == 0
    S = load('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\exp\Subjects.mat');
    Subjects = S.Subjects;
end

%% Extract Kinematics
clc;
IsStoFile = ~[KinematicsDir.isdir];
Log = zeros(1, length(KinematicsDir));
for i = 1:length(KinematicsDir)
    Log(i) = IsStoFile(i) == 1 &&  contains(KinematicsDir(i).name, 'Kinematics');
end

Files = find(Log);
L = length(Files);

Kinematics(L).Name = [];
Kinematics(L).Data = [];
NumSubj = length(Subjects);
% Subjects(NumSubj).Trials(:).Kinematics = [];

j = 1;
for i = Files
    
    % get filename
    Kinematics(j).Name = KinematicsDir(i).name;
    
    % load data
    Data = importdata(KinematicsDir(i).name);
    Data.filename = KinematicsDir(i).name;
    
    % get subject and trial
    Str = strsplit(Kinematics(j).Name, '_');
    Subj = Str{1};
    SubjInd = contains({Subjects.name}, Subj);
    Trial = Str{2};
    TrialInd = contains({Subjects(SubjInd).Trials.name}, Trial);
    Side = Str{3};
    
    TSData = Subjects(SubjInd).Trials(TrialInd).TSData; % identify temporal spatial data
    Kinematics(j).Data = FilterAndParseData(Data, TSData); % filter and parse to gait cycles
    
    ActData = Kinematics(j).Data;
    Fields2Del = {'data','textdata','Interp','Fdata','Filter','Fdata2','Filter2'};
    ActData = rmfield(ActData, Fields2Del);
    
    % match Kinematics to subject and trial structure
    if strcmp(Side, 'Left') % assing to subjects structure
        Subjects(SubjInd).Trials(TrialInd).Kinematics.Left =  ActData;
    elseif strcmp(Side, 'Right')
        Subjects(SubjInd).Trials(TrialInd).Kinematics.Right =  ActData;
    end
    
    j = j+1;
end

clearvars i j IsStoFile fs Wn1 Cutoff1 b a L j Str Subj SubjInd Trial TrialInd Side...
    TSData Data ResultsFolder KinematicsDir KinematicsFolder S ActData Log Dields2Del Files

%% Plot  Kinematics
PlotKinematics = 'Yes';
if strcmp(PlotKinematics, 'Yes')
    % LW = 2;
    SubjColors = colormap(jet(NumSubj));
    close all;
    clc
    
    H = figure('Position',[10 10 1500 900]);
    
    % define columns for left/right
        COL_left = { 'lumbar_extension', 'lumbar_bending', 'lumbar_rotation', 'pelvis_tilt', 'pelvis_list', 'pelvis_rotation',...
            'hip_flexion_l', 'hip_adduction_l', 'hip_rotation_l', 'knee_angle_l', 'ankle_angle_l'};
        COL_right = {'lumbar_extension', 'lumbar_bending', 'lumbar_rotation', 'pelvis_tilt', 'pelvis_list', 'pelvis_rotation',...
            'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r', 'knee_angle_r', 'ankle_angle_r'};

    
    % loop through subjects
    for subj = 1:NumSubj
        for trial = 1:length(Subjects(subj).Trials)
            if strcmp(Subjects(subj).Trials(trial).type, 'static')
                continue
            end
            
            % plot activation data
            for j = 1:length(COL_left)
                if j == 2 || j == 3 || j == 5 || j == 6
                    subplot(4, 3, j); hold on;
                    Ind = contains([Subjects(subj).Trials(trial).Kinematics.Left.colheaders], COL_left{j}); % find columns
                    plot(Subjects(subj).Trials(trial).Kinematics.Left.Parsed(:,Ind), '-', 'Color', SubjColors(subj, :));
                    Ind = contains([Subjects(subj).Trials(trial).Kinematics.Right.colheaders], COL_right{j}); % find columns
                    plot(-Subjects(subj).Trials(trial).Kinematics.Right.Parsed(:,Ind), '--',  'Color', SubjColors(subj, :));
                    Str = strrep(COL_left{j}, '_l',' ');
                    title(strrep(Str, '_',' '));
                    ylim([-30 30]);
                    
                else
                    subplot(4, 3, j); hold on;
                    Ind = contains([Subjects(subj).Trials(trial).Kinematics.Left.colheaders], COL_left{j}); % find columns
                    plot(Subjects(subj).Trials(trial).Kinematics.Left.Parsed(:,Ind), '-', 'Color', SubjColors(subj, :));
                    Ind = contains([Subjects(subj).Trials(trial).Kinematics.Right.colheaders], COL_right{j}); % find columns
                    plot(Subjects(subj).Trials(trial).Kinematics.Right.Parsed(:,Ind), '--',  'Color', SubjColors(subj, :));
                    Str = strrep(COL_left{j}, '_l',' ');
                    title(strrep(Str, '_',' '));
                    ylim([-30 30]);
                end
            end
        end
        
        text(150, -70 + subj*8, strcat('\bf ', Subjects(subj).name),  'Color', SubjColors(subj, :));
    end
    
    subplot(4,3,7);
    ylim([-30 60]);
    
    subplot(4,3,10);
    ylim([-80 20]);
    
    subplot(4,3,11);
    ylim([-60 30]);
    
     text(170, -20, strcat('\bf ', '-Left    --Right')); 
    
    
    % label axes
    for j = 1:8
        subplot(4, 3, j);
        ax = gca;
        ax.XTick = [0 20 40 60 80 100];
        ax.XTickLabel = [];
    end
    for j = 9:11
             subplot(4, 3, j);
        xlabel('% Gait Cycle');
    end
    
    for j = [1 4 7 10]
        subplot(4, 3, j);
%         ax = gca;
%         ax.YTick = [0 0.5 1];
%         ax.YTickLabel = {'0', '50', '100'};
        ylabel('Degrees');
    end

%     
    % save figure
    subplotsqueeze(H, 1.15);
    saveas(H, 'Kinematics.png');
    
end



%% Compare to literature?

%% Export

end