function [Kinematics, Subjects] = CheckKinematics(ResultsFolder, Subjects, PlotSubj, PlotCond)

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
% clc;
% IsStoFile = ~[KinematicsDir.isdir];
% Log = zeros(1, length(KinematicsDir));
% for i = 1:length(KinematicsDir)
%     Log(i) = IsStoFile(i) == 1 &&  contains(KinematicsDir(i).name, 'Kinematics');
% end
% 
% Files = find(Log);
% L = length(Files);
% 
% Kinematics(L).Name = [];
% Kinematics(L).Data = [];
% NumSubj = length(Subjects);
% % Subjects(NumSubj).Trials(:).Kinematics = [];
% 
% j = 1;
% for i = Files
%     
%     % get filename
%     Kinematics(j).Name = KinematicsDir(i).name;
%     
%     % load data
%     Data = importdata(KinematicsDir(i).name);
%     Data.filename = KinematicsDir(i).name;
%     
%     % get subject and trial
%     Str = strsplit(Kinematics(j).Name, '_');
%     Subj = Str{1};
%     SubjInd = contains({Subjects.name}, Subj);
%     Trial = Str{2};
%     TrialInd = contains({Subjects(SubjInd).Trials.name}, Trial);
%     Side = Str{3};
%     
%     TSData = Subjects(SubjInd).Trials(TrialInd).TSData; % identify temporal spatial data
%     Kinematics(j).Data = FilterAndParseData(Data, TSData); % filter and parse to gait cycles
%     
%     ActData = Kinematics(j).Data;
%     Fields2Del = {'data','textdata','Interp','Fdata','Filter','Fdata2','Filter2'};
%     ActData = rmfield(ActData, Fields2Del);
%     
%     % match Kinematics to subject and trial structure
%     if strcmp(Side, 'Left') % assing to subjects structure
%         Subjects(SubjInd).Trials(TrialInd).Kinematics.Left =  ActData;
%     elseif strcmp(Side, 'Right')
%         Subjects(SubjInd).Trials(TrialInd).Kinematics.Right =  ActData;
%     end
%     
%     j = j+1;
% end

NumSubj = length(Subjects);

[Kinematics, Subjects] = ExtractData(ResultsFolder, Subjects, 'kinematics');

clearvars i j IsStoFile fs Wn1 Cutoff1 b a L j Str Subj SubjInd Trial TrialInd Side...
    TSData Data KinematicsDir KinematicsFolder S ActData Log Dields2Del Files

%% Plot Kinematics
% PlotKinematics = 'Yes';
if strcmp(PlotSubj, 'Yes')
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

%% Combine left and right trials for the same subject into structure
clc;
for i = 1:length(Subjects) % loop through subjects
    
    for j = 1:length(Subjects(i).Trials) % loop through trials
        
        if strcmp(Subjects(i).Trials(j).type, 'static') == 0
            
            subj = contains({Kinematics.Subject}, Subjects(i).name);
            Str = strsplit(Subjects(i).Trials(j).name, '_');
            trial = contains({Kinematics.Trial}, Str{1});
            left =  contains({Kinematics.Side}, 'Left');
            right = contains({Kinematics.Side}, 'Right');
            
            SubjTrialLeft = subj + trial + left == 3;
            SubjTrialRight = subj + trial + right == 3;
            
            if sum(SubjTrialLeft) > 0
                Subjects(i).Trials(j).Kinematics.Left = Kinematics(SubjTrialLeft);
            end
            if sum(SubjTrialRight) > 0
                Subjects(i).Trials(j).Kinematics.Right = Kinematics(SubjTrialRight);
            end
            
            [Subjects(i).Trials(j).Kinematics.Left.Muscles, ~] = GetMuscles(Subjects(i).Trials(j).Kinematics.Left.Data.Parsed, ...
                Subjects(i).Trials(j).Kinematics.Left.Data.colheaders, 'notmet');
            
            [Subjects(i).Trials(j).Kinematics.Right.Muscles, ~] = GetMuscles(Subjects(i).Trials(j).Kinematics.Right.Data.Parsed, ...
                Subjects(i).Trials(j).Kinematics.Right.Data.colheaders, 'notmet');
            
            % average sides and combine together
            for k = 1:length(Subjects(i).Trials(j).Kinematics.Left.Muscles)
                Subjects(i).Trials(j).Kinematics.Muscles(k).name = Subjects(i).Trials(j).Kinematics.Left.Muscles(k).name;
                Subjects(i).Trials(j).Kinematics.Muscles(k).data = mean([...
                    Subjects(i).Trials(j).Kinematics.Left.Muscles(k).data.left_parsed, ...
                    Subjects(i).Trials(j).Kinematics.Right.Muscles(k).data.right_parsed],2);
                
                %             figure; hold on;
                %             plot(Subjects(i).Trials(j).Kinematics.Muscles(k).data, 'k', 'LineWidth', LW);
                %             plot(Subjects(i).Trials(j).Kinematics.Left.Muscles(k).data.left_parsed, 'b');
                %             plot(Subjects(i).Trials(j).Kinematics.Right.Muscles(k).data.right_parsed, 'r');
                
                Subjects(i).Trials(j).Kinematics.data(:,k) = Subjects(i).Trials(j).Kinematics.Muscles(k).data;
                Subjects(i).Trials(j).Kinematics.colheaders{k} = Subjects(i).Trials(j).Kinematics.Muscles(k).name;
            end
            
            % define new structure CondKinematics
            Conditions = {'Fm40','Fm20','NormF','Fp20','Fp40'};
            k = contains(Conditions, Str{1});
            CondKinematics(k).Condition = Conditions{k};
            CondKinematics(k).colheaders = Subjects(i).Trials(j).Kinematics.colheaders;
            CondKinematics(k).Data(:,:,i) = Subjects(i).Trials(j).Kinematics.data;
            
            
            clearvars SubjTrialLeft SubjTrialRight left right Str subj trial Array
        end
    end
end

clearvars MetabolicsDir MetabolicsFolder subjectPath TSData i j Nfiles Z ans

%% plot avg muscle forces by condition
if strcmp(PlotCond,'Yes')
clc; close all;
LW = 1.5; % set line width
% TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray');  rgb('DarkSlateGray'); rgb('Black')];
TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
FntSz = 9;

KinematicsFig = figure('Position',[50 50 1000 700]);
Stride = 1:100;
yPk = 1;
for i = [2 4 1 5 3]

CondKinematics(i).Avg = mean(CondKinematics(i).Data, 3); 
    
    % hip flexion
    subplot(2,3,1); hold on;
    Ind = contains([CondKinematics(i).colheaders], 'hip_flexion');
        p = sum( CondKinematics(i).Avg(:,Ind), 2) / 2;
    Pl(i).C = plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Hip Flexion');
       ylabel('Degrees'); 
    ylim([ -20 30]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTick = [0 0.5 1];
%     ax.YTickLabel = {'0', '50', '100'};
    ax.FontSize = FntSz;
    
       % hip adduciton
    subplot(2,3,2); hold on;
    Ind = contains([CondKinematics(i).colheaders], 'hip_adduction');
        p = sum( CondKinematics(i).Avg(:,Ind), 2) / 2;
   plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Hip Adduction');
%            ylabel('Degrees'); 
    ylim([ -10 10]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTick = [0 0.5 1];
%     ax.YTickLabel = {'0', '50', '100'};
    ax.FontSize = FntSz;
    
          % hip Rotation
    subplot(2,3,3); hold on;
    Ind = contains([CondKinematics(i).colheaders], 'hip_rotation');
        p = sum( CondKinematics(i).Avg(:,Ind), 2) / 2;
   plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Hip Rotation');
%             ylabel('Degrees'); 
    ylim([ -10 10]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTick = [0 0.5 1];
%     ax.YTickLabel = {'0', '50', '100'};
    ax.FontSize = FntSz;
    
    
          % knee flexion
    subplot(2,3,4); hold on;
    Ind = contains([CondKinematics(i).colheaders], 'knee_angle');
        p = sum( CondKinematics(i).Avg(:,Ind), 2) / 2;
   plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Knee Flexion');
          ylabel('Degrees'); 
    ylim([ -40 10]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTick = [0 0.5 1];
%     ax.YTickLabel = {'0', '50', '100'};
    ax.FontSize = FntSz;
    
    % ankle flexion
    subplot(2,3,5); hold on;
    Ind = contains([CondKinematics(i).colheaders], 'ankle_angle');
        p = sum( CondKinematics(i).Avg(:,Ind), 2) / 2;
   plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
    title('Ankle Flexion');
%         ylabel('Degrees'); 
    ylim([ -15 15]);
    ax = gca;
    ax.XTick = [0 25 50 75 100];
    ax.XTickLabel = [];
%     ax.YTick = [0 0.5 1];
%     ax.YTickLabel = {'0', '50', '100'};
    ax.FontSize = FntSz;
    
end

Conditions = {'+40', '+20','Norm','-20','-40'}; 
% legend([Pl(4).C, Pl(3).C, Pl(5).C, Pl(1).C, Pl(3).C], Conditions, 'Location','Best');
legend([Pl(5).C, Pl(4).C, Pl(3).C, Pl(2).C, Pl(1).C], Conditions, 'Location','Best', ...
    'FontSize',6.5);

subplotsqueeze(KinematicsFig, 1.12);
supertitle({'KINEMATICS'; ' '; ' '; ' '}, 'FontSize',16); 
saveas(KinematicsFig, 'KinematicsCond.png'); 
% saveas(MuscleForces, 'MuscleForces.pdf'); 
end

%% Compare to literature?

%% Export

end