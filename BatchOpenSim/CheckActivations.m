function [Activations, Subjects] = CheckActivations(ResultsFolder)

% Ensure Muscle Activations are within the threshold


%% Identify files
if exist('ResultsFolder','var') == 0
    ResultsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\MuscleControls';
end

ActivationsFolder = ResultsFolder;
ActivationsDir = dir(ActivationsFolder);
addpath(genpath(ActivationsFolder));
addpath(genpath('CodeLibrary'));
addpath(genpath('BatchOpenSim'));

% add subject data
if exist('Subjects', 'var') == 0
    S = load('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\exp\Subjects.mat');
    Subjects = S.Subjects;
end

%% Extract Activation data
clc;
IsStoFile = ~[ActivationsDir.isdir];
Log = zeros(1, length(ActivationsDir));
for i = 1:length(ActivationsDir)
    Log(i) = IsStoFile(i) == 1 &&  contains(ActivationsDir(i).name, 'controls');
end

Files = find(Log);
L = length(Files);

Activations(L).Name = [];
Activations(L).Data = [];
NumSubj = length(Subjects);
% Subjects(NumSubj).Trials(:).Activations = [];

j = 1;
for i = Files
    
    % get filename
    Activations(j).Name = ActivationsDir(i).name;
    
    % load data
    Data = importdata(ActivationsDir(i).name);
    Data.filename = ActivationsDir(i).name;
    
    % get subject and trial
    Str = strsplit(Activations(j).Name, '_');
    Subj = Str{1};
    SubjInd = contains({Subjects.name}, Subj);
    Trial = Str{2};
    TrialInd = contains({Subjects(SubjInd).Trials.name}, Trial);
    Side = Str{3};
    
    TSData = Subjects(SubjInd).Trials(TrialInd).TSData; % identify temporal spatial data
    Activations(j).Data = FilterAndParseData(Data, TSData); % filter and parse to gait cycles
    
    ActData = Activations(j).Data;
    Fields2Del = {'data','textdata','Interp','Fdata','Filter','Fdata2','Filter2'};
    ActData = rmfield(ActData, Fields2Del);
    
    % match activations to subject and trial structure
    if strcmp(Side, 'Left') % assing to subjects structure
        Subjects(SubjInd).Trials(TrialInd).Activations.Left =  ActData;
    elseif strcmp(Side, 'Right')
        Subjects(SubjInd).Trials(TrialInd).Activations.Right =  ActData;
    end
    
    j = j+1;
end

clearvars i j IsStoFile fs Wn1 Cutoff1 b a L j Str Subj SubjInd Trial TrialInd Side...
    TSData Data ResultsFolder ActivationsDir ActivationsFolder S ActData Log Dields2Del Files

%% Plot muscle activations
PlotActivations = 'Yes';
if strcmp(PlotActivations, 'Yes')
    % LW = 2;
    SubjColors = colormap(jet(NumSubj));
    close all;
    clc
    
    H = figure('Position',[10 10 1500 900]);
    
    % define columns and rename for left/right
    COLS = {'glut_med1','glut_med2','glut_med3','glut_max1','glut_max2','glut_max3', 'glut_min1','glut_min2','glut_min3',...
        'psoas','iliacus','tfl','semimem','semiten','bifemlh','bifemsh','grac','rect_fem','vas_med','vas_int','vas_lat',...
        'soleus','lat_gas','med_gas','tib_ant','tib_post','flex_dig'};
    COL_left = COLS;
    COL_right = COLS;
    for j = 1:length(COLS)
        COL_left{j}(end+1:end+2) = '_l';
        COL_right{j}(end+1:end+2) = '_r';
    end
    
    % loop through subjects
    for subj = 1:NumSubj
        for trial = 1:length(Subjects(subj).Trials)
            if strcmp(Subjects(subj).Trials(trial).type, 'static')
                continue
            end
            
            % plot activation data
            for j = 1:length(COLS)
                subplot(9, 3, j); hold on;
                Ind = contains([Subjects(subj).Trials(trial).Activations.Left.colheaders], COL_left{j}); % find columns
                plot(Subjects(subj).Trials(trial).Activations.Left.Parsed(:,Ind), '-', 'Color', SubjColors(subj, :));
                 Ind = contains([Subjects(subj).Trials(trial).Activations.Right.colheaders], COL_right{j}); % find columns
                plot(Subjects(subj).Trials(trial).Activations.Right.Parsed(:,Ind), '--',  'Color', SubjColors(subj, :));
                title(strrep(COLS{j}, '_',' '));
                ax = gca;
                ax.XTick = [];
                ax.YTick = [];
                ylim([0 1]);
            end
        end
        
        text(105, subj*0.8, strcat('\bf ', Subjects(subj).name),  'Color', SubjColors(subj, :));
    end
    
    % label axes
    for j = [1 4 7 10 13 16 19 22 25]
        subplot(9, 3, j);
        ax = gca;
        ax.YTick = [0 0.5 1];
        ax.YTickLabel = {'0', '50', '100'};
        ylabel('Activation (%)');
    end
    for j = [25 26 27]
        subplot(9, 3, j);
        ax = gca;
        ax.XTick = [0 20 40 60 80 100];
        ax.XTickLabel = {'0', '20', '40', '60', '80','100'};
        xlabel('% Gait Cycle');
    end
    
    % save figure
    subplotsqueeze(H, 1.15);
    saveas(H, 'MuscleActivations.png');
    
end

%% Compare to literature?

%% Export

end