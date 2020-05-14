function [Actuators, Subjects] = CheckActuators(ResultsFolder, Subjects, PlotGlobalFMs, PlotMuscleForces)

% Check Actuators forces and moments for muscles and joint/global residuals
% pass along whether these are within the threshold


% Ricky Pimentel
% Applied Biomechanics Lab


%% Settings and Inputs

if exist('ResultsFolder','var') == 0
    ResultsFolder = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\CMC_Results\ActuatorForces';
end

ActuatorsFolder = ResultsFolder;
ActuatorsDir = dir(ActuatorsFolder);
addpath(genpath(ActuatorsFolder));
addpath(genpath('CodeLibrary'));
addpath(genpath('BatchOpenSim'));

% add subject data
if exist('Subjects', 'var') == 0
    S = load('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\exp\Subjects.mat');
    Subjects = S.Subjects;
end

if exist('PlotGlobalFMs','var')
    PlotGlobalFMs = 'No';
end
if exist('PlotMuscleForces','var')
    PlotMuscleForces = 'No';
end

%% Set Actuators Tolerance Thresholds
% Threshold values from Hicks et al. 2015. Is My Model Good Enough? Best Practices for Verification and Validation of Musculoskeletal Models and Simulations of Movement
% DOI 10.1115/1.4029304
Thresh.Forces = 25;
Thresh.Moments = 75;

Thresh.AvgForces = 50;
Thresh.AvgMoments = 150;

%% Extract Actuators data
clc;
IsStoFile = ~[ActuatorsDir.isdir];
Log = zeros(1, length(ActuatorsDir));
for i = 1:length(ActuatorsDir)
    Log(i) = IsStoFile(i) == 1 &&  contains(ActuatorsDir(i).name, 'Actuation_force');
end

Files = find(Log);
L = length(Files);

Actuators(L).Name = [];
Actuators(L).Data = [];
NumSubj = length(Subjects);
for i = 1:NumSubj
    Subjects(i).Trials(6).Actuators = [];
end

j = 1;
for i = Files
    
    % get filename
    Actuators(j).Name = ActuatorsDir(i).name;
    
    % load data
    Data = importdata(ActuatorsDir(i).name);
    Data.filename = ActuatorsDir(i).name;
    
    % get subject and trial
    Str = strsplit(Actuators(j).Name, '_');
    Subj = Str{1};
    SubjInd = contains({Subjects.name}, Subj);
    Trial = Str{2};
    TrialInd = contains({Subjects(SubjInd).Trials.name}, Trial);
    Side = Str{3};
    
    TSData = Subjects(SubjInd).Trials(TrialInd).TSData; % identify temporal spatial data
    Actuators(j).Data = FilterAndParseData(Data, TSData); % filter and parse to gait cycles
    
    ActData = Actuators(j).Data;
    Fields2Del = {'data','textdata','Interp','Fdata','Filter','Fdata2','Filter2'};
    ActData = rmfield(ActData, Fields2Del);
    
    % match Actuators to subject and trial structure
    if strcmp(Side, 'Left') % assing to subjects structure
        Subjects(SubjInd).Trials(TrialInd).Actuators.Left =  ActData;
    elseif strcmp(Side, 'Right')
        Subjects(SubjInd).Trials(TrialInd).Actuators.Right =  ActData;
    end
    
    j = j+1;
end

clearvars i j IsStoFile fs Wn1 Cutoff1 b a L j Str Subj SubjInd Trial TrialInd Side...
    TSData Data ResultsFolder ActuatorsDir ActuatorsFolder S ActData Log Dields2Del Files

%% Check Overall Forces and Moments
COLS = {'FX','FY','FZ','MX','MY','MZ'};

for i = 1:length(Actuators)
    cols = contains({Actuators(i).Data.colheaders{:}}, COLS);
    Column = find(cols);
    
    
    % Check Actuators forces
    for j = 1:3
        Actuators(i).MaxForces(j) = max(abs(Actuators(i).Data.Parsed(:,Column(j))));
        Actuators(i).MeanForces(j) = mean(abs(Actuators(i).Data.Parsed(:,Column(j))));
    end
    if max(Actuators(i).MaxForces) > Thresh.Forces
        Actuators(i).Force_QC = 'fail';
    else
        Actuators(i).Force_QC = 'pass';
    end
    
    
    
    % Check Actuators moments
    for j = 4:6
        Actuators(i).MaxMoments(j-3) = max(abs(Actuators(i).Data.Parsed(:,Column(j))));
        Actuators(i).MeanMoments(j-3) = mean(abs(Actuators(i).Data.Parsed(:,Column(j))));
    end
    if max(Actuators(i).MaxMoments) > Thresh.Moments
        Actuators(i).Moment_QC = 'fail';
    else
        Actuators(i).Moment_QC = 'pass';
    end
    
    % determine
    if strcmp(Actuators(i).Force_QC, 'fail') || strcmp(Actuators(i).Moment_QC, 'fail')
        Actuators(i).QC_Fail = 1;
    else
        Actuators(i).QC_Fail = 0;
    end
    
    if max(Actuators(i).MeanForces) > Thresh.AvgForces || max(Actuators(i).MeanMoments) > Thresh.AvgMoments
        Actuators(i).Exclude = 1;
    elseif max(Actuators(i).MaxForces) > 75 || max(Actuators(i).MaxMoments) > 150
        Actuators(i).Exclude = 1;
    else
        Actuators(i).Exclude = 0;
    end
end

%% Plot global forces and moments
PlotGlobalFMs = 'Yes';
if strcmp(PlotGlobalFMs, 'Yes')
    SubjColors = colormap(jet(NumSubj));
    close all;
    clc;
    
    GlobalForcesFig = figure('Position',[10 10 1200 900]);
    
    COLS = {'FX','FY','FZ','MX','MY','MZ'};
    
    % loop through subjects
    for subj = 1:NumSubj
        for trial = 1:length(Subjects(subj).Trials)
            if strcmp(Subjects(subj).Trials(trial).type, 'static')
                continue
            end
            % find columns
            cols = contains([Subjects(subj).Trials(trial).Actuators.Left.colheaders], COLS);
            ColumnL = find(cols);
            cols = contains([Subjects(subj).Trials(trial).Actuators.Right.colheaders], COLS);
            ColumnR = find(cols);
            
            % plot Actuators data
            for j = 1:length(COLS)
                subplot(6, 1, j); hold on;
                ActL = Subjects(subj).Trials(trial).Actuators.Left.Parsed(:,ColumnL(j));
                ActR = Subjects(subj).Trials(trial).Actuators.Right.Parsed(:,ColumnR(j));
                plot(ActL, '-',  'Color', SubjColors(subj, :));
                plot(ActR, '--',  'Color', SubjColors(subj, :));
                ax = gca;
                ax.XTick = [0 20 40 60 80 100];
                ax.YTick = [];
            end
        end
        
        text(100, subj*120, strcat('\bf Subject', Subjects(subj).name(end-1:end)),  'Color', SubjColors(subj, :));
    end
    
    % add on plot titles, thresholds, and other details
    LW = 2;
    for j = 1:length(COLS)
        subplot(6,1,j); hold on;
        ylabel(strcat(COLS{j}, ' (N)'));
        if j < 4
            T = Thresh.Forces;
        else
            T = Thresh.Moments;
        end
        h = hline(T, '--k');
        h.LineWidth = LW;
        h = hline(-T, '--k');
        h.LineWidth = LW;
        xlim([1 100]);
        ax = gca;
        ax.FontSize = 10;
        ax.XTickLabel = {'0','20','40','60','80','100'}; 
        
        if j < 6
            ax.XTickLabels = [];
        end
    end
    
    for j = 1:3
          subplot(6,1,j); hold on;
          ax = gca; 
          ax.YTick = [-Thresh.Forces Thresh.Forces];
    end
    
    for j = 4:6
          subplot(6,1,j); hold on;
          ax = gca; 
          ax.YTick = [-Thresh.Moments Thresh.Moments];
    end
    
    subplot(6,1,6);
    xlabel('% of Gait Cycle');
    text(100, 50, '\bf - - Right');
    text(100, 0, '\bf - Left');
    
    subplotsqueeze(GlobalForcesFig, 1.1);
    saveas(GlobalForcesFig, 'GlobalForces.png');
end

%% Plot joint reserve moments
PlotGlobalFMs = 'Yes';
if strcmp(PlotGlobalFMs, 'Yes')
    SubjColors = colormap(jet(NumSubj));
    close all;
    clc;
    
    JtReservesFig = figure('Position',[10 10 1200 900]);
    
    COLS = {'lumbar_extension_reserve','lumbar_bending_reserve','lumbar_rotation_reserve', ...
        'hip_flexion_l_reserve',  'hip_flexion_r_reserve',...
        'hip_adduction_l_reserve','hip_adduction_r_reserve',  ...
        'hip_rotation_l_reserve','hip_rotation_r_reserve',  ...
        'knee_angle_l_reserve','knee_angle_r_reserve',  ...
        'ankle_angle_l_reserve', 'ankle_angle_r_reserve'};
    
    % loop through subjects
    for subj = 1:NumSubj
        for trial = 1:length(Subjects(subj).Trials)
            if strcmp(Subjects(subj).Trials(trial).type, 'static')
                continue
            end
            % find columns
            cols = contains([Subjects(subj).Trials(trial).Actuators.Left.colheaders], COLS);
            ColumnL = find(cols);
            cols = contains([Subjects(subj).Trials(trial).Actuators.Right.colheaders], COLS);
            ColumnR = find(cols);
            
            % plot Actuators data
            for j = 1:length(COLS)
                subplot(7, 2, j+1); hold on;
                ActL = Subjects(subj).Trials(trial).Actuators.Left.Parsed(:,ColumnL(j));
                ActR = Subjects(subj).Trials(trial).Actuators.Right.Parsed(:,ColumnR(j));
                plot(ActL, '-',  'Color', SubjColors(subj, :));
                plot(ActR, '--',  'Color', SubjColors(subj, :));
                ax = gca;
                ax.XTick = [0 20 40 60 80 100];
%                 ax.YTick = [];
            end
        end
        
        text(100, subj*120, strcat('\bf ', Subjects(subj).name),  'Color', SubjColors(subj, :));
    end
    
    % add on plot titles, thresholds, and other details
    LW = 2;
    for j = 1:length(COLS)
        subplot(7,2,j+1); hold on;
         Str = strrep(COLS{j}, '_reserve',''); 
         if strcmp(Str(end-1:end), '_r')
             Str = strcat(Str(1:end-2), ' Right'); 
         elseif strcmp(Str(end-1:end), '_l')
             Str = strcat(Str(1:end-2), ' Left');
         end

         Str = strrep(Str, 'lumbar', 'Lumbar');
         Str = strrep(Str, 'hip', 'Hip');
         Str = strrep(Str, 'knee', 'Knee');
         Str = strrep(Str, 'ankle', 'Ankle');
         Str = strrep(Str, 'extension', 'Extension');
         Str = strrep(Str, 'bending', 'Bending');
         Str = strrep(Str, 'rotation', 'Rotation');
         Str = strrep(Str, 'adduction', 'Adduction');
         Str = strrep(Str, 'flexion', 'Flexion');
         Str = strrep(Str, 'angle', 'Flexion');
         
        title(strrep(Str,'_',' ')); 
        ylabel('Nm');
        
        ax = gca;
        ax.FontSize = 8;
        ax.XTickLabels = [];
    end

    subplot(721);
    for subj = 1:length(SubjColors)
        t = text(subj/10, 0.1, strcat('\bf Subject ', Subjects(subj).name(3:4)),  'Color', SubjColors(subj, :));
        t.Rotation = 90;
    end
    text(0.866, 2.3, '\bf - - Right', 'HorizontalAlignment','center');
    text(0.433, 2.3, '\bf - Left',  'HorizontalAlignment','center');
    xlim([0 1.3]);
    ylim([0 2.5]);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    


for j = 13:14
    subplot(7,2,j);
    xlabel('% of Gait Cycle');
    ax = gca; 
    ax.XTick = [0 20 40 60 80 100]; 
    ax.XTickLabel = {'0','20','40','60','80','100'};
end


    
    subplotsqueeze(JtReservesFig, 1.1);
    saveas(JtReservesFig, 'JointReserves.png');
end

%% Plot Muscle Actuators by subject
PlotMuscleForces = 'Yes';
if strcmp(PlotMuscleForces, 'Yes')
    
    SubjColors = colormap(jet(NumSubj));
    close all;
    clc;
    
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
            
            % plot
            for j = 1:length(COLS)
                subplot(9, 3, j); hold on;
                Ind = contains([Subjects(subj).Trials(trial).Actuators.Left.colheaders], COL_left{j}); % find columns
                ActL = Subjects(subj).Trials(trial).Actuators.Left.Parsed(:,Ind) / Subjects(subj).Demo.mass; % normalize to body mass
                Ind = contains([Subjects(subj).Trials(trial).Actuators.Right.colheaders], COL_right{j}); % find columns
                ActR = Subjects(subj).Trials(trial).Actuators.Right.Parsed(:,Ind) / Subjects(subj).Demo.mass; % normalize to body mass
                plot(ActL, '-',  'Color', SubjColors(subj, :));
                plot(ActR, '--',  'Color', SubjColors(subj, :));
                title(strrep(COLS{j}, '_',' '));
                ax = gca;
                ax.XTick = [];
                %                 ax.YTick = [];
                %                 ylim([0 1]);
            end
        end
        
        text(105, subj, strcat('\bf ', Subjects(subj).name),  'Color', SubjColors(subj, :));
    end
    
    % label axes
    for j = [1 4 7 10 13 16 19 22 25]
        subplot(9, 3, j);
        %         ax = gca;
        %         ax.YTick = [0 0.5 1];
        %         ax.YTickLabel = {'0', '50', '100'};
        ylabel('N / kg');
    end
    for j = [25 26 27]
        subplot(9, 3, j);
        ax = gca;
        ax.XTick = [0 20 40 60 80 100];
        ax.XTickLabel = {'0', '20', '40', '60', '80','100'};
        xlabel('% Gait Cycle');
    end
    text(100, 0.5, '\bf - - Right');
    text(100, 0, '\bf - Left');
    
    % save figure
    subplotsqueeze(H, 1.15);
    saveas(H, 'MuscleForces.png');
end

%     %% Average forces bilaterally
%
%
%     Conditions = {'Fm40','Fm20','Norm','Fp20','Fp40'};
%     Cond = length(Conditions);
%     ActuatorCond(Cond).Name = [];
%     for i = 1:length(Conditions)
%         ActuatorCond(i).Name = Conditions{i};
%     end
%
% %     j = 1;
%     % loop through subjects to average across conditions
%     for subj = 1:length(Subjects)
%         for trial = 1:Cond
%
%             % pull left and right columns from data
%             ActuatorCond(trial).Left.ColHeaders = Subjects(subj).Trials(trial).Actuators.Left.colheaders;
%             ActuatorCond(trial).Right.ColHeaders = Subjects(subj).Trials(trial).Actuators.Right.colheaders;
%
%
%             % combine muscles with multiple lines of action
%             ActuatorCond(trial).Left.ColHeaders(:,end+1:end+8) =  ...
%                 {'glut_max_r', 'glut_med_r', 'glut_min_r', 'add_mag_r', 'glut_max_l', 'glut_med_l', 'glut_min_l', 'add_mag_l'};
%             ActuatorCond(trial).Right.ColHeaders(:,end+1:end+8) =  ...
%                 {'glut_max_r', 'glut_med_r', 'glut_min_r', 'add_mag_r', 'glut_max_l', 'glut_med_l', 'glut_min_l', 'add_mag_l'};
%
%             str = 'glut_max';
%             lcol = 117; rcol = 113;
%             Ind = contains(ActuatorCond(trial).Left.ColHeaders, {strcat(str, '1_l'), strcat(str, '2_l'), strcat(str, '3_l')});
%             ActuatorCond(trial).Left.Parsed(:,lcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Left.Parsed(:,Ind), 2);
%             Ind = contains(ActuatorCond(trial).Right.ColHeaders, {strcat(str, '1_r'), strcat(str, '2_r'), strcat(str, '3_r')});
%             ActuatorCond(trial).Right.Parsed(:,rcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Right.Parsed(:,Ind), 2);
%
%             str = 'glut_med';
%             lcol = 118; rcol = 114;
%             Ind = contains(ActuatorCond(trial).Left.ColHeaders, {strcat(str, '1_l'), strcat(str, '2_l'), strcat(str, '3_l')});
%             ActuatorCond(trial).Left.Parsed(:,lcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Left.Parsed(:,Ind), 2);
%             Ind = contains(ActuatorCond(trial).Left.ColHeaders, {strcat(str, '1_r'), strcat(str, '2_r'), strcat(str, '3_r')});
%             ActuatorCond(trial).Right.Parsed(:,rcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Right.Parsed(:,Ind), 2);
%
%             str = 'glut_min';
%             lcol = 119; rcol = 115;
%             Ind = contains(ActuatorCond(trial).Left.ColHeaders, {strcat(str, '1_l'), strcat(str, '2_l'), strcat(str, '3_l')});
%             ActuatorCond(trial).Left.Parsed(:,lcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Left.Parsed(:,Ind), 2);
%             Ind = contains(ActuatorCond(trial).Right.ColHeaders, {strcat(str, '1_r'), strcat(str, '2_r'), strcat(str, '3_r')});
%             ActuatorCond(trial).Right.Parsed(:,rcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Right.Parsed(:,Ind), 2);
%
%             str = 'add_mag';
%             lcol = 120; rcol = 116;
%             Ind = contains(ActuatorCond(trial).Left.ColHeaders, {strcat(str, '1_l'), strcat(str, '2_l'), strcat(str, '3_l')});
%             ActuatorCond(trial).Left.Parsed(:,lcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Left.Parsed(:,Ind), 2);
%             Ind = contains(ActuatorCond(trial).Right.ColHeaders, {strcat(str, '1_r'), strcat(str, '2_r'), strcat(str, '3_r')});
%             ActuatorCond(trial).Right.Parsed(:,rcol, subj) = sum(Subjects(subj).Trials(trial).Actuators.Right.Parsed(:,Ind), 2);
%
%              % copy over actuations
%             ActuatorCond(trial).Left.Parsed(:,1:112, subj) = Subjects(subj).Trials(trial).Actuators.Left.Parsed;
%             ActuatorCond(trial).Right.Parsed(:,1:112, subj) = Subjects(subj).Trials(trial).Actuators.Right.Parsed;
%
%             Left_Ind = contains(ActuatorCond(trial).Left.ColHeaders, '_l');
%             Left_Ind([13 32 41]) = 0; % remove muscles with _lat or _long
%             Left_Ind(110:112) = 1; % add L reserves
%             Right_Ind = contains(ActuatorCond(trial).Right.ColHeaders, '_r');
%             Right_Ind(105:109) = 0; % get rid of L reserves
%             sum(Right_Ind);
%             sum(Left_Ind);
%
%             a(:,:,1) = ActuatorCond(trial).Left.Parsed(:,Left_Ind, subj);
%             a(:,:,2) = ActuatorCond(trial).Right.Parsed(:,Right_Ind, subj);
%             ActuatorCond(trial).Data(:,:,subj) = mean(a, 3);
%
%             ActuatorCond(trial).ColHeaders = strrep(ActuatorCond(trial).Left.ColHeaders(Left_Ind), '_l', '');
%             ActuatorCond(trial).ColHeaders{12} = 'add_long';
%             ActuatorCond(trial).ColHeaders{31} = 'vas_lat';
%             ActuatorCond(trial).ColHeaders{40} = 'per_long';
%             %             ActuatorCond(trial).ColHeaders(13) = 'add_long';
%
%             % plot to verify averaging
% %             if subj == 1
% %                 figure; hold on;
% %                 Ind = contains(ActuatorCond(trial).Left.ColHeaders, 'soleus_l');
% %                 plot(ActuatorCond(trial).Left.Parsed(:,Ind,subj), '--r');
% %                 Ind = contains(ActuatorCond(trial).Left.ColHeaders, 'soleus_r');
% %                 plot(ActuatorCond(trial).Right.Parsed(:,Ind,subj), '--b');
% %                  Ind = contains(ActuatorCond(trial).ColHeaders, 'soleus');
% %                 plot(ActuatorCond(trial).Data(:,Ind, subj), '-k');
% %             end
%
%         end
%     end



%% Plot muscle forces by condition
MuscleActuatorsCond = figure('Position', [100 100 800 800]);
Colors = [rgb('PowderBlue'); rgb('DeepSkyBlue'); rgb('RoyalBlue');  rgb('Navy'); rgb('MidnightBlue')];
ColorMatch = [2 1 4 5 3];
LW = 1;


j =1;
for i = 1:length(Actuators)
    % define subject and trial from Actuators filename
    Str = strsplit(Actuators(i).Name, '_');
    Subj = Str{1};
    SubjInd = contains({Subjects.name}, Subj);
    Trial = Str{2};
    TrialInd = contains({Subjects(SubjInd).Trials.name}, Trial);
    
    cMatch = ColorMatch(TrialInd);
    
    % get subject mass from data structure
    Mass = Subjects(SubjInd).Demo.mass;
    
    if contains(Actuators(i).Name, 'Left')
        
        % ankle
        subplot(411); hold on;
        str = strcat( 'med_gas', '_l');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        H(cMatch) = plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        LegendName{cMatch} = Trial;
        
        subplot(412); hold on;
        str = strcat( 'lat_gas', '_l');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        
        subplot(413); hold on;
        str = strcat( 'soleus', '_l');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        
        subplot(414); hold on;
        str = strcat( 'tib_ant', '_l');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        
    elseif contains(Actuators(i).Name, 'Right')
        subplot(411); hold on;
        str = strcat( 'med_gas', '_r');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        R(cMatch) = plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        LegendNameR{cMatch} = Trial;
        
        subplot(412); hold on;
        str = strcat( 'lat_gas', '_r');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        
        subplot(413); hold on;
        str = strcat( 'soleus', '_r');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        
        subplot(414); hold on;
        str = strcat( 'tib_ant', '_r');
        Column = contains({Actuators(i).Data.colheaders{:}}, str);
        plot(Actuators(i).Data.Parsed(:,Column) / Mass, '-', 'Color', Colors(cMatch, :), ...
            'LineWidth', LW);
        
    end
    
end

subplot(411);
title('Med Gas Force');
ylabel('N / kg');
ylim([0 40]);

subplot(412);
title('Lat Gas Force');
ylabel('N / kg');
ylim([0 40]);

subplot(413);
title('Soleus Force');
ylabel('N / kg');
ylim([0 40]);

subplot(414);
title('Tib Ant Force');
ylabel('N / kg');
xlabel('% Gait Cycle');
ylim([0 40]);


if exist('LegendName', 'var') == 0
    LegendName = LegendNameR;
    H = R;
end

if length(Actuators) == 1
    LegendName = {LegendName{cMatch}};
    H = H(cMatch);
end

legend(H, LegendName);

subplotsqueeze(MuscleActuatorsCond, 1.1);

end

%%


% end