function [Actuators, Subjects] = CheckActuators(ResultsFolder, Subjects, PlotGlobalFMs, PlotMuscleForces, PlotCond)

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
if exist('PlotCond','var')
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
NumSubj = length(Subjects);

[Actuators, Subjects] = ExtractData(ResultsFolder, Subjects, 'forces');

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
        %         ylabel(strcat(COLS{j}, ' (N)'));
        ylabel(' N');
        title(COLS{j});
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
    
    
    subplot(6,1,1);
    title('Anterior-Posterior Forces'); 
    
    subplot(6,1,2);
    title('Vertical Forces'); 
    
    subplot(6,1,3);
    title('Medial-Lateral Forces'); 
    
    subplot(6,1,4);
    title('Frontal Plane Moments'); 
    
    subplot(6,1,5);
    title('Transverse Plane Moments'); 
    
    subplot(6,1,6);
    title('Sagittal Plane Moments'); 
    xlabel('% of Gait Cycle');
    text(100, 50, '\bf - - Right');
    text(100, 0, '\bf - Left');
    
    subplotsqueeze(GlobalForcesFig, 1.1);
    saveas(GlobalForcesFig, 'GlobalForces.png');
end

%% Plot joint reserve moments
% PlotGlobalFMs = 'Yes';
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
%     saveas(JtReservesFig, 'JointReserves.png');
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
%     saveas(H, 'MuscleForces.png');
end

%% Combine left and right trials for the same subject into structure
clc;
for i = 1:length(Subjects) % loop through subjects
    
    for j = 1:length(Subjects(i).Trials) % loop through trials
        
        if strcmp(Subjects(i).Trials(j).type, 'static') == 0
            
            subj = contains({Actuators.Subject}, Subjects(i).name);
            Str = strsplit(Subjects(i).Trials(j).name, '_');
            trial = contains({Actuators.Trial}, Str{1});
            left =  contains({Actuators.Side}, 'Left');
            right = contains({Actuators.Side}, 'Right');
            
            SubjTrialLeft = subj + trial + left == 3;
            SubjTrialRight = subj + trial + right == 3;
            
            if sum(SubjTrialLeft) > 0
                Subjects(i).Trials(j).Actuators.Left = Actuators(SubjTrialLeft);
            end
            if sum(SubjTrialRight) > 0
                Subjects(i).Trials(j).Actuators.Right = Actuators(SubjTrialRight);
            end
            
            [Subjects(i).Trials(j).Actuators.Left.Muscles, ~] = GetMuscles(Subjects(i).Trials(j).Actuators.Left.Data.Parsed, ...
                Subjects(i).Trials(j).Actuators.Left.Data.colheaders, 'notmet');
            
            [Subjects(i).Trials(j).Actuators.Right.Muscles, ~] = GetMuscles(Subjects(i).Trials(j).Actuators.Right.Data.Parsed, ...
                Subjects(i).Trials(j).Actuators.Right.Data.colheaders, 'notmet');
            
            % average sides and combine together
            for k = 1:length(Subjects(i).Trials(j).Actuators.Left.Muscles)
                Subjects(i).Trials(j).Actuators.Muscles(k).name = Subjects(i).Trials(j).Actuators.Left.Muscles(k).name;
                Subjects(i).Trials(j).Actuators.Muscles(k).data = mean([...
                    Subjects(i).Trials(j).Actuators.Left.Muscles(k).data.left_parsed, ...
                    Subjects(i).Trials(j).Actuators.Right.Muscles(k).data.right_parsed],2);
                
                %             figure; hold on;
                %             plot(Subjects(i).Trials(j).Actuators.Muscles(k).data, 'k', 'LineWidth', LW);
                %             plot(Subjects(i).Trials(j).Actuators.Left.Muscles(k).data.left_parsed, 'b');
                %             plot(Subjects(i).Trials(j).Actuators.Right.Muscles(k).data.right_parsed, 'r');
                
                Subjects(i).Trials(j).Actuators.data(:,k) = Subjects(i).Trials(j).Actuators.Muscles(k).data;
                Subjects(i).Trials(j).Actuators.colheaders{k} = Subjects(i).Trials(j).Actuators.Muscles(k).name;
            end
            
            % define new structure MuscleActuators
            Conditions = {'Fm40','Fm20','NormF','Fp20','Fp40'};
            k = contains(Conditions, Str{1});
            MuscleActuators(k).Condition = Conditions{k};
            MuscleActuators(k).colheaders = Subjects(i).Trials(j).Actuators.colheaders;
            MuscleActuators(k).Data(:,:,i) = Subjects(i).Trials(j).Actuators.data;
            
            
            clearvars SubjTrialLeft SubjTrialRight left right Str subj trial Array
        end
    end
end

clearvars MetabolicsDir MetabolicsFolder subjectPath TSData i j Nfiles Z ans


%% plot avg muscle forces by condition
PlotCond = 'Yes';
if strcmp(PlotCond, 'Yes')
    clc; close all;
    LW = 1.5; % set line width
    % TotalColors = [rgb('LightGray'); rgb('DarkGray'); rgb('Gray');  rgb('DarkSlateGray'); rgb('Black')];
    TotalColors = [rgb('Blue'); rgb('CornflowerBlue'); rgb('Black');  rgb('Coral'); rgb('OrangeRed')];
    FntSz = 9;
    
    MuscleForces = figure('Position',[50 50 1200 1000]);
    Stride = 1:100;
%     yPk = 1000;
    for i = [2 4 1 5 3]
        %     % erector spinae
        %     subplot(7,3,1); hold on;
        %     Ind = contains([MuscleActuators(i).colheaders], 'ercspn');
        %     plot(Stride, MuscleActuators(i).Avg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
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
        
        MuscleActuators(i).Avg = mean(MuscleActuators(i).Data, 3);
        
        % psoas
        subplot(7,3,1); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'psoas');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Psoas');
        ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % iliacus
        subplot(7,3,2); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'iliacus');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        Pl(i).C = plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Iliacus');
        %     ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % rec_fem
        subplot(7,3,3); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'rect_fem');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Rectus Femoris');
        %     ylabel('N');
        %     ylim([ 0 3]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % glut min
        subplot(7,3,4); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'glut_min');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Gluteus Minimus');
            ylabel('N');
        %     ylim([ 0 3]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % glut med
        subplot(7,3,5); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'glut_med');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Gluteus Medius');
        %     ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % glut max
        subplot(7,3,6); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'glut_max');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Gluteus Maximus');
        %     ylabel('N');
        %     ylim([ 0 3]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % sar
        subplot(7,3,7); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'sar');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Sartorius');
        ylabel('N');
        %     ylim([ 0 2]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % add_mag
        subplot(7,3,8); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'add_mag');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Adductor Magnus');
        %         ylabel('N'); ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % add_long
        subplot(7,3,9); hold on;
                Ind = contains([MuscleActuators(i).colheaders], 'add_long');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Adductor Longus');
        %         ylabel('N'); ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;

        % vas_lat
        subplot(7,3,10); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'vas_lat');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Vastus Lateralis');
            ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % vas_med
        subplot(7,3,11); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'vas_med');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Vastus Medialis');
%         ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % vas_int
        subplot(7,3,12); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'vas_int');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Vastus Intermedialis');
        %     ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        
        % bifemlh
        subplot(7,3,13); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'bifemlh');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Biceps Femoris - Long Head');
        ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % semimem
        subplot(7,3,14); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'semimem');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Semimembranosus');
%         ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % semiten
        subplot(7,3,15); hold on;
%         Ind = contains([MuscleActuators(i).colheaders], 'semiten');
%         plot(Stride, MuscleActuators(i).Avg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
         Ind = contains([MuscleActuators(i).colheaders], 'semiten');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Semitendinosus');
%         ylabel('N'); ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        
        % bifemsh
        subplot(7,3,16); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'bifemsh');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Biceps Femoris - Short Head');
            ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
     
         % tib_ant
        subplot(7,3,17); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'tib_ant');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Tibialis Anterior');
        %     ylabel('N');
        %     ylim([ 0 yPk]);
        %     ax = gca;
        %     ax.XTick = [0 25 50 75 100];
        %     ax.XTickLabel = {'0' '25' '50' '75' '100'};
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        %     xlabel('% Gait Cycle');
        %     ax.YTickLabel = [];
        %     ax.FontSize = FntSz;
        
        % ext_dig
        subplot(7,3,18); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'ext_dig');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Extensor Digitorum');
        %     ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        %     ax = gca;
        %     ax.XTick = [0 25 50 75 100];
        %     ax.XTickLabel = {'0' '25' '50' '75' '100'};
        %     ax.YTickLabel = [];
        %     xlabel('% Gait Cycle');
        %     ax.FontSize = FntSz;
        
    % lat_gas
        subplot(7,3,19); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'lat_gas');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Lateral Gastrocnemius');
        ylabel('N');
        %     ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = {'0' '25' '50' '75' '100'};
        xlabel('% Gait Cycle');
        ax.FontSize = FntSz;
        
        
        % med_gas
        subplot(7,3,20); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'med_gas');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Medial Gastrocnemius');
        %     ylabel('N');
        %     ylim([ 0 yPk]);
        %     ax = gca;
        %     ax.XTick = [0 25 50 75 100];
        %     ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        %     ax.FontSize = FntSz;
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = {'0' '25' '50' '75' '100'};
        %     ax.YTickLabel = [];
        xlabel('% Gait Cycle');
        ax.FontSize = FntSz;
        
        
           % soleus
        subplot(7,3,21); hold on;
        Ind = contains([MuscleActuators(i).colheaders], 'soleus');
        p = sum( MuscleActuators(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Soleus');
        %     ylabel('N');
        %     ylim([ 0 6]);
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
        
    end
    
    subplot(732); hold on; 
%     Conditions = {'-40', '-20','Norm','+20','+40'};
%     legend(Conditions, 'Location','Best');
    
    Conditions = {'+40', '+20','Norm','-20','-40'};
    % legend([Pl(4).C, Pl(3).C, Pl(5).C, Pl(1).C, Pl(3).C], Conditions, 'Location','Best');
    legend([Pl(5).C, Pl(4).C, Pl(3).C, Pl(2).C, Pl(1).C], Conditions, 'Location','Best', ...
        'FontSize',6.5);
    
    subplotsqueeze(MuscleForces, 1.12);
    supertitle({'MUSCLE FORCES'; ' '}, 'FontSize',16);
    saveas(MuscleForces, 'MuscleForces.png');
    % saveas(MuscleForces, 'MuscleForces.pdf');
    
end

end
