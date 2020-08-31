function [Activations, Subjects] = CheckActivations(ResultsFolder, Subjects, PlotActivations, PlotCond)

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
% clc;
% IsStoFile = ~[ActivationsDir.isdir];
% Log = zeros(1, length(ActivationsDir));
% for i = 1:length(ActivationsDir)
%     Log(i) = IsStoFile(i) == 1 &&  contains(ActivationsDir(i).name, 'controls');
% end
%
% Files = find(Log);
% L = length(Files);
%
% Activations(L).Name = [];
% Activations(L).Data = [];
% NumSubj = length(Subjects);
% % Subjects(NumSubj).Trials(:).Activations = [];
%
% j = 1;
% for i = Files
%
%     % get filename
%     Activations(j).Name = ActivationsDir(i).name;
%
%     % load data
%     Data = importdata(ActivationsDir(i).name);
%     Data.filename = ActivationsDir(i).name;
%
%     % get subject and trial
%     Str = strsplit(Activations(j).Name, '_');
%     Subj = Str{1};
%     SubjInd = contains({Subjects.name}, Subj);
%     Trial = Str{2};
%     TrialInd = contains({Subjects(SubjInd).Trials.name}, Trial);
%     Side = Str{3};
%
%     TSData = Subjects(SubjInd).Trials(TrialInd).TSData; % identify temporal spatial data
%     Activations(j).Data = FilterAndParseData(Data, TSData); % filter and parse to gait cycles
%
%     ActData = Activations(j).Data;
%     Fields2Del = {'data','textdata','Interp','Fdata','Filter','Fdata2','Filter2'};
%     ActData = rmfield(ActData, Fields2Del);
%
%     % match activations to subject and trial structure
%     if strcmp(Side, 'Left') % assing to subjects structure
%         Subjects(SubjInd).Trials(TrialInd).Activations.Left =  ActData;
%     elseif strcmp(Side, 'Right')
%         Subjects(SubjInd).Trials(TrialInd).Activations.Right =  ActData;
%     end
%
%     j = j+1;
% end

NumSubj = length(Subjects);

[Activations, Subjects] = ExtractData(ResultsFolder, Subjects, 'activations');

clearvars i j IsStoFile fs Wn1 Cutoff1 b a L j Str Subj SubjInd Trial TrialInd Side...
    TSData Data ActivationsDir ActivationsFolder S ActData Log Dields2Del Files

%% Plot muscle activations
% PlotActivations = 'Yes';
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

%% Combine left and right trials for the same subject into structure
clc;
for i = 1:length(Subjects) % loop through subjects
    
    for j = 1:length(Subjects(i).Trials) % loop through trials
        
        if strcmp(Subjects(i).Trials(j).type, 'static') == 0
            
            subj = contains({Activations.Subject}, Subjects(i).name);
            Str = strsplit(Subjects(i).Trials(j).name, '_');
            trial = contains({Activations.Trial}, Str{1});
            left =  contains({Activations.Side}, 'Left');
            right = contains({Activations.Side}, 'Right');
            
            SubjTrialLeft = subj + trial + left == 3;
            SubjTrialRight = subj + trial + right == 3;
            
            if sum(SubjTrialLeft) > 0
                Subjects(i).Trials(j).Activations.Left = Activations(SubjTrialLeft);
            end
            if sum(SubjTrialRight) > 0
                Subjects(i).Trials(j).Activations.Right = Activations(SubjTrialRight);
            end
            
            [Subjects(i).Trials(j).Activations.Left.Muscles, ~] = GetMuscles(Subjects(i).Trials(j).Activations.Left.Data.Parsed, ...
                Subjects(i).Trials(j).Activations.Left.Data.colheaders, 'notmet');
            
            [Subjects(i).Trials(j).Activations.Right.Muscles, ~] = GetMuscles(Subjects(i).Trials(j).Activations.Right.Data.Parsed, ...
                Subjects(i).Trials(j).Activations.Right.Data.colheaders, 'notmet');
            
            % average sides and combine together
            for k = 1:length(Subjects(i).Trials(j).Activations.Left.Muscles)
                Subjects(i).Trials(j).Activations.Muscles(k).name = Subjects(i).Trials(j).Activations.Left.Muscles(k).name;
                Subjects(i).Trials(j).Activations.Muscles(k).data = mean([...
                    Subjects(i).Trials(j).Activations.Left.Muscles(k).data.left_parsed, ...
                    Subjects(i).Trials(j).Activations.Right.Muscles(k).data.right_parsed],2);
                
                %             figure; hold on;
                %             plot(Subjects(i).Trials(j).Activations.Muscles(k).data, 'k', 'LineWidth', LW);
                %             plot(Subjects(i).Trials(j).Activations.Left.Muscles(k).data.left_parsed, 'b');
                %             plot(Subjects(i).Trials(j).Activations.Right.Muscles(k).data.right_parsed, 'r');
                
                Subjects(i).Trials(j).Activations.data(:,k) = Subjects(i).Trials(j).Activations.Muscles(k).data;
                Subjects(i).Trials(j).Activations.colheaders{k} = Subjects(i).Trials(j).Activations.Muscles(k).name;
            end
            
            % define new structure MuscleActivations
            Conditions = {'Fm40','Fm20','NormF','Fp20','Fp40'};
            k = contains(Conditions, Str{1});
            MuscleActivations(k).Condition = Conditions{k};
            MuscleActivations(k).colheaders = Subjects(i).Trials(j).Activations.colheaders;
            MuscleActivations(k).Data(:,:,i) = Subjects(i).Trials(j).Activations.data;
            
            
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
    
    MuscleActivationsFig = figure('Position',[50 50 1200 1000]);
    Stride = 1:100;
    yPk = 1;
    for i = [2 4 1 5 3]
        
        %     % erector spinae
        %     subplot(6,3,1); hold on;
        %     Ind = contains([MuscleActivations(i).colheaders], 'ercspn');
        %     plot(Stride, MuscleActivations(i).Avg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
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
        
        MuscleActivations(i).Avg = mean(MuscleActivations(i).Data, 3);
        
  % psoas
        subplot(7,3,1); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'psoas');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Psoas');
        ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
            ax.YTickLabel = [0 50 100];
        ax.FontSize = FntSz;
        
        % iliacus
        subplot(7,3,2); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'iliacus');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        Pl(i).C = plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Iliacus');
        %     ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % rec_fem
        subplot(7,3,3); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'rect_fem');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Rectus Femoris');
        %     ylabel('%');
         ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % glut min
        subplot(7,3,4); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'glut_min');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Gluteus Minimus');
            ylabel('%');
         ylim([ 0 yPk]);
        ax = gca;
          ax.YTickLabel = [0 50 100];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % glut med
        subplot(7,3,5); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'glut_med');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Gluteus Medius');
        %     ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % glut max
        subplot(7,3,6); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'glut_max');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Gluteus Maximus');
        %     ylabel('%');
         ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % sar
        subplot(7,3,7); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'sar');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Sartorius');
        ylabel('%');
         ylim([ 0 yPk]);
        ax = gca;
          ax.YTickLabel = [0 50 100];
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % add_mag
        subplot(7,3,8); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'add_mag');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Adductor Magnus');
         ylim([ 0 yPk]);
        %         ylabel('%'); ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % add_long
        subplot(7,3,9); hold on;
                Ind = contains([MuscleActivations(i).colheaders], 'add_long');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Adductor Longus');
        %         ylabel('%'); ylim([ 0 yPk]);
         ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        %     ax.YTick = [0 1 2 3];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;

        % vas_lat
        subplot(7,3,10); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'vas_lat');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Vastus Lateralis');
            ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
          ax.YTickLabel = [0 50 100];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % vas_med
        subplot(7,3,11); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'vas_med');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Vastus Medialis');
%         ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % vas_int
        subplot(7,3,12); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'vas_int');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Vastus Intermedialis');
        %     ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        
        % bifemlh
        subplot(7,3,13); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'bifemlh');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Biceps Femoris - Long Head');
        ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
          ax.YTickLabel = [0 50 100];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        % semimem
        subplot(7,3,14); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'semimem');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Semimembranosus');
%         ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
        % semiten
        subplot(7,3,15); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'semiten');
        plot(Stride, MuscleActivations(i).Avg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Semitendinosus');
%         ylabel('%'); 
        ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        ax.FontSize = FntSz;
        
        
        % bifemsh
        subplot(7,3,16); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'bifemsh');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Biceps Femoris - Short Head');
            ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
          ax.YTickLabel = [0 50 100];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        
     
         % tib_ant
        subplot(7,3,17); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'tib_ant');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Tibialis Anterior');
        %     ylabel('%');
            ylim([ 0 yPk]);
        %     ax = gca;
        %     ax.XTick = [0 25 50 75 100];
        %     ax.XTickLabel = {'0' '25' '50' '75' '100'};
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        ax.FontSize = FntSz;
        %     xlabel('% Gait Cycle');
        %     ax.YTickLabel = [];
        %     ax.FontSize = FntSz;
        
        % ext_dig
        subplot(7,3,18); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'ext_dig');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Extensor Digitorum');
        %     ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
        ax.YTickLabel = [];
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
        Ind = contains([MuscleActivations(i).colheaders], 'lat_gas');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Lateral Gastrocnemius');
        ylabel('%');
            ylim([ 0 yPk]);
        ax = gca;
          ax.YTickLabel = [0 50 100];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = {'0' '25' '50' '75' '100'};
        xlabel('% Gait Cycle');
        ax.FontSize = FntSz;
        
        
        % med_gas
        subplot(7,3,20); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'med_gas');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Medial Gastrocnemius');
        %     ylabel('%');
            ylim([ 0 yPk]);
        %     ax = gca;
        %     ax.XTick = [0 25 50 75 100];
        %     ax.XTickLabel = [];
        %     ax.YTickLabel = [];
        %     ax.FontSize = FntSz;
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = {'0' '25' '50' '75' '100'};
        %     ax.YTickLabel = [];
        xlabel('% Gait Cycle');
        ax.FontSize = FntSz;
        
        
           % soleus
        subplot(7,3,21); hold on;
        Ind = contains([MuscleActivations(i).colheaders], 'soleus');
        p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
        plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
        title('Soleus');
        %     ylabel('%');
        ylim([ 0 yPk]);
        %     ax = gca;
        %     ax.XTick = [0 25 50 75 100];
        %     ax.XTickLabel = [];
        %     ax.FontSize = FntSz;
        ax = gca;
        ax.YTickLabel = [];
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = {'0' '25' '50' '75' '100'};
        %     ax.YTickLabel = [];
        xlabel('% Gait Cycle');
        ax.FontSize = FntSz;
        
    end
    
    subplot(732); hold on; 
        
%         % psoas
%         subplot(6,3,1); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'psoas');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Psoas');
%         ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
% %         ax.YTick = [0 0.5 1];
% %         ax.YTickLabel = {'0', '50', '100'};
%         ax.FontSize = FntSz;
%         
%         % iliacus
%         subplot(6,3,2); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'iliacus');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         Pl(i).C = plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Iliacus');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         % glut max
%         subplot(6,3,6); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'glut_max');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Gluteus Maximus');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         % glut med
%         subplot(6,3,5); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'glut_med');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Gluteus Medius');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         %         % add_mag
%         %     subplot(6,3,6); hold on;
%         %     Ind = contains([MuscleActivations(i).colheaders], 'add_mag');
%         %     plot(Stride, MuscleActivations(i).Avg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
%         %     title('Adductor Mag');
%         %     ylabel('%'); ylim([ 0 yPk]);
%         
%         % rec_fem
%         subplot(6,3,3); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'rect_fem');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Rectus Femoris');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         %     ax.YTick = [0 1 2 3];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         
%         % sar
%         subplot(6,3,4); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'sar');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Sartorius');
%         ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         %     ax.YTick = [0 1 2 3];
%         ax.XTickLabel = [];
%         %     ax.YTickLabel = [];
%         ax.YTick = [0 0.5 1];
%         ax.YTickLabel = {'0', '50', '100'};
%         ax.FontSize = FntSz;
%         
%         % vas_med
%         subplot(6,3,7); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'vas_med');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Vastus Medialis');
%         ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTick = [0 0.5 1];
%         ax.YTickLabel = {'0', '50', '100'};
%         ax.FontSize = FntSz;
%         
%         % vas_int
%         subplot(6,3,8); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'vas_int');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Vastus Intermedialis');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         % vas_lat
%         subplot(6,3,9); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'vas_lat');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Vastus Lateralis');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         % bifemlh
%         subplot(6,3,10); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'bifemlh');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Biceps Femoris - Long Head');
%         ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTick = [0 0.5 1];
%         ax.YTickLabel = {'0', '50', '100'};
%         ax.FontSize = FntSz;
%         
%         % bifemsh
%         subplot(6,3,11); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'bifemsh');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Biceps Femoris - Short Head');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         % semimem
%         subplot(6,3,13); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'semimem');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Semimembranosus');
%         ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         %     ax.YTickLabel = [];
%         ax.YTick = [0 0.5 1];
%         ax.YTickLabel = {'0', '50', '100'};
%         ax.FontSize = FntSz;
%         
%         %     % semiten
%         %     subplot(6,3,13); hold on;
%         %     Ind = contains([MuscleActivations(i).colheaders], 'semiten');
%         %     plot(Stride, MuscleActivations(i).Avg(:,Ind), 'LineWidth', LW, 'Color', TotalColors(i,:));
%         %     title('Semitendinosus');
%         %     ylabel('%'); ylim([ 0 yPk]);
%         %     ax = gca;
%         %     ax.XTick = [0 25 50 75 100];
%         %     ax.XTickLabel = [];
%         %     ax.FontSize = FntSz;
%         %     % SPM shading
%         %     if i == 5
%         %         SPMshade('semiten', MetSPM, 2)
%         %     end
%         
%         % add_long
%         subplot(6,3,12); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'add_long');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Adductor Longus');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         
%         % soleus
%         subplot(6,3,18); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'soleus');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Soleus');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = {'0' '25' '50' '75' '100'};
%         ax.YTickLabel = [];
%         xlabel('% Gait Cycle');
%         ax.FontSize = FntSz;
%         
%         % med_gas
%         subplot(6,3,17); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'med_gas');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Medial Gastrocnemius');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = {'0' '25' '50' '75' '100'};
%         ax.YTickLabel = [];
%         xlabel('% Gait Cycle');
%         ax.FontSize = FntSz;
%         
%         % lat_gas
%         subplot(6,3,16); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'lat_gas');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Lateral Gastrocnemius');
%         ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = {'0' '25' '50' '75' '100'};
%         ax.YTick = [0 0.5 1];
%         ax.YTickLabel = {'0', '50', '100'};
%         xlabel('% Gait Cycle');
%         ax.FontSize = FntSz;
%         
%         % tib_ant
%         subplot(6,3,14); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'tib_ant');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Tibialis Anterior');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         %     ax = gca;
%         %     ax.XTick = [0 25 50 75 100];
%         %     ax.XTickLabel = {'0' '25' '50' '75' '100'};
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         %     xlabel('% Gait Cycle');
%         %     ax.YTickLabel = [];
%         %     ax.FontSize = FntSz;
%         
%         % ext_dig
%         subplot(6,3,15); hold on;
%         Ind = contains([MuscleActivations(i).colheaders], 'ext_dig');
%         p = sum( MuscleActivations(i).Avg(:,Ind), 2) / 2;
%         plot(Stride, p, 'LineWidth', LW, 'Color', TotalColors(i,:));
%         title('Extensor Digitorum');
%         %     ylabel('%');
%         ylim([ 0 yPk]);
%         ax = gca;
%         ax.XTick = [0 25 50 75 100];
%         ax.XTickLabel = [];
%         ax.YTickLabel = [];
%         ax.FontSize = FntSz;
%         %     ax = gca;
%         %     ax.XTick = [0 25 50 75 100];
%         %     ax.XTickLabel = {'0' '25' '50' '75' '100'};
%         %     ax.YTickLabel = [];
%         %     xlabel('% Gait Cycle');
%         %     ax.FontSize = FntSz;
%         
%     end
%     
%     subplot(632);
%     Conditions = {'-40', '-20','Norm','+20','+40'};
%     legend(Conditions, 'Location','Best');
    
    Conditions = {'+40', '+20','Norm','-20','-40'};
    % legend([Pl(4).C, Pl(3).C, Pl(5).C, Pl(1).C, Pl(3).C], Conditions, 'Location','Best');
    legend([Pl(5).C, Pl(4).C, Pl(3).C, Pl(2).C, Pl(1).C], Conditions, 'Location','Best', ...
        'FontSize',6.5);
    
    subplotsqueeze(MuscleActivationsFig, 1.12);
    supertitle({'MUSCLE ACTIVATIONS'; ' '}, 'FontSize',16);
    saveas(MuscleActivationsFig, 'MuscleActivationsCond.png');
    % saveas(MuscleForces, 'MuscleForces.pdf');
end
%% Compare to literature?

%% Export

end