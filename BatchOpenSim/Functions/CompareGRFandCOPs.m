function [Forces, GRF] = CompareGRFandCOPs(ForcesFileName, GRF, ToPlot)

%% Settings
if exist('ToPlot', 'var') == 0
    ToPlot = 'No';
end

SaveAsPNG = 'No'; % Set to yes to save comparison fig as PNG file


%% load both files and identify columns with vertical ground reaction force
% Load forces
if exist('ForcesFileName', 'var') == 0
    [Forces.name, Forces.path] = uigetfile('.forces','Select FORCES file','Select FORCES file to load');
    addpath(genpath(Forces.path));
else
    Forces.name = ForcesFileName;
    addpath(genpath(ForcesFileName));
end
Forces.data = importdata(Forces.name);
Forces.data.colheaders =  {'Sample','FX1','FY1','FZ1','X1','Y1','Z1','MZ1','FX2','FY2','FZ2','X2','Y2','Z2','MZ2'};

% load GRF data
if exist('GRF', 'var') == 0
    GRF = LoadGRF(0);
end

MainTrialName = strsplit(Forces.name, '_');

%% rename data
Forces.Right.X = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FX1'));
Forces.Right.Y = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FY1'));
Forces.Right.Z = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FZ1'));
Forces.Left.X = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FX2'));
Forces.Left.Y = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FY2'));
Forces.Left.Z = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FZ2'));

GRF.Right.X = GRF.AllData(:,strcmp(GRF.ColHeaders, 'ground_force_vx'));
GRF.Right.Y = GRF.AllData(:,strcmp(GRF.ColHeaders, 'ground_force_vz'));
GRF.Right.Z = GRF.AllData(:,strcmp(GRF.ColHeaders, 'ground_force_vy'));
GRF.Left.X = GRF.AllData(:,strcmp(GRF.ColHeaders, '1_ground_force_vx'));
GRF.Left.Y = GRF.AllData(:,strcmp(GRF.ColHeaders, '1_ground_force_vz'));
GRF.Left.Z = GRF.AllData(:,strcmp(GRF.ColHeaders, '1_ground_force_vy'));


%% Plot
if strcmp(ToPlot, 'Yes')
    LW = 1;
    MkrSz = 10;
    Ind = 1:5000;
    ForceComp = figure('Position',[100 100 1200 800]);
    % Left side
    subplot(321); hold on;
    OpenSimX = plot(GRF.Left.X(Ind), '-r', 'LineWidth',LW);
    ForcesX = plot(Forces.Left.X(Ind), '--b',  'LineWidth',LW);
    title('Left X Forces');
    xlabel('Frame'); ylabel('N');
    legend({'OpenSim', 'Lab'});
    
    subplot(323); hold on;
    OpenSimX = plot(-GRF.Left.Y(Ind), '-r', 'LineWidth',LW);
    ForcesX = plot(Forces.Left.Y(Ind), '--b',  'LineWidth',LW);
    title('Left Y Forces');
    xlabel('Frame'); ylabel('N');
    
    subplot(325); hold on;
    OpenSimZ = plot(GRF.Left.Z(Ind), '-r', 'LineWidth',LW);
    ForcesZ = plot(Forces.Left.Z(Ind), '--b',  'LineWidth',LW);
    title('Left Z Forces');
    xlabel('Frame'); ylabel('N');
    
    % Right side
    subplot(322); hold on;
    OpenSimX = plot(GRF.Right.X(Ind), '-r', 'LineWidth',LW);
    ForcesX = plot(Forces.Right.X(Ind), '--b',  'LineWidth',LW);
    title('Right X Forces');
    xlabel('Frame'); ylabel('N');
    legend({'OpenSim', 'Lab'});
    
    subplot(324); hold on;
    OpenSimX = plot(-GRF.Right.Y(Ind), '-r', 'LineWidth',LW);
    ForcesX = plot(Forces.Right.Y(Ind), '--b',  'LineWidth',LW);
    title('Right Y Forces');
    xlabel('Frame'); ylabel('N');
    
    subplot(326); hold on;
    OpenSimZ = plot(GRF.Right.Z(Ind), '-r', 'LineWidth',LW);
    ForcesZ = plot(Forces.Right.Z(Ind), '--b',  'LineWidth',LW);
    title('Right Z Forces');
    xlabel('Frame'); ylabel('N');
    
    ForceName = strrep(Forces.name,'_', ' ');
    GRFName = strrep(GRF.FileName,'_', ' ');
    Title = {strcat(ForceName, ' and OpenSimGRF.mot comparison'), ' '};
    supertitle(Title);
    
    % Save figure?
    if strcmp(SaveAsPNG, 'Yes')
        saveas(ForceComp, strcat(MainTrialName{1}, '_ForceComp.png'));
    end
end

end


