function[GRFdata] = LoadGRF(Files, ToPlot)

% Load OpenSim ground reaction force data

% INPUTS: 

% files         files to load (optional)

% ToPlot      if 'Yes', will plot GRF comparisons

% Ricky Pimentel - Applied Biomechanics Lab - UNC Chapel Hill
% January 2020

%% Input Parameters?
if exist('ToPlot','var') == 0
    ToPlot = 0; % default plot setting = off
end

%% Select Files for input
if Files == 0 % if no files defined
    [Files, path] = uigetfile('.mot','Select GRF .mot files to input', 'Multiselect','On'); % select GRF .mot files to load
    addpath(genpath(path)); % add file path to matlab path
else % if files already defined
    if iscell(Files) % if multiple files defined
        NumTrials = length(Files); % determine # of files load
        for i = 1:NumTrials 
            path = fileparts(Files(i));
            addpath(genpath(path)); % add file path to matlab path
        end
    else
        path = fileparts(Files); 
        addpath(genpath(path)); % add file path to matlab path
    end
end


%% load data
if iscell(Files)
    NumTrials = length(Files);
    GRFdata(NumTrials).FileName = 'placeholder'; % initialize variables
    GRFdata(NumTrials).AllData = 'placeholder';
    GRFdata(NumTrials).HeaderInfo = 'placeholder';
    for i = 1:NumTrials
        data = importdata(Files{i});
        GRFdata(i).FileName = Files{i};
        GRFdata(i).AllData = data.data; 
        GRFdata(i).HeaderInfo = data.textdata;
    end
    clearvars data
else
    NumTrials = 1;
    data = importdata(Files);
    GRFdata.FileName = Files;
    GRFdata.AllData = data.data; 
    GRFdata.HeaderInfo = data.textdata;
end

% quality assurance for trial loading
% for i = 1:NumTrials
%     if strfind(GRFdata(i).textdata(7,1), 'ground')
%         error(['No ground reaction force data detected. Please select ground reaction force data. Error on trial ' num2str(i)]);
%     end
% end

%% set up variables prior to plotting

% redefine headers
Headers = {'time',	'ground_force_vx',	'ground_force_vy',	'ground_force_vz',...
    'ground_force_px',	'ground_force_py',	'ground_force_pz',...
    '1_ground_force_vx',	'1_ground_force_vy',	'1_ground_force_vz',...
    '1_ground_force_px',	'1_ground_force_py',	'1_ground_force_pz',...
    'ground_torque_x',	'ground_torque_y',	'ground_torque_z',...
    '1_ground_torque_x',	'1_ground_torque_y',	'1_ground_torque_z'};
for i = 1:NumTrials
    GRFdata(i).ColHeaders = Headers; 
end

LW = 2; % set line width


%% plot ground reaction forces
if ToPlot == 1
    for i = 1:NumTrials
        
        GRFdata(i).Time = GRFdata(i).AllData(:,1); % define timeseries as first column
        
        FPData = figure('Position',[50 50 1200 800]); % create figure with set dimensions
        
        subplot(321); % plot GRFs
        plot(GRFdata(i).Time, GRFdata(i).AllData(:,[8 9 10]), 'LineWidth',LW);
        title('Plate 2 (Left) GRFs');
        legend({'X','Y','Z'});
        ylabel('N');
        
        subplot(322); % plot GRFs
        plot(GRFdata(i).Time, GRFdata(i).AllData(:,[2 3 4]), 'LineWidth',LW);
        title('Plate 1 (Right) GRFs');
        legend({'X','Y','Z'});
        
        subplot(323); % plot COPs
        plot(GRFdata(i).Time, GRFdata(i).AllData(:,[11 12 13]), 'LineWidth',LW);
        title('Plate 2 (Left) COPs');
        ylabel('m');
        
        subplot(324); % plot COPs
        plot(GRFdata(i).Time, GRFdata(i).AllData(:,[5 6 7]), 'LineWidth',LW);
        title('Plate 1 (Right) COPs');
        
        subplot(325); % plot GRFs
        plot(GRFdata(i).Time, GRFdata(i).AllData(:,[17 18 19]), 'LineWidth',LW);
        title('Plate 2 (Left) Torques');
        xlabel('s')
        ylabel('Nm');
        
        subplot(326); % plot Torques
        plot(GRFdata(i).Time, GRFdata(i).AllData(:,[14 15 16]), 'LineWidth',LW);
        title('Plate 1 (Right) Torques');
        xlabel('s')
        
        % set range to plot 
%         if exist('Range2Plot','var') == 0
%             Range2Plot = [GRFdata(i).Time(1) GRFdata(i).Time(end)]; % default range to plot
%         end
%         for j = 1:6
%             subplot(3,2,j);
%             xlim([Range2Plot(1) Range2Plot(2)]);
%         end
        
        % subplot squeeze
        % save figure?
    end
end

%% Display trials loaded
for i = 1:NumTrials
    disp([GRFdata(i).FileName ' loaded and saved in structure']); 
end

end

