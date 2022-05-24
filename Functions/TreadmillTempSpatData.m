function[TSData] = TreadmillTempSpatData(GRF) % , SubjDemos

% calculate gait events from treadmill force plate GRF data
% this script assumes no crossover steps on the treadmill (left foot on
% plate 2, right foot on plate 1)

%% Input Settings
L = length(GRF.AllData); % get length
TSData.FileName = GRF.FileName; % save filename
% save column headers for L/R_Strike and L/R_Off timing matricies. 
% Analog_Index = the time of gait event in the analog capture frequency (usually 1000 Hz)
% Event_Time = the time of gait event in seconds
% Digital_Index = the time of gait event in the digital capture frequency (usually 100 Hz)
TSData.StrikeOffColHeaders = {'Analog_Index','Event_Time', 'Digital_Index'}; 

% Locate correct columns for left and right sides and time
if sum(contains(GRF.ColHeaders, 'ground_force')) > 1
    RightVertGRF = contains(GRF.ColHeaders, 'ground_force_vy');
    LeftVertGRF = contains(GRF.ColHeaders, '1_ground_force_vy');  
    TimeCol = contains(GRF.ColHeaders, 'time');
elseif sum(contains(GRF.ColHeaders, 'FP')) > 1
    RightVertGRF = contains(GRF.ColHeaders, 'FP1_vy');
    LeftVertGRF = contains(GRF.ColHeaders, 'FP2_vy');  
    TimeCol = contains(GRF.ColHeaders, 'Header');
end

% get sampling frequency
Time = GRF.AllData(:, TimeCol);
TSData.AnalogSampFreq = 1 / (Time(2) - Time(1)); %

%% Find Step Times

% Left Side (Belt/Plate 2)
ZeroTimes = abs(GRF.AllData(:,LeftVertGRF)) < 20;  % find zero times
NumStrikes = 0; % counter for foot strikes
NumOffs = 0;  % counter for foot offs
for i = 2:L-1 % loop through all times 
    if ZeroTimes(i) == 0 && ZeroTimes(i-1) == 1 % find when foot hits treadmill and save as strike
        NumStrikes = NumStrikes + 1; 
        TSData.L_Strike(NumStrikes, :) = [i, Time(i)]; % save index and time
    elseif ZeroTimes(i-1) == 0 && ZeroTimes(i) == 1 % find when foot leaves treadmill and save as off
        NumOffs = NumOffs + 1;
        TSData.L_Off(NumOffs, :) = [i, Time(i)]; % save index and time
    end
end

TSData.L_NumStrikes = NumStrikes; % save number of strikes
TSData.L_NumOffs = NumOffs; % save number of offs
clearvars ZeroTimes NumStrikes NumOffs

% Right Side (Belt/Plate 1)
ZeroTimes = abs(GRF.AllData(:,RightVertGRF)) < 20;  % find zero times
NumStrikes = 0; % counter for foot strikes
NumOffs = 0;  % counter for foot offs
for i = 2:L-1 % loop through all times 
     if ZeroTimes(i) == 0 && ZeroTimes(i-1) == 1 % find when foot hits treadmill and save as strike
        NumStrikes = NumStrikes + 1; 
        TSData.R_Strike(NumStrikes, :) = [i, Time(i)]; % save index and time
    elseif ZeroTimes(i-1) == 0 && ZeroTimes(i) == 1 % find when foot leaves treadmill and save as off
        NumOffs = NumOffs + 1;
        TSData.R_Off(NumOffs, :) = [i, Time(i)]; % save index and time
    end
end

TSData.R_NumStrikes = NumStrikes; % save number of strikes
TSData.R_NumOffs = NumOffs; % save number of offs
clearvars ZeroTimes NumStrikes NumOffs


%% Make sure trials start with a strike and end with an off
% left
if TSData.L_Strike(end,1) > TSData.L_Off(end,1) % if trial ends with a foot strike -> delete last strike
    TSData.L_Strike(end,:) = []; 
    TSData.L_NumStrikes = TSData.L_NumStrikes - 1;
end
if TSData.L_Off(1,1) < TSData.L_Strike(1,1) % if trial starts with a foot off -> delete first off
    TSData.L_Off(1,:) = []; 
    TSData.L_NumOffs = TSData.L_NumOffs - 1;
end
% right
if TSData.R_Strike(end,1) > TSData.R_Off(end,1) % if trial ends with a foot strike -> delete last strike
    TSData.R_Strike(end,:) = []; 
    TSData.R_NumStrikes = TSData.R_NumStrikes - 1;
end
if TSData.R_Off(1,1) < TSData.R_Strike(1,1) % if trial starts with a foot off -> delete first off
    TSData.R_Off(1,:) = []; 
    TSData.R_NumOffs = TSData.R_NumOffs - 1;
end
clearvars ZeroTimes NumStrikes NumOffs


%% plot if desired
Plot = 'No'; 
if strcmp(Plot, 'Yes')
    figure;
    subplot(211); hold on;  % left side
    plot(GRF.AllData(:,LeftVertGRF), '-k');
    vline(TSData.L_Strike(:,1),'g');
    vline(TSData.L_Off(:,1), 'r');
    title('Left Gait Events');
    ylabel('vGRF'); xlabel('Frame');
    xlim([0 10000]);
    
    subplot(212); hold on;  % left side
    plot(GRF.AllData(:,RightVertGRF), '-k');
    vline(TSData.R_Strike(:,1),'g');
    vline(TSData.R_Off(:,1), 'r');
    title('Right Gait Events');
    ylabel('vGRF'); xlabel('Frame');
    xlim([0 10000]);
end

%% Ensure gait cycles are valid
% exclude step if time < 0.5 s and peak vGRF < 50 N
TimeThresh = 0.5; % set gait cycle timing threshold (in seconds) to identify questionable steps
VertThresh = 50; % vGRF thresh to delete invalid steps

% LEFT
ToDel = zeros(1, TSData.L_NumOffs); 
for i = 1:TSData.L_NumOffs - 1
    if TSData.L_Strike(i+1,2) - TSData.L_Strike(i,2) < TimeThresh
        frames = TSData.L_Strike(i,1):TSData.L_Strike(i+1,1);
        Step = GRF.AllData(frames,LeftVertGRF);
        if max(Step) < VertThresh
            ToDel(i) = 1; 
        end
        % plot 
%         time = TSData.L_Strike(i,1)-100:TSData.L_Strike(i+1,1)+100;
%         plot(GRF.AllData(time,LeftVertGRF), '-k');
%         title('Left vGRF'); 
        clearvars Step frames time 
    end
end
ToDelete = logical(ToDel); 
TSData.L_Strike(ToDelete, :) = []; 
TSData.L_Off(ToDelete, :) = []; 
TSData.L_NumStrikes = length(TSData.L_Strike); 
TSData.L_NumOffs = length(TSData.L_Off); 
clearvars ToDelete ToDel; 

% RIGHT
ToDel = zeros(1, TSData.R_NumOffs); 
for i = 1:TSData.R_NumOffs - 1
    if TSData.R_Strike(i+1,2) - TSData.R_Strike(i,2) < TimeThresh
        frames = TSData.R_Strike(i,1):TSData.R_Strike(i+1,1);
        Step = GRF.AllData(frames,RightVertGRF);
        if max(Step) < VertThresh
            ToDel(i) = 1; 
        end
        % plot 
%         time = TSData.R_Strike(i,1)-100:TSData.R_Strike(i+1,1)+100;
%         plot(GRF.AllData(time,RightVertGRF), '-k');
%         title('Right vGRF'); 
        clearvars Step frames time 
    end
end
ToDelete = logical(ToDel); 
TSData.R_Strike(ToDelete, :) = []; 
TSData.R_Off(ToDelete, :) = []; 
TSData.R_NumStrikes = length(TSData.R_Strike); 
TSData.R_NumOffs = length(TSData.R_Off); 
clearvars ToDelete ToDel; 


%% Calculate Temporal Spatial Paramaters from step timing info
% add later if necessary?

%%
TSData.NumGCs = min([TSData.L_NumStrikes TSData.R_NumStrikes]);

end