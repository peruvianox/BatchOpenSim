function[TSData] = TreadmillTempSpatData(GRF) % , SubjDemos

% calculate gait events from treadmill force plate GRF data
% this script assumes no crossover steps on the treadmill (left foot on
% plate 2, right foot on plate 1)

%% Input Settings
% SampFreq = 100; % sampling frequency in Hz
L = length(GRF.AllData); 

TSData.FileName = GRF.FileName; 

%% Locate correct columns
for i = 1:length(GRF.ColHeaders)
    if strcmp(GRF.ColHeaders{i}, 'ground_force_vy')
        RightVertGRF = i; 
    elseif strcmp(GRF.ColHeaders{i}, '1_ground_force_vy')
        LeftVertGRF = i;
    end
end

%% Find Step Times

% Left Side (Belt/Plate 2)
ZeroTimes = GRF.AllData(:,LeftVertGRF) == 0;  % find zero times
NumStrikes = 0; % counter for foot strikes
NumOffs = 0;  % counter for foot offs
for i = 2:L-1 % loop through all times 
    if ZeroTimes(i) == 1 && ZeroTimes(i+1) == 0 % find when foot hits treadmill and save as strike
        NumStrikes = NumStrikes + 1; 
        TSData.L_Strike(NumStrikes, :) = [i, GRF.AllData(i+1,1), round(i/10)]; % save index and time
    elseif ZeroTimes(i) == 0 && ZeroTimes(i+1) == 1 % find when foot leaves treadmill and save as off
        NumOffs = NumOffs + 1;
        TSData.L_Off(NumOffs, :) = [i, GRF.AllData(i+1,1), round(i/10)]; % save index and time
    end
end
if TSData.L_Strike(1,3) == 0
    TSData.L_Strike(1,:) = [];
    NumStrikes = NumStrikes - 1; 
    NumOffs = NumOffs - 1;
end
TSData.L_NumStrikes = NumStrikes; % save number of strikes
TSData.L_NumOffs = NumOffs; % save number of offs
clearvars ZeroTimes NumStrikes NumOffs

% Right Side (Belt/Plate 1)
ZeroTimes = GRF.AllData(:,RightVertGRF) == 0;  % find zero times
NumStrikes = 0; % counter for foot strikes
NumOffs = 0;  % counter for foot offs
for i = 2:L-1 % loop through all times 
    if ZeroTimes(i) == 1 && ZeroTimes(i+1) == 0 % find when foot hits treadmill and save as strike
        NumStrikes = NumStrikes + 1; 
        TSData.R_Strike(NumStrikes, :) = [i, GRF.AllData(i+1,1), round(i/10)]; % save index and time
    elseif ZeroTimes(i) == 0 && ZeroTimes(i+1) == 1% find when foot leaves treadmill and save as off
        NumOffs = NumOffs + 1;
        TSData.R_Off(NumOffs, :) = [i, GRF.AllData(i+1,1), round(i/10)]; % save index and time
    end
end
if TSData.L_Strike(1,3) == 0
    TSData.L_Strike(1,:) = [];
    NumStrikes = NumStrikes - 1; 
    NumOffs = NumOffs - 1;
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

%% Calculate Temporal Spatial Paramaters from step timing info
% add later if necessary


end