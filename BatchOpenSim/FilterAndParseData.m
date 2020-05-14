function[Data] = FilterAndParseData(Data, TSData, Settings)

% Filter and Parse data matrix with Time as the first column

%% Define Settings
Settings.PlotFiltering = 'No';

%% Interpolate to standard frequency
% check whether data frequency is an issue
% TimeDiff = Data.data(2:end,1) - Data.data(1:end-1,1);
% figure;
% plot(TimeDiff);

StartTime = Data.data(1,1);
EndTime = Data.data(end,1);
NewSampFreq = 0.001; % define new sampling frequency

NewTime = StartTime:NewSampFreq:EndTime;

% identiy all unique times, aka get rid of doubled times
% necessary for interpolating
[~, UniqueLog, ~] = unique(Data.data(:,1));
Data.data = Data.data(UniqueLog, :); 
% try
    NewSamples = interp1(Data.data(:,1), Data.data(:,2:end), NewTime);
% catch
%     F = griddedInterpolant(Data.data(:,1), Data.data(:,2:end));
% end
% keyboard;
% interp1(Data.data(4:end,1), Data.data(4:end,2:end), NewTime);

Data.Interp = horzcat(NewTime', NewSamples);

if strcmp(Settings.PlotFiltering, 'Yes')
    figure; hold on;
    plot(Data.data(:,1), Data.data(:,2), '.k');
    plot(NewTime, NewSamples(:,1), '*b');
    legend({'Original','Interpolated'});
end

%% Filter signals
close all;
TimeSamp = NewTime(2) - NewTime(1);
SampFreq = round(1/TimeSamp);

Cutoff1 = 20; % cutoff freq in Hz
Cutoff2 = 50;

fs = SampFreq / 2;
Wn1 = Cutoff1 / fs;
Wn2 = Cutoff2 / fs;

[b, a] = butter(4, Wn1, 'low'); % design filter
Data.Fdata = filtfilt(b,a,Data.Interp);
Data.Fdata(:,1) = NewTime'; % replace timing to undo filtering of timeseries
Data.Filter = Cutoff1;

[b, a] = butter(4, Wn2, 'low'); % design filter
Data.Fdata2 = filtfilt(b,a, Data.Interp);
Data.Fdata2(:,1) = NewTime'; % replace timing to undo filtering of timeseries
Data.Filter2 = Cutoff2;

if strcmp(Settings.PlotFiltering, 'Yes')
    % plot filtering
    figure('Position',[100 100 800 800]);
    subplot(511); hold on;
    plot(Data.Interp(:,2), '-k');
    plot(Data.Fdata(:,2), '-r');
    plot(Data.Fdata2(:,2), '-b');
    title('Total');
    legend({'raw', '20Hz Cutoff','50Hz Cutoff'});
    
    subplot(512); hold on;
    plot(Data.Interp(:,23), '-k');
    plot(Data.Fdata(:,23), '-r');
    plot(Data.Fdata2(:,23), '-b');
    title('Glut Max1');
    
    subplot(513); hold on;
    plot(Data.Interp(:,31), '-k');
    plot(Data.Fdata(:,31), '-r');
    plot(Data.Fdata2(:,31), '-b');
    title('Rec Fem');
    
    subplot(514); hold on;
    plot(Data.Interp(:,13), '-k');
    plot(Data.Fdata(:,13), '-r');
    plot(Data.Fdata2(:,13), '-b');
    title('Biceps Fem');
    
    subplot(515); hold on;
    plot(Data.Interp(:,37), '-k');
    plot(Data.Fdata(:,37), '-r');
    plot(Data.Fdata2(:,37), '-b');
    title('Soleus');
end

clearvars d1 fs Wn1 Cutoff1 Wn2 Cutoff2

%% Get Temporal Spatial Data and parse to gait cycles?
% [GRF] = LoadGRF(GRFfile, 0);
% [TSData] = TreadmillTempSpatData(GRF);

%% Identify gait cycles between trial start and end times
Time = Data.Fdata(:,strcmp(Data.colheaders, 'time'));
MetStart = Time(1);
MetEnd = Time(end);

GaitCycles = zeros(1, max([TSData.L_NumStrikes TSData.R_NumStrikes])-1); % initialize # of gait cycles
file = Data.filename; % define filename

if contains(file, 'Left') == 1
    % LEFT
    for i = 1:TSData.L_NumStrikes - 1
        if TSData.L_Strike(i,2) >=  MetStart && TSData.L_Strike(i+1,2) <= MetEnd
            GaitCycles(i) = 1;
            GaitCycleTime = [TSData.L_Strike(i,2) TSData.L_Strike(i+1,2)];
        end
    end
    
elseif contains(file, 'Right') == 1
    % RIGHT
    for i = 1:TSData.R_NumStrikes - 1
        if TSData.R_Strike(i,2) >=  MetStart && TSData.R_Strike(i+1,2) <= MetEnd
            GaitCycles(i) = 1;
            GaitCycleTime = [TSData.R_Strike(i,2) TSData.R_Strike(i+1,2)];
        end
    end
end


%% Parse to % of gait cycle
j = 1;
for i = find(GaitCycles)
    % identify gait cycles
    [~,Start] = min(abs(Time - GaitCycleTime(1)));
    [~,End] = min(abs(Time - GaitCycleTime(2)));
    
    % parse gait cycle to 100 points
    Data.Parsed(:,:,j) = interp1(Time(Start:End), Data.Fdata(Start:End,:), ...
        linspace(Time(Start), Time(End)));
    
    j = j+1;
end

end