function [RefData] = GenerateTypcialTorsoMotion
% TypicalTorsoMotion

% use walking data from trials with torso data to generalize torso motion
% while walking. This script will generate torso motion norms parsed to the
% gait cycle, so that it can be added to walking data without torso data.

% Ricky Pimentel - December 2019
% Applied Biomechanics Lab, UNC-Chapel Hill

%% Settings
addpath(genpath('Functions')); 
Plot = 'No'; % plotting variable Yes/No
StaticType = 'Each'; % each walking trial has its own static trial to reference
% StaticType = 'Every'; % "Every" static type will reference all loaded
% walking trials to one static trial - aka multiple conditions from the
% same walker

%% Load Example TRC files
[TRCfile, path] = uigetfile('.trc', 'Select TRC files to use as typical motion', 'MultiSelect','On');
addpath(genpath(path));

if iscell(TRCfile) % if multipe trials selected
    NumTrials = length(TRCfile);
    for i = 1:NumTrials
        Load = importdata(TRCfile{i}, '\t',5);
        data(i).file = TRCfile{i};
        if isnan(Load.data(1,1))
            Load.data(1,:) = [];
        end
        data(i).data = Load.data;
        data(i).textdata = Load.textdata;
        [NewData] = ClarifyHeaders(data(i).file, Load);
        data(i).colheaders = NewData.colheaders;
        clearvars load NewData
    end
else
    NumTrials = 1;
    Load = importdata(TRCfile, '\t',5);
    data.file = TRCfile;
    if isnan(Load.data(1,1))
        Load.data(1,:) = [];
    end
    data.data = Load.data;
    data.textdata = Load.textdata;
    [NewData] = ClarifyHeaders(data.file, Load);
    data.colheaders = NewData.colheaders;
    clearvars load NewData
end

%% Identify Forces files to load
% we will use force files to identify gait cycles
for i = 1:NumTrials
    % translate TRC file names to FORCES file names
    filename = strcat([data(i).file(1:end-4),'GRF.mot']);
    % load force files
    [data(i).GRF] = LoadGRF(filename, 0);
    if contains(data(i).file, 'Static')
        data(i).type = 'static'; 
    else
        data(i).type = 'walking';
    % get gait cycle times
    [data(i).TSData] = TreadmillTempSpatData(data(i).GRF);
    end
end
WalkingTrials = strcmp({data.type},'walking');  
StaticTrials = ~WalkingTrials; 

%% Parse Gait cycles
for i = find(WalkingTrials)
    % normalize gait cycles to left side
    data(i).TSData.NumGCs = min([data(i).TSData.L_NumStrikes data(i).TSData.L_NumOffs])-1;
    data(i).GCHeaders = {'Start Time','Start Frame','End Time','End Frame'};
    for j = 1: data(i).TSData.NumGCs
        data(i).GCTimes(j,1) = data(i).TSData.L_Strike(j, 2); % start time
        data(i).GCTimes(j,2) = data(i).TSData.L_Strike(j, 3); % start frame
        data(i).GCTimes(j,3) = data(i).TSData.L_Strike(j+1, 2); % end time
        data(i).GCTimes(j,4) = data(i).TSData.L_Strike(j+1, 3); % end frame
    end
end

%% Get desired marker trajectories
clc;
MissingClav = zeros(1,NumTrials); 
Markers = {'Sternum', 'L.Acr', 'R.Acr', 'Clav', 'S2', 'L.ASIS', 'R.ASIS', 'L.Knee','L.Ankle','R.Knee','R.Ankle'};
for i = find(WalkingTrials)
    [data(i).MarkerData] = GetMarkerTrajectories(data(i).data, data(i).colheaders, Markers);
    [data(i).AvgLegLength] = GetLegLength(data(i).MarkerData, Markers, data(i).TSData.L_Strike(:,3), data(i).TSData.R_Strike(:,3)); 
    %     VisualizeMarkers(data(i), Markers);
    if isempty(data(i).MarkerData(4).Name) % find trials with missing clav marker
        MissingClav(i) = 1;
    end
end
 
%% Replace Any Missing Markers
Markers = {'Sternum', 'L.Acr', 'R.Acr', 'Clav', 'S2', 'L.ASIS', 'R.ASIS'};
for i = 1:length(MissingClav)
    if MissingClav(i) == 1 % replace clav marker
        load('VirtualTorsoSternum.mat');
        
        % get Sternum-based torso orientation
        TorsoMatrix(:,:,1) = data(i).MarkerData(strcmp(Markers,'Sternum')).Trajectories;
        TorsoMatrix(:,:,2) = data(i).MarkerData(strcmp(Markers,'L.Acr')).Trajectories;
        TorsoMatrix(:,:,3) = data(i).MarkerData(strcmp(Markers,'R.Acr')).Trajectories;
        
        [UV, ~, ~] = GetOrientation(TorsoMatrix);
        GlobalVector = [1 1 1] ./ norm([1 1 1]);
        
        % generalize markers relative to torso orientation
        for j = 1:length(UV)
            % calculate cross and dot products
            C = cross(UV(j,:), GlobalVector);
            D = dot(UV(j,:), GlobalVector);
            NP0 = norm(UV(j,:)) ; % used for scaling
            if ~all(C==0) % check for colinearity
                Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ;
                Rot(:,:,j) = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % calculate rotation matrix
            else
                Rot(:,:,j) = sign(D) * (norm(UV) / NP0) ; % orientation and scaling
            end
            Clav(j,:) = mean(TorsoMatrix(j,:,:),3) + data(i).AvgLegLength .* SternumRelVector.Clav * Rot(:,:,j)';
        end
        data(i).MarkerData(4).Name = 'Clav'; % add clav data to structure
        data(i).MarkerData(4).Trajectories = Clav; 
    end
    % plot markers
    if strcmp(Plot, 'Yes')
        VisualizeMarkers(data(i), Markers, 'No');
    end
end

%% Generalize Motion to 2 adjacent Gait Cycles
CycLength = 200; % set cycle length for 2 gait cycles
Ind = horzcat(101:149, 50:100); 
for i = find(WalkingTrials)
    for j = 1:data(i).TSData.NumGCs-1
        % define gait cycles
        Start = data(i).GCTimes(j,2);
        End = data(i).GCTimes(j+1,4);
        CycDur =  End - Start + 1;
        % extract torso and pelvis marker trajectories
        data(i).AllCycles(j).Torso = [data(i).MarkerData(strcmp(Markers, 'Clav')).Trajectories(Start:End, :),...
            data(i).MarkerData(strcmp(Markers, 'L.Acr')).Trajectories(Start:End, :),...
            data(i).MarkerData(strcmp(Markers, 'R.Acr')).Trajectories(Start:End, :)];
        data(i).AllCycles(j).Pelvis = [data(i).MarkerData(strcmp(Markers, 'S2')).Trajectories(Start:End, :),...
            data(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(Start:End, :),...
            data(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(Start:End, :)];
        % do resampling
        data(i).Cycles2.Torso(:,:,j) = resample(data(i).AllCycles(j).Torso, CycLength, CycDur);
        data(i).Cycles2.Pelvis(:,:,j) = resample(data(i).AllCycles(j).Pelvis, CycLength, CycDur);
        
        % reduce double gait cycles down to one cycle
        % use the 2nd half of cycle 1 and first half of cycle 2 for continuous measure
        data(i).Cycles.Torso(:,:,j) =  data(i).Cycles2.Torso(Ind,:,j);
        data(i).Cycles.Pelvis(:,:,j) = data(i).Cycles2.Pelvis(Ind,:,j);
    end
        
    % mean and st dev of torso and pelvis data
    data(i).Cycles.TorsoAvg = mean(data(i).Cycles.Torso, 3);
    data(i).Cycles.TorsoStd = std(data(i).Cycles.Torso, 0, 3);
    data(i).Cycles.PelvisAvg = mean(data(i).Cycles.Pelvis, 3);
    data(i).Cycles.PelvisStd = std(data(i).Cycles.Pelvis, 0, 3);  
end
clearvars Start End CycDur i j k Ind

%% Define segment orientation
NewCycLength = 100; 
ChildAvg = zeros(NewCycLength, 3, 3);
ParentAvg = zeros(NewCycLength, 3, 3);
ChildStd = zeros(NewCycLength, 3, 3);
ParentStd = zeros(NewCycLength, 3, 3);
% create vector for organization
OrgVector = [1 2 3; 4 5 6; 7 8 9];
% Rows = time
% Columns = X,Y,Z coordinates
% Sheet = Markers
for i = find(WalkingTrials)
    % define child and parent segments in a nx3x3 matrix according to the
    % above definitions
    for j = 1:3
        ChildAvg(:,:,j) = data(i).Cycles.TorsoAvg(:,OrgVector(j,:));
        ParentAvg(:,:,j) = data(i).Cycles.PelvisAvg(:,OrgVector(j,:));
        
        ChildStd(:,:,j) = data(i).Cycles.TorsoStd(:,OrgVector(j,:));
        ParentStd(:,:,j) = data(i).Cycles.PelvisStd(:,OrgVector(j,:));
    end
    
    for j = 1:length(ParentAvg)
        % Average Curves
        P = ParentAvg(j,:,:);
        C = ChildAvg(j,:,:);
        [Orientation(j).Euler, Orientation(j).Rot, Orientation(j).Quat] = SegmentOrientation(P, C, 'No');
        data(i).Angles(j,:) = Orientation(j).Euler;
        
        % Std Curves Plus
        P = ParentAvg(j,:,:) + ParentStd(j,:,:);
        C = ChildAvg(j,:,:) + ChildStd(j,:,:);
        [StdOrientation(j).Euler, StdOrientation(j).Rot, StdOrientation(j).Quat] = SegmentOrientation(P, C, 'No');
        data(i).StdAnglesP(j,:) = StdOrientation(j).Euler;
        
        % Std Curves Minus
        P = ParentAvg(j,:,:) - ParentStd(j,:,:);
        C = ChildAvg(j,:,:) - ChildStd(j,:,:);
        [StdOrientation(j).Euler, StdOrientation(j).Rot, StdOrientation(j).Quat] = SegmentOrientation(P, C, 'No');
        data(i).StdAnglesM(j,:) = StdOrientation(j).Euler;
    end
end

clearvars P C i j ParentAvg ParentStd ChildAvg ChildStd

%% Get initial offsets from static trial
% StaticTrials = {'static_1.trc','Static_1_OpenSim.trc', 'Static_1_OpenSim.trc'};
% StaticTrials = uigetfile('.trc','Select static trials for the loaded trials','Multiselect','On'); 
% check for static trials in directory
% if contains(data(i).file,'Static')

%     StaticTrial = data(i).file;
if strcmp(StaticType, 'Every')
    StaticData = LoadStaticTRC(data(StaticTrials).file, 'No');
    
    for i = find(WalkingTrials)
        data(i).StaticData = StaticData;
    end
    
elseif strcmp(StaticType, 'Each')
    % Call for various trials from different subjects with their own static
    % trial
    for i = find(WalkingTrials)
        % define static trial name for walking trial ___
        data(i).StaticTrial = strcat([data(i).file(1:end-4),'_StaticTrial.trc']);
        if exist( data(i).StaticTrial ) ~= 2
            StaticTrials = uigetfile('.trc','Select static trials for the loaded trials','Multiselect','On');
        end
        % load static trial TRC data
          [data(i).StaticData] = LoadStaticTRC(data(i).StaticTrial, 'No');
    end
%     % load static data
%     for i = 1:NumTrials
%         if iscell(StaticTrials)
%           
%         else
%             [data(i).StaticData] = LoadStaticTRC(StaticTrials);
%         end
%     end
end

%% subtract offsets from static trial
SubMatrix = [-1, -1, 1] .* ones(size(data(1).Angles(:,:)));
for i = find(WalkingTrials)
    Math = SubMatrix .* data(i).StaticData.Euler;
    Offset.Angles = data(i).Angles - Math;
    Offset.StdAnglesP = data(i).StdAnglesP - Math;
    Offset.StdAnglesM = data(i).StdAnglesM - Math;
    data(i).Offset = Offset;
    clearvars Offset Math
end

%% Plot Generalized Torso Norms
% one line for each norm trial
close all; clc; 
Plot = 'Yes';
if strcmp(Plot, 'Yes')
    LW = 1;
    figure; hold on; grid on; 
    for i = find(WalkingTrials)
        % visualize orientation
        X =  plot(data(i).Offset.Angles(:,1), 'r', 'LineWidth',LW);
        plot(data(i).Offset.StdAnglesP(:,1),'r--', 'LineWidth',LW/2);
        plot(data(i).Offset.StdAnglesM(:,1),'r--', 'LineWidth',LW/2);
        Y =  plot(data(i).Offset.Angles(:,2), 'g', 'LineWidth',LW);
        plot(data(i).Offset.StdAnglesP(:,2),'g--', 'LineWidth',LW/2);
        plot(data(i).Offset.StdAnglesM(:,2),'g--', 'LineWidth',LW/2);
        Z = plot(data(i).Offset.Angles(:,3), 'b', 'LineWidth',LW);
        plot(data(i).Offset.StdAnglesP(:,3),'b--', 'LineWidth',LW/2);
        plot(data(i).Offset.StdAnglesM(:,3),'b--', 'LineWidth',LW/2);
    end
    legend([X, Y, Z], {'Sagittal','Frontal','Transverse'});
    xlabel('% Gait Cycle');
    ylabel('Degrees');
    title('Torso Motion Relative to Pelvis');
    
    
    clearvars X Y Z LW
end

%% Scale Offsets from pelvis and leg dimensions
% Get Markers
Markers = {'L.ASIS', 'R.ASIS','L.Knee','R.Knee','L.Ankle','R.Ankle'};
for i = find(WalkingTrials)
    % get Marker data
    [Scale.MarkerData] = GetMarkerTrajectories(data(i).data, data(i).colheaders, Markers);
    % strike times
    L_Strikes = data(i).TSData.L_Strike(:,3);
    R_Strikes = data(i).TSData.R_Strike(:,3);
    % get leg length
%     [Scale.AvgLegLength, Scale.Left, Scale.Right] = GetLegLength(Scale.MarkerData, Markers, L_Strikes, R_Strikes);
    % save scaling
    data(i).Scale = Scale; 
    clearvars Scale L_Strikes R_Strikes
end


%% Save Structure of Data
clc; 
RefData = data; 
RefData(StaticTrials) = []; 
SavedRefData = 'TypicalTorsoNew.mat';
save(SavedRefData, 'RefData'); 
TrialsWritten = {data.file};
disp('The following trials: ')
disp(TrialsWritten');
disp(strcat('Written to: ', SavedRefData)); 

 end


