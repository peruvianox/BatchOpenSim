function [RefData] = TypcialTorsoMotion
% TypicalTorsoMotion

% use walking data from trials with torso data to generalize torso motion
% while walking. This script will generate torso motion norms parsed to the
% gait cycle, so that it can be added to walking data without torso data.

% Ricky Pimentel - December 2019
% Applied Biomechanics Lab, UNC-Chapel Hill

%% Load Example TRC files
addpath(genpath('Functions')); 
[TRCfile, path] = uigetfile('.trc', 'Select TRC files to use as typical motion', 'MultiSelect','On');
addpath(genpath(path));
Plot = 'Yes'; % plotting variable Yes/No

if iscell(TRCfile) % if multipe trials selected
    NumTrials = length(TRCfile);
    for i = 1:NumTrials
        load = importdata(TRCfile{i}, '\t',5);
        data(i).file = TRCfile{i};
        if isnan(load.data(1,1))
            load.data(1,:) = [];
        end
        data(i).data = load.data;
        data(i).textdata = load.textdata;
        [NewData] = ClarifyHeaders(data(i).file, load);
        data(i).colheaders = NewData.colheaders;
        clearvars load NewData
    end
else
    NumTrials = 1;
    load = importdata(TRCfile, '\t',5);
    data.file = TRCfile;
    if isnan(load.data(1,1))
        load.data(1,:) = [];
    end
    data.data = load.data;
    data.textdata = load.textdata;
    [NewData] = ClarifyHeaders(data.file, load);
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
    % get gait cycle times
    [data(i).TSData] = TreadmillTempSpatData(data(i).GRF);
end

%% Parse Gait cycles
for i = 1:NumTrials
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
for i = 1:NumTrials
    [data(i).MarkerData] = GetMarkerTrajectories(data(i).data, data(i).colheaders, Markers);
    [data(i).AvgLegLength] = GetLegLength(data(i).MarkerData, Markers, data(i).TSData.L_Strike(:,3), data(i).TSData.R_Strike(:,3)); 
    %     VisualizeMarkers(data(i), Markers);
    if isempty(data(i).MarkerData(4).Name) % find trials with missing clav marker
        MissingClav(i) = 1;
    end
end
 
%% Replace Any Missing Markers
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
        VisualizeMarkers(data(i), Markers);
    end
end

%% Generalize Motion to 2 adjacent Gait Cycles
CycLength = 200; % set cycle length for 2 gait cycles
Ind = horzcat(101:149, 50:100); 
for i = 1:NumTrials
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

%% plot data to make sure its right
if strcmp(Plot, 'Yes')
    i = 1; % trial data to plot
    Xdata = [1 4 7 1];
    Ydata = [2 5 8 2];
    Zdata = [3 6 9 3];
    figure;
    for j = 1:100 % loop through time?
        hold on;
        plot3(data(i).Cycles.TorsoAvg(j,Xdata), data(i).Cycles.TorsoAvg(j,Ydata), data(i).Cycles.TorsoAvg(j,Zdata), '-k');
        plot3(data(i).Cycles.PelvisAvg(j,Xdata), data(i).Cycles.PelvisAvg(j,Ydata), data(i).Cycles.PelvisAvg(j,Zdata), '-k');
        for k = 1:4
            % torso stds
            [x,y,z] = ellipsoid(data(i).Cycles.TorsoAvg(j,Xdata(k)), data(i).Cycles.TorsoAvg(j,Ydata(k)), data(i).Cycles.TorsoAvg(j,Zdata(k)),...
                data(i).Cycles.TorsoStd(j,Xdata(k)), data(i).Cycles.TorsoStd(j,Ydata(k)), data(i).Cycles.TorsoStd(j,Zdata(k)));
            surf(x,y,z);
            % pelvis stds
            [x,y,z] = ellipsoid(data(i).Cycles.PelvisAvg(j,Xdata(k)), data(i).Cycles.PelvisAvg(j,Ydata(k)), data(i).Cycles.PelvisAvg(j,Zdata(k)),...
                data(i).Cycles.PelvisStd(j,Xdata(k)), data(i).Cycles.PelvisStd(j,Ydata(k)), data(i).Cycles.PelvisStd(j,Zdata(k)));
            surf(x,y,z);
        end
        az = -180;
        el =  6;
        view(az, el);
        axis equal;
        hold off;
        pause(0.1);
        clf;
        axis equal;
    end
end

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
for i = 1:NumTrials
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
for i = 1:NumTrials
    StaticTrials{i} = strcat([data(i).file(1:end-4),'_StaticTrial.trc']);
    if exist(StaticTrials{i}) ~= 2
        StaticTrials = uigetfile('.trc','Select static trials for the loaded trials','Multiselect','On'); 
        break
    end
end
    
% load static data
for i = 1:NumTrials 
    [data(i).StaticData] = LoadStaticTRC(StaticTrials{i});
end

%% subtract offsets from static trial
SubMatrix = [-1, -1, 1].*ones(size(data(i).Angles(:,:)));
LW = 1;
figure; hold on;
for i = 1:NumTrials
    Offset(i).Angles = data(i).Angles - SubMatrix.*data(i).StaticData.Euler;
    Offset(i).StdAnglesP = data(i).StdAnglesP - SubMatrix.*data(i).StaticData.Euler;
    Offset(i).StdAnglesM = data(i).StdAnglesM - SubMatrix.*data(i).StaticData.Euler;
    data(i).Offset = Offset(i); 
    
    % visualize orientation
    X =  plot(Offset(i).Angles(:,1), 'r', 'LineWidth',LW);
    plot(Offset(i).StdAnglesP(:,1),'r--', 'LineWidth',LW/2);
    plot(Offset(i).StdAnglesM(:,1),'r--', 'LineWidth',LW/2);
    Y =  plot(Offset(i).Angles(:,2), 'g', 'LineWidth',LW);
    plot(Offset(i).StdAnglesP(:,2),'g--', 'LineWidth',LW/2);
    plot(Offset(i).StdAnglesM(:,2),'g--', 'LineWidth',LW/2);
    Z = plot(Offset(i).Angles(:,3), 'b', 'LineWidth',LW);
    plot(Offset(i).StdAnglesP(:,3),'b--', 'LineWidth',LW/2);
    plot(Offset(i).StdAnglesM(:,3),'b--', 'LineWidth',LW/2);
end
legend([X, Y, Z], {'X','Y','Z'});
xlabel('% Gait Cycle');
ylabel('Degrees');
title('Torso Motion Average Movements');

clearvars X Y Z LW 

%% Scale Offsets from pelvis and leg dimensions
% Get Markers
Markers = {'L.ASIS', 'R.ASIS','L.Ankle','R.Ankle','S2', 'L.Acr','R.Acr','Sternum'};
for i = 1:NumTrials
    [Scale(i).MarkerData] = GetMarkerTrajectories(data(i).data, data(i).colheaders, Markers);
    
    % strike times
    L_Strikes = data(i).TSData.L_Strike(:,3);
    R_Strikes = data(i).TSData.R_Strike(:,3);
    
    % Compute average leg length at heel strikes
    % estimate using distance between ASIS and Ankle Markers
    Scale(i).LASIS = Scale(i).MarkerData(1).Trajectories(L_Strikes,:);
    Scale(i).LAnk = Scale(i).MarkerData(3).Trajectories(L_Strikes,:);
    Scale(i).L_LegLength = norm(mean(Scale(i).LASIS-Scale(i).LAnk));
    Scale(i).RASIS = Scale(i).MarkerData(2).Trajectories(R_Strikes,:);
    Scale(i).RAnk = Scale(i).MarkerData(4).Trajectories(R_Strikes,:);
    Scale(i).R_LegLength = norm(mean(Scale(i).RASIS-Scale(i).RAnk));
    % average left and right sides
    Scale(i).AvgLegLength = mean([Scale(i).L_LegLength Scale(i).R_LegLength]);
    clearvars LASIS RASIS LAnk RAnk
    
    % Scale pelvis
%     LASIS = Scale(i).MarkerData(1).Trajectories;
%     RASIS = Scale(i).MarkerData(2).Trajectories;
%     S2 = Scale(i).MarkerData(5).Trajectories;
%     InterASIS = mean(nansum(abs(LASIS-RASIS),2));
%     LASIS_S2 = mean(nansum(abs(LASIS-S2),2));
%     RASIS_S2 = mean(nansum(abs(RASIS-S2),2));
%     Scale(i).PelvisCircum = nansum([LASIS_S2 RASIS_S2 InterASIS]);
    
    % Scale location of torso and torso dimensions to pelvis size
%     LAcr = Scale(i).MarkerData(6).Trajectories;
%     RAcr = Scale(i).MarkerData(7).Trajectories;
%     Sternum = Scale(i).MarkerData(8).Trajectories;
%     InterAcr = mean(nansum(abs(LAcr-RAcr),2));
%     LAcr_Str = mean(nansum(abs(LAcr-Sternum),2));
%     RAcr_Str = mean(nansum(abs(RAcr-Sternum),2));
%     Scale(i).Torso = nansum([LAcr_Str RAcr_Str InterAcr]);
%     
%     % Determine overall scaling
%     Scale(i).Scale = Scale(i).Torso / sum([Scale(i).AvgLegLength Scale(i).PelvisCircum]);
    
    data(i).Scale = Scale(i); 
    
    clearvars LASIS RASIS S2 InterASIS LASIS_S2 RASIS_S2 LAcr RAcr Sternum...
        InterAcr LAcr_Str RAcr_Str L_LegLength L_Strikes R_LegLength R_Strikes
end

%% Save Structure of Data
RefData = data; 

end

%% Plot parent orientations to check defined axes
% Viz = 'Yes';
% if strcmp(Viz, 'Yes')
%     figure;
%     Tri = [1 2 3 1];
%     VizFactor = 100;
%     pLW = 3;
%     cLW = 1;
%     for i = 1:length(ParentAvg)
%         clf;
%         hold on;
%
%         % parent
%         plot3([0 ParentAxis1uv(i, 1)*VizFactor],...
%             [0 ParentAxis1uv(i, 2)*VizFactor],...
%             [0 ParentAxis1uv(i, 3)*VizFactor], 'r', 'LineWidth',pLW);
%         plot3([0 ParentAxis2uv(i, 1)*VizFactor],...
%             [0 ParentAxis2uv(i, 2)*VizFactor],...
%             [0 ParentAxis2uv(i, 3)*VizFactor], 'g', 'LineWidth',pLW);
%         plot3([0 ParentAxis3uv(i, 1)*VizFactor],...
%             [0 ParentAxis3uv(i, 2)*VizFactor],...
%             [0 ParentAxis3uv(i, 3)*VizFactor], 'b', 'LineWidth',pLW);
%
%         % child
%         plot3([0 ChildAxis1uv(i, 1)*VizFactor],...
%             [0 ChildAxis1uv(i, 2)*VizFactor],...
%             [0 ChildAxis1uv(i, 3)*VizFactor], 'r', 'LineWidth',cLW);
%         plot3([0 ChildAxis2uv(i, 1)*VizFactor],...
%             [0 ChildAxis2uv(i, 2)*VizFactor],...
%             [0 ChildAxis2uv(i, 3)*VizFactor], 'g', 'LineWidth',cLW);
%         plot3([0 ChildAxis3uv(i, 1)*VizFactor],...
%             [0 ChildAxis3uv(i, 2)*VizFactor],...
%             [0 ChildAxis3uv(i, 3)*VizFactor], 'b', 'LineWidth',cLW);
%
%         % child
%         plot3([0 ChildVectorRot(1,1,i)*VizFactor],...
%             [0 ChildVectorRot(1,2,i)*VizFactor],...
%             [0 ChildVectorRot(1,3,i)*VizFactor], '--r', 'LineWidth',cLW*5);
%         plot3([0 ChildVectorRot(2,1,i)*VizFactor],...
%             [0 ChildVectorRot(2,2,i)*VizFactor],...
%             [0 ChildVectorRot(2,3,i)*VizFactor], '--g', 'LineWidth',cLW*5);
%         plot3([0 ChildVectorRot(3,1,i)*VizFactor],...
%             [0 ChildVectorRot(3,2,i)*VizFactor],...
%             [0 ChildVectorRot(3,3,i)*VizFactor], '--b', 'LineWidth',cLW*5);
%
%         title({'Parent and Child positions at frame ' num2str(i)});
%         axis equal;
%         pause(0.1);
%     end
% end




