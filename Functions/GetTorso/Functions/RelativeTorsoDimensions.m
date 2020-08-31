% RelativeTorsoDimensions

% run this script to create a set of normalized torso vectors to use when
% analyzing or generating torso walking motion. This will compile all of
% the general torso markers (sternum, clavicle, L/R acromiums, C7, and T12)
% so (nearly) any combination of 3 of those markers can be used to track
% torso motion.

%% load static torso trial
[file, path] = uigetfile('.trc', 'Select Static TRC file to load');
addpath(genpath(path));

load = importdata(file, '\t',5);
data.file = file;
if isnan(load.data(1,1))
    load.data(1,:) = [];
end
data.data = load.data;
data.textdata = load.textdata;
data = ClarifyHeaders(data.file, load);
clearvars load NewData

%% get markers
Markers = {'Sternum', 'Clav', 'L.Acr','R.Acr', 'C7', 'T12','S2','L.ASIS', 'R.ASIS', 'L.Knee','R.Knee','L.Ankle','R.Ankle'};
i = 1; 
% get marker trajectories
[data(i).MarkerData] = GetMarkerTrajectories(data(i).data, data(i).colheaders, Markers);
% get average leg length for scaling
data(i).AvgLegLength = GetLegLength(data(i).MarkerData, Markers);

%% define pelvis and torso segments
% as a nx3*markers matrix where n = number of frames
i = 1;
% pelvis
MkrCol = data(i).MarkerData(strcmp(Markers,'S2')).Col;
data(i).Pelvis.S2 = data(i).data(:,MkrCol:MkrCol+2);
MkrCol = data(i).MarkerData(strcmp(Markers,'L.ASIS')).Col;
data(i).Pelvis.LAsis= data(i).data(:,MkrCol:MkrCol+2);
MkrCol = data(i).MarkerData(strcmp(Markers,'R.ASIS')).Col;
data(i).Pelvis.RAsis = data(i).data(:,MkrCol:MkrCol+2);

% torso
MkrCol = data(i).MarkerData(strcmp(Markers,'Sternum')).Col;
data(i).Torso.Sternum = data(i).data(:,MkrCol:MkrCol+2);
MkrCol = data(i).MarkerData(strcmp(Markers,'L.Acr')).Col;
data(i).Torso.LAcr = data(i).data(:,MkrCol:MkrCol+2);
MkrCol = data(i).MarkerData(strcmp(Markers,'R.Acr')).Col;
data(i).Torso.RAcr = data(i).data(:,MkrCol:MkrCol+2);
MkrCol = data(i).MarkerData(strcmp(Markers,'Clav')).Col;
data(i).Torso.Clav = data(i).data(:,MkrCol:MkrCol+2);
MkrCol = data(i).MarkerData(strcmp(Markers,'C7')).Col;
data(i).Torso.C7 = data(i).data(:,MkrCol:MkrCol+2);
MkrCol = data(i).MarkerData(strcmp(Markers,'T12')).Col;
data(i).Torso.T12 = data(i).data(:,MkrCol:MkrCol+2);

%% Clavicle based center
% Define virtual Markers from locations
% create Torso Center as center of Clav and Acromiums
% create 3D matrix
TorsoMatrix(:,:,1) = data(i).Torso.Clav;
TorsoMatrix(:,:,2) = data(i).Torso.LAcr;
TorsoMatrix(:,:,3) = data(i).Torso.RAcr;
% get center
data(i).Torso.Ctr = mean(TorsoMatrix, 3);
% get orientation
[UV, ~, ~] = GetOrientation(TorsoMatrix);

% get torso and pelvis orientation
GlobalVector = [1 1 1] ./ norm([1 1 1]);
Markers = {'Sternum', 'Clav', 'L.Acr','R.Acr', 'C7', 'T12'}; % only torso markers needed

% torso marker position vector - relative to torso center
OrigVector.Clav = data(i).Torso.Clav - data(i).Torso.Ctr;
OrigVector.Sternum = data(i).Torso.Sternum - data(i).Torso.Ctr;
OrigVector.C7 = data(i).Torso.C7 - data(i).Torso.Ctr;
OrigVector.T12= data(i).Torso.T12 - data(i).Torso.Ctr;
OrigVector.LAcr = data(i).Torso.LAcr - data(i).Torso.Ctr;
OrigVector.RAcr = data(i).Torso.RAcr - data(i).Torso.Ctr;

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
    
    % multiply original position vector by rotation matrix
    CPositionVector.Clav(j,:) = OrigVector.Clav(j,:)*Rot(:,:,j);
    CPositionVector.Sternum(j,:) = OrigVector.Sternum(j,:)*Rot(:,:,j);
    CPositionVector.C7(j,:) = OrigVector.C7(j,:)*Rot(:,:,j);
    CPositionVector.T12(j,:) = OrigVector.T12(j,:)*Rot(:,:,j);
    CPositionVector.LAcr(j,:) = OrigVector.LAcr(j,:)*Rot(:,:,j);
    CPositionVector.RAcr(j,:) = OrigVector.RAcr(j,:)*Rot(:,:,j);
end

% plot arrows to visualize
%         figure; hold on;
%         quiver3(zeros(100,1), zeros(100,1), zeros(100,1), UV(:, 1), UV(:, 2), UV(:, 3), 'b');
%         quiver3(0, 0, 0, 1, 1, 1, 'k', 'LineWidth', 2);
%         quiver3(zeros(100,1), zeros(100,1), zeros(100,1),...
%             TorsoMarkerUV(:, 1, k) , TorsoMarkerUV(:, 2, k), TorsoMarkerUV(:, 3, k) , 'r');
%         xlabel('X'); ylabel('Y'); zlabel('Z');

clearvars OrigVector TorsoMatrix UV

%% Sternum based center
% Define virtual Markers from locations
% create Torso Center as center of Sternum and Acromiums
% create 3D matrix
TorsoMatrix(:,:,1) = data(i).Torso.Sternum;
TorsoMatrix(:,:,2) = data(i).Torso.LAcr;
TorsoMatrix(:,:,3) = data(i).Torso.RAcr;
% get center
data(i).Torso.Ctr = mean(TorsoMatrix, 3);
% get orientation
[UV, ~, ~] = GetOrientation(TorsoMatrix);

% get torso and pelvis orientation
GlobalVector = [1 1 1] ./ norm([1 1 1]);
Markers = {'Sternum', 'Clav', 'L.Acr','R.Acr', 'C7', 'T12'}; % only torso markers needed

% torso marker position vector - relative to torso center
OrigVector.Clav = data(i).Torso.Clav - data(i).Torso.Ctr;
OrigVector.Sternum = data(i).Torso.Sternum - data(i).Torso.Ctr;
OrigVector.C7 = data(i).Torso.C7 - data(i).Torso.Ctr;
OrigVector.T12= data(i).Torso.T12 - data(i).Torso.Ctr;
OrigVector.LAcr = data(i).Torso.LAcr - data(i).Torso.Ctr;
OrigVector.RAcr = data(i).Torso.RAcr - data(i).Torso.Ctr;

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
    
    % multiply original position vector by rotation matrix
    SPositionVector.Clav(j,:) = OrigVector.Clav(j,:)*Rot(:,:,j);
    SPositionVector.Sternum(j,:) = OrigVector.Sternum(j,:)*Rot(:,:,j);
    SPositionVector.C7(j,:) = OrigVector.C7(j,:)*Rot(:,:,j);
    SPositionVector.T12(j,:) = OrigVector.T12(j,:)*Rot(:,:,j);
    SPositionVector.LAcr(j,:) = OrigVector.LAcr(j,:)*Rot(:,:,j);
    SPositionVector.RAcr(j,:) = OrigVector.RAcr(j,:)*Rot(:,:,j);
end

%% scale position vector by leg length
% clavicle
ClavRelVector.Clav = nanmean(CPositionVector.Clav ./  data(i).AvgLegLength);
ClavRelVector.Sternum = nanmean(CPositionVector.Sternum ./  data(i).AvgLegLength);
ClavRelVector.C7 = nanmean(CPositionVector.C7 ./  data(i).AvgLegLength);
ClavRelVector.T12 = nanmean(CPositionVector.T12 ./  data(i).AvgLegLength);
ClavRelVector.LAcr = nanmean(CPositionVector.LAcr ./  data(i).AvgLegLength);
ClavRelVector.RAcr = nanmean(CPositionVector.RAcr ./  data(i).AvgLegLength);

% sternum
SternumRelVector.Clav = nanmean(SPositionVector.Clav ./  data(i).AvgLegLength);
SternumRelVector.Sternum = nanmean(SPositionVector.Sternum ./  data(i).AvgLegLength);
SternumRelVector.C7 = nanmean(SPositionVector.C7 ./  data(i).AvgLegLength);
SternumRelVector.T12 = nanmean(SPositionVector.T12 ./  data(i).AvgLegLength);
SternumRelVector.LAcr = nanmean(SPositionVector.LAcr ./  data(i).AvgLegLength);
SternumRelVector.RAcr = nanmean(SPositionVector.RAcr ./  data(i).AvgLegLength);


%% save normalized vectors
save('VirtualTorsoClav.mat','ClavRelVector'); 
save('VirtualTorsoSternum.mat','SternumRelVector'); 

