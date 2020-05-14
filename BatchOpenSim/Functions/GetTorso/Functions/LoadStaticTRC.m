function[StaticData] = LoadStaticTRC(TRCfile, ToPlot)

% Load marker data (.trc format) from static trial

% INPUT
% filename
% segments

% OUTPUT
% StaticData   Structure containing


%% input data
if exist('TRCfile', 'var') == 0 % manual input if no file defined to load
    TRCfile = uigetfile('.trc', 'Select Static TRC file to load');
end
% load data
StaticData = importdata(TRCfile, '\t',5);
StaticData.file = TRCfile;
% get column headers
StaticData = ClarifyHeaders(StaticData.file, StaticData);

if isnan(StaticData.data(1,1)) % clear first row if NaN
    StaticData.data(1,:) = [];
end

%% Get initial segments
Markers = 'All';
[MarkerData] = GetMarkerTrajectories(StaticData.data, StaticData.colheaders, Markers);
Markers = {MarkerData.Name};
AvgLegLength = GetLegLength(MarkerData, Markers); 

%% Extract markers
Present = vertcat({MarkerData([1:end]).Name});

% identify markers that are present
% Torso
if sum(strcmp(Present,'R.Acr')) > 0
    Torso.RAcr = 1;
    Torso.RAcr_ind = find(strcmp(Present,'R.Acr'));
else
    Torso.RAcr = 0;
end
if sum(strcmp(Present,'L.Acr')) > 0
    Torso.LAcr = 1;
    Torso.LAcr_ind = find(strcmp(Present,'L.Acr'));
else
    Torso.LAcr = 0;
end
if sum(strcmp(Present,'Sternum')) > 0
    Torso.Ster = 1;
    Torso.Ster_ind = find(strcmp(Present,'Sternum'));
else
    Torso.Ster = 0;
end
if sum(strcmp(Present,'Clav')) > 0
    Torso.Clav = 1;
    Torso.Clav_ind = find(strcmp(Present,'Clav'));
else
    % estimate clav marker location
    load('VirtualTorsoSternum.mat');
    % get Sternum-based torso orientation
    TorsoMatrix(:,:,1) = MarkerData(strcmp(Markers,'Sternum')).Trajectories;
    TorsoMatrix(:,:,2) = MarkerData(strcmp(Markers,'L.Acr')).Trajectories;
    TorsoMatrix(:,:,3) = MarkerData(strcmp(Markers,'R.Acr')).Trajectories;
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
        Clav(j,:) = mean(TorsoMatrix(j,:,:),3) + AvgLegLength .* SternumRelVector.Clav * Rot(:,:,j)';
    end
    % save clav in structure
    NewLine = length(MarkerData) + 1;
    Torso.Clav_ind = NewLine;
    Torso.Clav = 1; 
    MarkerData(NewLine).Name = 'Clav'; % add clav data to structure
    MarkerData(NewLine).Trajectories = Clav;
end
if sum(strcmp(Present,'C7')) > 0
    Torso.C7 = 1;
    Torso.C7_ind = find(strcmp(Present,'C7'));
else
    Torso.C7 = 0;
end
if sum(strcmp(Present,'T12')) > 0
    Torso.T12 = 1;
    Torso.T12_ind = find(strcmp(Present,'T12'));
else
    Torso.T12 = 0;
end


% Pelvis
if sum(strcmp(Present,'R.ASIS')) > 0
    Pelvis.RASIS = 1;
    Pelvis.RASIS_ind = find(strcmp(Present,'R.ASIS'));
else
    Pelvis.RSIS = 0;
end
if sum(strcmp(Present,'L.ASIS')) > 0
    Pelvis.LASIS = 1;
    Pelvis.LASIS_ind = find(strcmp(Present,'L.ASIS'));
else
    Pelvis.LSIS = 0;
end
if sum(strcmp(Present,'S2')) > 0
    Pelvis.Sacr = 1;
    Pelvis.Sacr_ind = find(strcmp(Present,'S2'));
else
    Pelvis.Sacr= 0;
end

%% Get initial static position
% set time to average marker positions
[L,~] = size(StaticData.data);
if L > 10
    Time = 1:10; % Average over first 10 frames
else
    Time = 1:L;
end

%% Average Markers Positions
for i = 1:length(MarkerData)
    MarkerData(i).AvgLoc = mean(MarkerData(i).Trajectories(Time, :));
end

%% Get static orienatation between pelvis and torso segments
% define pelvis as parent and torso as child segments
Parent(:,:,1) = [MarkerData(Pelvis.Sacr_ind).AvgLoc];
Parent(:,:,2) = [MarkerData(Pelvis.LASIS_ind).AvgLoc];
Parent(:,:,3) = [MarkerData(Pelvis.RASIS_ind).AvgLoc];

Child(:,:,1) = [MarkerData(Torso.Clav_ind).AvgLoc];
Child(:,:,2) = [MarkerData(Torso.LAcr_ind).AvgLoc];
Child(:,:,3) = [MarkerData(Torso.RAcr_ind).AvgLoc];

% segment orientation
[StaticData.Euler, StaticData.Rot, StaticData.Quat] = SegmentOrientation(Parent, Child, ToPlot);

end
