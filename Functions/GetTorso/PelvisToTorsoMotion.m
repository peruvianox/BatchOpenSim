% PelvisToTorsoMotion

% use this script to add torso markers (Clav, L.Acr, and R.Acr) to walking
% data without torso markers from a TRC file. 



% Ricky Pimentel        January 2020
% Applied Biomechanics Laboratory - UNC Chapel Hill

%% Load TRC
if ~exist('TRCfile','var')
    [TRCfile, path] = uigetfile('.trc', 'Select TRC files to load and add Torso markers', 'MultiSelect','On');
end

addpath(genpath(path));
addpath(genpath('Functions'));
addpath(genpath('GetTorso'));

% intiialize counters
DynCounter = 1;
StatCounter = 1; 
StaticTrials = 0; 
DynamicTrials = 0; 

if iscell(TRCfile) % if multipe trials selected
    NumTrials = length(TRCfile);
    for i = 1:NumTrials
        if contains(TRCfile{i}, 'Static') || contains(TRCfile{i}, 'static')
            StaticTrials = i;
            StatCounter = StatCounter + 1; 
            data(i).Type = 'Static';
        else
            DynamicTrials(DynCounter) = i;
            DynCounter = DynCounter+1;
            data(i).Type = 'Dynamic';
        end
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
else % if only 1 trial selected 
    if contains(TRCfile, 'Static') || contains(TRCfile, 'static')
        StaticTrials = 1;
        data.Type = 'Static';
    else
        DynamicTrials = 1;
        data.Type = 'Dynamic';
    end
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

NumDynTrials = length(DynamicTrials);
clearvars DynCounter StatCounter

%% Load GRF files 
% we will use force files to identify gait cycles
for i = DynamicTrials
    % translate TRC file names to GRF file names
    filename = strcat([data(i).file(1:end-4),'GRF.mot']);
    [data(i).GRF] = LoadGRF(filename, 0); % load force files
    [data(i).TSData] = TreadmillTempSpatData(data(i).GRF); % get gait cycle times from force files
end

%% Parse Gait cycles
for i = DynamicTrials
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


%% Generalize torso marker locations relative to pelvis
clc;
% load reference data
if exist('TypicalTorso.mat', 'file') == 2 % if typical torso averages exist, use those
    load('TypicalTorso.mat');
else % otherwise select data to load, standard torso values included in Functions -> TorsoNorms folder
    [RefData] = TypcialTorsoMotion;
end

RefData(1) = [];

OrgVector = [1 2 3; 4 5 6; 7 8 9]; % vector for organizing columns
GlobalVector = [1 1 1] ./ norm([1 1 1]); 

%% get pelvis coordinate system from norm data
for i = 1:length(RefData)
    for k = 1:3
        % calculate pelvis center
        RefData(i).PelvisCtr(:,k) = mean(RefData(i).Cycles.PelvisAvg(:, OrgVector(:,k)),2);
        % get pelvis orientation matrix
        P(:,:,k) = RefData(i).Cycles.PelvisAvg(:,OrgVector(k,:));
    end
    % get pelvis orientation
    [UV,V, Coords] = GetOrientation(P);
    
    % apply rotation matrix of pelvis orientation to generalize markers
    % relative to pelvic orientation
    for k = 1:3 % select marker
        % torso marker position vector - relative to pelvis center
        TorsoMarker(:,:,k) = RefData(i).Cycles.TorsoAvg(:,OrgVector(k,:)) - RefData(i).PelvisCtr;
        TM = TorsoMarker(:,:,k);
        for j = 1:length(UV)
            % unit vector of the torso marker
%             TorsoMarkerUV(j,:, k) = TorsoMarker(j,:,k) ./ norm(TorsoMarker(j,:, k));
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
            RelativePositionVector(j,:,k) = TorsoMarker(j,:,k)*Rot(:,:,j);
        end
        % plot arrows to visualize
        %         figure; hold on;
        %         quiver3(zeros(100,1), zeros(100,1), zeros(100,1), UV(:, 1), UV(:, 2), UV(:, 3), 'b');
        %         quiver3(0, 0, 0, 1, 1, 1, 'k', 'LineWidth', 2);
        %         quiver3(zeros(100,1), zeros(100,1), zeros(100,1),...
        %             TorsoMarkerUV(:, 1, k) , TorsoMarkerUV(:, 2, k), TorsoMarkerUV(:, 3, k) , 'r');
        %         xlabel('X'); ylabel('Y'); zlabel('Z');
    end
    % scale position vector by leg length
    RefData(i).NormPV = RelativePositionVector ./  RefData(i).Scale.AvgLegLength;
end
clearvars TorsoMarker TorsoMarkerUV

% Average positions over multiple norm trials
for j = 1:3
    for i = 1:length(RefData)
        Data(:,:,i) = RefData(i).NormPV(:,:,j);
    end
    TorsoPositionVectors(:,:,j) = mean(Data, 3);
    AvgRot = mean(Rot,3);
end
clearvars Data C D NP0

%% plot to check reference data
close all;
figure;
for i = 1:length(RefData)
    subplot(131); hold on;
    Clav = plot(RefData(i).NormPV(:,1,1), 'r');
    LAcr = plot(RefData(i).NormPV(:,1,2), 'g');
    RAcr = plot(RefData(i).NormPV(:,1,3), 'b');
    title('X');
    subplot(132); hold on;
    plot(RefData(i).NormPV(:,2,1), 'r');
    plot(RefData(i).NormPV(:,2,2), 'g');
    plot(RefData(i).NormPV(:,2,3), 'b');
    title('Y');
    subplot(133); hold on;
    plot(RefData(i).NormPV(:,3,1), 'r');
    plot(RefData(i).NormPV(:,3,2), 'g');
    plot(RefData(i).NormPV(:,3,3), 'b');
    title('Z');
    
    legend([Clav, LAcr, RAcr], {'Clav','LAcr','RAcr'});
end
for i = 1:3
    subplot(131); plot(TorsoPositionVectors(:,1,i),'--k', 'LineWidth', 1);
    subplot(132); plot(TorsoPositionVectors(:,2,i),'--k', 'LineWidth', 1);
    subplot(133); plot(TorsoPositionVectors(:,3,i),'--k', 'LineWidth', 1);
end
clearvars Clav LAcr RAcr


%% get scaling for Trials to add torso data to
% assign markers to pull
Markers = {'S2','L.ASIS', 'R.ASIS','L.Knee','R.Knee', 'L.Ankle','R.Ankle'};
% Markers = {'L.PSIS', 'R.PSIS','L.ASIS', 'R.ASIS','L.Knee','R.Knee', 'L.Ankle','R.Ankle'};

for i = 1:NumTrials
    [data(i).MarkerData] = GetMarkerTrajectories(data(i).data, data(i).colheaders, Markers);
    
    if i == StaticTrials
        % use all time points for static trial
        L_Strikes = 1:data(i).data(end,1);
        R_Strikes = 1:data(i).data(end,1);
    else
        % get foot strike times from temporal spatial data
        L_Strikes = data(i).TSData.L_Strike(:,3);
        if L_Strikes(1) == 0
            L_Strikes(1) = [];
        end
        R_Strikes = data(i).TSData.R_Strike(:,3);
        if R_Strikes(1) == 0
            R_Strikes(1) = [];
        end
    end
    
    % Compute average leg length
    [Scale(i).AvgLegLength, Scale(i).Left, Scale(i).Right] = GetLegLength(data(i).MarkerData, Markers, L_Strikes, R_Strikes); 
    
end

%% Apply scales to Vectors
for i = 1:NumTrials
    % scale by leg length
    if i == StaticTrials % use average position vector for static trials
        AvgTorsoPosVector = mean(TorsoPositionVectors, 1); 
        NewTorso(i).Clav = AvgTorsoPosVector(:,:,1) .* Scale(i).AvgLegLength .* ones(data(i).data(end,1),3);
        NewTorso(i).LAcr = AvgTorsoPosVector(:,:,2) .* Scale(i).AvgLegLength .* ones(data(i).data(end,1),3);
        NewTorso(i).RAcr = AvgTorsoPosVector(:,:,3) .* Scale(i).AvgLegLength .* ones(data(i).data(end,1),3);
    else % use parsed position vector for dynamic trial
        NewTorso(i).Clav = TorsoPositionVectors(:,:,1) .* Scale(i).AvgLegLength;
        NewTorso(i).LAcr = TorsoPositionVectors(:,:,2) .* Scale(i).AvgLegLength;
        NewTorso(i).RAcr = TorsoPositionVectors(:,:,3) .* Scale(i).AvgLegLength;
    end
end

%% plot offsets
% figure;
% for i = 1:length(NewTorso)
%     subplot(131); hold on;
%     Ster = plot(NewTorso(i).Ster(:,1), 'r');
%     LAcr = plot(NewTorso(i).LAcr(:,1), 'g');
%     RAcr = plot(NewTorso(i).RAcr(:,1), 'b');
%     title('X');
%     subplot(132); hold on;
%    plot(NewTorso(i).Ster(:,2), 'r');
%     plot(NewTorso(i).LAcr(:,2), 'g');
%     plot(NewTorso(i).RAcr(:,2), 'b');
%     title('Y');
%     subplot(133); hold on;
%     plot(NewTorso(i).Ster(:,3), 'r');
%     plot(NewTorso(i).LAcr(:,3), 'g');
%    plot(NewTorso(i).RAcr(:,3), 'b');
%     title('Z');
%     legend([Ster, LAcr, RAcr], {'Ster','LAcr','RAcr'});
% end
% clearvars Ster LAcr RAcr

%% Generalize Motion According to Gait Cycles
CycLength = 200; % set cycle length for 2 gait cycles
Ind = horzcat(101:149, 50:100); 
for i = 1:NumTrials
    if i == StaticTrials
        
        % extract pelvis marker trajectories
        if sum(strcmp(Markers, 'S2')) == 0
            S2fromPSIS(:,:,1) = data(i).MarkerData(strcmp(Markers, 'L.PSIS')).Trajectories(:, :);
            S2fromPSIS(:,:,2) = data(i).MarkerData(strcmp(Markers, 'R.PSIS')).Trajectories(:, :);
            data(i).AllCycles.Pelvis = [mean(S2fromPSIS, 3),...
                data(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(:, :),...
                data(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(:, :)];
        else
            data(i).AllCycles.Pelvis = [data(i).MarkerData(strcmp(Markers, 'S2')).Trajectories(:, :),...
                data(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(:, :),...
                data(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(:, :)];
        end
        % compute pelvis origin
        for k = 1:3
            data(i).Cycles.PelvisOrigin(:,k,:) = mean(data(i).AllCycles.Pelvis(:,OrgVector(:,k),:),2);
        end
    else
        for j = 1: data(i).TSData.NumGCs-1
            % define gait cycles
            Start = data(i).GCTimes(j,2);
            End = data(i).GCTimes(j+1,4);
            CycDur =  End - Start + 1;
            % extract pelvis marker trajectories
            if sum(strcmp(Markers, 'S2')) == 0
                S2fromPSIS(:,:,1) = data(i).MarkerData(strcmp(Markers, 'L.PSIS')).Trajectories(Start:End, :);
                S2fromPSIS(:,:,2) = data(i).MarkerData(strcmp(Markers, 'R.PSIS')).Trajectories(Start:End, :);
                data(i).AllCycles.Pelvis = [mean(S2fromPSIS, 3),...
                    data(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(Start:End, :),...
                    data(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(Start:End, :)];
            else
                data(i).AllCycles(j).Pelvis = [data(i).MarkerData(strcmp(Markers, 'S2')).Trajectories(Start:End, :),...
                    data(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(Start:End, :),...
                    data(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(Start:End, :)];
            end
            % do resampling
            data(i).Cycles2.Pelvis(:,:,j) = resample(data(i).AllCycles(j).Pelvis, CycLength, CycDur);
            % reduce double gait cycles down to one cycle
            % use the 2nd half of cycle 1 and first half of cycle 2 for continuous measure
            data(i).Cycles.Pelvis(:,:,j) = data(i).Cycles2.Pelvis(Ind,:,j);
            
            % compute pelvis origin
            for k = 1:3
                data(i).Cycles.PelvisOrigin(:,k,j) = mean(data(i).Cycles.Pelvis(:,OrgVector(:,k),j),2);
            end
        end
    end
end
clearvars Start End CycDur i j k

%% combine vectors and average
clc; close all; 
OrgVector = [1 2 3; 4 5 6; 7 8 9];
for i = 1:NumTrials
     % put back into global coordinate system by multiplying by transposed rotation matrix
    if i == StaticTrials
        StaticClav = NewTorso(i).Clav*AvgRot' + data(i).Cycles.PelvisOrigin(:,:);
        StaticLAcr = NewTorso(i).LAcr*AvgRot' + data(i).Cycles.PelvisOrigin(:,:);
        StaticRAcr = NewTorso(i).RAcr*AvgRot' + data(i).Cycles.PelvisOrigin(:,:);
        NewTorso(i).Cycles = [StaticClav StaticLAcr StaticRAcr];
    else
        for j = 1:data(i).TSData.NumGCs-1 % number of gait cycle
            for k = 1:100 % percent of gait cycle
                NewClav(:,:,j) = NewTorso(i).Clav(k,:)*Rot(:,:,k)' + data(i).Cycles.PelvisOrigin(:,:,j);
                NewLAcr(:,:,j) = NewTorso(i).LAcr(k,:)*Rot(:,:,k)' + data(i).Cycles.PelvisOrigin(:,:,j);
                NewRAcr(:,:,j) = NewTorso(i).RAcr(k,:)*Rot(:,:,k)' + data(i).Cycles.PelvisOrigin(:,:,j);
            end
            NewTorso(i).Cycles(:,:,j) = [NewClav(:,:,j) NewLAcr(:,:,j) NewRAcr(:,:,j)];
        end
    end
end

% plot to show relative locations
for i = DynamicTrials
    figure; hold on;
    for j = 1:data(i).TSData.NumGCs-1
        Clav = plot3(NewClav(:,1,j), NewClav(:,2,j), NewClav(:,3,j), 'r');
        LAcr = plot3(NewLAcr(:,1,j), NewLAcr(:,2,j), NewLAcr(:,3,j), 'g');
        RAcr = plot3(NewRAcr(:,1,j), NewRAcr(:,2,j), NewRAcr(:,3,j), 'b');
        Pelvis = plot3(data(i).Cycles.Pelvis(:,[1 4 7],j), data(i).Cycles.Pelvis(:,[2 5 8],j), data(i).Cycles.Pelvis(:,[3 6 9],j), 'k');
        plot3(data(i).Cycles.PelvisOrigin(:,1,j), data(i).Cycles.Pelvis(:,2,j), data(i).Cycles.Pelvis(:,3,j), 'k')
    end
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend({'Clav','LAcr','RAcr', 'Pelvis'}); 
    title(strcat('Global positions for ',data(i).file));
end
clearvars Ster LAcr RAcr i j 

%% resample to time of each gait cycle
for i = 1:NumTrials
    if i == StaticTrials
        data(i).NewData = NewTorso(i).Cycles; % place zeros after array to make it the same length as original
        data(i).Export = horzcat(data(i).data, data(i).NewData);
    else
        for j = 1: data(i).TSData.NumGCs-1
            % define gait cycles
            Start = data(i).GCTimes(j,2);
            End = data(i).GCTimes(j,4);
            CycDur =  End - Start + 1;
            % resample back to original gait cycle times
            NewTorso(i).OrigCycles(j).Torso = resample(NewTorso(i).Cycles(:,:,j), CycDur, 100);
            data(i).NewData(Start:End, :) = NewTorso(i).OrigCycles(j).Torso;
        end
        % filter data to account for movements between cycles
        d = designfilt('lowpassiir','FilterOrder',4,  'HalfPowerFrequency',0.12,'DesignMethod','butter');
        data(i).FiltData = filtfilt(d,data(i).NewData);
        
        % plot 
        figure; subplot(311); hold on; title('X'); 
        plot(data(i).NewData(:,[1 4 7]),'--'); 
        plot(data(i).FiltData(:,[1 4 7]),'-'); 
         subplot(312); hold on; title('Y'); 
        plot(data(i).NewData(:,[2 5 8]),'--'); 
        plot(data(i).FiltData(:,[2 5 8]),'-');
         subplot(313); hold on; title('Z'); 
        plot(data(i).NewData(:,[3 6 9]),'--'); 
        plot(data(i).FiltData(:,[3 6 9]),'-');
        
        data(i).FiltData(End:length(data(i).data), :) = 0; % place zeros after array to make it the same length as original
        TimeTrim = data(i).NewData(:,1) == 0; % identify zeros for trimming to time (outside of left gait cycles)
        data(i).Export= horzcat(data(i).data, data(i).FiltData);
        data(i).Export(TimeTrim, :) = []; % trim new combine data to times between left gait cycles
        data(i).Export(1:25,:) = []; % trim again to make sure end points arent affected by filtering
        data(i).Export(end-25:end,:) = []; 
    clearvars Start End CycDur TimeTrim PreFilt
    end
end

%% plot orientations
% close all;
% Pelvis(:,:,1) = data(i).MarkerData(1).Trajectories;
% Pelvis(:,:,2) = data(i).MarkerData(2).Trajectories;
% Pelvis(:,:,3) = data(i).MarkerData(3).Trajectories;
%
% Torso(:,:,1) = data(i).NewData(:,1:3);
% Torso(:,:,2) = data(i).NewData(:,4:6);
% Torso(:,:,3) = data(i).NewData(:,7:9);
%
% [Euler, Rot, Quat] = SegmentOrientation(Pelvis, Torso, 'Yes');


%% Put New Cycles back into TRC file
clc;
for i = 1:NumTrials
    % Copy file to new name (_AddTorso.trc)
    data(i).NewFileName = strcat(path, data(i).file(1:end-4), '_AddTorso.trc');
    disp(strcat('Writing ', data(i).NewFileName));
    NewFilename = data(i).NewFileName;
    Src = strcat(path, data(i).file);
    copyfile(Src, NewFilename);
    % Correct Headers
    fid=fopen(NewFilename,'r'); %Read TRC files to determine correct headers
    FData = textscan(fid,'%q'); %Looks at headers
    FData = FData{1,1}; %Reorganizes
    FData{4,1}=data(i).file; %Renames file name
    fclose(fid);
    
    % find header locations in file
    LastCategory = find(strcmp(FData,'OrigNumFrames'), 1);
    ColumnHeaderStart = find(strcmp(FData,'Frame#'), 1);
    % get names of column headers
    [HeaderData] = ClarifyHeaders(Src, data(i));
    HeaderMarkers = HeaderData.colheaders;
    HeaderMarkers(2,:) = [];
    for j = 1:length(HeaderMarkers)
        if strcmp(HeaderMarkers{j},' ')
            ToDel(j) = 1;
        elseif isempty(HeaderMarkers{j})
            ToDel(j) = 1;
        else
            ToDel(j) = 0;
        end
    end
    HeaderMarkers(logical(ToDel)) = [];
    HeaderMarkers(1, end+1:end+3) = {'Clav','L.Acromium','R.Acromium'}; % Add labels for new torso markers
    
    % Write data to new TRC file
    fid=fopen(NewFilename,'w');
    
    % Copies over header info
    fprintf(fid,'%s\t', FData{1:5,1}); % Top line 
    fprintf(fid,'\n'); % new line
    fprintf(fid,'%s\t', FData{6:LastCategory,1}); % second line 
    fprintf(fid,'\n'); % new line
    fprintf(fid,'%s\t', FData{LastCategory+1:LastCategory+3,1});
    fprintf(fid,'%s\t', num2str(length(HeaderMarkers)));
    fprintf(fid,'%s\t', FData{LastCategory+5:ColumnHeaderStart-1,1});
    fprintf(fid,'\n'); % new line
    fprintf(fid,'%s\t', FData{ColumnHeaderStart:ColumnHeaderStart+1,1});
    
    % Copies over column headers with two tabs in between each marker
    NCol = length(HeaderMarkers);
    for m=1:NCol
        fprintf(fid,'%s\t',HeaderMarkers{m});
        fprintf(fid,'\t'); % tab
        fprintf(fid,'\t'); % tab
    end
    
    % new line and tab over twice to pass the Frame and Time columns
    fprintf(fid,'\n'); fprintf(fid,'\t'); fprintf(fid,'\t');
    
    % Labels x, y, and z columns
    d=0;
    for m=1:NCol
        d=d+1;
        fprintf(fid,'%s\t', ['X',num2str(d)]);
        fprintf(fid,'%s\t', ['Y',num2str(d)]);
        fprintf(fid,'%s\t', ['Z',num2str(d)]);
    end
    
    % Adds a space between xyz and data.
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
  % Inputs new/old marker x, y, z info
    [NumRow, ~] = size(data(i).Export);   
    for ii=1:NumRow
        fprintf(fid,'%3.4f\t',data(i).Export(ii,:)); fprintf(fid,'\n');
    end
    fclose(fid);
    fclose('all');
    
    % OpenSim won't recognized the file unless it is opened and saved again
    e=actxserver('excel.application');
    eW=e.Workbooks;
    eF=eW.Open(NewFilename); % open OutputTest.xls
    eS=eF.ActiveSheet; eF.Save;
    eF.Close; % close the file
    e.Quit; % close Excel entirely
    % export and clear variables
    disp(strcat(NewFilename, ' Exported!'));
    clearvars fid FData e eW eF eS ii Src NewFilename m j ii fid FData j d ans ToDel
end
close all;
disp('ALL TRIALS CONVERTED');


%% Plot torso and pelvis points for visualization
% X = 3;
% Y = 1;
% Z = 2;
% az = -180;
% el =  6;
%
% figure('Position', [50 50 1200 800]);
% for i = 1:length(Sacrum)
%     plot3([Sacrum(i,X), L.ASIS(i,X), R.ASIS(i,X), Sacrum(i,X)],...
%         [Sacrum(i,Y), L.ASIS(i,Y), R.ASIS(i,Y), Sacrum(i,Y)],...
%         [Sacrum(i,Z), L.ASIS(i,Z), R.ASIS(i,Z), Sacrum(i,Z)], 'b-');
%     hold on;
%     plot3([Sternum(i,X), L.Acromium(i,X), R.Acromium(i,X), Sternum(i,X)],...
%         [Sternum(i,Y), L.Acromium(i,Y), R.Acromium(i,Y), Sternum(i,Y)],...
%         [Sternum(i,Z), L.Acromium(i,Z), R.Acromium(i,Z), Sternum(i,Z)], 'r-');
%     plot3(Sacrum(i,X), Sacrum(i,Y), Sacrum(i,Z), 'b.','MarkerSize',20);
%     plot3(Sternum(i,X), Sternum(i,Y), Sternum(i,Z), 'r.','MarkerSize',20);
%     axis equal
%     title({'Pelvis (blue) and Torso (red) position at frame ' num2str(i)});
%     view(az, el);
%     pause(0.05);
%     hold off
% end


