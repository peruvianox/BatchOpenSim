function [Trials] = AddTorso (Trials, DynIndex, Settings)

% use this script to add torso markers (Clav, L.Acr, and R.Acr) to walking
% data without torso markers from a TRC file.


% Ricky Pimentel        January 2020
% Applied Biomechanics Laboratory - UNC Chapel Hill

%% Import TRC Data
for i = 1:length(DynIndex)
    Load = importdata(Trials(i).files.OpenSimTRC, '\t',5);
    if isnan(Load.data(1,1))
        Load.data(1,:) = [];
    end
    Trials(i).data = Load.data;
    Trials(i).textdata = Load.textdata;
    [NewData] = ClarifyHeaders(Trials(i).files.OpenSimTRC, Load);
    Trials(i).colheaders = NewData.colheaders;
    TimeCol = 2; % define time column
    Time = Trials(i).data(:, TimeCol);
    Trials(i).TSData.DigitalSampFreq = 1 / (Time(2) - Time(1)); % get digital sampling freq
    if DynIndex(i) == 1
        [~,n] = size(Trials(i).TSData.L_Strike); % get size of gait event matrix
        % 3rd column should contain gait event timing in digital samp freq
        if n == 2 % if 3rd column of gait events not present, calculate using digtal and
            Trials(i).TSData.L_Strike(:,3) = Trials(i).TSData.L_Strike(:,1) * ...
                Trials(i).TSData.DigitalSampFreq / Trials(i).TSData.AnalogSampFreq;
            Trials(i).TSData.R_Strike(:,3) = Trials(i).TSData.R_Strike(:,1) * ...
                Trials(i).TSData.DigitalSampFreq / Trials(i).TSData.AnalogSampFreq;
        end
    end
    clearvars Load NewData TimeCol n
end

%% Load GRF and temporal spatial data if not already present
% we will use force files to identify gait cycles
for i = find(DynIndex)
    if isfield(Trials, 'GRF')
        if isempty(Trials(i).GRF)
            [Trials(i).GRF] = LoadGRF(Trials(i).files.OpenSimGRF, 0); % load force files
        end
    end
    if isfield(Trials, 'TSData')
        if isempty(Trials(i).TSData)
            [Trials(i).TSData] = TreadmillTempSpatData(Trials(i).GRF); % get gait cycle times from force files
        end
    end
end

%% Parse Gait cycles
for i = find(DynIndex)
    % normalize gait cycles to left side
    Trials(i).TSData.NumGCs = min([Trials(i).TSData.L_NumStrikes Trials(i).TSData.L_NumOffs])-1;
    Trials(i).GCHeaders = {'Start Time','Start Frame','End Time','End Frame'};
    
    [~,n] = size(Trials(i).TSData.L_Strike);
    if n < 3
        Trials(i).Times.TRCSampleHz = 1 / (Trials(i).Times.TRC(2) - Trials(i).Times.TRC(1));
        Trials(i).Times.GRFSampleHz = 1 / (Trials(i).GRF.AllData(2,1) - Trials(i).GRF.AllData(1,1));
        Trials(i).TSData.L_Strike(:,3) = round( Trials(i).TSData.L_Strike(:,1) * Trials(i).Times.TRCSampleHz / Trials(i).Times.GRFSampleHz);
    end
    
    for j = 1: Trials(i).TSData.NumGCs
        Trials(i).GCTimes(j,1) = Trials(i).TSData.L_Strike(j, 2); % start time
        Trials(i).GCTimes(j,2) = Trials(i).TSData.L_Strike(j, 3); % start frame
        Trials(i).GCTimes(j,3) = Trials(i).TSData.L_Strike(j+1, 2); % end time
        Trials(i).GCTimes(j,4) = Trials(i).TSData.L_Strike(j+1, 3); % end frame
    end
end

%% Generalize torso marker locations relative to pelvis
% load reference data
if Settings.Billy == 1
    %       load('TypicalTorsoVariousSpeeds.mat');
    load('TypicalTorsoNew.mat', 'RefData');
    
else
    if exist('TypicalTorso.mat', 'file') == 2 % if typical torso averages exist, use those
       load('TypicalTorso3.mat', 'RefData');
    else % otherwise select data to load, standard torso values included in Functions -> TorsoNorms folder
        [RefData] = TypcialTorsoMotion;
    end
end

% Settings.Billy = 0; 

OrgVector = [1 2 3; 4 5 6; 7 8 9]; % vector for organizing columns
GlobalVector = [1 1 1] ./ norm([1 1 1]);

%% get pelvis coordinate system from norm data
if Settings.Billy == 1 % determine speeds of reference torso data
    for i = 1:length(RefData)
        Str = strsplit(RefData(i).file, '_');
        RefData(i).Speed = Str{1};
    end
end
% [m,n] = size(RefData(i).Cycles.PelvisAvg(:,OrgVector(:,:))); 
% P = zeros(m,n,3); 
for i = 1:length(RefData)
    for k = 1:3
        % calculate pelvis center
        RefData(i).PelvisCtr(:,k) = mean(RefData(i).Cycles.PelvisAvg(:, OrgVector(:,k)),2);
        % get pelvis orientation matrix
        P(:,:,k) = RefData(i).Cycles.PelvisAvg(:,OrgVector(k,:));
    end
    % get pelvis orientation
    [UV, ~, ~] = GetOrientation(P);
    
    % apply rotation matrix of pelvis orientation to generalize markers
    % relative to pelvic orientation
    
%     Rot = zeros(3, 3, length(UV)); % initialize variables
%     TorsoMarker = zeros(length(UV), 3, 3); 
    RelativePositionVector = zeros(length(UV), 3, 3); 
    for k = 1:3 % select marker
        % torso marker position vector - relative to pelvis center
        TorsoMarker(:,:,k) = RefData(i).Cycles.TorsoAvg(:,OrgVector(k,:)) - RefData(i).PelvisCtr;
%         TM = TorsoMarker(:,:,k);
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
    RefData(i).NormPV = RelativePositionVector ./  RefData(i).AvgLegLength;
end
clearvars TorsoMarker TorsoMarkerUV

%% Average positions over multiple norm trials
% [m,n,p] = size(RefData(1).NormPV); % initialize variables
% Data = zeros(m,n,p,length(RefData)); 

for i = 1:length(RefData)
    Data(:,:,:,i) = RefData(i).NormPV(:,:,:);
end
TorsoPositionVectors = mean(Data, 4);
clearvars Data C D NP0
AvgRot = mean(Rot,3); % average rotation matrix

%% plot to check reference data
if strcmp(Settings.PlotTorso, 'Yes')
    close all;
    LW = 2;
    
    figure('Position', [100 100 800 600]);
    subplot(131); hold on;
    for i = 1:length(RefData)
        Clav = plot(RefData(i).NormPV(:,1,1), '--r');
        LAcr = plot(RefData(i).NormPV(:,1,2), '--g');
        RAcr = plot(RefData(i).NormPV(:,1,3), '--b');
    end
    ClavAvg = plot(TorsoPositionVectors(:,1,1),'r', 'LineWidth', LW);
    LAcrAvg = plot(TorsoPositionVectors(:,1,2),'g', 'LineWidth', LW);
    RAcrAvg = plot(TorsoPositionVectors(:,1,3),'b', 'LineWidth', LW);
    title('X');
    xlabel('% of Gait Cycle');
    ylabel('Normalized Position Vector');
    legend([Clav, LAcr, RAcr, ClavAvg, LAcrAvg, RAcrAvg], ...
        {'Clav','LAcr','RAcr', 'Clav Avg','LAcr Avg','RAcr Avg'});
    
    subplot(132); hold on;
    for i = 1:length(RefData)
        plot(RefData(i).NormPV(:,2,1), '--r');
        plot(RefData(i).NormPV(:,2,2), '--g');
        plot(RefData(i).NormPV(:,2,3), '--b');
    end
    plot(TorsoPositionVectors(:,2,1),'r', 'LineWidth', LW);
    plot(TorsoPositionVectors(:,2,2),'g', 'LineWidth', LW);
    plot(TorsoPositionVectors(:,2,3),'b', 'LineWidth', LW);
    title('Y');
    xlabel('% of Gait Cycle');
    ylabel('Normalized Position Vector');
    
    subplot(133); hold on;
    for i = 1:length(RefData)
        plot(RefData(i).NormPV(:,3,1), '--r');
        plot(RefData(i).NormPV(:,3,2), '--g');
        plot(RefData(i).NormPV(:,3,3), '--b');
    end
    plot(TorsoPositionVectors(:,3,1),'r', 'LineWidth', LW);
    plot(TorsoPositionVectors(:,3,2),'g', 'LineWidth', LW);
    plot(TorsoPositionVectors(:,3,3),'b', 'LineWidth', LW);
    title('Z');
    xlabel('% of Gait Cycle');
    ylabel('Normalized Position Vector');
    
    clearvars Clav LAcr RAcr ClavAvg LAcrAvg RAcrAvg
end

%% get scaling for Trials to add torso data to
% assign markers to pull
Markers = {'S2','L.ASIS', 'R.ASIS','L.Knee','R.Knee', 'L.Ankle','R.Ankle', 'L_HJC','R_HJC',...
    'L.MKnee', 'R.MKnee','L.MAnkle','R.MAnkle','L.Heel','R.Heel','L.MT5','R.MT5','L.MT1','R.MT1'};
% Markers = {'L.PSIS', 'R.PSIS','L.ASIS', 'R.ASIS','L.Knee','R.Knee', 'L.Ankle','R.Ankle'};
Scale(length(DynIndex)).AvgLegLength = []; % initialize

for i = 1:length(DynIndex)
    [Trials(i).MarkerData] = GetMarkerTrajectories(Trials(i).data, Trials(i).colheaders, Markers);
    
    if i == find(DynIndex) % walking trials
        
        % get foot strike times from temporal spatial data
        L_Strikes = Trials(i).TSData.L_Strike(:,3);
        if L_Strikes(1) == 0
            L_Strikes(1) = [];
        end
        R_Strikes = Trials(i).TSData.R_Strike(:,3);
        if R_Strikes(1) == 0
            R_Strikes(1) = [];
        end
        
    else % static trials
        % use all time points for static trial
        L_Strikes = 1:Trials(i).data(end,1);
        R_Strikes = 1:Trials(i).data(end,1);
        if isnan(L_Strikes)
            L_Strikes = 1:Trials(i).data(end,2);
        end
        if isnan(R_Strikes)
            R_Strikes = 1:Trials(i).data(end,2);
        end
    end
    
    % make sure no strikes occur after length of trial
    Length = size(Trials(i).data, 1);
    L_Strikes(L_Strikes > Length) = [];
    R_Strikes(R_Strikes > Length) = [];
    
    % Compute average leg length
    [Scale(i).AvgLegLength] = GetLegLength(Trials(i).MarkerData, Markers, L_Strikes, R_Strikes);
    
end

%% Apply scales to Vectors
NewTorso(length(DynIndex)).Clav = []; 
NewTorso(length(DynIndex)).LAcr = []; 
NewTorso(length(DynIndex)).RAcr = []; 

for i = 1:length(DynIndex)
    % scale by leg length
    if i == find(~DynIndex) % use average position vector for static trials
        AvgTorsoPosVector = mean(TorsoPositionVectors, 1);
        NewTorso(i).Clav = AvgTorsoPosVector(:,:,1) .* Scale(i).AvgLegLength .* ones(length(Trials(i).data),3);
        NewTorso(i).LAcr = AvgTorsoPosVector(:,:,2) .* Scale(i).AvgLegLength .* ones(length(Trials(i).data),3);
        NewTorso(i).RAcr = AvgTorsoPosVector(:,:,3) .* Scale(i).AvgLegLength .* ones(length(Trials(i).data),3);
    else % use parsed position vector for dynamic trial
        if Settings.Billy == 1 % separate torso position vector by speed
            Str = strsplit(Trials(i).name, '_');
            Match = find(contains({RefData.Speed}, Str{1})); % match Trial to speed of Reference data
            NewTorso(i).Clav = RefData(Match).NormPV(:,:,1) .* Scale(i).AvgLegLength;
            NewTorso(i).LAcr = RefData(Match).NormPV(:,:,2) .* Scale(i).AvgLegLength;
            NewTorso(i).RAcr = RefData(Match).NormPV(:,:,3) .* Scale(i).AvgLegLength;
        else
            NewTorso(i).Clav = TorsoPositionVectors(:,:,1) .* Scale(i).AvgLegLength;
            NewTorso(i).LAcr = TorsoPositionVectors(:,:,2) .* Scale(i).AvgLegLength;
            NewTorso(i).RAcr = TorsoPositionVectors(:,:,3) .* Scale(i).AvgLegLength;
        end
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
for i = 1:length(DynIndex)
    if i == find(~DynIndex)
        
        % extract pelvis marker trajectories
        if sum(strcmp(Markers, 'S2')) == 0
            S2fromPSIS(:,:,1) = Trials(i).MarkerData(strcmp(Markers, 'L.PSIS')).Trajectories(:, :);
            S2fromPSIS(:,:,2) = Trials(i).MarkerData(strcmp(Markers, 'R.PSIS')).Trajectories(:, :);
            Trials(i).AllCycles.Pelvis = [mean(S2fromPSIS, 3),...
                Trials(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(:, :),...
                Trials(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(:, :)];
        else
            Trials(i).AllCycles.Pelvis = [Trials(i).MarkerData(strcmp(Markers, 'S2')).Trajectories(:, :),...
                Trials(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(:, :),...
                Trials(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(:, :)];
        end
        % compute pelvis origin
        for k = 1:3
            Trials(i).Cycles.PelvisOrigin(:,k,:) = mean(Trials(i).AllCycles.Pelvis(:,OrgVector(:,k),:),2);
        end
        
    else
        
        for j = 1: Trials(i).TSData.NumGCs-1
            % define gait cycles
            Start = round(Trials(i).GCTimes(j,2));
            if Start < 1
                Start = 1;
            end
            End = round(Trials(i).GCTimes(j+1,4));
            CycDur =  End - Start + 1; % define gait cycle duration
            
            % make sure Start and End times are within length of TRC data
            % gait cycles are pulled from GRF data potentially leading to
            % discrepancies
            TrialLength = size(Trials(i).Times.TRC, 1);
            if Start > TrialLength || End > TrialLength % if time points are longer than TRC data
                Trials(i).TSData.NumGCs = j; % log # of gait cycles as current loop iteration
                break % exit "for" loop through gait cycles and continue to next step
            end
            
            % extract pelvis marker trajectories
            if sum(strcmp(Markers, 'S2')) == 0
                S2fromPSIS(:,:,1) = Trials(i).MarkerData(strcmp(Markers, 'L.PSIS')).Trajectories(Start:End, :);
                S2fromPSIS(:,:,2) = Trials(i).MarkerData(strcmp(Markers, 'R.PSIS')).Trajectories(Start:End, :);
                Trials(i).AllCycles.Pelvis = [mean(S2fromPSIS, 3),...
                    Trials(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(Start:End, :),...
                    Trials(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(Start:End, :)];
            else
                Trials(i).AllCycles(j).Pelvis = [Trials(i).MarkerData(strcmp(Markers, 'S2')).Trajectories(Start:End, :),...
                    Trials(i).MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(Start:End, :),...
                    Trials(i).MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(Start:End, :)];
            end
            % do resampling
            Trials(i).Cycles2.Pelvis(:,:,j) = resample(Trials(i).AllCycles(j).Pelvis, CycLength, CycDur);
            % reduce double gait cycles down to one cycle
            % use the 2nd half of cycle 1 and first half of cycle 2 for continuous measure
            Trials(i).Cycles.Pelvis(:,:,j) = Trials(i).Cycles2.Pelvis(Ind,:,j);
            
            % compute pelvis origin
            for k = 1:3
                Trials(i).Cycles.PelvisOrigin(:,k,j) = mean(Trials(i).Cycles.Pelvis(:,OrgVector(:,k),j),2);
            end
        end
        
    end
end
clearvars Start End CycDur i j k

%% combine vectors and average

for i = 1:length(DynIndex)
    % put back into global coordinate system by multiplying by transposed rotation matrix
    if i == find(~DynIndex)
        StaticClav = NewTorso(i).Clav*AvgRot' + Trials(i).Cycles.PelvisOrigin(:,:);
        StaticLAcr = NewTorso(i).LAcr*AvgRot' + Trials(i).Cycles.PelvisOrigin(:,:);
        StaticRAcr = NewTorso(i).RAcr*AvgRot' + Trials(i).Cycles.PelvisOrigin(:,:);
        NewTorso(i).Cycles = [StaticClav StaticLAcr StaticRAcr];
    else
        [m,n] = size(Trials(i).Cycles.PelvisOrigin(:,:,1)); % initialize variables
        NewClav = zeros(m,n,Trials(i).TSData.NumGCs-1); 
        NewLAcr = zeros(m,n,Trials(i).TSData.NumGCs-1); 
        NewRAcr = zeros(m,n,Trials(i).TSData.NumGCs-1); 
        
        for j = 1:Trials(i).TSData.NumGCs-1 % number of gait cycle
            for k = 1:100 % percent of gait cycle
                NewClav(:,:,j) = NewTorso(i).Clav(k,:)*Rot(:,:,k)' + Trials(i).Cycles.PelvisOrigin(:,:,j);
                NewLAcr(:,:,j) = NewTorso(i).LAcr(k,:)*Rot(:,:,k)' + Trials(i).Cycles.PelvisOrigin(:,:,j);
                NewRAcr(:,:,j) = NewTorso(i).RAcr(k,:)*Rot(:,:,k)' + Trials(i).Cycles.PelvisOrigin(:,:,j);
            end
            NewTorso(i).Cycles(:,:,j) = [NewClav(:,:,j) NewLAcr(:,:,j) NewRAcr(:,:,j)];
        end

        % plot to show relative locations
        % for i = DynamicTrials
        %         figure; hold on;
        %         for j = 1:data(i).TSData.NumGCs-1
        %             Clav = plot3(NewClav(:,1,j), NewClav(:,2,j), NewClav(:,3,j), 'r');
        %             LAcr = plot3(NewLAcr(:,1,j), NewLAcr(:,2,j), NewLAcr(:,3,j), 'g');
        %             RAcr = plot3(NewRAcr(:,1,j), NewRAcr(:,2,j), NewRAcr(:,3,j), 'b');
        %             Pelvis = plot3(data(i).Cycles.Pelvis(:,[1 4 7],j), data(i).Cycles.Pelvis(:,[2 5 8],j), data(i).Cycles.Pelvis(:,[3 6 9],j), 'k');
        %             plot3(data(i).Cycles.PelvisOrigin(:,1,j), data(i).Cycles.Pelvis(:,2,j), data(i).Cycles.Pelvis(:,3,j), 'k')
        %         end
        %         axis equal;
        %         xlabel('X'); ylabel('Y'); zlabel('Z');
        %         legend({'Clav','LAcr','RAcr', 'Pelvis'});
        %         title(strcat('Global positions for ',data(i).name));
    end
end
clearvars Ster LAcr RAcr i j

%% resample to time of each gait cycle
for i = 1:length(DynIndex)
    if i == find(~DynIndex)
        % use whoelle time for static trials
        Trials(i).NewData = NewTorso(i).Cycles;
        Trials(i).Export = horzcat(Trials(i).data, Trials(i).NewData);
    else
        for j = 1: Trials(i).TSData.NumGCs-1
            % define gait cycles
            Start = round(Trials(i).GCTimes(j,2));
            if Start < 1
                Start = 1;
            end
            End = round(Trials(i).GCTimes(j,4));
            CycDur =  End - Start + 1;
            % resample back to original gait cycle times
            NewTorso(i).OrigCycles(j).Torso = resample(NewTorso(i).Cycles(:,:,j), CycDur, 100);
            Trials(i).NewData(Start:End, :) = NewTorso(i).OrigCycles(j).Torso;
        end
        
        % check for gaps or undefined data
        if sum(isnan(Trials(i).NewData)) > 0
            Msg = ['Trials(i).subject' ,' ', 'Trials(i).name', ' ', 'has gaps or invalid data'];
            error(Msg);
            % get rid of NANs?
        end
        
        % first apply median filter to remove spikes
        Trials(i).MedFiltData = medfilt1(Trials(i).NewData, 3);
        
        % then apply butterworth filter data to smooth
        % accounts for movements between cycles
        T = Trials(i).data(2,2) - Trials(i).data(1,2); % sample time
        cutoff = 6; % cutoff frequency in Hz
        [b, a] = butter(4,cutoff*T*2,'low'); % design filter
        Trials(i).FiltData = filtfilt(b,a,Trials(i).MedFiltData);
        
        % plot filtering
        figure; subplot(311); hold on; title('X');
        plot(Trials(i).NewData(:,[1 4 7]),'.');
        plot(Trials(i).MedFiltData(:,[1 4 7]),'--');
        plot(Trials(i).FiltData(:,[1 4 7]),'-');
        subplot(312); hold on; title('Y');
        plot(Trials(i).NewData(:,[2 5 8]),'.');
        plot(Trials(i).MedFiltData(:,[2 5 8]),'--');
        plot(Trials(i).FiltData(:,[2 5 8]),'-');
        subplot(313); hold on; title('Z');
        plot(Trials(i).NewData(:,[3 6 9]),'.');
        plot(Trials(i).MedFiltData(:,[3 6 9]),'--');
        plot(Trials(i).FiltData(:,[3 6 9]),'-');
        
        
        Trials(i).FiltData(End:length(Trials(i).data), :) = 0; % place zeros after array to make it the same length as original
        TimeTrim = Trials(i).NewData(:,1) == 0; % identify zeros for trimming to time (outside of left gait cycles)
        Trials(i).Export= horzcat(Trials(i).data, Trials(i).FiltData);
        Trials(i).Export(TimeTrim, :) = []; % trim new combine data to times between left gait cycles
        Trials(i).Export(1:25,:) = []; % trim again to make sure end points arent affected by filtering
        Trials(i).Export(end-25:end,:) = [];
        clearvars Start End CycDur TimeTrim PreFilt
    end
end

%% Put New Cycles back into TRC file
for i = 1:length(DynIndex)

    % Copy file to new name (AddTorso.trc)
    NewFilename = strcat(Trials(i).folder, '\OpenSim\', Trials(i).files.OpenSimTRC(1:end-4), 'AddTorso.trc');
    Trials(i).AddTorsoFullFileName = NewFilename;
    Trials(i).files.OpenSimAddTorso = strcat(Trials(i).files.OpenSimTRC(1:end-4), 'AddTorso.trc');
    disp(strcat('Writing ', Trials(i).files.OpenSimAddTorso));
    Src = strcat(Trials(i).folder, '\OpenSim\', Trials(i).files.OpenSimTRC);
    copyfile(Src, NewFilename);
    % Correct Headers
    fid=fopen(NewFilename,'r'); % Read TRC files to determine correct headers
    FData = textscan(fid,'%q'); % Looks at headers
    FData = FData{1,1}; % Reorganizes
    FData{4,1}=NewFilename; % Renames file name
    fclose(fid);
    
    % find header locations in
    StartMetrics = find(strcmp(FData,'DataRate'), 1);
    LastCategory = find(strcmp(FData,'OrigNumFrames'), 1);
    ColumnHeaderStart = find(strcmp(FData,'Frame#'), 1);
    
    % get names of column headers
    [HeaderData] = ClarifyHeaders(Src, Trials(i));
    HeaderMarkers = HeaderData.colheaders;
    HeaderMarkers(2,:) = [];
    ToDel = zeros(length(HeaderMarkers), 1); 
    for j = 1:length(HeaderMarkers)
        if strcmp(HeaderMarkers{j},' ')
            ToDel(j) = 1;
        elseif isempty(HeaderMarkers{j})
            ToDel(j) = 1;
        end
    end
    HeaderMarkers(logical(ToDel)) = [];
    HeaderMarkers(1, end+1:end+3) = {'Clav','L.Acromium','R.Acromium'}; % Add labels for new torso markers
    
    % Write data to new TRC file
    fid=fopen(NewFilename,'w');
    
    % Copies over header info
    fprintf(fid,'%s\t', FData{1:StartMetrics-1,1}); % Top line
    fprintf(fid,'\n'); % new line
    fprintf(fid,'%s\t', FData{StartMetrics:LastCategory,1}); % second line
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
    [NumRow, ~] = size(Trials(i).Export);
    for ii=1:NumRow
        fprintf(fid,'%3.4f\t',Trials(i).Export(ii,:)); fprintf(fid,'\n');
    end
    fclose(fid);
    fclose('all');
    
    % OpenSim won't recognized the file unless it is opened and saved again
    e = actxserver('excel.application');
    eW = e.Workbooks;
    eF = eW.Open(NewFilename); % open OutputTest.xls
    eS = eF.ActiveSheet; 
    eF.Save;
    eF.Close; % close the file
    e.Quit; % close Excel entirely
    
    % export and clear variables
%     disp(strcat(NewFilename, ' Exported!'));
    clearvars fid FData e eW eF eS ii Src NewFilename m j ii fid FData j d ans ToDel
 
end
close all;
disp('ALL TRIALS CONVERTED');


end




   
    
    %     else % shorten IK time to region of interest
    %
    %         if strcmp(Trials(i).type, 'static') % if static trial, use whole time
    %             Trials(i).ExportTrim = Trials(i).Export;
    %         else
    %             Times = Trials(i).Export(:,2); % time series
    %             [~, StartInd] = min(abs(Trials(i).Times.Start(1) - Times));
    %             [~, EndInd] = min(abs(Trials(i).Times.Start(1)+10 - Times));
    %             Trials(i).ExportTrim = Trials(i).Export(StartInd:EndInd, :);
    %         end
    %
    %         % Copy file to new name (_AddTorso.trc)
    %         Trials(i).NewFileName = strcat(Trials(i).folder, '\OpenSim\', Trials(i).files.OpenSimTRC(1:end-4), '_AddTorso.trc');
    %         disp(strcat('Writing ', Trials(i).NewFileName));
    %         NewFilename = Trials(i).NewFileName;
    %         Src = strcat(Trials(i).folder, '\OpenSim\', Trials(i).files.OpenSimTRC);
    %         copyfile(Src, NewFilename);
    %         % Correct Headers
    %         fid=fopen(NewFilename,'r'); % Read TRC files to determine correct headers
    %         FData = textscan(fid,'%q'); % Looks at headers
    %         FData = FData{1,1}; % Reorganizes
    %         FData{4,1}=NewFilename; % Renames file name
    %         fclose(fid);
    %
    %         % find header locations in file
    %         LastCategory = find(strcmp(FData,'OrigNumFrames'), 1);
    %         ColumnHeaderStart = find(strcmp(FData,'Frame#'), 1);
    %         % get names of column headers
    %         [HeaderData] = ClarifyHeaders(Src, Trials(i));
    %         HeaderMarkers = HeaderData.colheaders;
    %         HeaderMarkers(2,:) = [];
    %         for j = 1:length(HeaderMarkers)
    %             if strcmp(HeaderMarkers{j},' ')
    %                 ToDel(j) = 1;
    %             elseif isempty(HeaderMarkers{j})
    %                 ToDel(j) = 1;
    %             else
    %                 ToDel(j) = 0;
    %             end
    %         end
    %         HeaderMarkers(logical(ToDel)) = [];
    %         HeaderMarkers(1, end+1:end+3) = {'Clav','L.Acromium','R.Acromium'}; % Add labels for new torso markers
    %
    %         % Write data to new TRC file
    %         fid=fopen(NewFilename,'w');
    %
    %         % Copies over header info
    %         fprintf(fid,'%s\t', FData{1:5,1}); % Top line
    %         fprintf(fid,'\n'); % new line
    %         fprintf(fid,'%s\t', FData{6:LastCategory,1}); % second line
    %         fprintf(fid,'\n'); % new line
    %         fprintf(fid,'%s\t', FData{LastCategory+1:LastCategory+3,1});
    %         fprintf(fid,'%s\t', num2str(length(HeaderMarkers)));
    %         fprintf(fid,'%s\t', FData{LastCategory+5:ColumnHeaderStart-1,1});
    %         fprintf(fid,'\n'); % new line
    %         fprintf(fid,'%s\t', FData{ColumnHeaderStart:ColumnHeaderStart+1,1});
    %
    %         % Copies over column headers with two tabs in between each marker
    %         NCol = length(HeaderMarkers);
    %         for m=1:NCol
    %             fprintf(fid,'%s\t',HeaderMarkers{m});
    %             fprintf(fid,'\t'); % tab
    %             fprintf(fid,'\t'); % tab
    %         end
    %
    %         % new line and tab over twice to pass the Frame and Time columns
    %         fprintf(fid,'\n'); fprintf(fid,'\t'); fprintf(fid,'\t');
    %
    %         % Labels x, y, and z columns
    %         d=0;
    %         for m=1:NCol
    %             d=d+1;
    %             fprintf(fid,'%s\t', ['X',num2str(d)]);
    %             fprintf(fid,'%s\t', ['Y',num2str(d)]);
    %             fprintf(fid,'%s\t', ['Z',num2str(d)]);
    %         end
    %
    %         % Adds a space between xyz and data.
    %         fprintf(fid,'\n');
    %         fprintf(fid,'\n');
    %
    %         % Inputs new/old marker x, y, z info
    %         [NumRow, ~] = size(Trials(i).ExportTrim);
    %         for ii=1:NumRow
    %             fprintf(fid,'%3.4f\t',Trials(i).ExportTrim(ii,:)); fprintf(fid,'\n');
    %         end
    %         fclose(fid);
    %         fclose('all');
    %
    %         % OpenSim won't recognized the file unless it is opened and saved again
    %         e=actxserver('excel.application');
    %         eW=e.Workbooks;
    %         eF=eW.Open(NewFilename); % open OutputTest.xls
    %         eS=eF.ActiveSheet; eF.Save;
    %         eF.Close; % close the file
    %         e.Quit; % close Excel entirely
    %         % export and clear variables
    %         disp(strcat(NewFilename, ' Exported!'));
    %         clearvars fid FData e eW eF eS ii Src NewFilename m j ii fid FData j d ans ToDel
    %
    %     end
