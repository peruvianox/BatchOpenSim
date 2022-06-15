function [mass, FPused, files, directoryf, Analog] = convertFPdata_OpenSim(input_file, FPcal_data, directory, Settings)
% done=writeForcesFile(time,forcedata,directory,file)
% writes a motion analysis style forces file which includes the
% forces (Fx, Fy, Fz) ceneter of pressure (x, y, z) and moment (x, y, z)
% Inputs:
% input_file: file name to be converted to forces file
% FPcal_file: .cal file that includes force plate data
% zero_file: zero forces file to zero forceplates
% directory: Where the output file should be written

% edited by Ricky Pimentel October 2019
% added comments and split code into functional sections to improve readability


%% Check for FPcal file otherwise load one in selected directory
% if narg >= 2
    % Load FPcal_file
    %     [S,pos,origin,R] = load_fpcal(FPcal_file);
    
% define FPcal values from loaded structure
S = FPcal_data.S;
pos = FPcal_data.pos;
origin = FPcal_data.origin;
R = FPcal_data.R;

files = input_file; 

%% Fy threshold value to eliminate undefined COP Calulations
fythresh = 10;
if isstring(files) == 0
    loop_index = 1;
else
    loop_index = size(files,2);
    k = 0;
end

%% Loop through files
Analog(loop_index).voltdata = [];
Analog(loop_index).voltdata_filt = [];
Analog(loop_index).channel_names = [];

for w = 1:loop_index
    
    % index filenames for loading
    if loop_index == 1
        input_file = files(1:end);
    else
        k = k+1;
        input_file = cell2mat(files(k));
        input_file = input_file(1:end);
    end
    
    %% Load .anc file from input_file name
    %     [data, time, ~, channel_names, ~, voltdata] = load_anc(input_file);
    [~, time, ~, channel_names, ~, voltdata] = load_anc_fast(input_file);
    
    %     if zero == 'y'
    %         % Zero the voltdata if user desires
    %         voltdata=voltdata-ones(size(voltdata,1),1)*mean_zero;
    %     end
    
    % filter force data
    T = time(2)-time(1); % sample time
    cutoff = 12; % cutoff frequency in Hz for GRFs
    [b, a] = butter(3,cutoff*T*2,'low'); % design filter
    voltdata_filt = filtfilt(b,a,voltdata); % do filtering
    
    Analog(w).voltdata = voltdata;
    Analog(w).voltdata_filt = voltdata_filt;
    Analog(w).channel_names = channel_names;
    
    % frequency analysis ?
    %     NFFT=2^nextpow2(size(voltdata,1));
    %     Y=fft(voltdata(:,3),NFFT)/(size(voltdata,1));
    %     f=(1/T)/2*linspace(0,1,NFFT/2+1);
    %     figure(2); hold on; plot(f,2*abs(Y(1:NFFT/2+1)),'g'); axis([0 200 0 0.1]);
    
    %% Plot Voltage
    if strcmp(Settings.PlotVoltage, 'Yes')
        % select time window to plot
        s = 1; % start frame
        window = s:s+(3/T); % next 3 seconds
        LW = 2; % set line width
        %     time(window) % uncomment to display time window plotted in command window
        disp_name = strrep(input_file, '_', ' ');
        
        % plot filtering comparison
        figure('Position', [100 100 1200 800]); % FiltComp =
        
        subplot(231); hold on;
        plot(time(window),voltdata(window,strcmp(channel_names, 'F1X')),'-k', 'LineWidth', LW/2);
        plot(time(window),voltdata_filt(window,strcmp(channel_names, 'F1X')),'--b', 'LineWidth', LW);
        legend({'raw','filt'});
        xlabel('Time (s)'); ylabel('Analog Voltage (V)');
        xlabel(num2str(s))
        title(strcat(disp_name, 'FP # 1 X'));
        
        subplot(232); hold on;
        plot(time(window),voltdata(window,strcmp(channel_names, 'F1Y')),'-k', 'LineWidth', LW/2);
        plot(time(window),voltdata_filt(window,strcmp(channel_names, 'F1Y')),'--b', 'LineWidth', LW);
        legend({'raw','filt'});
        xlabel('Time (s)'); ylabel('Analog Voltage (V)');
        title(strcat(disp_name, 'FP # 1 Y'));
        
        subplot(233); hold on;
        plot(time(window),voltdata(window,strcmp(channel_names, 'F1Z')),'-k', 'LineWidth', LW/2);
        plot(time(window),voltdata_filt(window,strcmp(channel_names, 'F1Z')),'--b', 'LineWidth', LW);
        legend({'raw','filt'});
        xlabel('Time (s)'); ylabel('Analog Voltage (V)');
        title(strcat(disp_name, 'FP # 1 Z'));
        
        subplot(234); hold on;
        plot(time(window),voltdata(window,strcmp(channel_names, 'F2X')),'-k', 'LineWidth', LW/2);
        plot(time(window),voltdata_filt(window,strcmp(channel_names, 'F2X')),'--b', 'LineWidth', LW);
        legend({'raw','filt'});
        xlabel('Time (s)'); ylabel('Analog Voltage (V)');
        title(strcat(disp_name, 'FP # 2 X'));
        
        subplot(235); hold on;
        plot(time(window),voltdata(window,strcmp(channel_names, 'F2Y')),'-k', 'LineWidth', LW/2);
        plot(time(window),voltdata_filt(window,strcmp(channel_names, 'F2Y')),'--b', 'LineWidth', LW);
        legend({'raw','filt'});
        xlabel('Time (s)'); ylabel('Analog Voltage (V)');
        title(strcat(disp_name, 'FP # 2 Y'));
        
        subplot(236); hold on;
        plot(time(window),voltdata(window,strcmp(channel_names, 'F2Z')),'-k', 'LineWidth', LW/2);
        plot(time(window),voltdata_filt(window,strcmp(channel_names, 'F2Z')),'--b', 'LineWidth', LW);
        legend({'raw','filt'});
        xlabel('Time (s)'); ylabel('Analog Voltage (V)');
        title(strcat(disp_name, 'FP # 2 Z'));
    end
    
    %% Find which forceplates are being used
    %     mean_data = mean(data(:,1:size(data,2)));
    FPused=[];
    
    % New lab at UNC-CH MEJ only has two plates
    F1=sum(strcmp(channel_names,'F1Z'));
    F2=sum(strcmp(channel_names,'F2Z'));
    
    if F1==1 && F2==1 % set number of force plates used
        FPused=[1 2];
    end
    
    
    %% Find active Forceplates and create new voltdata for those forceplates
    FP_Loc = [];
    VDataNew = [];
    ChannelsNew=[];
    
    for i = 1:1:length(FPused) % loop through # of force plates
        FPnumber = FPused(i);
        FPnumber = int2str(FPnumber);
        FPname = ['F' FPnumber 'X'];
        FP_location = strmatch(FPname, channel_names(:),'exact');
        FP_Loc = [FP_Loc FP_location];
        VDataNew = [VDataNew voltdata_filt(:,FP_location:(FP_location+5))];
        ChannelsNew = [ChannelsNew channel_names(FP_location:(FP_location+5))];
    end
    
    %% Convert Voltdata to Force, Moment, and CoP Data
    Output = [];
    [g,~] = size(VDataNew);
    Xcur = zeros(g,1);
    Ycur = zeros(g,1);
    
    for j = 1:1:length(FPused) % loop through force plates used
        FPnumber = FPused(j);
        a = FPnumber;
        %For forceplates 1 through 3
        %     if FPnumber <= 3
        %     G=4000; %FP Amplifier Gain
        %     Vo=10;  %FP Excitation Voltage
        %     CF=10^6; %conversion factor from uV to V
        %     else
        % undefined for treadmill forceplates
        G = 1;
        Vo = 1;
        CF = 1;
        %     end
        
        FPname = ['F' int2str(a) 'X'];
        FP_location = find(strcmp(FPname,ChannelsNew(:)));    %Finds location of data based on channel names
        Trans = [R(:,:,a) zeros(size(R(:,:,a))); zeros(size(R(:,:,a))) R(:,:,a)]; %Transformation matrix from FPcal file
        ForMom = -VDataNew(:,FP_location:(FP_location+5))*S(:,:,a)'*(CF)*Trans/(Vo*G); % Negative to get force on person
        toffset = cross(ones(size(ForMom,1),1)*(-R(:,:,a)*origin(1:3,a))',ForMom(:,1:3)); %Calculate Torque due to transducer offset
        
        % Allow for zero value of z force by setting threshold and setting COP to 0 if z force is below threshold
        for jk = 1:1:g
            if  ForMom(jk,3) >= fythresh || ForMom(jk,3) <= -fythresh
                x = ((-ForMom(jk,5)-toffset(jk,2))./ForMom(jk,3)+pos(1,a))*1000; %x axis Center of pressure
                y = ((ForMom(jk,4)+toffset(jk,1))./ForMom(jk,3)+pos(2,a))*1000; %y axis Center of Pressure
                
                Xcur(jk) = x; %[Xcur; x];
                Ycur(jk) = y; %[Ycur; y];
            else
                x = 0;
                y = 0;
                Xcur(jk) = x;%[Xcur; x];
                Ycur(jk) = y;%[Ycur; y];
            end
        end
        
        % Calculate COPs
        % dont threshold COPs until after fine-tuning
%         Xcur_UnFilt = ((-ForMom(:,5)-toffset(:,2))./ForMom(:,3)+pos(1,a))*1000; %x axis Center of pressure
%         Ycur_UnFilt = ((ForMom(:,4)+toffset(:,1))./ForMom(:,3)+pos(2,a))*1000; %y axis Center of Pressure
%         T = time(2)-time(1); % sample time
%         cutoff = 5; % cutoff frequency in Hz for GRFs
%         [b, a] = butter(3,cutoff*T*2,'low'); % design filter
%         Xcur = filtfilt(b,a,Xcur_UnFilt); % do filtering
%         Ycur = filtfilt(b,a,Ycur_UnFilt); % do filtering
        
% plot to check CoPs
%         figure; 
%         subplot(211); hold on;
%         plot(Xcur_UnFilt, '--r');
%         plot(Xcur, 'r');
%         subplot(212); hold on;
%         plot(Ycur_UnFilt, '--b');
%         plot(Ycur, 'b');
        

        %          for jk = 1:1:g
        % %             if  ForMom(jk,3) >= fythresh || ForMom(jk,3) <= -fythresh
        %                 Xcur(jk) = ((-ForMom(jk,5)-toffset(jk,2))./ForMom(jk,3)+pos(1,a))*1000; %x axis Center of pressure
        %                 Ycur(jk) = ((ForMom(jk,4)+toffset(jk,1))./ForMom(jk,3)+pos(2,a))*1000; %y axis Center of Pressure
        %         end
        
        % z axis Center of pressure Origin (-) to transform to lab reference frame
        Zcur = (zeros(size(Xcur)) - origin(3,a)+pos(3,a)) * 1000;
        % Zcur=zeros(size(Xcur)); %z axis Center of Pressure
        
        % Calculation of the Torque Mz
        torque = cross((ones(size(Xcur,1),1)*pos(1:3,a)'-[Xcur Ycur Zcur]*.001),ForMom(:,1:3));
        tx = (ForMom(:,4)+torque(:,1)+toffset(:,1)); %Mx output in Nm
        ty = (ForMom(:,5)+torque(:,2)+toffset(:,2)); %My output in Nm
        tz = (ForMom(:,6)+torque(:,3)+toffset(:,3)); %Mz output in Nm
        
        % Gather all data and add to previous data for output
        Output = [Output ForMom(:,1:3) Xcur Ycur Zcur tx ty tz];
    end
    
    
    %% Standard Thresholding
    
    % set indexing threshold for applying conversions
    FC_threshold = 20; % Force threshold at 20 N - crude
    forces = zeros(size(Output));
    fc_temp1 = find(Output(:,3) >= FC_threshold);
    fc_temp2 = find(Output(:,12) >= FC_threshold);
    Analog.Output = Output;
    
    % Rotate Outputs while applying thresholding
    
    % right ground forces - plate 1
    forces(fc_temp1,1) = Output(fc_temp1,1);
    forces(fc_temp1,2) = Output(fc_temp1,3);
    forces(fc_temp1,3) = -Output(fc_temp1,2);
    % right CoPs - plate 1
    forces(fc_temp1,4) = Output(fc_temp1,4)/1000; % divide by 1000 to convert from m to mm
    forces(fc_temp1,5) = Output(fc_temp1,6)/1000;
    forces(fc_temp1,6) = -Output(fc_temp1,5)/1000;
    % left ground forces - plate 2
    forces(fc_temp2,7) = Output(fc_temp2,10);
    forces(fc_temp2,8) = Output(fc_temp2,12);
    forces(fc_temp2,9) = -Output(fc_temp2,11);
    % left CoPs
    forces(fc_temp2,10) = Output(fc_temp2,13)/1000; % divide by 1000 to convert from m to mm
    forces(fc_temp2,11) = Output(fc_temp2,15)/1000;
    forces(fc_temp2,12) = -Output(fc_temp2,14)/1000;
    % ground torques
    % right plate (#1)
    forces(fc_temp1,13) = Output(fc_temp1,7);
    forces(fc_temp1,14) = Output(fc_temp1,9);
    forces(fc_temp1,15) = -Output(fc_temp1,8);
    % left plate (#2)
    forces(fc_temp2,16) = Output(fc_temp2,16);
    forces(fc_temp2,17) = Output(fc_temp2,18);
    forces(fc_temp2,18) = -Output(fc_temp2,17);
    
    %% Fine Tune Thresholding
    % if static trial, use all points that are on the FPs
    strings = {'Static','static', 'Static_Marker_Reference','HJC','hjc'};
    if contains(input_file, strings) 
        
        % do nothing, GRFs dont need any additional thresholding for static trials
        
    else % fine tune thresholding for walking trials
        
        Forces = zeros(size(Output));
        NewThresh = 10; % threshold in N for steps
        BackupThresh = 10; % in unable to apply new threshold, use backup
        window = 30; % search window for fine-tuning step timing
        
        % First force plate
        vGRF = Output(:,3); % pull vertical ground reaction force
        Threshed = forces(:,2);
        Changes = find(ischange(double(Threshed >= FC_threshold)) == 1);
        NewChanges = zeros(length(Changes),1);
        
        for i = 1:length(Changes)
            % if heel strike, next point is above zero and previous point is 0
            if Threshed(Changes(i)) > 0 && Threshed(Changes(i)-1) == 0
                if Changes(i)-window <= 0 % make sure to not go past current index
                    Srch = 1;
                else
                    Srch = Changes(i)-window;
                end
                % search backwards to find first (last) positive vGRF
                Count = find(vGRF(Changes(i) : -1 : Srch) < NewThresh,1) - 1;
                if isempty(Count) % backup
                    Count = find(vGRF(Changes(i) : -1 : Srch) < BackupThresh,1) - 1;
                    if isempty(Count) % worst case, use original threshold
                        Count = find(vGRF(Changes(i) : -1 : Srch) < FC_threshold,1) - 1;
                    end
                end
                NewChanges(i) = Changes(i) - Count;
                
                if i == 1
                    StartEvent = 'Strike'; % label first event
                end
                
                % if toe off, current point is zero and previous point is above 0
            elseif Threshed(Changes(i)) == 0 && Threshed(Changes(i)-1) > 0
                if Changes(i)+window > length(vGRF) % make sure to not go past current index
                    Srch = length(vGRF);
                else
                    Srch = Changes(i)+window;
                end
                % search forwards to find last (first) positive vGRF
                Count = find(vGRF(Changes(i):Srch) < NewThresh,1) - 1;
                if isempty(Count) % backup
                    Count = find(vGRF(Changes(i):Srch) < BackupThresh,1) - 1;
                    if isempty(Count) % worst case, use original threshold
                        Count = find(vGRF(Changes(i):Srch) < FC_threshold,1) - 1;
                    end
                end
                NewChanges(i) = Count + Changes(i);
                
                if i == 1
                    StartEvent = 'Off'; % label first event
                end
                
            end
        end
        
        New_fc_temp1 = zeros(length(Output),1); % initialize
        
        % apply new times to new logical
        if strcmp(StartEvent, 'Strike')
            for i = 1:length(NewChanges)
                if mod(i,2) == 1 % odd event is strike
                    if i == length(NewChanges)
                        END = length(Output);
                    else
                        END = NewChanges(i+1);
                    end
                    New_fc_temp1(NewChanges(i):END,1) = 1;
                end
            end
        elseif strcmp(StartEvent, 'Off')
            for i = 1:length(NewChanges)
                if i == 1 % signal to start off pulling data
                    New_fc_temp1(1:NewChanges(i)) = 1;
                end
                if mod(i,2) == 0 % even event is strike
                    if i == length(NewChanges)
                        END = length(Output);
                    else
                        END = NewChanges(i+1);
                    end
                    New_fc_temp1(NewChanges(i):END,1) = 1;
                end
            end
        end
        New_fc_temp1 = logical(New_fc_temp1);
        
%         clearvars vGRF Threshed ThreshLog Changes NewChanges Count StartEvent
        
        % Second force plate
        vGRF = Output(:,12); % pull vertical ground reaction force
        Threshed = forces(:,8);
        ThreshLog = ischange(double(Threshed >= FC_threshold));
        Changes = find(ThreshLog == 1);
        NewChanges = zeros(length(Changes),1);
        
        for i = 1:length(Changes)
            % if heel strike, current point is above zero and previous point is 0
            if Threshed(Changes(i)) > 0 && Threshed(Changes(i)-1) == 0
                if Changes(i)-window <= 0 % make sure to not go past current index
                    Srch = 1;
                else
                    Srch = Changes(i)-window;
                end
                % search backwards to find first (last) positive vGRF
                Count = find(vGRF(Changes(i) : -1 : Srch) < NewThresh,1) - 1;
                if isempty(Count) % backup
                    Count = find(vGRF(Changes(i) : -1 : Srch) < BackupThresh,1) - 1;
                    if isempty(Count) % worst case, use original threshold
                        Count = find(vGRF(Changes(i) : -1 : Srch) < FC_threshold,1) - 1;
                    end
                end
                NewChanges(i) = Changes(i) - Count;
                
                if i == 1
                    StartEvent = 'Strike'; % label first event
                end
                
                % if toe off, current point is zero and previous point is above 0
            elseif Threshed(Changes(i)) == 0 && Threshed(Changes(i)-1) > 0
                if Changes(i)+window > length(vGRF) % make sure to not go past current index
                    Srch = length(vGRF);
                else
                    Srch = Changes(i)+window;
                end
                % search forwards to find last (first) positive vGRF
                Count = find(vGRF(Changes(i):Srch) < NewThresh,1) - 1;
                if isempty(Count)
                    Count = find(vGRF(Changes(i):Srch) < BackupThresh,1) - 1;
                    if isempty(Count) % worst case, use original threshold
                        Count = find(vGRF(Changes(i):Srch) < FC_threshold,1) - 1;
                    end
                end
                NewChanges(i) = Count + Changes(i);
                
                if i == 1
                    StartEvent = 'Off'; % label first event
                end
                
            end
        end
        
        New_fc_temp2 = zeros(length(Output),1); % initialize new timing logical
        
        % apply new times to new logical
        if strcmp(StartEvent, 'Strike')
            for i = 1:length(NewChanges)
                if mod(i,2) == 1 % odd event is strike
                    if i == length(NewChanges)
                        END = length(Output);
                    else
                        END = NewChanges(i+1)-1;
                    end
                    New_fc_temp2(NewChanges(i)+1:END,1) = 1;
                end
            end
        elseif strcmp(StartEvent, 'Off')
            for i = 1:length(NewChanges)
                if i == 1 % signal to start off pulling data
                    New_fc_temp2(1:NewChanges(i)) = 1;
                end
                if mod(i,2) == 0 % even event is strike
                    if i == length(NewChanges)
                        END = length(Output);
                    else
                        END = NewChanges(i+1)-1;
                    end
                    New_fc_temp2(NewChanges(i)+1:END,1) = 1;
                end
            end
        end
        New_fc_temp2 = logical(New_fc_temp2);
        
        %% rotate outputs while applying new thresholding
        % right ground forces - plate 1
        Forces(New_fc_temp1,1) = Output(New_fc_temp1,1);
        Forces(New_fc_temp1,2) = Output(New_fc_temp1,3);
        Forces(New_fc_temp1,3) = -Output(New_fc_temp1,2);
        % right CoPs - plate 1
        Forces(New_fc_temp1,4) = Output(New_fc_temp1,4) ./ 1000; % divide by 1000 to convert from m to mm
        Forces(New_fc_temp1,5) = Output(New_fc_temp1,6) ./ 1000;
        Forces(New_fc_temp1,6) = -Output(New_fc_temp1,5) ./ 1000;
        
        % left ground Forces - plate 2
        Forces(New_fc_temp2,7) = Output(New_fc_temp2,10);
        Forces(New_fc_temp2,8) = Output(New_fc_temp2,12);
        Forces(New_fc_temp2,9) = -Output(New_fc_temp2,11);
        % left CoPs
        Forces(New_fc_temp2,10) = Output(New_fc_temp2,13) ./ 1000; % divide by 1000 to convert from m to mm
        Forces(New_fc_temp2,11) = Output(New_fc_temp2,15) ./ 1000;
        Forces(New_fc_temp2,12) = -Output(New_fc_temp2,14) ./ 1000;
        
        % ground torques
        % right plate (#1)
        Forces(New_fc_temp1,13) = Output(New_fc_temp1,7);
        Forces(New_fc_temp1,14) = Output(New_fc_temp1,9);
        Forces(New_fc_temp1,15) = -Output(New_fc_temp1,8);
        % left plate (#2)
        Forces(New_fc_temp2,16) = Output(New_fc_temp2,16);
        Forces(New_fc_temp2,17) = Output(New_fc_temp2,18);
        Forces(New_fc_temp2,18) = -Output(New_fc_temp2,17);
        
        % overwrite old forces values
        forces = Forces;
        clearvars Forces;
    end
    
    %% plot updated thresholds
    %     s = 79707;
    %     ind = 1:length(forces);
    %     figure; hold on;
    %     plot(forces(ind, 8), 'b');
    %     plot(Output(ind, 12), '--r');
    % %     plot(vGRF(ind), '--k')
    %    xlim([ s-1000 s+1000]);
    %     plot(Forces(ind,8), 'g');
    %
    %     legend({'Original','Raw','New'});
    %
    %% if y and z are correct?
    %     forces(fc_temp1,1:3)=Output(fc_temp1,1:3);
    %     forces(fc_temp1,4:6)=Output(fc_temp1,4:6)/1000;
    %     forces(fc_temp2,7:9)=Output(fc_temp2,10:12);
    %     forces(fc_temp2,10:12)=Output(fc_temp2,13:15)/1000;
    %     forces(fc_temp1,13:15)=Output(fc_temp1,7:9)/1000;
    %     forces(fc_temp2,16:18)=Output(fc_temp2,16:18)/1000;
    
    %% Plot GRFs
    if strcmp(Settings.PlotGRFs, 'Yes')
        Forces.name = strrep(input_file, '.anc', '.forces');
        Forces.data = importdata(Forces.name);
        Forces.data.colheaders =  {'Sample','FX1','FY1','FZ1','X1','Y1','Z1','MZ1',...
            'FX2','FY2','FZ2','X2','Y2','Z2','MZ2'};
        MainTrialName = strsplit(Forces.name, '_');
        
        % pull columns with data
        Forces.Right.X = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FX1'));
        Forces.Right.Y = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FY1'));
        Forces.Right.Z = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FZ1'));
        Forces.Left.X = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FX2'));
        Forces.Left.Y = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FY2'));
        Forces.Left.Z = Forces.data.data(:,strcmp(Forces.data.colheaders, 'FZ2'));
        
        GRF.ColHeaders = {'ground_force_vx','ground_force_vy','ground_force_vz',...
            'ground_force_px','ground_force_py','ground_force_pz',...
            '1_ground_force_vx','1_ground_force_vy','1_ground_force_vz',...
            '1_ground_force_px','1_ground_force_py','1_ground_force_pz',...
            'ground_torque_x','ground_torque_y','ground_torque_z',...
            '1_ground_torque_x','1_ground_torque_y','1_ground_torque_z'};
        GRF.AllData = forces;
        GRF.FileName = input_file;
        GRF.Right.X = GRF.AllData(:,strcmp(GRF.ColHeaders, 'ground_force_vx'));
        GRF.Right.Y = GRF.AllData(:,strcmp(GRF.ColHeaders, 'ground_force_vz'));
        GRF.Right.Z = GRF.AllData(:,strcmp(GRF.ColHeaders, 'ground_force_vy'));
        GRF.Left.X = GRF.AllData(:,strcmp(GRF.ColHeaders, '1_ground_force_vx'));
        GRF.Left.Y = GRF.AllData(:,strcmp(GRF.ColHeaders, '1_ground_force_vz'));
        GRF.Left.Z = GRF.AllData(:,strcmp(GRF.ColHeaders, '1_ground_force_vy'));
        
        % Plot
        % close all;
        LW = 1;
        %         MkrSz = 10;
        Ind = 1:5000;
        ForceComp = figure('Position',[100 100 1200 800]);
        % Left side
        subplot(321); hold on;
        OpenSimX = plot(GRF.Left.X(Ind), '-r', 'LineWidth',LW);
        ForcesX = plot(Forces.Left.X(Ind), '--b',  'LineWidth',LW);
        title('Left X Forces');
        xlabel('Frame'); ylabel('N');
        legend({'OpenSim', 'Lab'});
        
        subplot(323); hold on;
        OpenSimX = plot(-GRF.Left.Y(Ind), '-r', 'LineWidth',LW);
        ForcesX = plot(Forces.Left.Y(Ind), '--b',  'LineWidth',LW);
        title('Left Y Forces');
        xlabel('Frame'); ylabel('N');
        
        subplot(325); hold on;
        OpenSimZ = plot(GRF.Left.Z(Ind), '-r', 'LineWidth',LW);
        ForcesZ = plot(Forces.Left.Z(Ind), '--b',  'LineWidth',LW);
        title('Left Z Forces');
        xlabel('Frame'); ylabel('N');
        
        % Right side
        subplot(322); hold on;
        OpenSimX = plot(GRF.Right.X(Ind), '-r', 'LineWidth',LW);
        ForcesX = plot(Forces.Right.X(Ind), '--b',  'LineWidth',LW);
        title('Right X Forces');
        xlabel('Frame'); ylabel('N');
        legend({'OpenSim', 'Lab'});
        
        subplot(324); hold on;
        OpenSimX = plot(-GRF.Right.Y(Ind), '-r', 'LineWidth',LW);
        ForcesX = plot(Forces.Right.Y(Ind), '--b',  'LineWidth',LW);
        title('Right Y Forces');
        xlabel('Frame'); ylabel('N');
        
        subplot(326); hold on;
        OpenSimZ = plot(GRF.Right.Z(Ind), '-r', 'LineWidth',LW);
        ForcesZ = plot(Forces.Right.Z(Ind), '--b',  'LineWidth',LW);
        title('Right Z Forces');
        xlabel('Frame'); ylabel('N');
        
        
        ForceName = strrep(Forces.name,'_', ' ');
        GRFName = strrep(GRF.FileName,'_', ' ');
        supertitle(strcat(MainTrialName{1}, ' -  .forces and GRF.mot file comparison')); 
        % saveas(ForceComp, strcat(MainTrialName{1}, '_ForceComp.png'));
    end
    
    %% Plot Moments
    if strcmp( Settings.PlotMoments, 'Yes')
        
        Moments.name = strrep(input_file, '.anc', '.forces');
        Moments.data = importdata(Moments.name);
        Moments.data.colheaders =  {'Sample','FX1','FY1','FZ1','X1','Y1','Z1','MZ1',...
            'FX2','FY2','FZ2','X2','Y2','Z2','MZ2'};
        MainTrialName = strsplit(Moments.name, '_');
        
        % pull columns with data
        Moments.Right.Z = Moments.data.data(:,strcmp(Moments.data.colheaders, 'MZ1'));
        Moments.Left.Z = Moments.data.data(:,strcmp(Moments.data.colheaders, 'MZ2'));
        
        GRF.ColHeaders = {'ground_force_vx','ground_force_vy','ground_force_vz',...
            'ground_force_px','ground_force_py','ground_force_pz',...
            '1_ground_force_vx','1_ground_force_vy','1_ground_force_vz',...
            '1_ground_force_px','1_ground_force_py','1_ground_force_pz',...
            'ground_torque_x','ground_torque_y','ground_torque_z',...
            '1_ground_torque_x','1_ground_torque_y','1_ground_torque_z'};
        GRF.AllData = forces;
        GRF.FileName = input_file;
        GRF.Right.Z = GRF.AllData(:,strcmp(GRF.ColHeaders, 'ground_torque_y'));
        GRF.Left.Z = GRF.AllData(:,strcmp(GRF.ColHeaders, '1_ground_torque_y'));
        
        % Plot
        % close all;
        LW = 1;
        %         MkrSz = 10;
        Ind = 1:2000;
        ForceComp = figure('Position',[100 100 1200 800]);
        % Left side
        subplot(121); hold on;
        OpenSimZ = plot(GRF.Left.Z(Ind), '-r', 'LineWidth',LW);
        MomentsZ = plot(Moments.Left.Z(Ind)./1000, '--b',  'LineWidth',LW);
        title('Left Z Moments');
        xlabel('Frame'); ylabel('Nm');
        legend({'OpenSim', 'Lab'});
        
        % Right side
        subplot(122); hold on;
        OpenSimZ = plot(GRF.Right.Z(Ind), '-r', 'LineWidth',LW);
        MomentsZ = plot(Moments.Right.Z(Ind)./1000, '--b',  'LineWidth',LW);
        title('Right Z Moments');
        xlabel('Frame'); ylabel('Nm');
        legend({'OpenSim', 'Lab'});
        supertitle(strcat(MainTrialName{1}, ' -  .forces and GRF.mot Moment comparison'));
    end
    
    %% rename files and export data
    output_file =  strrep(input_file, '.anc',  '_OpenSimGRF.mot'); % define output file name
    
    % define column headers
    ChanOutds = {'time','ground_force_vx','ground_force_vy','ground_force_vz',...
        'ground_force_px','ground_force_py','ground_force_pz',...
        '1_ground_force_vx','1_ground_force_vy','1_ground_force_vz',...
        '1_ground_force_px','1_ground_force_py','1_ground_force_pz',...
        'ground_torque_x','ground_torque_y','ground_torque_z',...
        '1_ground_torque_x','1_ground_torque_y','1_ground_torque_z'};
    
    % ensure same dimensions of output array
    [~, Ncols] = size(Output); 
    if Ncols < length(ChanOutds)
        Output(:,2:end+1) = Output(:,1:end);
        Output(:, 1) = time; 
    end
    Analog.Output = Output; % save output
    O = array2table(Output);
    O.Properties.VariableNames = ChanOutds;
    
    disp(['Writing   ' output_file]); 
    Osim.writeMOT(O, 'FilePath', ['OpenSim/' output_file]);
    
    
end


end









