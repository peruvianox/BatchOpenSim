function [GC] = RankGC(Trials)

% identify and rank gait cycles from GRF, kinematics, and kinetics
% also excludes crossover steps during ranking

close all;
% get directory path
SubjPath = [Trials.folder '\OpenSim\'];

% get times of each step
L_GC_t = Trials.TSData.L_Strike(:,2);
R_GC_t = Trials.TSData.R_Strike(:,2);

LGC = struct([]);  % initialize structures
RGC = struct([]);

%% GRFs
% load
GRF = Osim.readMOT([SubjPath Trials.files.OpenSimGRF]);
Time = GRF.Header; % get time series

% delete extra gait cycles if there any after 120 s
d = find(L_GC_t > 120); % left
if isempty(d) == 0
    L_GC_t(d) = [];
end
d = find(R_GC_t > 120); % right
if isempty(d) == 0
    R_GC_t(d) = [];
end


% LEFT
LgrfCols = 8:13; % identify grf and cop columns
for i = 1:length(L_GC_t)-1 % cut and resample for each gait cycle
    [~,A] = min(abs(Time - L_GC_t(i))); % identify gait cycles
    [~,B] = min(abs(Time - L_GC_t(i+1)));
    
    % parse gait cycle to 100 points
    LGC(i).GRFs = interp1(Time(A:B), table2array(GRF(A:B,LgrfCols)), linspace(Time(A), Time(B)));
    LGC(i).GRFnames = GRF.Properties.VariableNames(LgrfCols); % save variable names
end
% plot(LGC(1).GRFs(:,1:3))

% RIGHT
RgrfCols = 2:7;
for i = 1:length(R_GC_t)-1 % cut and resample for each gait cycle
    [~,A] = min(abs(Time - R_GC_t(i))); % identify gait cycles
    [~,B] = min(abs(Time - R_GC_t(i+1)));
    
    % parse gait cycle to 100 points
    RGC(i).GRFs = interp1(Time(A:B), table2array(GRF(A:B,RgrfCols)), linspace(Time(A), Time(B)));
    RGC(i).GRFnames = GRF.Properties.VariableNames(RgrfCols); % save variable names
end


%% kinematics
TName = Trials.name;
% load
IK = Osim.readMOT([SubjPath 'IK_Files\' TName '_IK.mot']);
Time = IK.Header; % get time series

% LEFT
LikCols = zeros(1, length(IK.Properties.VariableNames));
for n = 1:length(IK.Properties.VariableNames)
    if strcmp(IK.Properties.VariableNames{n}(end-1:end), '_r') == 0
        LikCols(n) = 1;
    end
end
LikCols = logical(LikCols);
for i = 1:length(L_GC_t)-1 % cut and resample for each gait cycle
    [~,A] = min(abs(Time - L_GC_t(i))); % identify gait cycles
    [~,B] = min(abs(Time - L_GC_t(i+1)));
    
    % parse gait cycle to 100 points
    LGC(i).IK = interp1(Time(A:B), table2array(IK(A:B,LikCols)), linspace(Time(A), Time(B)));
    LGC(i).IKnames = IK.Properties.VariableNames(LikCols); % save variable names
end


% RIGHT
RikCols = zeros(1, length(IK.Properties.VariableNames));
for n = 1:length(IK.Properties.VariableNames)
    if strcmp(IK.Properties.VariableNames{n}(end-1:end), '_l') == 0
        RikCols(n) = 1;
    end
end
RikCols = logical(RikCols);
for i = 1:length(R_GC_t)-1 % cut and resample for each gait cycle
    [~,A] = min(abs(Time - R_GC_t(i))); % identify gait cycles
    [~,B] = min(abs(Time - R_GC_t(i+1)));
    
    % parse gait cycle to 100 points
    RGC(i).IK = interp1(Time(A:B), table2array(IK(A:B,RikCols)), linspace(Time(A), Time(B)));
    RGC(i).IKnames = IK.Properties.VariableNames(RikCols); % save variable names
end

%% kinetics
IncludeKinetics = 'No';
if strcmp(IncludeKinetics, 'Yes')
    % load
    ID = Osim.readSTO([SubjPath 'ID_Files\' TName '_ID.sto']);
    Time = ID.Header; % get time series
    
    % LEFT
    LidCols = zeros(1, length(ID.Properties.VariableNames));
    for n = 1:length(ID.Properties.VariableNames)
        if contains(ID.Properties.VariableNames{n}, '_r_moment') == 0
            LidCols(n) = 1;
        end
    end
    LidCols = logical(LidCols);
    for i = 1:length(L_GC_t)-1 % cut and resample for each gait cycle
        [~,A] = min(abs(Time - L_GC_t(i))); % identify gait cycles
        [~,B] = min(abs(Time - L_GC_t(i+1)));
        
        % parse gait cycle to 100 points
        LGC(i).ID = interp1(Time(A:B), table2array(ID(A:B,LidCols)), linspace(Time(A), Time(B)));
        LGC(i).IDnames = ID.Properties.VariableNames(LidCols); % save variable names
    end
    
    % RIGHT
    RidCols = zeros(1, length(ID.Properties.VariableNames));
    for n = 1:length(ID.Properties.VariableNames)
        if contains(ID.Properties.VariableNames{n}, '_l_moment') == 0
            RidCols(n) = 1;
        end
    end
    RidCols = logical(RidCols);
    for i = 1:length(R_GC_t)-1 % cut and resample for each gait cycle
        [~,A] = min(abs(Time - R_GC_t(i))); % identify gait cycles
        [~,B] = min(abs(Time - R_GC_t(i+1)));
        
        % parse gait cycle to 100 points
        RGC(i).ID = interp1(Time(A:B), table2array(ID(A:B,RidCols)), linspace(Time(A), Time(B)));
        RGC(i).IDnames = ID.Properties.VariableNames(RidCols); % save variable names
    end
end

%% Create 3D matrix, average, and standard dev of all variables

% LEFT
for i = 1:length(LGC) % include all measures into one table for error measurement
    if strcmp(IncludeKinetics, 'Yes')
        GC.L.All(:,:,i) = [LGC(i).GRFs, LGC(i).IK, LGC(i).ID];
    else
        GC.L.All(:,:,i) = [LGC(i).GRFs, LGC(i).IK];
    end
end
GC.L.Avg = mean(GC.L.All, 3); % get average
GC.L.Std = std(GC.L.All, 0, 3);  % get standard deviation
GC.L.Diff = GC.L.All - GC.L.Avg;
GC.L.Diff2 = GC.L.Diff .^2;
GC.L.SumDcols = sum(GC.L.Diff2, 1); % sum of squares
GC.L.SumD = sum(GC.L.SumDcols, 2); % sum of squares
if strcmp(IncludeKinetics, 'Yes')
    GC.L.ColNames = [LGC(i).GRFnames, LGC(i).IKnames, LGC(i).IDnames];
else
    GC.L.ColNames = [LGC(i).GRFnames, LGC(i).IKnames];
end
[GC.L.SortVal, GC.L.SortInd] = sort(reshape(GC.L.SumD,[length(LGC),1])); % find strides with minimal error

% RIGHT
for i = 1:length(RGC) % include all measures into one table for error measurement
    if strcmp(IncludeKinetics, 'Yes')
        GC.R.All(:,:,i) = [RGC(i).GRFs, RGC(i).IK, RGC(i).ID];
    else
        GC.R.All(:,:,i) = [RGC(i).GRFs, RGC(i).IK];
    end
end
GC.R.Avg = mean(GC.R.All, 3); % get average
GC.R.Std = std(GC.R.All, 0, 3);  % get standard deviation
GC.R.Diff = GC.R.All - GC.R.Avg;
GC.R.Diff2 = GC.R.Diff .^2;
GC.R.SumDcols = sum(GC.R.Diff2, 1); % sum of squares
GC.R.SumD = sum(GC.R.SumDcols, 2); % sum of squares
if strcmp(IncludeKinetics, 'Yes')
    GC.R.ColNames = [RGC(i).GRFnames, RGC(i).IKnames, RGC(i).IDnames];
else
    GC.R.ColNames = [RGC(i).GRFnames, RGC(i).IKnames];
end
[GC.R.SortVal, GC.R.SortInd] = sort(reshape(GC.R.SumD,[length(RGC),1])); % find strides with minimal error

%% Exclude crossover steps
if Trials.Cross.L_Num>0
    ind = GC.L.SortInd == Trials.Cross.L_Stride;
    Ind = sum(ind, 2) > 0;
    GC.L.SortVal(Ind) = [];
    GC.L.SortInd(Ind) = [];
end
if Trials.Cross.R_Num>0
    ind = GC.R.SortInd == Trials.Cross.R_Stride;
    Ind = sum(ind, 2) > 0;
    GC.R.SortVal(Ind) = [];
    GC.R.SortInd(Ind) = [];
end

% save times of ranked GC into structure
GC.L.RankTimes = [Trials.TSData.L_Strike(GC.L.SortInd, 2), Trials.TSData.L_Strike(GC.L.SortInd+1, 2)];
GC.R.RankTimes = [Trials.TSData.R_Strike(GC.R.SortInd, 2), Trials.TSData.R_Strike(GC.R.SortInd+1, 2)];

%% Plot all gait cycles and highest ranked
PlotGCRanks = 'No';
if strcmp(PlotGCRanks, 'Yes')
    % a = 0.25;
    LW = 0.25;
    C = 2;
    
    % LEFT
    n = GC.L.SortInd;
    h = figure('Position',[100 100 1000 800]);
    % GRFs
    subplot(331); hold on;
    col = contains(GC.L.ColNames, 'ground_force_vx');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    ylabel('GRFs (N)');
    
    subplot(332); hold on;
    col = contains(GC.L.ColNames, 'ground_force_vy');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    
    subplot(333); hold on;
    col = contains(GC.L.ColNames, 'ground_force_vz');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    
    % Kinematics
    subplot(334); hold on;
    col = strcmp(GC.L.ColNames, 'hip_flexion_l');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    ylabel('Kinematics (deg)');
    
    subplot(335); hold on;
    col = strcmp(GC.L.ColNames, 'knee_angle_l');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    
    subplot(336); hold on;
    col = strcmp(GC.L.ColNames, 'ankle_angle_l');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    
    % Kinetics
    subplot(337); hold on;
    col = strcmp(GC.L.ColNames, 'hip_flexion_l_moment');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    ylabel('Kinetics (Nm)');
    
    subplot(338); hold on;
    col = strcmp(GC.L.ColNames, 'knee_angle_l_moment');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    
    subplot(339); hold on;
    col = strcmp(GC.L.ColNames, 'ankle_angle_l_moment');
    for i = 1:length(n)
        plot(GC.L.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.L.All(:,col, GC.L.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.L.All(:,col, GC.L.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.L.ColNames{col})
    
    saveas(h, [SubjPath 'Figures\' TName 'AllGCs_L.png']);
    
    
    % RIGHT
    n = GC.R.SortInd;
    h = figure('Position',[100 100 1000 800]);
    % GRFs
    subplot(331); hold on;
    col = contains(GC.R.ColNames, 'ground_force_vx');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    ylabel('GRFs (N)');
    
    subplot(332); hold on;
    col = contains(GC.R.ColNames, 'ground_force_vy');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    
    subplot(333); hold on;
    col = contains(GC.R.ColNames, 'ground_force_vz');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    
    % Kinematics
    subplot(334); hold on;
    col = strcmp(GC.R.ColNames, 'hip_flexion_r');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    ylabel('Kinematics (deg)');
    
    subplot(335); hold on;
    col = strcmp(GC.R.ColNames, 'knee_angle_r');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    
    subplot(336); hold on;
    col = strcmp(GC.R.ColNames, 'ankle_angle_r');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    
    % Kinetics
    subplot(337); hold on;
    col = strcmp(GC.R.ColNames, 'hip_flexion_r_moment');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    ylabel('Kinetics (Nm)');
    
    subplot(338); hold on;
    col = strcmp(GC.R.ColNames, 'knee_angle_r_moment');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    
    subplot(339); hold on;
    col = strcmp(GC.R.ColNames, 'ankle_angle_r_moment');
    for i = 1:length(n)
        plot(GC.R.All(:,col, n(i)), '-k', 'LineWidth', LW);
    end
    plot(GC.R.All(:,col, GC.R.SortInd(3)), '-r', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(2)), '-b', 'LineWidth', C);
    plot(GC.R.All(:,col, GC.R.SortInd(1)), '-g', 'LineWidth', C);
    title(GC.R.ColNames{col})
    
    saveas(h, [SubjPath 'Figures\' TName 'AllGCs_R.png']);
    close all;
end

%% save all data into a subject specific mat file


end

