function [Trials] = CrossoverAnalysis(Trials, Index, SubjMass, Settings, CrossoverFig)

% Identify Crossover steps by analyzing vGRFs from treadmill walking

% Excludes Crossover Steps that meet any of the following criteria: 

%       loading phase < 9*subject mass (lowered floor from gravity*subject mass
%       time loading peak > 3 standard deviations from the mean
%       rms error > 3 standard deviations from the mean vGRF curve
%       foot off timing > 3 standard deviations from mean time

% option to plot all crossover and clean steps (and save a PNG file)


% Ricky Pimentel        February 2020
% Applied Biomechanics Laboratory - UNC Chapel Hill

%% create crossover figure if it doesnt exist
% PlotCrossover = 'No'; 
% if strcmp(PlotCrossover, 'Yes')
    if exist('CrossoverFig', 'var') == 0
        CrossoverFig = figure('Position',[100 100 1200 800]);
    end
% end

%% get average peak vGRFs profile and identify outliers
for i = find(Index)
    clearvars LPeaks RPeaks Cross
    
    % load GRF data if the field doesnt exist in Trials structure
    % or if the structure is empty for Trial of interest
    if ~isfield(Trials(i), 'GRF')
        Trials(i).GRF = LoadGRF(Trials(i).files.OpenSimGRF, 0);
    elseif isempty(Trials(i).GRF)
        Trials(i).GRF = LoadGRF(Trials(i).files.OpenSimGRF, 0);
    end
    
    % get temporal spatial data
    Trials(i).TSData = TreadmillTempSpatData(Trials(i).GRF);
    
    VgrfCol = contains(Trials(i).GRF.ColHeaders, '_vy');% identify vGRF columns
    Trials(i).GRF.Vgrf = Trials(i).GRF.AllData(:,VgrfCol);
    
    %% plot vGRFs from all steps

    % right
    subplot(222); hold on;
    for j = 1:Trials(i).TSData.R_NumStrikes - 1
        Start = Trials(i).TSData.R_Strike(j,1);
        End = Trials(i).TSData.R_Strike(j+1,1);
        CycDur = End - Start + 1;
        Trials(i).GRF.R_vGRF_cyc(:,:,j) = resample(Trials(i).GRF.Vgrf(Start:End,1), 100, CycDur);
        plot(Trials(i).GRF.R_vGRF_cyc(:,:,j), '-k');
        title(strcat(Trials(i).name, 'All Right vGRFs'));
        [RPeaks.Val(j,:), RPeaks.Loc(j,:)] = findpeaks(Trials(i).GRF.R_vGRF_cyc(:,:,j),...
            'NPeaks', 2);
        
        
        plot(RPeaks.Loc(j,:), RPeaks.Val(j,:), 'og', 'MarkerSize', 12); % circle the first peak
    end
    hline(SubjMass * 9.81, '--k');
    hline(SubjMass * 8, '--b');
    
    % left
    subplot(221); hold on;
    for j = 1:Trials(i).TSData.L_NumStrikes - 1
        Start = Trials(i).TSData.L_Strike(j,1);
        End = Trials(i).TSData.L_Strike(j+1,1);
        CycDur = End - Start + 1;
        Trials(i).GRF.L_vGRF_cyc(:,:,j) = resample(Trials(i).GRF.Vgrf(Start:End,2), 100, CycDur);
        plot(Trials(i).GRF.L_vGRF_cyc(:,:,j), '-k');
        title(strcat(Trials(i).name, 'All Left vGRFs'));
        [LPeaks.Val(j,:), LPeaks.Loc(j,:)] = findpeaks(Trials(i).GRF.L_vGRF_cyc(:,:,j),...
            'NPeaks', 2);
        plot(LPeaks.Loc(j,:), LPeaks.Val(j,:), 'og', 'MarkerSize', 12); % circle the first peak
    end
    hline(SubjMass * 9.81, '--k');
    hline(SubjMass * 8, '--b');
    
    
    %% analyze vGRF peak values
    % loading vGRF should be greater than subject weight
    WtThresh1 = SubjMass * 9; % approximate using 9 rather than 9.81 m/ss
    % prop vGRF should be present (above 0.5 * weight)
    %      WtThresh2 = SubjMass * 5; % approximate using 5 rather than 9.81 m/ss
    RPeaks.BadLocs = zeros(5, length(RPeaks.Val));
    LPeaks.BadLocs = zeros(5, length(LPeaks.Val));
    for j = 1:length(RPeaks.Val)
        if RPeaks.Val(j,1) < WtThresh1
            RPeaks.BadLocs(5,j) = 1; % right peaks less that weight threshold
            RPeaks.Loc(j,:) = NaN; % remove that location from timing analysis
            RPeaks.Val(j,:) = NaN; % remove that location from timing analysis
        end
    end
    for j = 1:length(LPeaks.Val)
        if LPeaks.Val(j,1) < WtThresh1
            LPeaks.BadLocs(5,j) = 1; % right peaks less that weight threshold
            LPeaks.Loc(j,:) = NaN; % remove that location from timing analysis
            LPeaks.Val(j,:) = NaN; % remove that location from timing analysis
        end
    end
    
    % Peak timing analysis
    Var = 3; % 2.56; % set variance threshold for identifying off-timing of peak vGRF (in number of SDs)
    % RIGHT side
    % vGRF first peak location
    RPeaks.AvgLoc1 = nanmean(RPeaks.Loc(:,1)); % average peak location (in time)
    RPeaks.StdLoc1 = nanstd(RPeaks.Loc(:,1)); % std in peak location (in time)
    RPeaks.Loc_UB1 = RPeaks.AvgLoc1 + Var .* RPeaks.StdLoc1; % create upper bound
    RPeaks.Loc_LB1 = RPeaks.AvgLoc1 - Var .* RPeaks.StdLoc1; % create lower bound
    RPeaks.BadLocs(1,:) = RPeaks.Loc(:,1) > RPeaks.Loc_UB1;
    RPeaks.BadLocs(2,:) = RPeaks.Loc(:,1) < RPeaks.Loc_LB1;
    % vGRF second peak location
    %     RPeaks.AvgLoc2 = nanmean(RPeaks.Loc(:,2)); % average peak location (in time)
    %     RPeaks.StdLoc2 = nanstd(RPeaks.Loc(:,2)); % std in peak location (in time)
    %     RPeaks.Loc_UB2 = RPeaks.AvgLoc2 + Var .* RPeaks.StdLoc2; % create upper bound
    %     RPeaks.Loc_LB2 = RPeaks.AvgLoc2 - Var .* RPeaks.StdLoc2; % create lower bound
    %     RPeaks.BadLocs(3,:) = RPeaks.Loc(:,2) > RPeaks.Loc_UB2;
    %     RPeaks.BadLocs(4,:) = RPeaks.Loc(:,2) < RPeaks.Loc_LB2;
    
    RPeaks.BADLOCs = sum(RPeaks.BadLocs, 1) > 0; % determine overall crossover steps
    
    % LEFT side
    % vGRF first peak location
    LPeaks.AvgLoc1 = mean(LPeaks.Loc(:,1)); % average peak location (in time)
    LPeaks.StdLoc1 = std(LPeaks.Loc(:,1)); % std in peak location (in time)
    LPeaks.Loc_UB1 = LPeaks.AvgLoc1 + Var .* LPeaks.StdLoc1; % create upper bound
    LPeaks.Loc_LB1 = LPeaks.AvgLoc1 - Var .* LPeaks.StdLoc1; % create lower bound
    LPeaks.BadLocs(1,:) = LPeaks.Loc(:,1) > LPeaks.Loc_UB1;
    LPeaks.BadLocs(2,:) = LPeaks.Loc(:,1) < LPeaks.Loc_LB1;
    % vGRF second peak location
    %     LPeaks.AvgLoc2 = mean(LPeaks.Loc(:,2)); % average peak location (in time)
    %     LPeaks.StdLoc2 = std(LPeaks.Loc(:,2)); % std in peak location (in time)
    %     LPeaks.Loc_UB2 = LPeaks.AvgLoc2 + Var .* LPeaks.StdLoc2; % create upper bound
    %     LPeaks.Loc_LB2 = LPeaks.AvgLoc2 - Var .* LPeaks.StdLoc2; % create lower bound
    %     LPeaks.BadLocs(3,:) = LPeaks.Loc(:,2) > LPeaks.Loc_UB2;
    %     LPeaks.BadLocs(4,:) = LPeaks.Loc(:,2) < LPeaks.Loc_LB2;
    
    LPeaks.BADLOCs = sum(LPeaks.BadLocs, 1) > 0; % determine overall crossover steps
    
    % plot the crossover steps - identify with red * at first peak
    %     subplot(222); hold on; % right
    %     plot(RPeaks.Loc(RPeaks.BADLOCs), RPeaks.Val(RPeaks.BADLOCs), '*r', 'MarkerSize', 20);
    %     subplot(221); hold on; % left
    %     plot(LPeaks.Loc(LPeaks.BADLOCs), LPeaks.Val(LPeaks.BADLOCs), '*r', 'MarkerSize', 20);
    
    
    %% exclude using large rms error from average vGRF curve
    rmsThresh = 3; %2.56; % threshold for exclusion (in multiples of standard dev)
    % RIGHT side
    Trials(i).GRF.R_vGRF_avg = mean(...
        Trials(i).GRF.R_vGRF_cyc(:,:,~logical(RPeaks.BADLOCs)), 3);
    Trials(i).GRF.R_vGRF_std = std(...
        Trials(i).GRF.R_vGRF_cyc(:,:,~logical(RPeaks.BADLOCs)), 0, 3);
    Trials(i).GRF.R_vGRF_SD = mean(Trials(i).GRF.R_vGRF_std);
    for j = 1:length(RPeaks.BADLOCs)
        Trials(i).GRF.R_vGRF_rms(:,:,j) = rms(Trials(i).GRF.R_vGRF_avg - Trials(i).GRF.R_vGRF_cyc(:,:,j));
        if Trials(i).GRF.R_vGRF_rms(:,:,j) > rmsThresh * Trials(i).GRF.R_vGRF_SD
            RPeaks.BADLOCs(j) = 1;
        end
    end
    % LEFT side
    Trials(i).GRF.L_vGRF_avg = mean(...
        Trials(i).GRF.L_vGRF_cyc(:,:,~logical(LPeaks.BADLOCs)),3);
    Trials(i).GRF.L_vGRF_std = std(...
        Trials(i).GRF.L_vGRF_cyc(:,:,~logical(LPeaks.BADLOCs)), 0, 3);
    Trials(i).GRF.L_vGRF_SD = mean(Trials(i).GRF.L_vGRF_std);
    for j = 1:length(LPeaks.BADLOCs)
        Trials(i).GRF.L_vGRF_rms(:,:,j) = rms(Trials(i).GRF.L_vGRF_avg - Trials(i).GRF.L_vGRF_cyc(:,:,j));
        if Trials(i).GRF.L_vGRF_rms(:,:,j) > rmsThresh * Trials(i).GRF.L_vGRF_SD
            LPeaks.BADLOCs(j) = 1;
        end
    end
    
    %% analyze timing of foot off to identify any missed crossovers
    Var = 3; % 2.56;
    % RIGHT
    RPeaks.Off = NaN(length(RPeaks.BADLOCs),1);
    for j = 1:length(RPeaks.BADLOCs)
        val = find(flipud(Trials(i).GRF.R_vGRF_cyc(:,:,j))>2, 1) + 1;
        Trials(i).GRF.R_cyc(:,:,j) = zeros(100, 1);
        Trials(i).GRF.R_cyc(1:100-val,:,j) = Trials(i).GRF.R_vGRF_cyc(1:100-val, :, j);
        if RPeaks.BADLOCs(j) == 0
            RPeaks.Off(j) = 100 - val;   % save timing of foot off if not a crossover
        end
    end
    RPeaks.OffAvg = nanmean(RPeaks.Off); % get average foot off timing
    RPeaks.OffSD = nanstd(RPeaks.Off); % standard deviation of foot off
    RPeaks.OffLB = RPeaks.OffAvg - (Var .* RPeaks.OffSD); % create lower and upper bounds based on 3 * standard dev
    RPeaks.OffUB = RPeaks.OffAvg + (Var .* RPeaks.OffSD);
    RPeaks.BADLOCs(RPeaks.Off > RPeaks.OffUB) = 1; % exclude
    RPeaks.BADLOCs(RPeaks.Off < RPeaks.OffLB) = 1;
    clearvars Off cyc val Avg SD LB UB
    
    % LEFT
    LPeaks.Off = NaN(length(LPeaks.BADLOCs),1);
    for j = 1:length(LPeaks.BADLOCs)
        val = find(flipud(Trials(i).GRF.L_vGRF_cyc(:,:,j))>2, 1) + 1;
        Trials(i).GRF.L_cyc(:,:,j) = zeros(100, 1);
        Trials(i).GRF.L_cyc(1:100-val,:,j) = Trials(i).GRF.L_vGRF_cyc(1:100-val, :, j);
        if LPeaks.BADLOCs(j) == 0
            LPeaks.Off(j) = 100 - val;   % save timing of foot off if not a crossover
        end
    end
    LPeaks.OffAvg = nanmean(LPeaks.Off);
    LPeaks.OffSD = nanstd(LPeaks.Off);
    LPeaks.OffLB = LPeaks.OffAvg - (Var .* LPeaks.OffSD);
    LPeaks.OffUB = LPeaks.OffAvg + (Var .* LPeaks.OffSD);
    LPeaks.BADLOCs(LPeaks.Off > LPeaks.OffUB) = 1;
    LPeaks.BADLOCs(LPeaks.Off < LPeaks.OffLB) = 1;
    clearvars Off cyc val Avg SD LB UB
    
    %% save Bad locations as crossover steps
    LPeaks.Crossovers = LPeaks.BADLOCs;
    RPeaks.Crossovers = RPeaks.BADLOCs;
    
    %% plot all clean (non-crossover) steps
    % RIGHT
    for j = 1:length(~logical(RPeaks.Crossovers))
        subplot(224); hold on;
        if RPeaks.Crossovers(j) == 0
            grfs = Trials(i).GRF.R_vGRF_cyc(:,:,j);
            plot(grfs, '-k');
        elseif RPeaks.Crossovers(j) == 1
            subplot(222); hold on;
            grfs = Trials(i).GRF.R_vGRF_cyc(:,:,j);
            plot(grfs, '-r');
        end
    end
    hline(SubjMass * 9.81, '--k');
    title('Clean vGRFs - Right');
    ylabel('N'); xlabel('% of Gait Cycle');
    
    % LEFT
    for j = 1:length(~logical(LPeaks.Crossovers))
        subplot(223); hold on;
        if LPeaks.Crossovers(j) == 0
            grfs = Trials(i).GRF.L_vGRF_cyc(:,:,j);
            plot(grfs, '-k');
        elseif LPeaks.Crossovers(j) == 1
            subplot(221); hold on;
            grfs = Trials(i).GRF.L_vGRF_cyc(:,:,j);
            plot(grfs, '-r');
        end
    end
    hline(SubjMass * 9.81, '--k');
    title('Clean vGRFs - Left');
    ylabel('N'); xlabel('% of Gait Cycle');
    
    %% save figure as PNG
    if strcmp(Settings.SaveCrossoverPlot, 'Yes')
        supertitle(strcat(Trials(i).subject, '-', Trials(i).name)); 
        % define folder location to save
        FigFileName = strcat(Trials(i).folder, '\', Trials(i).subject, '-Crossovers-',Trials(i).name, '.png'); 
        saveas(CrossoverFig, FigFileName); 
    end
    clf;
    
    %% save crossover times and strides in Cross field within Trials structure
    % LEFT
    Cross.L_Num = sum(LPeaks.Crossovers);
    Cross.L_Stride = find(LPeaks.Crossovers);
    Cross.L_Times = [Trials(i).TSData.L_Strike(Cross.L_Stride, 2), Trials(i).TSData.L_Strike(Cross.L_Stride+1, 2)];
    if strcmp(Settings.PrintResults, 'Yes') % print results if desired
        disp(strcat('Indentified _', num2str(Cross.L_Num), ' LEFT crossover step(s) in trial _', Trials(i).name));
        disp('Crossover steps ocurred at times:');
        for k = 1:Cross.L_Num
            disp(Cross.L_Times(k,:));
        end
    end
    % RIGHT
    Cross.R_Num = sum(RPeaks.Crossovers);
    Cross.R_Stride = find(RPeaks.Crossovers);
    Cross.R_Times = [Trials(i).TSData.R_Strike(Cross.R_Stride, 2), Trials(i).TSData.R_Strike(Cross.R_Stride+1, 2)];
    if strcmp(Settings.PrintResults, 'Yes') % print results if desired
        disp(strcat('Indentified _', num2str(Cross.R_Num), ' RIGHT crossover step(s) in trial _', Trials(i).name));
        disp('Crossover steps ocurred at times:');
        for k = 1:Cross.R_Num
            disp(Cross.R_Times(k,:));
        end
    end
    
    % save crossover data in Trials structure
    Cross.LPeaks = LPeaks;
    Cross.RPeaks = RPeaks;
    Trials(i).Cross = Cross;
    
    clearvars Start End CycDur i j Var VgrfCol 
end


end


%% apply known crossovers to contralateral step before and after
% for every crossover step, identity the other side's
% preceeding and following step also as crossovers
%     NumSteps = min([Trials(i).TSData.L_NumStrikes Trials(i).TSData.R_NumStrikes]);
%     [~, ind] = min([Trials(i).TSData.L_Strike(1,1) Trials(i).TSData.R_Strike(1,1)]);
%     LPeaks.Crossovers = zeros(length(LPeaks.BADLOCs),1); % initialize crossover index
%     RPeaks.Crossovers = zeros(length(RPeaks.BADLOCs),1);

%     if ind == 1 % if left strike starts first
%         for j = 1:NumSteps % loop through all viable steps
%             if LPeaks.BADLOCs(j) == 1
%                 LPeaks.Crossovers(j) = 1;
% %                 RPeaks.Crossovers(j) = 1;
% %                 RPeaks.Crossovers(j-1) = 1;
%             end
%             if RPeaks.BADLOCs(j) == 1
%                 RPeaks.Crossovers(j) = 1;
% %                 LPeaks.Crossovers(j) = 1;
% %                 LPeaks.Crossovers(j+1) = 1;
%             end
%         end
%     elseif ind == 2 % if right strike starts first
%         for j = 1:NumSteps % loop through all viable steps
%             if LPeaks.BADLOCs(j) == 1
%                 LPeaks.Crossovers(j) = 1;
% %                 RPeaks.Crossovers(j) = 1;
% %                 RPeaks.Crossovers(j+1) = 1;
%             end
%             if RPeaks.BADLOCs(j) == 1
%                 RPeaks.Crossovers(j) = 1;
% %                 LPeaks.Crossovers(j) = 1;
% %                 LPeaks.Crossovers(j-1) = 1;
%             end
%         end
%     end

%% more exlusions if desired
%     reply = input('Exclude more crossover steps in this trial? Y/N [Y]:','s');
%     if isempty(reply)
%         reply = 'N';
%     end
%
%      rmsThresh = 2; % threshold for exclusion (in multiples of standard dev)
%     while strcmp(reply, 'Y') % if want to exclude more crossover steps
%         clf;
%         % exclude using large rms error from average vGRF curve
%         % RIGHT side
%         Trials(i).GRF.R_vGRF_avg = mean(...
%             Trials(i).GRF.R_vGRF_cyc(:,:,~logical(RPeaks.Crossovers)), 3);
%         Trials(i).GRF.R_vGRF_std = std(...
%             Trials(i).GRF.R_vGRF_cyc(:,:,~logical(RPeaks.Crossovers)), 0, 3);
%         Trials(i).GRF.R_vGRF_SD = mean(Trials(i).GRF.R_vGRF_std);
%         for j = 1:length(RPeaks.Crossovers)
%             Trials(i).GRF.R_vGRF_rms(:,:,j) = rms(Trials(i).GRF.R_vGRF_avg - Trials(i).GRF.R_vGRF_cyc(:,:,j));
%             if Trials(i).GRF.R_vGRF_rms(:,:,j) > rmsThresh * Trials(i).GRF.R_vGRF_SD
%                 RPeaks.Crossovers(j) = 1;
%             end
%         end
%         % LEFT side
%         Trials(i).GRF.L_vGRF_avg = mean(...
%             Trials(i).GRF.L_vGRF_cyc(:,:,~logical(LPeaks.Crossovers)),3);
%         Trials(i).GRF.L_vGRF_std = std(...
%             Trials(i).GRF.L_vGRF_cyc(:,:,~logical(LPeaks.Crossovers)), 0, 3);
%         Trials(i).GRF.L_vGRF_SD = mean(Trials(i).GRF.L_vGRF_std);
%         for j = 1:length(LPeaks.Crossovers)
%             Trials(i).GRF.L_vGRF_rms(:,:,j) = rms(Trials(i).GRF.L_vGRF_avg - Trials(i).GRF.L_vGRF_cyc(:,:,j));
%             if Trials(i).GRF.L_vGRF_rms(:,:,j) > rmsThresh * Trials(i).GRF.L_vGRF_SD
%                 LPeaks.Crossovers(j) = 1;
%             end
%         end
%
%         % re-plot
%         % RIGHT all vGRFs
%         subplot(222); hold on;
%         for j = 1:Trials(i).TSData.R_NumStrikes - 1
%             plot(Trials(i).GRF.R_vGRF_cyc(:,:,j), '-k');
%             title(strcat(Trials(i).name, 'All Right vGRFs'));
%         end
%         hline(SubjMass * 9.81, '--k');
%         ylabel('N'); xlabel('% of Gait Cycle');
%         % left
%         subplot(221); hold on;
%         for j = 1:Trials(i).TSData.L_NumStrikes - 1
%             plot(Trials(i).GRF.L_vGRF_cyc(:,:,j), '-k');
%             title(strcat(Trials(i).name, 'All Left vGRFs'));
%         end
%         hline(SubjMass * 9.81, '--k');
%         ylabel('N'); xlabel('% of Gait Cycle');
%
%         % plot all clean steps (non-crossover)
%         subplot(224); hold on; % right
%         for j = 1:length(~logical(RPeaks.Crossovers))
%             if RPeaks.Crossovers(j) == 0
%                 grfs = Trials(i).GRF.R_vGRF_cyc(:,:,j);
%                 plot(grfs, '-k');
%                 title('Clean vGRFs - Right');
%                 ylabel('N'); xlabel('% of Gait Cycle');
%             end
%         end
%
%         subplot(223); hold on; % left
%         for j = 1:length(~logical(LPeaks.Crossovers))
%             if LPeaks.Crossovers(j) == 0
%                 grfs = Trials(i).GRF.L_vGRF_cyc(:,:,j);
%                 plot(grfs, '-k');
%                 title('Clean vGRFs - Right');
%                 ylabel('N'); xlabel('% of Gait Cycle');
%             end
%         end
%
%         reply = input('Exclude more crossovers again? Y/N [Y]:','s');
%         if isempty(reply)
%             reply = 'N';
%         elseif strcmp(reply, 'Y')
%             rmsThresh = rmsThresh - 0.5; % reduce threshold for identifying crossover steps
%         end
%     end % end while loop if reply == N