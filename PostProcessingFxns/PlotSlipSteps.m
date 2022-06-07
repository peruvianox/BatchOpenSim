

%% Create figure
close all; 
figure('Position',[100 100 800 600]); 
LW = 2;
addpath(genpath('C:\Users\richa\Documents\Packages\MoCapTools'));

%% plot norm ankle moment
NI = 3; 
A = 0.2; 
% left side = subplot 1
subplot(211); hold on; 
ColNames = Subjects.Trials(NI).GC.L.ColNames; 
Ind = contains(ColNames, 'ankle_angle_l_moment'); 
step = Subjects.Trials(NI).GC.L.SortInd(1); 
D = -Subjects.Trials(NI).GC.L.All(:,Ind, step); 
T = D + Subjects.Trials(NI).GC.L.Std(:,Ind);
B = D - Subjects.Trials(NI).GC.L.Std(:,Ind);
x = 1:100;
fill([x fliplr(x)], [T' fliplr(B')], 'k', 'FaceAlpha', A, 'EdgeAlpha', 0); 
L1 = plot(D, '-k', 'LineWidth', LW); 
N1 = 'Norm';

% right side = subplot 2
subplot(212); hold on; 
ColNames = Subjects.Trials(NI).GC.R.ColNames; 
Ind = contains(ColNames, 'ankle_angle_r_moment'); 
D = -Subjects.Trials(NI).GC.R.Avg(:,Ind); 
T = D + Subjects.Trials(NI).GC.L.Std(:,Ind);
B = D - Subjects.Trials(NI).GC.L.Std(:,Ind);
x = 1:100;
fill([x fliplr(x)], [T' fliplr(B')], 'k', 'FaceAlpha', A, 'EdgeAlpha', 0); 
plot(D, '-k', 'LineWidth', LW); 


%% Load slip data
load('D:\UNC_ABL\SlipProcess\YA10SlipHS.mat'); % load slip times
SlipID = Osim.readSTO('D:\UNC_ABL\SlipProcess\YA10\OpenSim\ID_Files\Slip_1_ID.sto'); % load slip inverse dynamics trial

%% plot slip ankle moment

% get frame and convert to time by dividing by 100
SlipTimes.Left = SlipHS.Left / 100;
SlipTimes.Right = SlipHS.Right / 100;
SI = 4; 
LW2 = 1;
Slips2plot = 5; 
% left - loop through all slips
subplot(211); 
for i = 1:Slips2plot
    % get steps after slip pert
    TI = find(Subjects.Trials(SI).TSData.L_Strike(:,2) > SlipTimes.Left(i), 1);
    SlipStart = Subjects.Trials(SI).TSData.L_Strike(TI,2);
    SlipEnd = Subjects.Trials(SI).TSData.L_Strike(TI + 1,2);
    Slip2Start = Subjects.Trials(SI).TSData.L_Strike(TI + 1,2);
    Slip2End = Subjects.Trials(SI).TSData.L_Strike(TI + 2,2);
    
    OppTI = find(Subjects.Trials(SI).TSData.R_Strike(:,2) > SlipTimes.Left(i), 1);
    OppSlipStart = Subjects.Trials(SI).TSData.R_Strike(OppTI,2);
    OppSlipEnd = Subjects.Trials(SI).TSData.R_Strike(OppTI + 1,2);
    OppSlip2Start = Subjects.Trials(SI).TSData.R_Strike(OppTI + 1,2);
    OppSlip2End = Subjects.Trials(SI).TSData.R_Strike(OppTI + 2,2);
    
    % plot slip step ankle moment
    Slip1 = parseID(SlipID, SlipStart, SlipEnd);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_l_moment');
    D = Slip1(:,Ind);
    L2 = plot(-D, 'r-', 'LineWidth', LW2);
    N2 = 'Slip';
    
    % plot opposite slip step
    OppSlip1 = parseID(SlipID, OppSlipStart, OppSlipEnd);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_r_moment');
    D = OppSlip1(:,Ind);
    L3 = plot(-D, 'b--', 'LineWidth', LW2);
    N3 = 'Opp Slip';
    
    % plot recover step ankle moment
    Slip2 = parseID(SlipID, Slip2Start, Slip2End);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_l_moment');
    D = Slip2(:,Ind);
    L4 = plot(-D, 'm-', 'LineWidth', LW2);
    N4 = 'Recovery';
    
    % plot opposite recovery step
    OppSlip2 = parseID(SlipID, OppSlip2Start, OppSlip2End);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_r_moment');
    D = OppSlip2(:,Ind);
    L5 = plot(-D, 'c--', 'LineWidth', LW2);
    N5 = 'Opp Recovery';
end

%%
% right
subplot(212); 
for i = 1:Slips2plot
    % get steps after slip pert
    TI = find(Subjects.Trials(SI).TSData.R_Strike(:,2) > SlipTimes.Right(i), 1);
    SlipStart = Subjects.Trials(SI).TSData.R_Strike(TI,2);
    SlipEnd = Subjects.Trials(SI).TSData.R_Strike(TI + 1,2);
    Slip2Start = Subjects.Trials(SI).TSData.R_Strike(TI + 1,2);
    Slip2End = Subjects.Trials(SI).TSData.R_Strike(TI + 2,2);
    
    OppTI = find(Subjects.Trials(SI).TSData.L_Strike(:,2) > SlipTimes.Right(i), 1);
    OppSlipStart = Subjects.Trials(SI).TSData.L_Strike(OppTI,2);
    OppSlipEnd = Subjects.Trials(SI).TSData.L_Strike(OppTI + 1,2);
    OppSlip2Start = Subjects.Trials(SI).TSData.L_Strike(OppTI + 1,2);
    OppSlip2End = Subjects.Trials(SI).TSData.L_Strike(OppTI + 2,2);
    
    % plot slip step ankle moment
    Slip1 = parseID(SlipID, SlipStart, SlipEnd);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_r_moment');
    D = Slip1(:,Ind);
    L2 = plot(-D, 'r-', 'LineWidth', LW2);
    N2 = 'Slip';
    
    % plot opposite slip step
    OppSlip1 = parseID(SlipID, OppSlipStart, OppSlipEnd);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_l_moment');
    D = OppSlip1(:,Ind);
    L3 = plot(-D, 'b--', 'LineWidth', LW2);
    N3 = 'Opp Slip';
    
    % plot recover step ankle moment
    Slip2 = parseID(SlipID, Slip2Start, Slip2End);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_r_moment');
    D = Slip2(:,Ind);
    L4 = plot(-D, 'm-', 'LineWidth', LW2);
    N4 = 'Recovery';
    
    % plot opposite recovery step
    OppSlip2 = parseID(SlipID, OppSlip2Start, OppSlip2End);
    Ind = contains(SlipID.Properties.VariableNames, 'ankle_angle_l_moment');
    D = OppSlip2(:,Ind);
    L5 = plot(-D, 'c--', 'LineWidth', LW2);
    N5 = 'Opp Recovery';
end

%% format axes
subplot(211); 
title('Left Slip Ankle Moment - Young Adult')
legend([L1, L2, L3, L4, L5], {N1, N2, N3, N4, N5}); 
ylabel('Ankle Moment (N)'); 
ylim([-50 200]); 

subplot(212); 
title('Right Slip Ankle Moment - Young Adult')
xlabel('% Gait Cycle'); 
ylabel('Ankle Moment (N)'); 
ylim([-50 200]); 

%% define functions
