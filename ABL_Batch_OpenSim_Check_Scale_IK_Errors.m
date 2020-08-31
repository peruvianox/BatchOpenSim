% function[ ] = ABL_OpenSim_Scale(Settings)

% Batch process scaling 


%% Settings
% cd('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\BatchOpenSim');
% clear; clc; warning off;
% Settings.TrialWindows = 'Yes'; % set to 'No' to run IK for entire duration of each trial
% Settings.NumStrides = [1]; % number of strides to check for quality data if results don't pass quality control
% Settings.LockModel = 'Yes'; % set to 'Yes' to use models that have locked subtalar and MTP joints
% % Settings.LockModel = 'Full'; % set to Full to use models that have locked subtalar, MTP, and lumbar joints
% Settings.ModelStrength = 'Normal'; % set to 'Strong' to use model with 1.5x strength
% Settings.LowPassKinematicsFilter = 6; % define all kinematics filtering to a specified value in Hz
% 
% % Pull in the modeling classes straight from the OpenSim distribution
% import org.opensim.modeling.*
% import java.io.*
% 
% % add geometry files to path - MAY NEED TO SET UP USERS SPECIFIC PATH
% Geopath = 'C:\Users\richa\Documents\OpenSim 4.0\Geometry'; % set to location of geometry files
% ModelVisualizer.addDirToGeometrySearchPaths(Geopath); % add geometry files to path
% 
% % add BatchOpenSim Functions to path
% addpath(genpath('Functions'));

%% move to directory where this subject's files are kept
D = dir; % get current path
CurrFolder = D(1).folder;

GenericFilePath = strcat(CurrFolder, '\OpenSimProcessingFiles'); % set generic folder path
addpath(genpath(GenericFilePath)); % add generic folder to path
GenericDir = dir(GenericFilePath);

% select subject directory
subjectPath = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\exp';
% subjectPath = 'D:\Billys Study';
if ~exist(subjectPath, 'dir')
    subjectPath = uigetdir(CurrFolder, 'Select Folder Containing Subject Data');
end
addpath(genpath(subjectPath));
load('Subjects.mat'); % load subjects from conversion batch

%% load subject anthropometrics
% set your demographics file here
Demo.File = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Metabolics_Feedback_Demographics';
[~,~,Demo.Raw] = xlsread(Demo.File);
Demo.AgeCol = strcmp(Demo.Raw(1,:), 'Age');
Demo.HtCol = strcmp(Demo.Raw(1,:), 'Height (cm)');
Demo.GenderCol = strcmp(Demo.Raw(1,:), 'Gender');
% Demo.MassCol = strcmp(Demo.Raw(1,:), 'Mass'); % dont overwrite mass from static trial

%% loop through subjects
for subj = 1%:length(Subjects)
    
    % set names and demographics
    SubjRow = find(strcmp(Demo.Raw(:,1), Subjects(subj).name),1);
    Subjects(subj).Demo.age = Demo.Raw{SubjRow, Demo.AgeCol};
    Subjects(subj).Demo.height = Demo.Raw{SubjRow, Demo.HtCol} * 10; % height from cm to mm
    
    % identify trial types
    StaticTrials = logical(strcmp({Subjects(subj).Trials.type}, 'static'));
    DynamicTrials = ~StaticTrials;
    
    % get static trial with hip joint center
    if contains(Subjects(subj).Trials(StaticTrials).name, 'nohjc')
        D = dir(Subjects(subj).Trials(1).folder);
        A = contains({D.name}, {'static', 'Static'});
        B = contains({D.name}, 'trc');
        C =  ~contains({D.name}, 'nohjc');
       
       E = logical(A+B + C > 2);
       Subjects(subj).Trials(StaticTrials).name = D(E).name; 
    end
    
    if sum(StaticTrials) > 1
        error('Multiple static trials detected. Reduce to only one.');
    end
    
    clearvars A B C D
    
    %% Calculate virtual markers to use in scaling
    clc; disp('Adding Virtual Markers for scaling');
    % assign markers to pull
    Markers = {'S2','L.ASIS', 'R.ASIS','L.PSIS', 'R.PSIS','L.Knee','R.Knee', 'L.MKnee','R.MKnee', 'L_HJC','R_HJC',...
        'L.Ankle', 'R.Ankle','L.MAnkle','R.MAnkle','L.Heel','R.Heel','L.MT5','R.MT5','L.MT1','R.MT1'};
    
    [Subjects(subj).Trials(StaticTrials).MarkerData] = GetMarkerTrajectories(Subjects(subj).Trials(StaticTrials).TRC.data,...
        Subjects(subj).Trials(StaticTrials).TRC.colheaders, Markers);
    
    PlotVirtual = 'Yes';
    close all;
% markers loaded in subjects structure                
    % pelvis mid HJC
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L_HJC');
    L.HJC = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R_HJC');
    R.HJC = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = L.HJC; 
    A(:,:,2) = R.HJC; 
    Mid.HJC = mean(A, 3); 
    
    % mid ASIS
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.ASIS');
    L.ASIS = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.ASIS');
    R.ASIS = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = L.ASIS;
    A(:,:,2) = R.ASIS;
    Mid.ASIS = mean(A, 3);
    
    % mid PSIS
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.PSIS');
    L.PSIS = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.PSIS');
    R.PSIS = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = L.PSIS;
    A(:,:,2) = R.PSIS;
    Mid.PSIS = mean(A, 3);
    
    % mid pelvis
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'S2');
    Mid.S2 = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = Mid.ASIS;
    A(:,:,2) = Mid.PSIS;
    Mid.Pelvis = mean(A, 3);
    
    % knee joint center
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.Knee');
    L.Knee = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.MKnee');
    L.MKnee = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = L.Knee;
    A(:,:,2) = L.MKnee;
    L.KJC = mean(A, 3);
    
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.Knee');
    R.Knee = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.MKnee');
    R.MKnee = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = R.Knee;
    A(:,:,2) = R.MKnee;
    R.KJC = mean(A, 3);
    
    % ankle joint center
     Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.Ankle');
    L.Ankle = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.MAnkle');
    L.MAnkle = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = L.Ankle;
    A(:,:,2) = L.MAnkle;
    L.AJC = mean(A, 3);
    
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.Ankle');
    R.Ankle = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.MAnkle');
    R.MAnkle = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    A(:,:,1) = R.Ankle;
    A(:,:,2) = R.MAnkle;
    R.AJC = mean(A, 3);
    
    % AJC floor
    L.AJC_Floor = L.AJC;
    L.AJC_Floor(:,2) = 0;
    R.AJC_Floor = R.AJC;
    R.AJC_Floor(:,2) = 0;
    
    % heel floor
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.Heel');
    L.Heel = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    L.Heel_Floor = L.Heel;
    L.Heel_Floor(:,2) = 0;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.Heel');
    R.Heel = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    R.Heel_Floor = R.Heel;
    R.Heel_Floor(:,2) = 0;
    
    % MT1 floor
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.MT1');
    L.MT1 = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    L.MT1_Floor = L.MT1;
    L.MT1_Floor(:,2) = 0;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.MT1');
    R.MT1 = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    R.MT1_Floor = R.MT1;
    R.MT1_Floor(:,2) = 0;
    
    % MT5 floor
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'L.MT5');
    L.MT5 = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    L.MT5_Floor = L.MT5;
    L.MT5_Floor(:,2) = 0;
    Ind = contains({Subjects(subj).Trials(StaticTrials).MarkerData.Name}, 'R.MT5');
    R.MT5 = Subjects(subj).Trials(StaticTrials).MarkerData(Ind).Trajectories;
    R.MT5_Floor = R.MT5;
    R.MT5_Floor(:,2) = 0;
    
    % MidMT floor
    A(:,:,1) = L.MT1_Floor;
    A(:,:,2) = L.MT5_Floor;
    L.MidMT_Floor = mean(A, 3);
    A(:,:,1) = R.MT1_Floor;
    A(:,:,2) = R.MT5_Floor;
    R.MidMT_Floor = mean(A, 3);

    
    % plot all virtual markers to check accuracy
    if strcmp(PlotVirtual, 'Yes')
        a = 3;
        b = 2; 
        c = 1; 
        MkrSz = 20; 
        figure('Position',[50 50 1000 800]); 
        hold on; grid on; 
        
        % hip joint centers
        plot3(L.HJC(1,a), L.HJC(1,b), L.HJC(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(R.HJC(1,a), R.HJC(1,b), R.HJC(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(Mid.HJC(1,a), Mid.HJC(1,b), Mid.HJC(1,c), 'g.', 'MarkerSize', MkrSz);
        % pelvic points
        plot3(L.ASIS(1,a), L.ASIS(1,b), L.ASIS(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(R.ASIS(1,a), R.ASIS(1,b), R.ASIS(1,c), 'r.', 'MarkerSize', MkrSz);
         plot3(L.PSIS(1,a), L.PSIS(1,b), L.PSIS(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(R.PSIS(1,a), R.PSIS(1,b), R.PSIS(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(Mid.ASIS(1,a), Mid.ASIS(1,b), Mid.ASIS(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(Mid.PSIS(1,a), Mid.PSIS(1,b), Mid.PSIS(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(Mid.Pelvis(1,a), Mid.Pelvis(1,b), Mid.Pelvis(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(Mid.S2(1,a), Mid.S2(1,b), Mid.S2(1,c), 'k.', 'MarkerSize', MkrSz);
        % knee and ankle joint centers
        plot3(L.Knee(1,a), L.Knee(1,b), L.Knee(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(L.MKnee(1,a), L.MKnee(1,b), L.MKnee(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(L.KJC(1,a), L.KJC(1,b), L.KJC(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(R.Knee(1,a), R.Knee(1,b), R.Knee(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(R.MKnee(1,a), R.MKnee(1,b), R.MKnee(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(R.KJC(1,a), R.KJC(1,b), R.KJC(1,c), 'g.', 'MarkerSize', MkrSz);
        
        plot3(L.Ankle(1,a), L.Ankle(1,b), L.Ankle(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(L.MAnkle(1,a), L.MAnkle(1,b), L.MAnkle(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(L.AJC(1,a), L.AJC(1,b), L.AJC(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(R.Ankle(1,a), R.Ankle(1,b), R.Ankle(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(R.MAnkle(1,a), R.MAnkle(1,b), R.MAnkle(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(R.AJC(1,a), R.AJC(1,b), R.AJC(1,c), 'g.', 'MarkerSize', MkrSz);
        
        % floor virtual markers
        plot3(L.Heel(1,a), L.Heel(1,b), L.Heel(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(L.Heel_Floor(1,a), L.Heel_Floor(1,b), L.Heel_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(L.AJC_Floor(1,a), L.AJC_Floor(1,b), L.AJC_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(L.MT1(1,a), L.MT1(1,b), L.MT1(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(L.MT1_Floor(1,a), L.MT1_Floor(1,b), L.MT1_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(L.MT5(1,a), L.MT5(1,b), L.MT5(1,c), 'b.', 'MarkerSize', MkrSz);
        plot3(L.MT5_Floor(1,a), L.MT5_Floor(1,b), L.MT5_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(L.MidMT_Floor(1,a), L.MidMT_Floor(1,b), L.MidMT_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        
         plot3(R.Heel(1,a), R.Heel(1,b), R.Heel(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(R.Heel_Floor(1,a), R.Heel_Floor(1,b), R.Heel_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(R.AJC_Floor(1,a), R.AJC_Floor(1,b), R.AJC_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(R.MT1(1,a), R.MT1(1,b), R.MT1(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(R.MT1_Floor(1,a), R.MT1_Floor(1,b), R.MT1_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(R.MT5(1,a), R.MT5(1,b), R.MT5(1,c), 'r.', 'MarkerSize', MkrSz);
        plot3(R.MT5_Floor(1,a), R.MT5_Floor(1,b), R.MT5_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        plot3(R.MidMT_Floor(1,a), R.MidMT_Floor(1,b), R.MidMT_Floor(1,c), 'g.', 'MarkerSize', MkrSz);
        
        axis equal;
    end
    
    % Export static trial with virtual makers to new TRC file
    VirtualData = [Mid.HJC, Mid.ASIS, Mid.PSIS, Mid.Pelvis, R.KJC, L.KJC, R.AJC, L.AJC, R.AJC_Floor, L.AJC_Floor,...
        R.Heel_Floor, L.Heel_Floor, R.MT1_Floor, L.MT1_Floor, R.MT5_Floor, L.MT5_Floor, R.MidMT_Floor, L.MidMT_Floor];
    % load TRC data
    file = strcat(Subjects(subj).Trials(StaticTrials).folder, '\OpenSim\', ...
        Subjects(subj).Trials(StaticTrials).files.OpenSimAddTorso);
    file2write = strcat(file(1:end-4), 'Virtual.trc');
    data = importdata(file,'\t',5);
    TRC = data.data;
    
    if isnan(TRC(1,2))  %delete first row if necessary
        TRC(1,:) = [];
    end
    if isnan(TRC(2,1))  %delete first column if necessary
        TRC(:,1) = [];
    end
    [NumRow, NumColumn] = size(TRC);
    
    % Correct Headers
    fid = fopen(file,'r'); %Read TRC files to determine correct headers
    FData = textscan(fid,'%q'); % Looks at headers
    FData = FData{1,1}; % Reorganizes
    FData{4,1} = file; % Renames file name
    fclose(fid);
    
    % find header locations in file
    StartMetrics = find(strcmp(FData,'DataRate'), 1);
    LastCategory = find(strcmp(FData,'OrigNumFrames'), 1);
    ColumnHeaderStart = find(strcmp(FData,'Frame#'), 1);
    Col = str2double(FData{18,1});
    VirtualHeaders = {'Mid.HJC', 'Mid.ASIS', 'Mid.PSIS', 'Mid.Pelvis', 'R.KJC', 'L.KJC', 'R.AJC', 'L.AJC', 'R.AJC_Floor', 'L.AJC_Floor',...
        'R.Heel_Floor', 'L.Heel_Floor', 'R.MT1_Floor', 'L.MT1_Floor', 'R.MT5_Floor', 'L.MT5_Floor', 'R.MidMT_Floor', 'L.MidMT_Floor'};
    FData{18,1} = num2str(Col + length(VirtualHeaders)); 
    
    fid = fopen(file2write,'w'); % Open TRC files to write
    % Copies over header info
    fprintf(fid,'%s\t', FData{1:StartMetrics-1,1});
    fprintf(fid,'\n'); % new line
    fprintf(fid,'%s\t', FData{StartMetrics:LastCategory,1});
    fprintf(fid,'\n'); % new line
    fprintf(fid,'%s\t', FData{LastCategory+1:ColumnHeaderStart-1,1});
    fprintf(fid,'\n'); % new line
    fprintf(fid,'%s\t', FData{ColumnHeaderStart:ColumnHeaderStart+1,1});
    
    % Copies over column headers with two tabs in between each marker
    NCol = (NumColumn-2) / 3;
       % add virtual marker names to header
    NewHeaders = FData(ColumnHeaderStart+2:ColumnHeaderStart+NCol+1,1);
    Ind = contains(NewHeaders, 'L_HJC'); 
    NewHeaders{Ind} = 'L.HJC';
    Ind = contains(NewHeaders, 'R_HJC');
    NewHeaders{Ind} = 'R.HJC';
    
    NewHeaders(end+1:end+length(VirtualHeaders)) = VirtualHeaders;
    for m=1:NCol+length(VirtualHeaders)
        fprintf(fid,'%s\t', NewHeaders{m}); %FData{m+ColumnHeaderStart+1,1});
        fprintf(fid,'\t'); % tab
        fprintf(fid,'\t'); % tab
    end
    
    fprintf(fid,'\n');  % new line and
    fprintf(fid,'\t');  % tab over twice to pass the Frame and Time columns
    fprintf(fid,'\t');
    
    % Labels x, y, and z columns
    d = 0;
    for m = 1:NCol+length(VirtualHeaders)
        d = d+1;
        fprintf(fid,'%s\t', ['X',num2str(d)]);
        fprintf(fid,'%s\t', ['Y',num2str(d)]);
        fprintf(fid,'%s\t', ['Z',num2str(d)]);
    end
    
    % Adds a space between xyz and data.
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    % Inputs new/old marker x, y, z info
    TRCnew = [TRC, VirtualData];
    for ii = 1:NumRow
        fprintf(fid,'%3.4f\t',TRCnew(ii,:)); fprintf(fid,'\n');
    end
    
    fclose('all');
    
%     = strcat(Subjects(subj).Trials(StaticTrials).folder, '\OpenSim\ScaleFiles\', file2write); 
     Dst = strcat(Subjects(subj).Trials(StaticTrials).folder, '\OpenSim\ScaleFiles\', ...
        Subjects(subj).Trials(StaticTrials).files.OpenSimAddTorso(1:end-4), 'Virtual.trc');
    movefile(file2write, Dst); 
    Subjects(subj).Trials(StaticTrials).files.Virtual = Dst;
    
    % This seems pointless, but OpenSim won't recognize the file unless it is opened and saved again
    e = actxserver('excel.application');
    eW = e.Workbooks;
    eF = eW.Open(Dst); % open file
    eS = eF.ActiveSheet;
    eF.Save;
    eF.Close; % close the file
    e.Quit; % close Excel entirely
    
    clearvars data FData fid TRC TRCnew ii d m NCol NumRow NumColumn e eF eR eS eW Col MkrSz NewHeaders ...
        file file2write Ind L R Mid  LastCategory StartMetrics A ans ColumnHeaderStart Dst VirtualData VirtualHeaders a b c
    
    %% Scale Model
    clc; disp('Setting up Scaling');
    
    import org.opensim.modeling.*
    import java.io.*
    
    % create scale folder
    ScaleFolder = strcat(Subjects(subj).Folders.OpenSimFolder, '\ScaleFiles');
    mkdir(ScaleFolder);
    
    % move static trial into scale folder
    if strcmp(Settings.LockModel, 'Full')
        StaticMarkerFile = strcat(Subjects(subj).Folders.OpenSimFolder, '\', ...
            Subjects(subj).Trials(StaticTrials).files.OpenSimTRC);
        NewStaticMarkerFile = strcat(ScaleFolder, '\', Subjects(subj).Trials(StaticTrials).files.OpenSimTRC);
        copyfile(StaticMarkerFile, NewStaticMarkerFile);
        StaticMarkerFileGRF = strrep(StaticMarkerFile,'.trc','GRF.mot');
        NewStaticMarkerFileGRF = strrep(NewStaticMarkerFile,'.trc','GRF.mot');
        copyfile(StaticMarkerFileGRF, NewStaticMarkerFileGRF);
    elseif strcmp(Settings.LockModel, 'Yes')
        StaticMarkerFile = strcat(Subjects(subj).Folders.OpenSimFolder, '\', ...
            Subjects(subj).Trials(StaticTrials).files.OpenSimAddTorso);
        NewStaticMarkerFile = strcat(ScaleFolder, '\', Subjects(subj).Trials(StaticTrials).files.OpenSimAddTorso);
        copyfile(StaticMarkerFile, NewStaticMarkerFile);
        StaticMarkerFileGRF = strrep(StaticMarkerFile,'AddTorso.trc','GRF.mot');
        NewStaticMarkerFileGRF = strrep(NewStaticMarkerFile,'.trc','GRF.mot');
        copyfile(StaticMarkerFileGRF, NewStaticMarkerFileGRF);
    end
    % define and copy original model OSIM file
    if strcmp(Settings.LockModel, 'Yes')
        Orig.Model = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name},...
            'gait2392_Scale_ABLMarkerSet_UMBprobed_locked_STRONG.osim')).name);
        SubjModel = strcat(ScaleFolder, '\gait2392_Scale_ABLMarkerSet_UMBprobed_locked_STRONG.osim');
    elseif strcmp(Settings.LockModel, 'Full')
        Orig.Model = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name},...
            'gait2392_Scale_ABLMarkerSet_UMBprobed_lockedTorsoFeet.osim')).name);
        SubjModel = strcat(ScaleFolder, '\gait2392_Scale_ABLMarkerSet_UMBprobed_lockedTorsoFeet.osim');
    elseif strcmp(Settings.ModelStrength, 'Normal')
        Orig.Model = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name},...
            'gait2392_Scale_ABLMarkerSet_UMBprobed_locked.osim')).name);
        SubjModel = strcat(ScaleFolder, '\gait2392_Scale_ABLMarkerSet_UMBprobed_locked.osim');
        %     else
        %         Orig.Model = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, ...
        %             'gait2392_Scale_ABLMarkerSet_UMBprobed.osim')).name);
        %         SubjModel = strcat(ScaleFolder, '\gait2392_Scale_ABLMarkerSet_UMBprobed.osim');
    end
    copyfile(Orig.Model, SubjModel);
    
    % Load the model and initialize
    model = Model(fullfile(Orig.Model));
    model.initSystem();
    
    % copy over original markerset file
    Orig.MkrSetFile = strcat(GenericFilePath, '\',...
        GenericDir(contains({GenericDir.name}, 'gait2392_Scale_MarkerSet_ABL')).name);
    MkrSetFile = strcat(ScaleFolder, '\', Subjects(subj).name, '_MkrSet.xml');
    copyfile(Orig.MkrSetFile, MkrSetFile); % copy generic scale file to subject setup directory
    
    % identify setup XML file
%         Orig.ScaleSetupFile = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup_Scale')).name);
    Orig.ScaleSetupFile = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup__Scale')).name);
    [scaleXML, scaleRootName, ~] = xml_read(Orig.ScaleSetupFile);
    
    % Change attributes in structure
    scaleXML.ScaleTool.mass = Subjects(subj).Demo.mass;
    scaleXML.ScaleTool.height = Subjects(subj).Demo.height;
    scaleXML.ScaleTool.age = Subjects(subj).Demo.age;
    scaleXML.ScaleTool.ATTRIBUTE.name = Subjects(subj).name;
    
    % GenericModelMaker
    if strcmp(Settings.LockModel, 'Yes')
        scaleXML.ScaleTool.GenericModelMaker.model_file = 'gait2392_Scale_ABLMarkerSet_UMBprobed_locked_STRONG.osim';
    elseif strcmp(Settings.LockModel, 'Full')
        scaleXML.ScaleTool.GenericModelMaker.model_file = 'gait2392_Scale_ABLMarkerSet_UMBprobed_lockedTorsoFeet.osim';
    elseif strcmp(Settings.LockModel, 'Yes')
        scaleXML.ScaleTool.GenericModelMaker.model_file = 'gait2392_Scale_ABLMarkerSet_UMBprobed_locked.osim';
        %     else
        %         scaleXML.ScaleTool.GenericModelMaker.model_file = 'gait2392_Scale_ABLMarkerSet_UMBprobed.osim';
    end
    scaleXML.ScaleTool.GenericModelMaker.marker_set_file = strcat(Subjects(subj).name, '_MkrSet.xml');
    
    % only use first 3 frames from static trial
    StaticTrialNum = contains({Subjects(subj).Trials.type}, 'static');
    ScaleStart = Subjects(subj).Trials(StaticTrialNum).Times.TRC(1);
    ScaleEnd = Subjects(subj).Trials(StaticTrialNum).Times.TRC(3);
    
    % ModelScaler
    MarkerFile = strcat(Subjects(subj).Trials(StaticTrials).files.OpenSimAddTorso);
    scaleXML.ScaleTool.ModelScaler.marker_file = MarkerFile;
    OutputModelFile = strcat(Subjects(subj).name, '_Scaled.osim');
    scaleXML.ScaleTool.ModelScaler.output_model_file = OutputModelFile;
    scaleXML.ScaleTool.ModelScaler.time_range = [ScaleStart ScaleEnd];
    
    % MarkerPlacer
    scaleXML.ScaleTool.MarkerPlacer.marker_file = MarkerFile;
    scaleXML.ScaleTool.MarkerPlacer.output_model_file = OutputModelFile;
    scaleXML.ScaleTool.MarkerPlacer.time_range = [ScaleStart ScaleEnd];
    
    % export  XML to SetupDir (specific to each subject)
    SetupScale = strcat(ScaleFolder, '\', Subjects(subj).name, '_Setup_Scale.xml');
    xml_write(SetupScale, scaleXML, scaleRootName);
    
    scale = ScaleTool(SetupScale); % open scaling tool with new attributes
    % scale.ScaleTool(ModelScaler.getOutputModelFileName(OutputModelFile));
    
    %%
    fopen('ScaleResults.txt','w'); % clear old scale results
    diary('ScaleResults.txt');
    diary on
    scale.run(); % run scaling
    diary off
    
    
    FID = fopen('ScaleResults.txt');
    %        FID = fopen('ScaleResults - Copy.txt');
    TXT = textscan(FID, '%s');
    
    for i = 1:length(TXT{1})
        if contains(TXT{1}(i), {'total'}) && contains(TXT{1}(i+1), {'square'}) && ...
                contains(TXT{1}(i+2), {'error'}) && contains(TXT{1}(i+3), {'='})
            Ind = i;
            break
        end
    end
    
    Subjects(subj).Trials(StaticTrialNum).ScaleErr.TotalSqErr = str2double(TXT{1}(Ind+4));
    String = cell2mat(TXT{1}(Ind+7));
    Subjects(subj).Trials(StaticTrialNum).ScaleErr.RMSErr = str2double(String(5:end-1));
    String = cell2mat(TXT{1}(Ind+8));
    Subjects(subj).Trials(StaticTrialNum).ScaleErr.MaxErr = str2double(String(5:end));
    String = cell2mat(TXT{1}(Ind+9));
    Subjects(subj).Trials(StaticTrialNum).ScaleErr.MaxMkr = String(2:end-1);
    
    % scale strength and other muscle parameters to age, height, mass, etc?
    
    clearvars ScaleRootName Geopath StaticMarkerFile NewStaticMarkerFile...
        StaticMarkerFileGRF NewStaticMarkerFileGRF i col scale Ind String
    
    
    
end % End Subject loop

% save('SubjectsScaleIKCheck_UpdatedScaling2.mat', 'Subjects');
% save('SubjectsScaleIKCheck.mat', 'Subjects');


%%
% ----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information. OpenSim is developed at Stanford University       %
% and supported by the US National Institutes of Health (U54 GM072970,    %
% R24 HD065690) and by DARPA through the Warrior Web program.             %
%                                                                         %
% Copyright (c) 2005-2020 Stanford University and the Authors                      %
%                                                                         %
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file except in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         %
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
% ----------------------------------------------------------------------- %
% This batch processing code originated from Edith Arnold - setupAndRunIKBatchExample.m

% Author: Ricky Pimentel, December 2019
% Applied Biomechanics Lab, University of North Carolina at Chapel Hill

% To Set up OpenSim API for MATLAB, follow directions at:
% https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab

