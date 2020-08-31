%% ABL_BatchOpenSim (BOS)
% Batch process motion files in OpenSim
% See README for documentation 

% This batch processing code originated from Edith Arnold - setupAndRunIKBatchExample.m

% To Set up OpenSim API for MATLAB, follow directions at:
% https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab

% Author: Ricky Pimentel, December 2019
% Applied Biomechanics Lab, University of North Carolina at Chapel Hill

% All code within 


%% Settings
cd('C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\BatchOpenSim');
clear; clc; warning off; dbstop if error;

Settings.SetupBatch = 'No'; % Set up batch processing by converting TRC
Settings.TrialWindows = 'Yes'; % set to 'No' to run IK for entire duration of each trial
Settings.NumStrides = [1]; % number of strides to check for quality data if results don't pass quality control
Settings.LockModel = 'Yes'; % set to 'Yes' to use models that have locked subtalar and MTP joints
% Settings.LockModel = 'Full'; % set to Full to use models that have locked subtalar, MTP, and lumbar joints
Settings.ModelStrength = 'Strong'; % set to 'Strong' to use model with 1.5x strength; otherwise 'Normal'
Settings.LowPassKinematicsFilter = 6; % define all kinematics filtering to a specified value in Hz

Settings.Scale = 'Yes';
Settings.IK = 'Yes';
Settings.RRA = 'Yes';
Settings.SO = 'No';
Settings.CMC = 'Yes';
Settings.MA = 'No';
Settings.ID = 'No';

% plot settings
Settings.PlotIKErrors = 'No';

% Pull in the modeling classes straight from the OpenSim distribution
import org.opensim.modeling.*
import java.io.*

% MAY NEED TO SET UP USERS SPECIFIC PATH TO MODEL GEOMETRY
Geopath = 'C:\Users\richa\Documents\OpenSim 4.0\Geometry'; % set to location of geometry files
ModelVisualizer.addDirToGeometrySearchPaths(Geopath); % add geometry files to path

addpath(genpath('Functions')); % add BatchOpenSim Functions to path

%% move to study directory
D = dir; % get current path
CurrFolder = D(1).folder; % assign folder
GenericFilePath = strcat(CurrFolder, '\OpenSimProcessingFiles'); % set generic folder path
addpath(genpath(GenericFilePath)); % add generic folder to path
GenericDir = dir(GenericFilePath);
Settings.GenericPath = GenericFilePath;
Settings.GenericDir = GenericDir;

%% load subject anthropometrics
% set your demographics file here
Demo.File = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Metabolics_Feedback_Demographics';
[~,~,Demo.Raw] = xlsread(Demo.File);
Demo.AgeCol = strcmp(Demo.Raw(1,:), 'Age');
Demo.HtCol = strcmp(Demo.Raw(1,:), 'Height (cm)');
Demo.GenderCol = strcmp(Demo.Raw(1,:), 'Gender');
% Demo.MassCol = strcmp(Demo.Raw(1,:), 'Mass'); % overwrite mass from static trial

%% Setup Batch Process
if strcmp(Settings.SetupBatch, 'Yes')
    % copy
    [Subjects, subjectPath] = ABL_OpenSim_Setup_Batch;
    addpath(genpath(subjectPath));
else
    % hard code subject directory
    subjectPath = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\exp';
    
    % manually select subject directory (comment above)
    if ~exist(subjectPath, 'dir')
        subjectPath = uigetdir(CurrFolder, 'Select Folder Containing Subject Data');
    end
    addpath(genpath(subjectPath));
    load('Subjects.mat'); % load subjects from conversion batch
end


%% Save Demographics and Setup Folders for IK, RRA, and CMC
for subj = 1:length(Subjects)
    
    % set names and demographics
    SubjRow = find(strcmp(Demo.Raw(:,1), Subjects(subj).name),1);
    Subjects(subj).Demo.age = Demo.Raw{SubjRow, Demo.AgeCol};
    Subjects(subj).Demo.height = Demo.Raw{SubjRow, Demo.HtCol} * 10; % height from cm to mm
    
    % identify trial types
    StaticTrials = logical(strcmp({Subjects(subj).Trials.type}, 'static'));
    DynamicTrials = ~StaticTrials;
    
    if sum(StaticTrials) > 1
        error('Multiple static trials detected. Reduce to only one.');
    end
    
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
    
    clearvars A B C D E
    
end
clearvars New Old scaleRootName SetupScale col

%% Scale Model
[Subjects] = ABL_Scale(Settings, StaticTrials, Subjects);
[scaleError] = ComputeScaleError(Subjects);

ScaleErrTbl = struct2table(scaleError); 
avg = mean(ScaleErrTbl.RMSError); 
SD = std(ScaleErrTbl.RMSError); 
disp(strcat('Average RMS Scale Error = ', num2str(avg))); 
disp(strcat('Standard Dev RMS Scale Error = ', num2str(SD))); 
[m,i] = max(ScaleErrTbl.MaxError); 
disp(strcat('Max scale error = ', num2str(m))); 
disp(strcat('Max scale error Marker = ', ScaleErrTbl.MaxMkr(i)));
disp('All errors in m'); 


%% Main Processing Loop
for subj = 9%1:length(Subjects) % Subject Loop
    
    for trial = 4%find(DynamicTrials) % Trial Loop
        
        %% Create Processing Folders
        % create scale folder
        ScaleFolder = strcat(Subjects(subj).Folders.OpenSimFolder, '\ScaleFiles');
        mkdir(ScaleFolder);
        
        % create IK folder
        IKFolder = strcat(Subjects(subj).Folders.OpenSimFolder , '\IK_Files');
        mkdir(IKFolder);
        
        % copy generic setup XML file
        Orig.IK_Setup = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup_IK.xml')).name);
        IKfile.Setup = strcat(IKFolder, '\', Subjects(subj).name, '_Setup_IK.xml');
        copyfile(Orig.IK_Setup, IKfile.Setup);
        
        % create RRA folders
        RRAFolder = strcat(Subjects(subj).Folders.OpenSimFolder , '\RRA_Files');
        mkdir(RRAFolder);
        RRAFolder2 = strcat(Subjects(subj).Folders.OpenSimFolder , '\RRA_Files2');
        mkdir(RRAFolder2);
        
        % CMC folder
        CMCFolder = strcat(Subjects(subj).Folders.OpenSimFolder , '\CMC_Files');
        mkdir(CMCFolder);
        
        %% Initalize Settings for Trial loop
        tic;
        import org.opensim.modeling.*
        import java.io.*
        
        % define trial names and labels
        if strcmp(Settings.LockModel, 'Full')
            MarkerFile = Subjects(subj).Trials(trial).files.OpenSimTRC;    % Get the trial file name
            
        elseif strcmp(Settings.LockModel, 'Yes')
            MarkerFile = Subjects(subj).Trials(trial).files.OpenSimAddTorso;    % Get the trial file name
        end
        TrialNameLong = regexprep(MarkerFile,'.trc','');  % Create name of trial from .trc file name
        Str = strsplit(TrialNameLong, '_');
        TrialName = Str{1};
        SubjName = Subjects(subj).name;
        fullpath = fullfile(Subjects(subj).Folders.OpenSimFolder, MarkerFile);
        markerData = MarkerData(fullpath);  % Get trc data to determine time range
        
        GRFtrialname = strcat(Subjects(subj).Folders.OpenSimFolder, '\', Subjects(subj).Trials(trial).files.OpenSimGRF);
        
        % Use trial windows, or process whole trial
        if strcmp(Settings.TrialWindows, 'Yes')
            % Get start and end time of each stride for RRA and CMC
            
            % add small time window on to the start and end to offset CMC delay of 0.03 s
            Window = 0.05;
            
            % LEFT STRIDE
            L_Strikes = Subjects(subj).Trials(trial).TSData.L_Strike(:,2);
            L_Inds = find(L_Strikes >  Subjects(subj).Trials(trial).Times.StartWindow, Settings.NumStrides);
            % make sure time starts after 1 second into the trial
            if L_Strikes(L_Inds(1)) < 1
                L_Inds = L_Inds + 1;
            end
            % make sure strides arent crossover steps
            while sum(L_Inds == Subjects(subj).Trials(trial).Cross.L_Stride) > 0
                L_Inds = L_Inds + 1;
            end
            Subjects(subj).Trials(trial).Times.Start_Left = L_Strikes(L_Inds) - Window;
            Subjects(subj).Trials(trial).Times.End_Left = L_Strikes(L_Inds+1) + Window;
            
            % RIGHT STRIDE
            R_Strikes = Subjects(subj).Trials(trial).TSData.R_Strike(:,2);
            R_Inds = find(R_Strikes >  Subjects(subj).Trials(trial).Times.StartWindow, Settings.NumStrides);
            % make sure time starts after 1 second into the trial
            if R_Strikes(R_Inds(1)) < 1
                R_Inds = R_Inds + 1;
            end
            % make sure strides arent crossover steps
            while sum(R_Inds == Subjects(subj).Trials(trial).Cross.R_Stride) > 0
                R_Inds = R_Inds + 1;
            end
            Subjects(subj).Trials(trial).Times.Start_Right = R_Strikes(R_Inds) - Window;
            Subjects(subj).Trials(trial).Times.End_Right = R_Strikes(R_Inds+1) + Window;
            
            % define overall start and end time for IK (all strides)
            Subjects(subj).Trials(trial).Times.Start_IK = ...
                min([Subjects(subj).Trials(trial).Times.Start_Left(1), Subjects(subj).Trials(trial).Times.Start_Right(1)]);
            Subjects(subj).Trials(trial).Times.End_IK = ...
                max([Subjects(subj).Trials(trial).Times.End_Left(end), Subjects(subj).Trials(trial).Times.End_Right(end)]);
            
        else % use full trial time
            Subjects(subj).Trials(trial).initial_time = markerData.getStartFrameTime();
            Subjects(subj).Trials(trial).final_time = markerData.getLastFrameTime();
        end
        
        %% IK
        import org.opensim.modeling.*
        import java.io.*
        
        if strcmp(Settings.IK, 'Yes')
            
            % create IK folder
            IKFolder = strcat(Subjects(subj).Folders.OpenSimFolder , '\IK_Files');
            mkdir(IKFolder);
            
            % copy generic setup XML file
            Orig.IK_Setup = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup_IK.xml')).name);
            IKfile.Setup = strcat(IKFolder, '\', Subjects(subj).name, '_Setup_IK.xml');
            copyfile(Orig.IK_Setup, IKfile.Setup);
            
            ScaleFolder = strcat(Subjects(subj).Folders.OpenSimFolder, '\ScaleFiles');
            OutputModelFile = strcat(Subjects(subj).name, '_Scaled.osim');
            ikTool = InverseKinematicsTool(IKfile.Setup); % initialize IK tool
            model = Model(strcat(ScaleFolder, '\', OutputModelFile)); % load scaled model for IK
            model.initSystem(); % initialize model
            ikTool.setModel(model);  % Tell Tool to use the loaded model
            
            % Setup ikTool
            ikTool.setName(TrialName);
            ikTool.setMarkerDataFileName(fullpath);
            
            % Run IK for full time of multiple strides
            ikTool.setStartTime( Subjects(subj).Trials(trial).Times.Start_IK);
            ikTool.setEndTime( Subjects(subj).Trials(trial).Times.End_IK);
            IKfile.Output = fullfile(IKFolder, [TrialName '_ik.mot']);
            ikTool.setOutputMotionFileName(IKfile.Output);
            
            % marker weighting already set in generic file
            
            % Save the Setup XML settings in a setup file
            outfile = ['Setup_IK_' TrialName '.xml'];
            ikTool.print(fullfile(IKFolder, outfile));
            ikTool.run(); % Run IK
            
            
            %% Process IK Error files and interpret
            clearvars ind
            D = dir;
            ind = contains({D.name}, 'errors.sto');
            ErrorFileName = D(ind).name;
            movefile(ErrorFileName, IKFolder);
            ErrorData = importdata(strcat(IKFolder, '\', ErrorFileName));
            Subjects(subj).Trials(trial).IKErrors.Average = mean(ErrorData.data(:,3)); % mean of the mean
            Subjects(subj).Trials(trial).IKErrors.Max = max(ErrorData.data(:,4)); % max of the max
            
            % define appropriate error thresholds
            Threshold.RMSError = 0.025;
            Threshold.MaxError = 0.05;
            
            % document if errors are low enough for inverse kinematics
            if Subjects(subj).Trials(trial).IKErrors.Average > Threshold.RMSError
                Subjects(subj).Trials(trial).Pass_IK = 'No';
            elseif  Subjects(subj).Trials(trial).IKErrors.Max > Threshold.MaxError
                Subjects(subj).Trials(trial).Pass_IK = 'No';
            else
                Subjects(subj).Trials(trial).Pass_IK = 'Yes';
            end
            
            % move IK marker locations file
            ind = contains({D.name}, 'model_marker_locations.sto');
            MarkerLocFileName = D(ind).name;
            movefile(MarkerLocFileName, IKFolder);
            
            %         clearvars ErrorData D ind ErrorFileName
            
            %% Plot IK marker errors
            % load model marker data
            if strcmp(Settings.PlotIKErrors, 'Yes')
                ModelMkrData = importdata(strcat(IKFolder, '\', MarkerLocFileName));
                
                Ind = contains(ModelMkrData.colheaders, '_tx');
                ModelMkrData.Xdata = ModelMkrData.data(:,Ind);
                Ind = contains(ModelMkrData.colheaders, '_ty');
                ModelMkrData.Ydata = ModelMkrData.data(:,Ind);
                Ind = contains(ModelMkrData.colheaders, '_tz');
                ModelMkrData.Zdata = ModelMkrData.data(:,Ind);
                
                Markers = strrep([ModelMkrData.colheaders(Ind)], '_tz','');
                Times = Subjects(subj).Trials(trial).TRC.data(:,2);
                
                [StartVal, StartInd] = min(abs(Times - Subjects(subj).Trials(trial).Times.Start_IK));
                [EndVal, EndInd] = min(abs(Times - Subjects(subj).Trials(trial).Times.End_IK));
                
                % get IK marker data
                for i = 1:length(Markers)
                    Ind = find(strcmp(Subjects(subj).Trials(trial).TRC.colheaders(1,:), Markers{i}));
                    MarkerData(i).name = Markers{i};
                    MarkerData(i).x = Subjects(subj).Trials(trial).TRC.data(StartInd:EndInd,Ind);
                    MarkerData(i).y = Subjects(subj).Trials(trial).TRC.data(StartInd:EndInd,Ind+1);
                    MarkerData(i).z = Subjects(subj).Trials(trial).TRC.data(StartInd:EndInd,Ind+2);
                end
                
                % plot
                MkrSz = 10;
                close all; clc;
                H = figure('Position',[100 100 900 800]);
                L = min([length(ModelMkrData.Zdata), length(MarkerData(1).z)]);
                
                filename = strcat(Subjects(subj).name,'_', Subjects(subj).Trials(trial).name, '_IKMkrErr.gif');
                
                for r = 1:100%L
                    mdl_z = 1000*ModelMkrData.Zdata(r,:);
                    mdl_x = 1000*ModelMkrData.Xdata(r,:);
                    mdl_y = 1000*ModelMkrData.Ydata(r,:);
                    
                    for i = 1:length(MarkerData)
                        if isempty(MarkerData(i).x)
                            continue
                        end
                        orig_x(i) = MarkerData(i).x(r);
                        orig_y(i) = MarkerData(i).y(r);
                        orig_z(i) = MarkerData(i).z(r);
                    end
                    subplot(221);
                    m1 = plot3(mdl_z, mdl_x, mdl_y, '.b', 'MarkerSize', MkrSz);
                    hold on; grid on;
                    m2 = plot3(orig_z, orig_x, orig_y, '.r', 'MarkerSize', MkrSz);
                    line([mdl_z(1:28); orig_z], [mdl_x(1:28); orig_x], [mdl_y(1:28); orig_y],'Color','k');
                    title(['Frame ' num2str(r)]);
                    view([-1400, 1200, 1200]);
                    
                    axis equal;
                    xlim([-1000 0]);
                    ylim([0 1000]);
                    zlim([0 1300]);
                    %                 pause(0.01);
                    hold off;
                    
                    subplot(222);
                    text(0.2, 0.8, '\bf Model','Color', 'b', 'FontSize',20);
                    text(0.2, 0.6, '\bf Original','Color', 'r', 'FontSize',20);
                    
                    subplot(223);
                    plot(mdl_x, mdl_y, '.b', 'MarkerSize', MkrSz);
                    hold on; grid on;
                    plot(orig_x, orig_y, '.r', 'MarkerSize', MkrSz);
                    line([mdl_x(1:28); orig_x], [mdl_y(1:28); orig_y],'Color','k');
                    hold off; axis equal;
                    title(['Time = ' num2str(ModelMkrData.data(r,1))]);
                    
                    subplot(224);
                    plot(mdl_z, mdl_y, '.b', 'MarkerSize', MkrSz);
                    hold on; grid on;
                    plot(orig_z, orig_y, '.r', 'MarkerSize', MkrSz);
                    line([mdl_z(1:28); orig_z], [mdl_y(1:28); orig_y],'Color','k');
                    hold off; axis equal;
                    title(['Frame ' num2str(r)]);
                    
                    % Capture the plot as an image
                    frame = getframe(H);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    % Write to the GIF File
                    if r == 1
                        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                    else
                        imwrite(imind,cm,filename,'gif','WriteMode','append');
                    end
                    
                end
            end
            
        end
        
        
        %% Update External Loads file
        Orig.ExternalLoads = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'ExternalLoads')).name);
        RRAfile.ExternalLoads = strcat(RRAFolder, '\', SubjName, '_ExternalLoads.xml');
        copyfile(Orig.ExternalLoads, RRAfile.ExternalLoads);  % external loads file
        
        [ExtLdsXML, ExtLdsRootName, ~] = xml_read(RRAfile.ExternalLoads); % read XML file
        xls = strrep(MarkerFile, 'OpenSim.trc', 'grf.xls');
        ExtLdsXML.ExternalLoads.objects.ExternalForce(1).applied_to_body = 'calcn_r';
        ExtLdsXML.ExternalLoads.objects.ExternalForce(1).data_source_name = xls;
        ExtLdsXML.ExternalLoads.objects.ExternalForce(1).torque_identifier = 'ground_torque';
        ExtLdsXML.ExternalLoads.objects.ExternalForce(2).applied_to_body = 'calcn_l';
        ExtLdsXML.ExternalLoads.objects.ExternalForce(2).data_source_name = xls;
        ExtLdsXML.ExternalLoads.objects.ExternalForce(2).torque_identifier = '1_ground_torque';
        ExtLdsXML.ExternalLoads.datafile = GRFtrialname;
        ExtLdsXML.ExternalLoads.external_loads_model_kinematics_file = IKfile.Output;
        ExtLdsXML.ExternalLoads.lowpass_cutoff_frequency_for_load_kinematics = ...
            Settings.LowPassKinematicsFilter; % update low pass filter for kinematics
        ExtLdsTrial = strcat(RRAfile.ExternalLoads(1:end-4), '_', TrialName, '.xml'); % name trial-specific external loads XML file
        xml_write(ExtLdsTrial, ExtLdsXML, ExtLdsRootName); % write external loads XML file for specific trial
        
        
        %% Copy Over Files for RRA and CMC
        
        % copy over RRA setup file from OpenSimProcessingFiles folder
        Orig.RRA_Setup = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup_RRA')).name);
        RRAfile.Setup = strcat(RRAFolder, '\', SubjName, '_Setup_RRA.xml');
        copyfile(Orig.RRA_Setup , RRAfile.Setup);  % RRA setup file
        
        % copy scaled model file into RRA folder
        Old = strcat(ScaleFolder, '\', OutputModelFile);
        New = strcat(RRAFolder, '\', OutputModelFile);
        copyfile(Old, New);
        
        % define RRA Actuator and Tasks files
        if strcmp(Settings.LockModel, 'Yes')
            Orig.RRA_Actuators = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'RRA_Actuators_Locked.xml')).name);
            Orig.RRA_Tasks = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'RRA_Tasks_Locked.xml')).name);
        elseif strcmp(Settings.LockModel, 'Full')
            Orig.RRA_Actuators = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'RRA_Actuators_LockedTorsoFeet.xml')).name);
            Orig.RRA_Tasks = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'RRA_Tasks_LockedTorsoFeet.xml')).name);
        else
            Orig.RRA_Actuators = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'RRA_Actuators.xml')).name);
            Orig.RRA_Tasks = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'RRA_Tasks.xml')).name);
        end
        
        % copy RRA actuators file
        RRAfile.Actuators = strcat(RRAFolder, '\', SubjName, '_RRA_Actuators.xml');
        RRAfile.ActuatorFile = strcat(SubjName, '_RRA_Actuators.xml');
        copyfile(Orig.RRA_Actuators, RRAfile.Actuators);
        RRAfile.Actuators2 = strcat(RRAFolder2, '\', SubjName, '_RRA_Actuators.xml');
        RRAfile.ActuatorFile2 = strcat(SubjName, '_RRA_Actuators.xml');
        copyfile(Orig.RRA_Actuators, RRAfile.Actuators2);
        
        % copy RRA tasks file
        RRAfile.Tasks = strcat(RRAFolder, '\', SubjName, '_RRA_Tasks.xml');
        copyfile(Orig.RRA_Tasks, RRAfile.Tasks);
        RRAfile.Tasks2 = strcat(RRAFolder2, '\', SubjName, '_RRA_Tasks.xml');
        copyfile(Orig.RRA_Tasks, RRAfile.Tasks2);
        
        
        %% RRA and CMC Loops for each side and stride
        import org.opensim.modeling.*
        import java.io.*
        
        Side = {'Left','Right'};
        for side = [2]
            for Stride = Settings.NumStrides
                
                SubjTrialSideStride = strcat(SubjName, '_', TrialName, '_', Side{side}, '_', num2str(Stride));
                
                %% Setup and run RRA 1
                % make trial-specific results directory
                RRATrialFolder = strcat(RRAFolder, '\', TrialName, '\', Side{side}, '_', num2str(Stride));
                mkdir(RRATrialFolder);
                
                % read in new XML file preferences
                [rraXML, rraRootName, ~] = xml_read(RRAfile.Setup);
                
                % update general RRA attributes
                OutputModelFile = strcat(RRAFolder, '\', Subjects(subj).name, '_Scaled.osim');
                rraXML.RRATool.model_file = OutputModelFile;  % use model output from IK
                rraXML.RRATool.force_set_files = RRAfile.ActuatorFile;
                rraXML.RRATool.results_directory = RRATrialFolder;
                rraXML.RRATool.task_set_file = RRAfile.Tasks;
                rraXML.RRATool.external_loads_file = ExtLdsTrial;
                
                % update RRA settings for specific trial
                % set gait cycle times
                if strcmp(Side{side}, 'Left')
                    Start = Subjects(subj).Trials(trial).Times.Start_Left(Stride);
                    End = Subjects(subj).Trials(trial).Times.End_Left(Stride);
                elseif strcmp(Side{side}, 'Right')
                    Start = Subjects(subj).Trials(trial).Times.Start_Right(Stride);
                    End = Subjects(subj).Trials(trial).Times.End_Right(Stride);
                end
                rraXML.RRATool.initial_time = Start; % start time
                rraXML.RRATool.final_time = End; % end time
                rraXML.RRATool.desired_kinematics_file = IKfile.Output;
                rraXML.RRATool.output_model_file = strcat(Subjects(subj).name, '_RRA_Adj.osim'); % output model
                rraXML.RRATool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                rraXML.RRATool.ATTRIBUTE = strcat(Subjects(subj).name, '_RRA');
                RRAfile.Trial = strcat(RRAfile.Setup(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific RRA setup XML file
                xml_write(RRAfile.Trial, rraXML, rraRootName);  % write new RRA xml specific to trial
                
                % setup and run RRA 1
                rraTool = RRATool(RRAfile.Trial); %   initialize RRA tool
                rraTool.setAdjustCOMToReduceResiduals(1);
                
                % create text file for mass changes
                RRAfile.MassChanges = strcat(RRATrialFolder, '\', SubjTrialSideStride, '_MassChanges.txt');
                
                fopen(RRAfile.MassChanges,'w'); % clear old scale results
                diary(RRAfile.MassChanges);
                diary on;
                rraTool.run(); % run RRA
                diary off;
                
                %% Update Segment Masses from RRA1 results
                % read in segment masses from txt file to change
                fid = fopen(RRAfile.MassChanges, 'r');
                TxtData1 = textscan(fid,'%q');
                TxtDataString = TxtData1{1,1};
                [StarWars,~] = find(contains(TxtDataString, '*****'));
                fclose(fid);
                
                TxtData = TxtDataString(StarWars(1):end);
                CellDif = 8; % number of cells between segment name and new mass value
                Segments = {'pelvis','femur_r','tibia_r','talus_r','calcn_r','toes_r', ...
                    'femur_l','tibia_l','talus_l','calcn_l','toes_l', 'torso'};
                
                SegInOutput = find(contains(TxtData, Segments));
                SegInOutput([1, 14, 15, 16]) = []; % remove mentions of pelvis that are not for segment masses
                
                for i = 1:length(SegInOutput) % save mass changes into Subjects structure
                    Subjects(subj).Trials(trial).RRAmassChanges(i).segment = Segments{i};
                    Subjects(subj).Trials(trial).RRAmassChanges(i).value = ...
                        str2double(cell2mat(TxtData(SegInOutput(i)+CellDif)));
                end
                
                % Update model segment masses based on RRA results
                newModel = Model(model); % create new model
                for i = 1:length(SegInOutput)
                    Seg = Segments{i};
                    Segment = newModel.getBodySet().get(Seg);
                    Segment.setMass(Subjects(subj).Trials(trial).RRAmassChanges(i).value);
                end
                
                % save model as osim file
                RRAUpdatedModel = strcat(RRAFolder2, '\', Subjects(subj).name, '_RRA_Post1.osim');
                newModel.print(RRAUpdatedModel);
                
                clearvars TxtData CellDif R2D2 StarWars TxtData1 fid
                
                %% repeat RRA (RRA2)
                import org.opensim.modeling.*
                import java.io.*
                
                % copy RRA1 setup file
                RRAfile.Setup2 = strcat(RRAFolder, '\', SubjName, '_Setup_RRA2.xml');
                copyfile(RRAfile.Setup, RRAfile.Setup2);
                
                % load RRA2 xml and change attributes
                [rraXML, rraRootName, ~] = xml_read(RRAfile.Setup2);
                
                % make stride and side specific results directory within folder 2
                RRATrialFolder = strcat(RRAFolder2, '\', TrialName, '\', Side{side}, '_', num2str(Stride));
                mkdir(RRATrialFolder);
                
                % call for new model and results directory
                rraXML.RRATool.model_file = strcat(RRAUpdatedModel);  % use model output from IK
                rraXML.RRATool.results_directory = RRATrialFolder;
                % same actuators, tasks, and external loads files
                rraXML.RRATool.force_set_files = RRAfile.ActuatorFile;
                rraXML.RRATool.task_set_file = RRAfile.Tasks;
                rraXML.RRATool.external_loads_file = ExtLdsTrial;
                
                % update RRA settings for specific trial
                rraXML.RRATool.initial_time = Start; % start time
                rraXML.RRATool.final_time = End; % end time
                rraXML.RRATool.desired_kinematics_file = IKfile.Output;
                rraXML.RRATool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                
                OutputModelFile = strcat(Subjects(subj).name, '_', Subjects(subj).Trials(trial).name, '_RRA_AdjPost2.osim');
                rraXML.RRATool.output_model_file = OutputModelFile; % output model
                rraXML.RRATool.ATTRIBUTE = strcat(Subjects(subj).name, '_RRA2');
                RRAfile.Setup2 = strcat(RRAFolder2, '\', Subjects(subj).name, '_Setup_RRA2.xml');
                RRAfile.Trial2 = strcat(RRAfile.Setup2(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific RRA setup XML file
                xml_write(RRAfile.Trial2, rraXML, rraRootName);  % write new RRA xml specific to trial
                
                rraTool = RRATool(RRAfile.Trial2); % specify new setup file
                rraTool.setAdjustCOMToReduceResiduals(1);
                
                RRAfile.MassChanges2 = strcat(RRATrialFolder, '\', SubjTrialSideStride, '_MassChanges.txt');
                fopen(RRAfile.MassChanges2,'w'); % clear old scale results
                diary(RRAfile.MassChanges2);
                diary on;
                rraTool.run(); % run RRA2
                diary off;
                
                % make sure mass changes are < 0.5 kg
                % read in segment masses from txt file to change
                %                 fid = fopen(RRAfile.MassChanges2, 'r');
                %                 TxtData1 = textscan(fid,'%q');
                %                 TxtDataString = TxtData1{1,1};
                %                 [StarWars,~] = find(contains(TxtDataString, '*****'));
                %                 fclose(fid);
                %
                %                 TxtData = TxtDataString(StarWars(1):end);
                %                 CellDif = 8; % number of cells between segment name and new mass value
                
                %% Add metabolic probes to RRA-updated model for CMC metabolics
                if strcmp(Settings.CMC, 'Yes')
                    clc;
                    import org.opensim.modeling.*
                    import java.io.*
                    
                    % define model for CMC
                    CMCmodel = Model(RRAUpdatedModel);
                    CMCmodel.initSystem();
                    
                    % get slow twitch data from text file
                    fn = 'C:\Users\richa\Documents\OpenSim 4.0\Metabolix\Scripts\metabolicsSlowTwitchRatios_Gait2392.txt';
                    fdata = importdata(fn);
                    muscleName = fdata.textdata;
                    twitchRatio = fdata.data;
                    
                    % The following booleans are constructor arguments for the Umberger probe.
                    % These settings are used for all probes.
                    activationMaintenanceRateOn = true;
                    maintenanceRateOn = true;
                    shorteningRateOn = true;
                    basalRateOn = false;
                    mechanicalWorkRateOn = true;
                    reportTotalMetabolicsOnly = false;
                    
                    % ADD UMBERGER METABOLIC PROBES
                    % The mass of each muscle will be calculated using data from the model:
                    % muscleMass = (maxIsometricForce / sigma) * rho * optimalFiberLength
                    % where sigma = 0.25e6 is the specific tension of mammalian muscle (in
                    % Pascals) and rho = 1059.7 is the density of mammalian muscle (in kg/m^3).
                    
                    % The slow-twitch ratio used for muscles that either do not appear in the
                    % file, or appear but whose proportion of slow-twitch fibers is unknown.
                    defaultTwitchRatio = 0.5;
                    % Get the slow-twitch ratio as the default value.
                    slowTwitchRatio = ones(1,length(twitchRatio)) * defaultTwitchRatio;
                    
                    % Define a whole-body probe that will report the total metabolic energy
                    % consumption over the simulation.
                    UMBProbe = Umberger2010MuscleMetabolicsProbe(...
                        activationMaintenanceRateOn,...
                        shorteningRateOn,...
                        basalRateOn,...
                        mechanicalWorkRateOn);
                    UMBProbe.setOperation("value");
                    UMBProbe.set_report_total_metabolics_only(reportTotalMetabolicsOnly);
                    
                    % Add the probe to the model and provide a name.
                    CMCmodel.addProbe(UMBProbe);
                    UMBProbe.setName("metabolics");
                    
                    % ADD BHARGAVA PROBES
                    % from Bhargava et al. 2004. DOI: 10.1016/s0021-9290(03)00239-2
                    % https://www.ncbi.nlm.nih.gov/pubmed/14672571
                    
                    BHARProbe = Bhargava2004MuscleMetabolicsProbe(...
                        activationMaintenanceRateOn,...
                        maintenanceRateOn,...
                        shorteningRateOn,...
                        basalRateOn,...
                        mechanicalWorkRateOn);
                    BHARProbe.setOperation("value");
                    BHARProbe.set_report_total_metabolics_only(reportTotalMetabolicsOnly);
                    
                    % Add the probe to the model and provide a name.
                    CMCmodel.addProbe(BHARProbe);
                    BHARProbe.setName("metabolics");
                    
                    % add muscle twitch parameters
                    ratio_slow_twitch_fibers = 0.5;
                    activation_constant_slow_twitch = 40;
                    activation_constant_fast_twitch = 133;
                    maintenance_constant_slow_twitch = 74;
                    maintenance_constant_fast_twitch = 111;
                    
                    % The default setting is false, in which case, the muscle mass is calculated from the following formula:
                    % m = (Fmax/specific_tension)*density*Lm_opt,
                    % where specific_tension and density are properties defined above
                    % (note that their default values are set based on mammalian muscle, 0.25e6 N/m^2 and 1059.7 kg/m^3, respectively);
                    % Fmax and Lm_opt are the maximum isometric force and optimal fiber length, respectively, of the muscle.
                    use_provided_muscle_mass = false;
                    
                    % Loop through all muscles, adding parameters for each into the whole-body probe.
                    for j = 1:CMCmodel.getMuscles().getSize()
                        
                        % offset for python/C++ looping in OpenSim API
                        iMuscle = j -1;
                        % save muscle name in text file via diary function
                        diary('temp.txt'); % open diary in command window
                        thisMuscle = CMCmodel.getMuscles().get(iMuscle) % allow printing of muscle name into command window
                        diary off
                        
                        % get muscle name and save
                        FID = fopen('temp.txt');
                        TXT = textscan(FID, '%s');
                        if strcmp(TXT{1}(end), 'off')
                            TXT{1}(end) = [];
                        end
                        if strcmp(TXT{1}(end), 'diary')
                            TXT{1}(end) = [];
                        end
                        txt = TXT{1}(end);
                        Muscle = txt{1}(1:end-2);
                        
                        % UMBERGER PROBE
                        % Set the slow-twitch ratio to the physiological value, if it is known.
                        %     for i = 1:length(twitchRatio) % in twitchRatios.items():
                        ind = find(contains({muscleName{:}}, Muscle));
                        if strcmp(Muscle, muscleName{ind}) && twitchRatio(ind) ~= -1
                            slowTwitchRatio(ind) = twitchRatio(ind);
                        end
                        
                        % Add this muscle to the whole-body probe. The arguments are muscle
                        % name, slow-twitch ratio, and muscle mass. Note that the muscle mass
                        % is ignored unless we set useProvidedMass to True.
                        UMBProbe.addMuscle(thisMuscle.getName(),...
                            slowTwitchRatio(ind));
                        
                        % BHARGAVA PROBE
                        BHARProbe.addMuscle(thisMuscle.getName(),...
                            ratio_slow_twitch_fibers,...
                            activation_constant_slow_twitch,...
                            activation_constant_fast_twitch,...
                            maintenance_constant_slow_twitch,...
                            maintenance_constant_fast_twitch);
                    end
                    
                    % Save the new model to a file
                    CMCModelPath = strcat(CMCFolder, '\', SubjTrialSideStride, '_CMC.osim');
                    CMCmodel.print(CMCModelPath);
                    
                    
                    %% CMC
                    clc;
                    import org.opensim.modeling.*
                    import java.io.*
                    
                    % copy over generic files into general CMC folder
                    Orig.CMC_Setup = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'CMC_Setup')).name);
                    CMCfile.Setup = strcat(CMCFolder, '\', Subjects(subj).name, '_Setup_CMC.xml');
                    copyfile(Orig.CMC_Setup , CMCfile.Setup);  % CMC setup file
                    Orig.CMC_Actuators = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'CMC_Actuators_Locked.xml')).name);
                    CMCfile.Actuators = strcat(CMCFolder, '\', Subjects(subj).name, '_CMC_Actuators_Locked.xml');
                    CMCfile.ActuatorFile = strcat(Subjects(subj).name, '_CMC_Actuators_Locked.xml');
                    copyfile(Orig.CMC_Actuators, CMCfile.Actuators);  % CMC actuators file
                    Orig.CMC_Tasks = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'CMC_Tasks_locked_deleted.xml')).name);
                    CMCfile.Tasks = strcat(CMCFolder, '\', Subjects(subj).name, '_CMC_Tasks.xml');
                    copyfile(Orig.CMC_Tasks, CMCfile.Tasks);  % CMC tasks file
                    Orig.CMC_Constraints = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'CMC_ControlConstraints_locked_deleted.xml')).name);
                    CMCfile.Constraints = strcat(CMCFolder, '\', Subjects(subj).name, '_CMC_Constraints.xml');
                    copyfile(Orig.CMC_Constraints, CMCfile.Constraints);  % CMC control contraints
                    
                    RRA2Dir = dir(RRATrialFolder); % get RRA2 results directory
                    
                    % create trial specific CMC folder
                    CMCTrialFolder = strcat(CMCFolder, '\', TrialName, '\', Side{side}, '_', num2str(Stride));
                    mkdir(CMCTrialFolder);
                    
                    % clear CMCTrialFolderContents to prevent extra copies of old data
                    CMCTrialFolderDir = dir(CMCTrialFolder);
                    ToDel = find(~[CMCTrialFolderDir.isdir]);
                    %                 ToDel(contains({CMCTrialFolderDir.name}, 'states')) == 0 % dont delete states
                    for i = ToDel
                        delete(strcat(CMCTrialFolder, '\', CMCTrialFolderDir(i).name));
                    end
                    
                    % % load CMC xml and change attributes
                    [cmcXML, cmcRootName, ~] = xml_read(CMCfile.Setup);
                    
                    % call for new model and results directory
                    cmcXML.CMCTool.model_file = CMCModelPath; %CMCModel;
                    cmcXML.CMCTool.results_directory = CMCTrialFolder;
                    % same actuators, tasks, and external loads files
                    cmcXML.CMCTool.force_set_files = CMCfile.ActuatorFile;
                    cmcXML.CMCTool.task_set_file = CMCfile.Tasks;
                    cmcXML.CMCTool.external_loads_file = ExtLdsTrial;
                    
                    % % update CMC settings for specific trial
                    % set gait cycle times
                    if strcmp(Side{side}, 'Left')
                        Start = Subjects(subj).Trials(trial).Times.Start_Left(Stride);
                        End = Subjects(subj).Trials(trial).Times.End_Left(Stride);
                    elseif strcmp(Side{side}, 'Right')
                        Start = Subjects(subj).Trials(trial).Times.Start_Right(Stride);
                        End = Subjects(subj).Trials(trial).Times.End_Right(Stride);
                    end
                    cmcXML.CMCTool.initial_time = Start; % start time
                    cmcXML.CMCTool.final_time = End;  % end time
                    cmcXML.CMCTool.AnalysisSet.objects.Kinematics.start_time = Start;
                    cmcXML.CMCTool.AnalysisSet.objects.Kinematics.end_time = End;
                    cmcXML.CMCTool.AnalysisSet.objects.Actuation.start_time = Start;
                    cmcXML.CMCTool.AnalysisSet.objects.Actuation.end_time = End;
                    RRA2KinematicsQ = strcat(RRATrialFolder, '\', RRA2Dir(contains({RRA2Dir.name}, 'Kinematics_q.sto')).name);
                    cmcXML.CMCTool.desired_kinematics_file = RRA2KinematicsQ;
                    cmcXML.CMCTool.constraints_file = CMCfile.Constraints;
                    % no filter on Q-file kinematics
                    %     cmcXML.CMCTool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                    
                    OutputModelFile = strcat(SubjTrialSideStride,  '_CMC_Adj.osim'); % output model
                    cmcXML.CMCTool.output_model_file = OutputModelFile;
                    cmcXML.CMCTool.ATTRIBUTE = strcat(Subjects(subj).name, '_CMC');
                    %                 CMCfile.Setup = strcat(CMCFolder, '\', Subjects(subj).name, '_Setup_CMC.xml');
                    CMCfile.Trial = strcat(CMCFolder, '\',  SubjTrialSideStride, '_Setup_CMC.xml'); % name trial-specific CMC setup XML file
                    xml_write(CMCfile.Trial, cmcXML, cmcRootName);
                    
                    cmcTool = CMCTool(CMCfile.Trial); % load CMC tool with new setup file
                    
                    CMCResults = strcat(CMCTrialFolder, '\', SubjTrialSideStride, '_CMC_Results.txt');
                    fopen(CMCResults,'w'); % clear old scale results
                    diary(CMCResults);
                    diary on;
                    cmcTool.run(); % run CMC
                    diary off;
                    fclose('all');
                    
                    % rename CMC results files
                    CMCdir = dir(CMCTrialFolder);
                    Output = {CMCdir(~[CMCdir.isdir]).name};
                    
                    % look for outputs
                    for j = 1:length(Output)
                        if contains(Output{j}, 'states')
                            continue
                        else
                            Src = strcat(CMCdir(1).folder, '\', Output{j});
                            Dst = strcat(CMCdir(1).folder, '\', SubjTrialSideStride, '_', Output{j});
                            movefile(Src, Dst);
                        end
                    end
                    clearvars CMCdir Name Src Dst Start End
                    
                    %%  SO
                    import org.opensim.modeling.*
                    import java.io.*
                    
                    if strcmp(Settings.SO, 'Yes')
                        % copy over generic files into general SO folder
                        Orig.SO_Setup = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'SO_Setup')).name);
                        SOfile.Setup = strcat(SOFolder, '\', Subjects(subj).name, '_Setup_SO.xml');
                        copyfile(Orig.SO_Setup , SOfile.Setup);  % SO setup file
                        Orig.SO_Actuators = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'CMC_Actuators_NoTorsoX')).name);
                        SOfile.Actuators = strcat(SOFolder, '\', Subjects(subj).name, '_CMC_Actuators_NoTorsoX.xml');
                        SOfile.ActuatorFile = strcat(Subjects(subj).name, '_CMC_Actuators_NoTorsoX.xml');
                        copyfile(Orig.SO_Actuators, SOfile.Actuators);  % SO actuators file
                        %     Orig.SO_Tasks = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'SO_Tasks_locked_deleted.xml')).name);
                        %     SOfile.Tasks = strcat(SOFolder, '\', Subjects(subj).name, '_SO_Tasks.xml');
                        %     copyfile(Orig.SO_Tasks, SOfile.Tasks);  % SO tasks file
                        
                        % SO control contraints
                        
                        if strcmp(Settings.LeftLeg, 'Yes')
                            % create trial specific SO folder
                            SOTrialFolder = strcat(SOFolder, '\', TrialName, '\Left_1');
                            mkdir(SOTrialFolder);
                            
                            % no filter on Q-file kinematics
                            
                            % % load SO2 xml and change attributes
                            [soXML, soRootName, ~] = xml_read(SOfile.Setup);
                            
                            RRA2Dir = dir(RRATrialFolder);
                            % call for new model and results directory
                            soXML.AnalyzeTool.model_file = strcat(RRAFolder2, '\', Subjects(subj).name, '_RRA_AdjPost2.osim');  % use model output from IK
                            soXML.AnalyzeTool.results_directory = SOTrialFolder;
                            % same actuators, tasks, and external loads files
                            soXML.AnalyzeTool.force_set_files = SOfile.ActuatorFile;
                            %     soXML.AnalyzeTool.task_set_file = SOfile.Tasks;
                            soXML.AnalyzeTool.external_loads_file = ExtLdsTrial;
                            % % update SO settings for specific trial
                            % set for first left gait cycle
                            Ind = find(Subjects(subj).Trials(trial).TSData.L_Strike(:,2) > Subjects(subj).Trials(trial).initial_time, 3);
                            L_Start = Subjects(subj).Trials(trial).TSData.L_Strike(Ind(1), 2) - 0.05;
                            L_End = Subjects(subj).Trials(trial).TSData.L_Strike(Ind(1)+1, 2) + 0.05;
                            soXML.AnalyzeTool.initial_time = L_Start; % start time
                            soXML.AnalyzeTool.final_time = L_End; % end time
                            RRA2KinematicsQ = strcat(RRATrialFolder, '\', RRA2Dir(contains({RRA2Dir.name}, 'Kinematics_q.sto')).name);
                            soXML.AnalyzeTool.coordinates_file = RRA2KinematicsQ;
                            soXML.AnalyzeTool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                            
                            soXML.AnalyzeTool.output_model_file = strcat(Subjects(subj).name, '_SO_Adj.osim'); % output model
                            soXML.AnalyzeTool.ATTRIBUTE = strcat(Subjects(subj).name, '_SO');
                            SOfile.Setup = strcat(SOFolder, '\', Subjects(subj).name, '_Setup_SO.xml');
                            SOfile.Trial = strcat(SOfile.Setup(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific SO setup XML file
                            xml_write(SOfile.Trial, soXML, soRootName);
                            
                            % scale optimal force or max force to body mass?
                            % default settings time 1.5 body weight?
                            
                            staticoptimization = AnalyzeTool(SOfile.Trial); % specify new
                            staticoptimization.run(); % run SO
                        end
                        
                        % RIGHT Stride SO
                        
                        % create trial specific SO folder
                        SOTrialFolder = strcat(SOFolder, '\', TrialName, '\Right_1');
                        mkdir(SOTrialFolder);
                        
                        % no filter on Q-file kinematics
                        
                        % % load SO2 xml and change attributes
                        [soXML, soRootName, ~] = xml_read(SOfile.Setup);
                        
                        RRA2Dir = dir(RRATrialFolder);
                        % call for new model and results directory
                        soXML.AnalyzeTool.model_file = strcat(RRAFolder2, '\', Subjects(subj).name, '_RRA_AdjPost2.osim');  % use model output from IK
                        soXML.AnalyzeTool.results_directory = SOTrialFolder;
                        % same actuators, tasks, and external loads files
                        soXML.AnalyzeTool.force_set_files = SOfile.ActuatorFile;
                        %     soXML.AnalyzeTool.task_set_file = SOfile.Tasks;
                        soXML.AnalyzeTool.external_loads_file = ExtLdsTrial;
                        % % update SO settings for specific trial
                        % set for first left gait cycle
                        Ind = find(Subjects(subj).Trials(trial).TSData.R_Strike(:,2) > Subjects(subj).Trials(trial).initial_time, 3);
                        R_Start = Subjects(subj).Trials(trial).TSData.R_Strike(Ind(1), 2) - 0.05;
                        R_End = Subjects(subj).Trials(trial).TSData.R_Strike(Ind(1)+1, 2) + 0.05;
                        soXML.AnalyzeTool.initial_time = R_Start; % start time
                        soXML.AnalyzeTool.final_time = R_End; % end time
                        RRA2KinematicsQ = strcat(RRATrialFolder, '\', RRA2Dir(contains({RRA2Dir.name}, 'Kinematics_q.sto')).name);
                        soXML.AnalyzeTool.coordinates_file = RRA2KinematicsQ;
                        soXML.AnalyzeTool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                        
                        soXML.AnalyzeTool.output_model_file = strcat(Subjects(subj).name, '_SO_Adj.osim'); % output model
                        soXML.AnalyzeTool.ATTRIBUTE = strcat(Subjects(subj).name, '_SO');
                        SOfile.Setup = strcat(SOFolder, '\', Subjects(subj).name, '_Setup_SO.xml');
                        SOfile.Trial = strcat(SOfile.Setup(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific SO setup XML file
                        xml_write(SOfile.Trial, soXML, soRootName);
                        
                        % scale optimal force or max force to body mass?
                        % default settings time 1.5 body weight?
                        
                        staticoptimization = AnalyzeTool(SOfile.Trial); % specify new
                        try
                            staticoptimization.run(); % run SO
                        catch
                            FailedSO(subj,trial) =1;
                        end
                    else
                        SOFolder = strcat(Subjects(subj).Folders.OpenSimFolder , '\SO_Files');
                        SOTrialFolder = strcat(SOFolder, '\', TrialName, '\Right_1');
                    end
                    
                    %% MuscleAnalysis
                    if strcmp(Settings.MA, 'Yes')
                        % copy over generic files into general MA folder
                        Orig.MA_Setup = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'MA_Setup')).name);
                        MAfile.Setup = strcat(MAFolder, '\', Subjects(subj).name, '_Setup_MA.xml');
                        copyfile(Orig.MA_Setup , MAfile.Setup);  % MA setup file
                        Orig.MA_Actuators = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'CMC_Actuators_NoTorsoX')).name);
                        MAfile.Actuators = strcat(MAFolder, '\', Subjects(subj).name, '_CMC_Actuators_NoTorsoX.xml');
                        MAfile.ActuatorFile = strcat(Subjects(subj).name, '_CMC_Actuators_NoTorsoX.xml');
                        copyfile(Orig.MA_Actuators, MAfile.Actuators);  % MA actuators file
                        %     Orig.MA_Tasks = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'MA_Tasks_locked_deleted.xml')).name);
                        %     MAfile.Tasks = strcat(MAFolder, '\', Subjects(subj).name, '_MA_Tasks.xml');
                        %     copyfile(Orig.MA_Tasks, MAfile.Tasks);  % MA tasks file
                        
                        % MA control contraints
                        
                        if strcmp(Settings.LeftLeg, 'Yes')
                            % create trial specific MA folder
                            MATrialFolder = strcat(MAFolder, '\', TrialName, '\Left_1');
                            mkdir(MATrialFolder);
                            
                            % no filter on Q-file kinematics
                            
                            % % load MA2 xml and change attributes
                            [maXML, maRootName, ~] = xml_read(MAfile.Setup);
                            
                            RRA2Dir = dir(RRATrialFolder);
                            % call for new model and results directory
                            maXML.AnalyzeTool.model_file = strcat(RRAFolder2, '\', Subjects(subj).name, '_RRA_AdjPost2.osim');  % use model output from IK
                            maXML.AnalyzeTool.results_directory = MATrialFolder;
                            % same actuators, tasks, and external loads files
                            maXML.AnalyzeTool.force_set_files = MAfile.ActuatorFile;
                            %     maXML.AnalyzeTool.task_set_file = MAfile.Tasks;
                            maXML.AnalyzeTool.external_loads_file = ExtLdsTrial;
                            % % update MA settings for specific trial
                            % set for first left gait cycle
                            Ind = find(Subjects(subj).Trials(trial).TSData.L_Strike(:,2) > Subjects(subj).Trials(trial).initial_time, 3);
                            L_Start = Subjects(subj).Trials(trial).TSData.L_Strike(Ind(1), 2) - 0.05;
                            L_End = Subjects(subj).Trials(trial).TSData.L_Strike(Ind(1)+1, 2) + 0.05;
                            maXML.AnalyzeTool.initial_time = L_Start; % start time
                            maXML.AnalyzeTool.final_time = L_End; % end time
                            RRA2KinematicsQ = strcat(RRATrialFolder, '\', RRA2Dir(contains({RRA2Dir.name}, 'Kinematics_q.sto')).name);
                            maXML.AnalyzeTool.coordinates_file = RRA2KinematicsQ;
                            
                            SODir = dir(SOTrialFolder);
                            SOControlsFiles = strcat(SOTrialFolder, '\', SODir(contains({SODir.name}, 'StaticOptimization_controls.xml')).name);
                            maXML.AnalyzeTool.controls_file = SOControlsFiles;
                            maXML.AnalyzeTool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                            
                            maXML.AnalyzeTool.output_model_file = strcat(Subjects(subj).name, '_MA_Adj.osim'); % output model
                            maXML.AnalyzeTool.ATTRIBUTE = strcat(Subjects(subj).name, '_MA');
                            MAfile.Setup = strcat(MAFolder, '\', Subjects(subj).name, '_Setup_MA.xml');
                            MAfile.Trial = strcat(MAfile.Setup(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific MA setup XML file
                            xml_write(MAfile.Trial, maXML, maRootName);
                            
                            % scale optimal force or max force to body mass?
                            % default settings time 1.5 body weight?
                            
                            muscleanalysis = AnalyzeTool(MAfile.Trial); % specify new
                            muscleanalysis.run(); % run MA
                        end
                        
                        % RIGHT Stride MA
                        
                        % create trial specific MA folder
                        MATrialFolder = strcat(MAFolder, '\', TrialName, '\Right_1');
                        mkdir(MATrialFolder);
                        
                        % no filter on Q-file kinematics
                        
                        % % load MA2 xml and change attributes
                        [maXML, maRootName, ~] = xml_read(MAfile.Setup);
                        
                        RRA2Dir = dir(RRATrialFolder);
                        % call for new model and results directory
                        SODir = dir(SOTrialFolder);
                        SOControlsFiles = strcat(SOTrialFolder, '\', SODir(contains({SODir.name}, 'StaticOptimization_controls.xml')).name);
                        maXML.Controllers.controls_file = SOControlsFiles;
                        maXML.AnalyzeTool.model_file = strcat(RRAFolder2, '\', Subjects(subj).name, '_RRA_AdjPost2.osim');  % use model output from IK
                        maXML.AnalyzeTool.results_directory = MATrialFolder;
                        % same actuators, tasks, and external loads files
                        maXML.AnalyzeTool.force_set_files = MAfile.ActuatorFile;
                        %     maXML.AnalyzeTool.task_set_file = MAfile.Tasks;
                        maXML.AnalyzeTool.external_loads_file = ExtLdsTrial;
                        % % update MA settings for specific trial
                        % set for first left gait cycle
                        Ind = find(Subjects(subj).Trials(trial).TSData.R_Strike(:,2) > Subjects(subj).Trials(trial).initial_time, 3);
                        R_Start = Subjects(subj).Trials(trial).TSData.R_Strike(Ind(1), 2) - 0.05;
                        R_End = Subjects(subj).Trials(trial).TSData.R_Strike(Ind(1)+1, 2) + 0.05;
                        maXML.AnalyzeTool.initial_time = R_Start; % start time
                        maXML.AnalyzeTool.final_time = R_End; % end time
                        
                        maXML.AnalyzeTool.start_time = R_Start; % start time
                        maXML.AnalyzeTool.end_time = R_End; % end time
                        
                        RRA2KinematicsQ = strcat(RRATrialFolder, '\', RRA2Dir(contains({RRA2Dir.name}, 'Kinematics_q.sto')).name);
                        maXML.AnalyzeTool.coordinates_file = RRA2KinematicsQ;
                        maXML.AnalyzeTool.lowpass_cutoff_frequency_for_coordinates = Settings.LowPassKinematicsFilter;
                        
                        %             maXML.AnalyzeTool.output_model_file = strcat(Subjects(subj).name, '_MA_Adj.osim'); % output model
                        maXML.AnalyzeTool.ATTRIBUTE = strcat(Subjects(subj).name, '_MA');
                        MAfile.Setup = strcat(MAFolder, '\', Subjects(subj).name, '_Setup_MA.xml');
                        MAfile.Trial = strcat(MAfile.Setup(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific MA setup XML file
                        xml_write(MAfile.Trial, maXML, maRootName);
                        %             keyboard
                        % scale optimal force or max force to body mass?
                        % default settings time 1.5 body weight?
                        
                        muscleanalysis = AnalyzeTool(MAfile.Trial); % specify new
                        try
                            muscleanalysis.run(); % run MA
                        catch
                            FailedMA(subj,trial) = 1;
                        end
                    end
                    
                    %% Inverse Dynamics
                    if strcmp(Settings.ID, 'Yes')
                        % copy over generic files into general SO folder
                        Orig.ID_Setup = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'ID_Setup')).name);
                        IDfile.Setup = strcat(IDFolder, '\', Subjects(subj).name, '_Setup_ID.xml');
                        copyfile(Orig.ID_Setup , IDfile.Setup);  % ID setup file
                        %             Orig.ID_Actuators = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'CMC_Actuators_NoTorsoX')).name);
                        %             IDfile.Actuators = strcat(IDFolder, '\', Subjects(subj).name, '_CMC_Actuators_NoTorsoX.xml');
                        %             IDfile.ActuatorFile = strcat(Subjects(subj).name, '_CMC_Actuators_NoTorsoX.xml');
                        %             copyfile(Orig.ID_Actuators, IDfile.Actuators);  % ID actuators file
                        %     Orig.ID_Tasks = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'ID_Tasks_locked_deleted.xml')).name);
                        %     IDfile.Tasks = strcat(IDFolder, '\', Subjects(subj).name, '_ID_Tasks.xml');
                        %     copyfile(Orig.ID_Tasks, IDfile.Tasks);  % ID tasks file
                        
                        % ID control contraints
                        
                        if strcmp(Settings.LeftLeg, 'Yes')
                            % create trial specific ID folder
                            IDTrialFolder = strcat(IDFolder, '\', TrialName, '\Left_1');
                            mkdir(IDTrialFolder);
                            
                            % no filter on Q-file kinematics
                            
                            % % load ID2 xml and change attributes
                            [idXML, idRootName, ~] = xml_read(IDfile.Setup);
                            
                            RRA2Dir = dir(RRATrialFolder);
                            % call for new model and results directory
                            idXML.InverseDynamicsTool.model_file = strcat(RRAFolder2, '\', Subjects(subj).name, '_RRA_AdjPost2.osim');  % use model output from IK
                            idXML.InverseDynamicsTool.results_directory = IDTrialFolder;
                            % same actuators, tasks, and external loads files
                            %                 idXML.InverseDynamicsTool.force_set_files = IDfile.ActuatorFile;
                            %     idXML.InverseDynamicsTool.task_set_file = IDfile.Tasks;
                            idXML.InverseDynamicsTool.external_loads_file = ExtLdsTrial;
                            % % update ID settings for specific trial
                            % set for first left gait cycle
                            Ind = find(Subjects(subj).Trials(trial).TSData.L_Strike(:,2) > Subjects(subj).Trials(trial).initial_time, 3);
                            L_Start = Subjects(subj).Trials(trial).TSData.L_Strike(Ind(1), 2) - 0.05;
                            L_End = Subjects(subj).Trials(trial).TSData.L_Strike(Ind(1)+1, 2) + 0.05;
                            idXML.InverseDynamicsTool.time_range = [L_Start,L_End]; % start time
                            %                 idXML.InverseDynamicsTool.initial_time = L_Start; % start time
                            %                 idXML.InverseDynamicsTool.final_time = L_End; % end time
                            RRA2KinematicsQ = strcat(RRATrialFolder, '\', RRA2Dir(contains({RRA2Dir.name}, 'Kinematics_q.sto')).name);
                            idXML.InverseDynamicsTool.coordinates_file = RRA2KinematicsQ;
                            idXML.InverseDynamicsTool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                            
                            %                 idXML.InverseDynamicsTool.output_model_file = strcat(Subjects(subj).name, '_ID_Adj.osim'); % output model
                            idXML.InverseDynamicsTool.ATTRIBUTE = strcat(Subjects(subj).name, '_ID');
                            IDfile.Setup = strcat(IDFolder, '\', Subjects(subj).name, '_Setup_ID.xml');
                            IDfile.Trial = strcat(IDfile.Setup(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific ID setup XML file
                            xml_write(IDfile.Trial, idXML, idRootName);
                            
                            % scale optimal force or max force to body mass?
                            % default settings time 1.5 body weight?
                            
                            inversedynamics = InverseDynamicsTool(IDfile.Trial); % specify new
                            inversedynamics.run(); % run ID
                        end
                        
                        % RIGHT Stride ID
                        
                        % create trial specific ID folder
                        IDTrialFolder = strcat(IDFolder, '\', TrialName, '\Right_1');
                        mkdir(IDTrialFolder);
                        
                        % no filter on Q-file kinematics
                        
                        % % load ID2 xml and change attributes
                        [idXML, idRootName, ~] = xml_read(IDfile.Setup);
                        
                        RRA2Dir = dir(RRATrialFolder);
                        % call for new model and results directory
                        idXML.InverseDynamicsTool.model_file = strcat(RRAFolder2, '\', Subjects(subj).name, '_RRA_AdjPost2.osim');  % use model output from IK
                        idXML.InverseDynamicsTool.results_directory = IDTrialFolder;
                        % same actuators, tasks, and external loads files
                        %             idXML.InverseDynamicsTool.force_set_files = IDfile.ActuatorFile;
                        %     idXML.InverseDynamicsTool.task_set_file = IDfile.Tasks;
                        idXML.InverseDynamicsTool.external_loads_file = ExtLdsTrial;
                        % % update ID settings for specific trial
                        % set for first left gait cycle
                        Ind = find(Subjects(subj).Trials(trial).TSData.R_Strike(:,2) > Subjects(subj).Trials(trial).initial_time, 3);
                        R_Start = Subjects(subj).Trials(trial).TSData.R_Strike(Ind(1), 2) - 0.05;
                        R_End = Subjects(subj).Trials(trial).TSData.R_Strike(Ind(1)+1, 2) + 0.05;
                        idXML.InverseDynamicsTool.time_range = [R_Start, R_End]; % start time
                        %             idXML.InverseDynamicsTool.initial_time = R_Start; % start time
                        %             idXML.InverseDynamicsTool.final_time = R_End; % end time
%                         RRA2KinematicsQ = strcat(RRATrialFolder, '\', RRA2Dir(contains({RRA2Dir.name}, 'Kinematics_q.sto')).name);
%                         idXML.InverseDynamicsTool.coordinates_file = RRA2KinematicsQ;
                        idXML.InverseDynamicsTool.lowpass_cutoff_frequency = Settings.LowPassKinematicsFilter;
                        
                        %             idXML.InverseDynamicsTool.output_model_file = strcat(Subjects(subj).name, '_ID_Adj.osim'); % output model
                        idXML.InverseDynamicsTool.ATTRIBUTE = strcat(Subjects(subj).name, '_ID');
                        IDfile.Setup = strcat(IDFolder, '\', Subjects(subj).name, '_Setup_ID.xml');
                        IDfile.Trial = strcat(IDfile.Setup(1:end-4), '_', MarkerFile(1:end-12), '.xml'); % name trial-specific ID setup XML file
                        xml_write(IDfile.Trial, idXML, idRootName);
                        
                        % scale optimal force or max force to body mass?
                        % default settings time 1.5 body weight?
                        
                        inversedynamics = InverseDynamicsTool(IDfile.Trial); % specify new
                        try
                            inversedynamics.run(); % run ID
                        catch
                            FailedID(subj,trial) =1;
                        end
                    end
                    
                end
                
            end % end Stride Loop
            
        end % end Side Loop
        
    end % End Trial Loop
    
end % End Subject loop

%% save matlab metadata to file?

