function[Subjects] = ABL_Scale(Settings, Subjects)

% Batch process OpenSim Scaling

%% Settings
import org.opensim.modeling.*
import java.io.*

GenericFilePath = Settings.GenericPath;
GenericDir = Settings.GenericDir;

%% Subject Loop
clc;
for subj = 1:length(Subjects)

    %% Create scale folder and copy over TRC and GRF files
    % create scale folder
    ScaleFolder = strcat(Subjects(subj).Folders.OpenSimFolder, '\ScaleFiles');
    mkdir(ScaleFolder);
    StaticTrials = strcmp({Subjects(subj).Trials.type}, 'static');
    disp(strcat('Scaling ', Subjects(subj).name));
    cd(ScaleFolder);

    % make virtual markers and copy into scale folder
    useHJC = 0;
    if useHJC == 1
        StaticMarkerFile = strcat(Subjects(subj).Folders.path, '\HJC\', ...
            Subjects(subj).Trials(StaticTrials).files.TRC(1:end-4),'_hjc.trc');
    else
        StaticMarkerFile = strcat(Subjects(subj).Folders.OpenSimFolder, '\', ...
            Subjects(subj).Trials(StaticTrials).files.OpenSimTRC);
    end
    StaticVirtualFile = [ScaleFolder '\' Subjects(subj).Trials(StaticTrials).files.OpenSimTRC(1:end-4) '_Virtual.trc'];
    Osim.MakeVirtualMkr(StaticMarkerFile, StaticVirtualFile);

    % add mass if not defined
    if isempty(Subjects(subj).Demo)
        fn = [Subjects(subj).name '_Info.csv'];
        T = readtable(fn);
        r = contains(T{:,1}, 'Mass');
        Subjects(subj).Demo.mass = str2double(char(T{r, 2}));
    end


    %% Prep To Scale - Define inputs & outputs
    % define and copy original model OSIM file
    if strcmp(Settings.Model, 'exo')
        modelName = 'gait2392_ABL_probed_locked_exo.osim';
    elseif strcmp(Settings.Model, 'generic')
        modelName = 'gait2392_Scale_ABLMarkerSet_UMBprobed_locked.osim';
    end
    if strcmp(Settings.CustomModel, 'Yes') % use customized model
        SubjModel = strcat(ScaleFolder, '\', modelName);
    else % or copy from generic
        Orig.Model = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name},modelName)).name);
        SubjModel = strcat(ScaleFolder, '\', modelName);
        copyfile(Orig.Model, SubjModel);
    end

    % Load the model and initialize
    model = Model(fullfile(SubjModel));
    model.initSystem();

    % copy over original markerset file
    if strcmp(Settings.CustomModel, 'Yes') % use customized model
        Orig.MkrSetFile = strcat(GenericFilePath, '\',...
            GenericDir(contains({GenericDir.name}, 'gait2392_Scale_MarkerSet_GT_Virtual.xml')).name);
        MkrSetFile = strcat(ScaleFolder, '\', Subjects(subj).name, '_MkrSet.xml');
    else
        Orig.MkrSetFile = strcat(GenericFilePath, '\',...
            GenericDir(contains({GenericDir.name}, 'gait2392_Scale_MarkerSet_ABL_Virtual.xml')).name);
        MkrSetFile = strcat(ScaleFolder, '\', Subjects(subj).name, '_MkrSet.xml');
    end
    copyfile(Orig.MkrSetFile, MkrSetFile); % copy generic scale file to subject setup directory

    % identify setup XML file
    if strcmp(Settings.Model, 'exo')
        Orig.ScaleSetupFile = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup_Scale_GT_Exo.xml')).name);
    elseif strcmp(Settings.Model, 'generic')
        Orig.ScaleSetupFile = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup_Scale_GT.xml')).name);
    end
    [scaleXML, scaleRootName, ~] = xml_read(Orig.ScaleSetupFile);

    % Change attributes in structure
    scaleXML.ScaleTool.mass = Subjects(subj).Demo.mass;
    scaleXML.ScaleTool.ATTRIBUTE.name = Subjects(subj).name;

    % GenericModelMaker
    scaleXML.ScaleTool.GenericModelMaker.model_file = SubjModel;
    scaleXML.ScaleTool.GenericModelMaker.marker_set_file = strcat(Subjects(subj).name, '_MkrSet.xml');

    % ModelScaler
    scaleXML.ScaleTool.ModelScaler.marker_file = [Subjects(subj).Trials(StaticTrials).files.OpenSimTRC(1:end-4) '_Virtual.trc'];
    OutputModelFile = strcat(Subjects(subj).name, '_Scaled.osim');
    scaleXML.ScaleTool.ModelScaler.output_model_file = OutputModelFile;

    % MarkerPlacer
    scaleXML.ScaleTool.MarkerPlacer.marker_file = [Subjects(subj).Trials(StaticTrials).files.OpenSimTRC(1:end-4) '_Virtual.trc'];
    scaleXML.ScaleTool.MarkerPlacer.output_model_file = OutputModelFile;
    scaleXML.ScaleTool.MarkerPlacer.output_motion_file = strcat(Subjects(subj).name, '_Scale_Motion.sto');
    scaleXML.ScaleTool.MarkerPlacer.output_marker_file = strcat(Subjects(subj).name, '_Scale_Markers.xml');
    scaleXML.ScaleTool.MarkerPlacer.output_scale_file = strcat(Subjects(subj).name, '_ScaleSet.xml');

    % export  XML to SetupDir (specific to each subject)
    SetupScale = strcat(ScaleFolder, '\', Subjects(subj).name, '_Setup_Scale.xml');
    xml_write(SetupScale, scaleXML, scaleRootName);


    %% Run Scaling
    command = ['opensim-cmd run-tool ' SetupScale];

    delete ScaleLog.txt % clear scale log file
    diary ScaleLog.txt; % start recording outputs
    diary on;
    system(command); % run scaling
    diary off;

    %% save printed results in opensim log file
    disp('Analyzing Scale Error')
    FID = fopen('ScaleLog.txt');
    TXT = textscan(FID, '%s');
    fclose(FID);
    ScaleFactors = struct();
    j = 1;
    % look for scaling factors - for OpenSim version 4.4
    for i = 1:length(TXT{1})
        if contains(TXT{1}(i), {'overall'}) &&  ...
                contains(TXT{1}(i+1), {'scale'}) && contains(TXT{1}(i+2), {'factor'})

            Segment = char(TXT{1}(i-13));
            Factor = char(TXT{1}(i+4));
            disp([Segment ' scaled by ' Factor])
            ScaleFactors(j).Segment = Segment;
            ScaleFactors(j).ScaleFactor = str2double(Factor);
            j = j + 1;
        end
    end

    % save scale factors in matlab structure
    Subjects(subj).Scale = ScaleFactors;

    clearvars ScaleRootName Geopath StaticMarkerFile NewStaticMarkerFile...
        StaticMarkerFileGRF NewStaticMarkerFileGRF i col scale Ind String

end % end subject loop

close all;

end