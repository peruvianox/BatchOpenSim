function[Subjects] = ABL_Scale(Settings, Subjects)

% Batch process OpenSim Scaling

%% Settings
GenericFilePath = Settings.GenericPath;
GenericDir = Settings.GenericDir;

%% Subject Loop
for subj = 1:length(Subjects)
    
    %% Create scale folder and copy over TRC and GRF files
    % create scale folder
    ScaleFolder = strcat(Subjects(subj).Folders.OpenSimFolder, '\ScaleFiles');
    mkdir(ScaleFolder);
    StaticTrials = strcmp({Subjects(subj).Trials.type}, 'static');
    
    % make virtual markers and copy into scale folder
    StaticMarkerFile = strcat(Subjects(subj).Folders.OpenSimFolder, '\', ...
        Subjects(subj).Trials(StaticTrials).files.OpenSimTRC);
    StaticVirtualFile = [ScaleFolder '\' Subjects(subj).Trials(StaticTrials).files.OpenSimTRC(1:end-4) '_Virtual.trc'];
    Osim.MakeVirtualMkr(StaticMarkerFile, StaticVirtualFile);
    %      copyfile([StaticMarkerFile(1:end-4) '_Virtual.trc'], );
    
    disp(strcat('Scaling ', Subjects(subj).name));
    
    
    %% Prep To Scale - Define inputs & outputs
    import org.opensim.modeling.*
    import java.io.*
    
    % define and copy original model OSIM file
    Orig.Model = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name},...
        'gait2392_Scale_ABLMarkerSet_UMBprobed_locked.osim')).name);
    SubjModel = strcat(ScaleFolder, '\gait2392_Scale_ABLMarkerSet_UMBprobed_locked.osim');
    copyfile(Orig.Model, SubjModel);
    
    % Load the model and initialize
    model = Model(fullfile(Orig.Model));
    model.initSystem();
    
    % copy over original markerset file
    Orig.MkrSetFile = strcat(GenericFilePath, '\',...
        GenericDir(contains({GenericDir.name}, 'gait2392_Scale_MarkerSet_ABL_Virtual.xml')).name);
    MkrSetFile = strcat(ScaleFolder, '\', Subjects(subj).name, '_MkrSet.xml');
    copyfile(Orig.MkrSetFile, MkrSetFile); % copy generic scale file to subject setup directory
    
    % identify setup XML file
    Orig.ScaleSetupFile = strcat(GenericFilePath, '\', GenericDir(contains({GenericDir.name}, 'Setup_Scale')).name);
    [scaleXML, scaleRootName, ~] = xml_read(Orig.ScaleSetupFile);
    
    % Change attributes in structure
    scaleXML.ScaleTool.mass = Subjects(subj).Demo.mass;
    scaleXML.ScaleTool.ATTRIBUTE.name = Subjects(subj).name;
    
    % GenericModelMaker
    scaleXML.ScaleTool.GenericModelMaker.model_file = 'gait2392_Scale_ABLMarkerSet_UMBprobed_locked.osim';
    scaleXML.ScaleTool.GenericModelMaker.marker_set_file = strcat(Subjects(subj).name, '_MkrSet.xml');
    
    % only use first 3 frames from static trial
    %     StaticTrialNum = contains({Subjects(subj).Trials.type}, 'static');
    %     ScaleStart = Subjects(subj).Trials(StaticTrialNum).Times.TRC(1);
    %     ScaleEnd = Subjects(subj).Trials(StaticTrialNum).Times.TRC(3);
    
    % ModelScaler
    scaleXML.ScaleTool.ModelScaler.marker_file = StaticVirtualFile;
    OutputModelFile = strcat(Subjects(subj).name, '_Scaled.osim');
    scaleXML.ScaleTool.ModelScaler.output_model_file = OutputModelFile;
    %     scaleXML.ScaleTool.ModelScaler.time_range = [ScaleStart ScaleEnd];
    
    % MarkerPlacer
    scaleXML.ScaleTool.MarkerPlacer.marker_file = StaticVirtualFile;
    scaleXML.ScaleTool.MarkerPlacer.output_model_file = OutputModelFile;
    %     scaleXML.ScaleTool.MarkerPlacer.time_range = [ScaleStart ScaleEnd];
    scaleXML.ScaleTool.MarkerPlacer.output_motion_file = strcat(Subjects(subj).name, '_Scale_Motion.sto');
    scaleXML.ScaleTool.MarkerPlacer.output_marker_file = strcat(Subjects(subj).name, '_Scale_Markers.xml');
    scaleXML.ScaleTool.MarkerPlacer.output_scale_file = strcat(Subjects(subj).name, '_ScaleSet.xml');
    
    % export  XML to SetupDir (specific to each subject)
    SetupScale = strcat(ScaleFolder, '\', Subjects(subj).name, '_Setup_Scale.xml');
    xml_write(SetupScale, scaleXML, scaleRootName);
    
    scale = ScaleTool(SetupScale); % open scaling tool with new attributes
    
    % Run Scaling and 
    scale.run(); % run scaling
    
    % save printed results in opensim log file
    FID = fopen('opensim.log');
    TXT = textscan(FID, '%s');
    % look for each subject's model scaling
    for i = 1:length(TXT{1})
        if contains(TXT{1}(i), {'Deleted'}) &&  ...
                contains(TXT{1}(i+2), {'unused'}) && contains(TXT{1}(i+3), {'markers'})...
                && contains(TXT{1}(i+6), {char(Subjects(subj).name)})
            ind = i;
            break
        end
    end
    % extract marker error locations
    for i = ind:ind+50
        if contains(TXT{1}(i), {'total'}) && contains(TXT{1}(i+1), {'square'}) && ...
                contains(TXT{1}(i+2), {'error'}) && contains(TXT{1}(i+3), {'='})
            Ind = i;
            break
        end
    end
    
    % save scale sets
    
    % save marker errors
    Subjects(subj).Trials(StaticTrials).ScaleErr.TotalSqErr = str2double(TXT{1}(Ind+4));
    Subjects(subj).Trials(StaticTrials).ScaleErr.RMSErr = str2double(TXT{1}(Ind+9));
    Subjects(subj).Trials(StaticTrials).ScaleErr.MaxErr = str2double(TXT{1}(Ind+12));
    Subjects(subj).Trials(StaticTrials).ScaleErr.MaxMkr = TXT{1}(Ind+13);
    
    disp(['Total Square Error = ' TXT{1}(Ind+4)]);
    disp(' ');
    
    clearvars ScaleRootName Geopath StaticMarkerFile NewStaticMarkerFile...
        StaticMarkerFileGRF NewStaticMarkerFileGRF i col scale Ind String
    
end % end subject loop


end