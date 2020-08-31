function [Struct, Subjects] = ExtractData(ResultsFolder, Subjects, type)

% Extract data from files within folder

%% initial settings
% add paths
addpath(genpath('Scripts'));
addpath(genpath('CodeLibrary'));

%% Extract Metabolic data
clc;
Dir = dir(ResultsFolder);
IsStoFile = ~[Dir.isdir];
dbstop if error
j = 1;
Nfiles =  sum([Dir.isdir]==0);
Struct(Nfiles).Subject = [];

for i = 1:length(Dir)
    if IsStoFile(i) == 1
        
        filename = Dir(i).name;
        
        % get stride times from GRF for each trial
        Str = strsplit(filename, '_');
        Struct(j).Subject = Str{1};
        
        TrialName = Str{2};
        Subj = contains({Subjects.name}, Str{1});
        Trial = contains({Subjects(Subj).Trials.name}, Str{2});
        TSData = Subjects(Subj).Trials(Trial).TSData;
        
        % create structure
        Struct(j).Trial = TrialName;
        % label as left or right side
        if contains(Dir(i).name, 'Left')
            Struct(j).Side = 'Left';
        else
            Struct(j).Side = 'Right';
        end
        
        %% Select STO file to input
        % if exist('file', 'var') == 0
        %     [file, Path] = uigetfile('.sto','Select .sto file of Muscle Metabolic Data to load');
        % end
        % addpath(genpath(Path));
        
        % load data
        Data = importdata(filename);
        
        % define filename in data structure
        Data.filename = filename;
        
        %% Filter and Parse Raw Data
        [Data] = FilterAndParseData(Data, TSData);
        
        %% Specific processing depending on data type
        if strcmp(type, 'met')
            %% Check to see what probes are present
            % define duplicate columns of metabolic reports as various types of
            % metabolic models
            UMB2ind = contains([Data.colheaders], '_0');
            BHARind = contains([Data.colheaders], '_1');
            
            % Umberger
            % check to make sure the first and second set of metabolics are the same
            % sum(sum(Data.Parsed(:,UMB1ind) - Data.Parsed(:,UMB2ind)))
            
            % only save one copy of UMB data
            Data.UMB_Parsed = Data.Parsed(:,UMB2ind);
            Data.UMB_Interp = Data.Interp(:,UMB2ind);
            
            % Bhar
            Data.BHAR_Parsed = Data.Parsed(:,BHARind);
            Data.BHAR_Interp = Data.Interp(:,BHARind);
            
            % umberger metabolics
            [UMB.ParsedMuscles, UMB.ParsedJoints] = GetMuscles(Data.UMB_Parsed, Data.colheaders(UMB2ind), type);
            [UMB.InterpMuscles, UMB.InterpJoints] = GetMuscles(Data.UMB_Interp, Data.colheaders(UMB2ind), type);
            
            % bhargava metabolics
            [BHAR.ParsedMuscles, BHAR.ParsedJoints] = GetMuscles(Data.BHAR_Parsed, Data.colheaders(BHARind), type);
            [BHAR.InterpMuscles, BHAR.InterpJoints] = GetMuscles(Data.BHAR_Interp, Data.colheaders(BHARind), type);
            
            
            %% Define variables on left and right sides
            % and define muscles spanning hip, knee, and ankle joints
            Struct(j).Data = Data;
            % get cycle times for computing integrals later on
            Struct(j).Time = Struct(j).Data.Interp(:,1);
            
            % extract gait events for each side
            if strcmp(Struct(j).Side, 'Left')
                [~, Ind] = min(abs(Struct(j).Time(1) - TSData.L_Strike(:,2)));
                Struct(j).GC_Start = TSData.L_Strike(Ind,2);
                [~, Struct(j).GC_StartInd] = min(abs(Struct(j).Time - TSData.L_Strike(Ind,2)));
                Struct(j).GC_End = TSData.L_Strike(Ind+1,2);
                [~, Struct(j).GC_EndInd] = min(abs(Struct(j).Time - TSData.L_Strike(Ind+1,2)));
                Struct(j).GC_Off = TSData.L_Off(Ind,2);
                [~, Struct(j).GC_OffInd] = min(abs(Struct(j).Time - TSData.L_Off(Ind,2)));
                
                OppInd = find(TSData.R_Strike(:, 2) > TSData.L_Strike(Ind, 2), 1);
                Struct(j).GC_2DS = TSData.R_Strike(OppInd, 2);
                [~,Struct(j).GC_2DSInd] =  min(abs(Struct(j).Time - TSData.R_Strike(OppInd,2)));
                Struct(j).GC_1DS = TSData.R_Off(OppInd-1, 2);
                [~, Struct(j).GC_1DSInd] =  min(abs(Struct(j).Time - TSData.R_Off(OppInd-1,2)));
                
            else
                % for right side
                [~, Ind] = min(abs(Struct(j).Time(1) - TSData.R_Strike(:,2)));
                Struct(j).GC_Start = TSData.R_Strike(Ind,2);
                [~, Struct(j).GC_StartInd] = min(abs(Struct(j).Time - TSData.R_Strike(Ind,2)));
                Struct(j).GC_End = TSData.R_Strike(Ind+1,2);
                [~, Struct(j).GC_EndInd] = min(abs(Struct(j).Time - TSData.R_Strike(Ind+1,2)));
                Struct(j).GC_Off = TSData.R_Off(Ind,2);
                [~, Struct(j).GC_OffInd] = min(abs(Struct(j).Time - TSData.R_Off(Ind,2)));
                
                OppInd = find(TSData.L_Strike(:, 2) > TSData.R_Strike(Ind, 2), 1);
                Struct(j).GC_2DS = TSData.L_Strike(OppInd, 2);
                [~,Struct(j).GC_2DSInd] =  min(abs(Struct(j).Time - TSData.L_Strike(OppInd,2)));
                Struct(j).GC_1DS = TSData.L_Off(OppInd-1, 2);
                [~, Struct(j).GC_1DSInd] =  min(abs(Struct(j).Time - TSData.L_Off(OppInd-1,2)));
            end
            % create logical of TS times for each file
            Z = zeros(length(Struct(j).Time),1);                       % stance
            Z(Struct(j).GC_StartInd: Struct(j).GC_OffInd-1) = 1;
            Struct(j).LogTimes.Stance = logical(Z);
            Z = zeros(length(Struct(j).Time),1);                       % 1DS
            Z(Struct(j).GC_StartInd: Struct(j).GC_1DSInd-1) = 1;
            Struct(j).LogTimes.DS1 = logical(Z);
            Z = zeros(length(Struct(j).Time),1);                       % single support
            Z(Struct(j).GC_1DSInd: Struct(j).GC_2DSInd-1) = 1;
            Struct(j).LogTimes.SingSup = logical(Z);
            Z = zeros(length(Struct(j).Time),1);                       % 2DS
            Z(Struct(j).GC_2DSInd: Struct(j).GC_OffInd-1) = 1;
            Struct(j).LogTimes.DS2 = logical(Z);
            Z = zeros(length(Struct(j).Time),1);                       % swing
            Z(Struct(j).GC_OffInd: Struct(j).GC_EndInd) = 1;
            Struct(j).LogTimes.Swing = logical(Z);
            Z = zeros(length(Struct(j).Time),1);                       % stride
            Z(Struct(j).GC_StartInd: Struct(j).GC_EndInd) = 1;
            Struct(j).LogTimes.Stride = logical(Z);
            
            % save data structures
            Struct(j).UMB = UMB;
            Struct(j).BHAR = BHAR;
            
            if strcmp(Struct(j).Side, 'Left')
                Subjects(Subj).Trials(Trial).Met.Left = Struct(j);
            elseif strcmp(Struct(j).Side, 'Right')
                Subjects(Subj).Trials(Trial).Met.Right = Struct(j);
            end
            
            clearvars Zeros Ind OppInd
                    
        elseif strcmp(type, 'forces')
            
            Struct(j).Data = Data;
            if strcmp(Struct(j).Side, 'Left') % assing to subjects structure
                Subjects(Subj).Trials(Trial).Actuators.Left =  Struct(j).Data;
            elseif strcmp(Struct(j).Side, 'Right')
                Subjects(Subj).Trials(Trial).Actuators.Right =  Struct(j).Data;
            end
            
            

        elseif strcmp(type, 'activations')
            Struct(j).Data = Data;
            if strcmp(Struct(j).Side, 'Left') % assing to subjects structure
                Subjects(Subj).Trials(Trial).Activations.Left =  Struct(j).Data;
            elseif strcmp(Struct(j).Side, 'Right')
                Subjects(Subj).Trials(Trial).Activations.Right =  Struct(j).Data;
            end
            
        elseif strcmp(type, 'kinematics')
            Struct(j).Data = Data;
            if strcmp(Struct(j).Side, 'Left') % assing to subjects structure
                Subjects(Subj).Trials(Trial).Kinematics.Left =  Struct(j).Data;
            elseif strcmp(Struct(j).Side, 'Right')
                Subjects(Subj).Trials(Trial).Kinematics.Right =  Struct(j).Data;
            end
        end
            

        j = j + 1;
        
            

    end
end

% clearvars i j k Str Subj Trial TrialName Z



end
