function[UMB, BHAR, Data] = GetMuscleMetabolics(file, TSData, Path)

% obtain muscle metabolincs from .sto files

%% initial settings
% add paths
addpath(genpath('Scripts'));
addpath(genpath('CodeLibrary'));

%% Select STO file to input
if exist('file', 'var') == 0
    [file, Path] = uigetfile('.sto','Select .sto file of Muscle Metabolic Data to load');
end
addpath(genpath(Path));

% load data
Data = importdata(file);

% define filename in data structure
Data.filename = file; 

close all;

%% Filter and Parse Raw Data
[Data] = FilterAndParseData(Data, TSData);


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

%% Define variables on left and right sides
% and define muscles spanning hip, knee, and ankle joints

% umberger metabolics
[UMB.ParsedMuscles, UMB.ParsedJoints] = GetMuscles(Data.UMB_Parsed, Data.colheaders(UMB2ind));
[UMB.InterpMuscles, UMB.InterpJoints] = GetMuscles(Data.UMB_Interp, Data.colheaders(UMB2ind));

% bhargava metabolics
[BHAR.ParsedMuscles, BHAR.ParsedJoints] = GetMuscles(Data.BHAR_Parsed, Data.colheaders(BHARind));
[BHAR.InterpMuscles, BHAR.InterpJoints] = GetMuscles(Data.BHAR_Interp, Data.colheaders(BHARind));

end

