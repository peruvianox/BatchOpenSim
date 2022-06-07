function[] = findHJC(staticFile, Dir)

%       Using a calibration full hip range-of-motion trial, findHJC.m locates the
%       hip joint center between a pelvis and thigh and writes the hip joint
%       center location to one or many selected *.trc files.  

% OTHER M.Files REQUIRED:
%       soder, load_A, calcHJC
%
% INPUTS:
%       Inputs are selected from the prompts displayed while running finHJC.m.
%
%     'Select Static Calibration Trial':
%         Pick a static trial to set a starting reference location for
%         finding marker locations throughout calibration trials and for a
%         reference to determine HJC locations when writting to files
%
%     'Select Right Leg HJC Calibration Trial':
%         Find and open the calibration trial for finding the HJC relative
%         to the anatomical pelvis frame.
%
%     'Select Left Leg HJC Calibration Trial':
%         Find and open the calibration trial for finding the HJC relative
%         to the anatomical pelvis frame.
%
%     'Select Pelvis Markers for HJC Locating':
%         Of the markers from the selected calibration trial, select the
%         pelvis markers you want to use in HJC calculation. Multiple
%         markers can be selected.
%
%     'Select Thigh Markers for HJC Locating':
%         Of the markers from the selected calibration trial, select the
%         right thigh markers you want to use in HJC calculation.  Multiple
%         markers may be selected.
%
%     'Selct files to Write HJC Locations':
%         Select the files that you would like to write the HJC locations to.
%         A new folder 'TRC Files with Hip Joint Centers Added' will be
%         created in the same directory as the selected files.  Multiple
%         files may be selected
%
%
% OUTPUTS:
%       The selected files to write the HJC locations to are rewritten in the
%       new folder with the HJC's added.  Also, a text file with the HJC's
%       relative to the pelvic mid-ASIS frame, the mean difference calculated
%       between HJC locations relative to the pelvis and the thigh, and the
%       standard deviation of the means; is written as /HJCstatistics.txt
%
% AUTHORS: Joseph Farron, NMBL, University of Wisconsin-Madison
% DATE: November 2, 2005
% UPDATES:
% Amy Silder, June 5, 2006
% Ricky Pimentel, September 2021

%% Set Defaults
dbstop if error;

% Default Marker Names for pelvis and thighs
PelvMarks={'L.ASIS';'R.ASIS';'L.PSIS';'R.PSIS';'S2'}; % pelvis
RTMarks={'R.TH1';'R.TH2';'R.TH3';'R.Knee'}; % right thigh
LTMarks={'L.TH1';'L.TH2';'L.TH3';'L.TH4';'L.Knee'}; % left thigh

if length(RTMarks)<3
    error('not enough right thigh markers selected for HJC analysis'); 
end
if length(LTMarks)<3
    error('not enough left thigh markers selected for HJC analysis'); 
end

addpath(genpath('C:\Users\richa\Documents\Packages\MoCapTools'))

%% FIRST STEP
% Use a static trial to average all marker locations.  These will be used to:
%   1.  Find a starting reference from which movements can be related
%   2.  Note marker locations with reference to each other, so that a
%       reference frame can be determined in HJC-written files in which
%       markers normally used for reference are missing.

if exist('staticFile', 'var') == 0
    %User selects bilateral static trial if not specified
    [staticFile,Dir]=uigetfile('*.trc','Select Static Calibration Trial');
end

if exist('Dir', 'var') == 0
    Dir = fileparts(which(staticFile));
end

addpath(genpath(Dir)); 
cd(Dir);

% create static marker reference file as matlab table
trc = Osim.readTRC(staticFile);
Markers = trc.Properties.VariableNames;
trcData = table2array(trc); 
meanTRC = mean(trcData); % average static marker data

% AllMarkers = [PelvMarks; RTMarks; LTMarks]'; % extract pelvis and thigh marker data only
% Inds = contains(Markers, AllMarkers);

PelvRef = meanTRC(contains(Markers, PelvMarks)); 
RTRef = meanTRC(contains(Markers, RTMarks)); 
LTRef = meanTRC(contains(Markers, LTMarks)); 

% create table of static reference data and save
% StaticRef = array2table(meanTRC(:,Inds)); 
% StaticRef.Properties.VariableNames = string(Markers(Inds));


%% SECOND STEP
% Find the hip joint center location using a least squares method to find where
% pelvis and thigh share a common point.  It finds the point relative to thigh and pelvic frames and averages the distance between the two.  It
% outputs the averaged HJC in the mid-asis pelvic frame, along with the average distance between the two HJC's calculated and the standard
% deviation of all locations.

D = dir(Dir(1:end-1));
Dt = struct2table(D);
Dc = table2cell(Dt);
Rind = strcmp('RHJC_1.trc', Dc(:,1));
RhjcFile = Dc{Rind};
Lind = strcmp('LHJC_1.trc', Dc(:,1));
LhjcFile = Dc{Lind};

% calculate HJCs
[R_HJC, R_HJC_avg, R_HJC_std] = Osim.calcHJC(RhjcFile, PelvMarks, RTMarks, PelvRef, RTRef);
[L_HJC, L_HJC_avg, L_HJC_std] = Osim.calcHJC(LhjcFile, PelvMarks, LTMarks, PelvRef, LTRef);


%% THIRD STEP
%Select files to write HJC's to, and calculate the HJC's
%
%User is prompted to select files to write HJC's.  Then, each file is
%analyzed, and a pelvic reference frame is determined based on available
%markers.  The HJC is then transformed into this reference frame, and from
%the marker data in the file, the HJC is calculated in the global frame and
%written into the file.

%Select files to which the HJC locations should be added (multiple must be
%selected).  The same directory will be used to create a file for new files
%with HJC locations.
close all;
disp('Select Files to Add Hip Joint Center Locations')
[files,directory] = uigetfile('*.trc','Select Files to Add HJC Locations','multiselect','on');

% Locate pelvic markers in each data file, find the pelvic center, locate
% the HJC's and marker data and HJC data into a new file
for j=1:length(files)
    
    trc = Osim.readTRC(files{j});
    Markers = trc.Properties.VariableNames;
    trcData = table2array(trc); 
    meanTRC = mean(trcData); % average static marker data
    sampFreq = 1 / (trc.Time(2) - trc.Time(1));
    
    % Find the center of the markers and define the coordinate
    % system of the marker set as identical to the global frame.
    PelvRef = meanTRC(contains(Markers, PelvMarks)); 
    D_ref = [mean(reshape(PelvRef, 3, length(PelvMarks))')]';
    
    % For static data, locate the rotation and location of the
    % coordinate system that the known HJC is in.
    rasis = trcData(1, contains(Markers, 'R.ASIS'))';
    lasis = trcData(1, contains(Markers, 'L.ASIS'))';
    sacral = trcData(1, contains(Markers, 'S2'))';
    
    % create pelvis coodinate system
    midasis = (lasis+rasis)/2;
    y = lasis-rasis;
    y = y/sqrt(sum(y.*y));
    z = cross((sacral-lasis),(rasis-lasis));
    z = z/sqrt(sum(z.*z));
    X = cross(y,z);
    R = [X y z];
    
    %Find the Transformation Matrix from the HJC system to the marker-set system.
    D = midasis - D_ref;
    T = [R D;0 0 0 1];
    % Find the HJC in the m-s system.
    rhjcms = T * R_HJC;
    lhjcms = T * L_HJC;
    
    % Locate the markers in the marker set throughout the trial being written to
    % then use soder to find the transformations to each time set of markers.
    center = zeros(length(trcData), 3);
    r_hjc = zeros(length(trcData), 3);
    l_hjc = zeros(length(trcData), 3);
    time = zeros(length(trcData), 1);

    % get all pelvis markers
    marks = trcData(:, contains(Markers, PelvMarks));
    
    for i = 1:length(trcData)

        % get pelvis orientation for each frame
        [T_pelv,~] = Osim.soder([PelvRef; marks(i,:)]);
        
        % From these T, find the HJC in the global frame.
        center(i,:) = mean(reshape(marks(i,:),3,length(marks(1,:))/3)');
        
        rr = [center(i,:)]'+[T_pelv(1:3,1:3)*rhjcms(1:3,1)];
        r_hjc(i,(1:3)) = [rr(1:3)]';

        ll = [center(i,:)]'+[T_pelv(1:3,1:3)*lhjcms(1:3,1)];
        l_hjc(i,(1:3)) = [ll(1:3)]';

        time(i,1) = i/sampFreq - 1/sampFreq;
    end

    Markers{1} = 'Header';
    trcTable = array2table([trcData r_hjc l_hjc]);
    trcTable.Properties.VariableNames = cellstr([Markers, ...
        {'R.HJC_x'}, {'R.HJC_y'}, {'R.HJC_z'},...
        {'L.HJC_x'}, {'L.HJC_y'}, {'L.HJC_z'}]);
    
    % write HJC data to new TRC file
    HJCfile = [files{j}(1:end-4) '_hjc.trc'];
    Osim.writeTRC(trcTable, 'FilePath', HJCfile);
    
    % move non-hjc files to new folder
    mkdir('HJC')
    movefile(HJCfile, strcat('HJC/', HJCfile));

end


for j = 1:length(files)
     Osim.plotMarkers(strcat('HJC/', HJCfile))
end


%% Write a text file that gives statistical information about HJC location
fid=fopen('HJCstatistics.txt','w');
fprintf(fid,['Statistics of Calibration of HJC (locations relative to pelvis)      \n']);
fprintf(fid,['                                                                     \n']);
fprintf(fid,['Right Leg HJC                      Left Leg HJC                      \n']);
fprintf(fid,['X          Y           Z           X          Y          Z           \n']);
fprintf(fid,'%-f',R_HJC(1)); fprintf(fid,['   ']); 
fprintf(fid,'%-f',R_HJC(2)); fprintf(fid,['   ']); 
fprintf(fid,'%-f',R_HJC(3)); fprintf(fid,['   ']);
fprintf(fid,'%-f',L_HJC(1)); fprintf(fid,['   ']); 
fprintf(fid,'%-f',L_HJC(2)); fprintf(fid,['   ']); 
fprintf(fid,'%-f',L_HJC(3)); fprintf(fid,['   \n']);
fprintf(fid,['                                                                     \n']);
fprintf(fid,['Mean         Std. Dev.             Mean         Std. Dev.            \n']);
fprintf(fid,'%-f',R_Ave);fprintf(fid,['     ']);
fprintf(fid,'%-f',R_Std_Dev);fprintf(fid,['              ']);
fprintf(fid,'%-f',L_Ave);fprintf(fid,['     ']);
fprintf(fid,'%-f\n',L_Std_Dev);
fclose(fid);

disp(' '); 
disp('All files written!'); 

end


