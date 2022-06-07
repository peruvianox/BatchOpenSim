function [S,pos,origin,R]=load_fpcal(infile)
%   [S,pos,origin,R]=load_fpcal(infile)
%   LOAD_FPCAL is used to load a forceplate calibration file (forcepla.cal)
%   that is in a format specified in the Motion Analysis EvaRT software
%   
%   Inputs:
%       infile - forceplate calibration file to be loaded
%                if unspecified, the 'forcepla.cal' file in the current
%                directory is loaded
%
%   Outputs:
%       S       forceplate calibration matrices, is 6x6xn, where n is the #
%               of forceplates
%       pos     position vectors (specified in m) to the top center of the forceplates in
%               the motion analysis reference frame
%       origin  origin of the transducer relative to the top center of the
%               forceplate (specified in mm)
%       R       transformation matrix to convert from the forceplate
%               reference frame to the motion capture reference frame
%
%   Created: July 11, 2008
%
%   MATLAB Version 7.1

n = nargin;
if (n==0);
    infile='forcepla_new_bertec.cal';
end
fid=fopen(infile,'r');

if (fid==-1);
    disp('Forceplate calibration file not found');
    S = [];
    pos = [];
    origin = [];
    R = [];
    return;
end

disp(['Loading file...' infile] );

i=0;
s=zeros(6,6);
r=zeros(3,3);
while feof(fid)==0
%     keyboard
    i=i+1;
    % First line is the forceplate #
    line=fgetl(fid);
    fpn=sscanf(line,'%d');
    % Next line is the dimensions of the forceplate - skip this
    line=fgetl(fid);
    % Next 6 lines are the calibration matrix
    for j=1:6
        line=fgetl(fid);
%         keyboard
        s(j,1:6) = sscanf(line,'%f %f %f %f %f %f');
    end
    S(:,:,i)=s;
    % Next line is the origin of the forceplate
    line=fgetl(fid);
    origin(:,i) = sscanf(line,'%f %f %f');
    % Next line is the top center of the forceplate
    line=fgetl(fid);
    pos(1:3,i) = sscanf(line,'%f %f %f');
    % Next 3 lines are the rotation matrix
    for j=1:3
        line=fgetl(fid);
        r(j,1:3) = sscanf(line,'%f %f %f');
    end
    R(:,:,i)=r;
end

% Return the position and origin data in m
pos = 0.01*pos;
origin = 0.01*origin;
