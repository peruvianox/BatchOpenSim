function [anc, time, samp_rate, channel_names, range, voltdata]=load_anc(infile,inpath)

%% LOAD_ANC is used to open analog data files from Motion Analysis Realtime
%   output (*.anc).  The number of channels, samp rate and range are determined
%   from the file header.
%
% Last modified:    Amy Silder; October 7, 2005
%                   Darryl Thelen, July 11, 2008 to also output data in volts
%                   Ricky Pimentel, October 2019 to organize and improve comments

%% Set input parameters
narg = nargin;
% keyboard
if (narg==0)
    [infile, inpath]=uigetfile('*.anc','Select input file');
    if infile==0
        disp('No file selected');
        return;
    end
    fid = fopen([inpath infile],'r');
elseif (narg>0)
    infile=char(infile);
    if (infile((length(infile)-3):length(infile))~='.anc') 
        infile = [infile(1:length(infile)) '.anc'];
    end
    if (narg==2)
        fid = fopen([inpath infile],'r');
    else 
        fid=fopen(infile,'r');
    end
end
file = infile(1:length(infile));

if (fid==-1) % error message if file not found
    disp('File not found');
    anc=[];
    return
end

%% disregard header info
for h=1:2
    crap=fgetl(fid);
end
crap=fscanf(fid,'%s',5);
duration=fscanf(fid,'%f',1);
crap=fscanf(fid,'%s',1);
num_channels=fscanf(fid,'%f',1);
for h=1:5
    crap=fgetl(fid);
end
crap=fscanf(fid,'%s',1);

%% Define channels
channel_name='';
channel_name=fgetl(fid);
j=1;
jl=length(channel_name);

for i=1:(num_channels)
    name=sscanf(channel_name(j:jl),'%s',1);
    ii=findstr(channel_name(j:jl),name);
    j=j+ii(1)+length(name);
    channel_names(i,1)=cellstr(name);
end

%% Load channel data
crap=fscanf(fid,'%s',1);
for i=1:num_channels
    samp_rate(i)=fscanf(fid,'%f',1);
end
crap=fscanf(fid,'%s',1);
for i=1:num_channels
    range(i)=fscanf(fid,'%f',1);
end

npts=round(duration*samp_rate(1)+1);
time=zeros(npts,1);
data=zeros(npts,num_channels);
i=1;
linetemp=(fscanf(fid,'%f',num_channels+1))';

while ((feof(fid)==0) && (length(linetemp)>0))
    time(i,1)=linetemp(1);
    data(i,1:num_channels)=linetemp(2:(num_channels+1));
    i=i+1;
    linetemp=(fscanf(fid,'%f',num_channels+1))';
end

%% Convert to volts
% voltdata=(ones(size(data,1),1)*range).*(data/2048)*0.001;
voltdata=(ones(size(data,1),1)*range).*(data/32768)*0.001;%JRF 16 bit

%% Package data for export
if (nargout>1)
    % return trc as the matrix of kinematic data
    anc=data;
else
    % return all information in a single structure
    anc.data=voltdata;
    anc.addata=data;
    anc.time=time;
    anc.freq=samp_rate;
    anc.nframes=size(data,1);
    anc.nchan=num_channels;
    anc.channel_names=channel_names;
    anc.range=range;
    anc.file=file;
    anc.inpath=inpath;
end

end

