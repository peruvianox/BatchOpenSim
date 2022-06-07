function [samp_rate, channel_names, range, time, data, fid, voltdata]=load_zero(infile, inpath)
%LOAD_ANC is used to open analog data files from Motion Analysis Realtime
%   output (*.anc).  The number of channels, samp rate and range are determined
%   from the file header.
%
% Last modified:    Amy Silder; October 7, 2005


file=[infile(1:length(infile)) '.anc'];
%inpath=inpath;
fid=fopen(infile,'r');
fileroot = file(1:length(file)-4);
disp(['Opening file...' file]);

%disregard header info
for h=1:2;
    crap=fgetl(fid);
end
crap=fscanf(fid,'%s',5);
duration=fscanf(fid,'%f',1);
crap=fscanf(fid,'%s',1);
num_channels=fscanf(fid,'%f',1);
for h=1:5;
    crap=fgetl(fid);
end
crap=fscanf(fid,'%s',1);

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

crap=fscanf(fid,'%s',1);
for i=1:num_channels;
    samp_rate(i)=fscanf(fid,'%f',1);
end
crap=fscanf(fid,'%s',1);
for i=1:num_channels;
    range(i)=fscanf(fid,'%f',1);
end
time=zeros(round(duration*samp_rate(1)),1);
data=zeros(round(duration*samp_rate(1)),num_channels);
i=1;
while feof(fid)==0;
    line_temp=fscanf(fid,'%f',num_channels+1);
   for j=1:length(line_temp);
   	data(i,j)=line_temp(j);
   end
   i=i+1;
end;
time=data(:,1);
data=data(:,2:size(data,2));
voltdata=(ones(size(data,1),1)*range).*(data/2048)*0.001;
