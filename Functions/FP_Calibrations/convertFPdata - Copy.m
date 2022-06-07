function [mass,FPused, files, directoryf] = convertFPdata(input_file, FPcal_file, directory, zero_file)
% done=writeForcesFile(time,forcedata,directory,file) 
% writes a motion analysis style forces file which includes the 
% forces (Fx, Fy, Fz) ceneter of pressure (x, y, z) and free moment for
% each force plate
% Inputs:
% input_file: file name to be converted to forces file  
% FPcal_file: .cal file that includes force plate data  
% zero_file: zero forces file to zero forceplates 
% directory: Where the output file should be written

narg = nargin;
if narg < 1
    display('Select Files for force plate.');
    [files,directoryf]=uigetfile('*.anc','Select Files for force plate','multiselect','on');
else
    files=input_file
    directoryf=directory;
end

% Check for FPcal file otherwise load one in selected directory
if narg >= 2
    % Load FPcal_file
    [S,pos,origin,R]=load_fpcal(FPcal_file);
    % S       forceplate calibration matrices, is 6x6xn, where n is the #
    %         of forceplates
    % pos     position vectors (specified in m) to the top center of the forceplates in
    %         the motion analysis reference frame
    % origin  origin of the transducer relative to the top center of the
    %         forceplate (specified in m)
    % R       transformation matrix to convert from the forceplate
    %         reference frame to the motion capture reference frame
else
    % Loads .cal file in current directory
    [S,pos,origin,R]=load_fpcal('forcepla.cal');
end

% Check for directory input otherwise ask for one
% if narg >= 3
%    
% else
%     %Ask user for the directory
%     directory = input('Enter directory where file will be created ', 's');
% end

% Check for zero file input otherwise check if one is needed
if narg >= 4
    % Load zero file
    [samp_ratez, channel_namesz, rangez, timez, dataz, fidz, voltdataz]=load_zero(zero_file);
    mean_zero = mean(voltdataz(:,1:size(voltdataz,2)));
else
    %Ask user if they want to use zero file
%     zero = input('Do you want to use a zero file (y/n)', 's');
    zero='n';
%     gait = input('Are you analyzing running using 2 belts (y/n)', 's');
    if zero == 'y'
        % Load zero file
        [samp_ratez, channel_namesz, rangez, timez, dataz, voltdataz]=load_anc();
        mean_zero = mean(voltdataz(:,1:size(voltdataz,2)));
    end
end

% Fy threshold value to eliminate undefined COP Calulations
fythresh = 20.0;

if iscellstr(files)==0
    loop_index=1;
else 
    loop_index=size(files,2);
    k=0;
end

for w = 1:loop_index
    
    if loop_index==1
        input_file=files(1:end);%JRF removed '-4'
    else 
      k=k+1;
      input_file=cell2mat(files(k));
      input_file=input_file(1:end);%JRF removed '-4'
    end
    
    % Load .anc file from input_file name

    [data,time,samp_rate, channel_names, range, voltdata]=load_anc(input_file,directory);%;input_file, directoryf);
    if zero == 'y'
    % Zero the voltdata if user desires
    voltdata=voltdata-ones(size(voltdata,1),1)*mean_zero;  
    end
   
    %filtering the data
    T=time(2)-time(1);
    cutoff=100;
    [b a]=butter(3,cutoff*T*2,'low');
%     [b1 a1]=butter(3,80*T*2,'low');
    voltdata=filter(b,a,voltdata);
%     voltdata_filt1=filter(b1,a1,voltdata);
    %frequency analysis
   
%     NFFT=2^nextpow2(size(voltdata,1));
%     Y=fft(voltdata(:,3),NFFT)/(size(voltdata,1));
%     f=(1/T)/2*linspace(0,1,NFFT/2+1);
%     figure(2); hold on; plot(f,2*abs(Y(1:NFFT/2+1)),'g'); axis([0 200 0 0.1]);

     %figure(1); hold on; plot(time,voltdata(:,2),'k'); plot(time,voltdata_filt(:,2),'c:'); %plot(time,voltdata_filt1(:,1),'r:'); %plot(time,voltdata_filt1(:,3),'c');
%      %figure(2); hold on; plot(time,voltdata(:,2),'g'); plot(time,voltdata_filt(:,2),'k'); %plot(time,voltdata_filt1(:,2),'c');
%       plot(time,voltdata(:,1),'g'); plot(time,voltdata_filt(:,1),'k'); %plot(time,voltdata_filt1(:,1),'c');
% Find which forceplates are being used 
    mean_data=mean(data(:,1:size(data,2)));
    FPused=[];

    
    %UNC only has two plates
F1=sum(strcmp(channel_names,'F1Z'));
F2=sum(strcmp(channel_names,'F2Z'));
F3=sum(strcmp(channel_names,'F3Z'));
F4=sum(strcmp(channel_names,'F4Z'));
% F5=sum(strcmp(channel_names,'F5Z'));

% %   if F1==1 && F2==1 && F3==1 && F4==1 && F5==1
% %     FPused=[1 2 3 4 5];
 %elseif F1==1 && F2==1 && F3==1  

% % if F1==1 && F2==1 && F3==1 
% %     FPused=[1 2 3];
%   elseif F4==1 
%     FPused=4;
%   elseif F5==1
%     FPused=5;
%   end
% 
  if F1==1 && F2==1 && F3==1 && F4==1
    FPused=[1 2 3 4];
    
  elseif  F1==1 && F2==1
      FPused=[1 2 3 4];
    
   end



%   %JRF Added for CNNINM study (no F3)
%   if F1==1 && F2==1 && F4==1 && F5==1
%       FPused=[1 2 4 5];
%   end
  
 
      
% for j=3:6:size(data,2)
%     if mean_data(j)>=8
%         s = channel_names(j);
%         for t=1:1:6
%             r=['F' int2str(t) 'Z'];
%             w=strmatch(r, s);
%             if w>=1
%                 FPused=[FPused t];
%             end
%         end
%     end
% end


%Find active Forceplates and create new voltdata for those forceplates
% FP_Loc=zeros(length(FPused),2);VDataNew=zeros(length(FPused),2);ChannelsNew=zeros(length(FPused),2);
FP_Loc=[];VDataNew=[];ChannelsNew=[];
  for i=1:1:length(FPused)
    FPnumber=FPused(i);
    FPnumber=int2str(FPnumber);
    FPname=['F' FPnumber 'X'];
    FP_location=strmatch(FPname,channel_names(:),'exact');
    FP_Loc=[FP_Loc FP_location];
    VDataNew=[VDataNew voltdata(:,FP_location:(FP_location+5))];
    ChannelsNew=[ChannelsNew channel_names(FP_location:(FP_location+5))];
  end

%Convert Voltdata to Force and Moment Data
Output=[];
ChanOut=[];
% X=[];
% Y=[];
% Z=[];
[g,h]=size(VDataNew);
Xcur=zeros(g,1);Ycur=zeros(g,1);
  for j=1:1:length(FPused)
    FPnumber=FPused(j);
    a=FPnumber;
    %For forceplates 1 through 3
%     if FPnumber <= 3 
%     G=4000; %FP Amplifier Gain
%     Vo=10;  %FP Excitation Voltage
%     CF=10^6; %conversion factor from uV to V
%     else
    %undefined for treadmill forceplates
    G=1;
    Vo=1;
    CF=1;
%     end
    FPname=['F' int2str(a) 'X'];
    FP_location=strmatch(FPname,ChannelsNew(:),'exact');    %Finds location of data based on channel names
    Trans=[R(:,:,a) zeros(size(R(:,:,a))); zeros(size(R(:,:,a))) R(:,:,a)]; %Transformation matrix from FPcal file
    ForMom=-VDataNew(:,FP_location:(FP_location+5))*S(:,:,a)'*(CF)*Trans/(Vo*G); % Negative to get force on person
    toffset=cross(ones(size(ForMom,1),1)*(-R(:,:,a)*origin(1:3,a))',ForMom(:,1:3)); %Calculate Torque due to transducer offset
    % Allow for zero value of z force by setting threshold and setting COP
    % to 0 if z force is below threshold
    
    for jk=1:1:g
      if  ForMom(jk,3)>=fythresh %|| ForMom(jk,3)<=-fythresh
    x=((-ForMom(jk,5)-toffset(jk,2))./ForMom(jk,3)+pos(1,a))*1000; %x axis Center of pressure 
    y=((ForMom(jk,4)+toffset(jk,1))./ForMom(jk,3)+pos(2,a))*1000; %y axis Center of Pressure

    Xcur(jk)=x; %[Xcur; x];
    Ycur(jk)=y; %[Ycur; y];
      else
    x=0;
    y=0;
    Xcur(jk)=x;%[Xcur; x];
    Ycur(jk)=y;%[Ycur; y];
      end
    end
    %z axis Center of pressure Origin (-) to transform to lab reference
    %frame
    Zcur=(zeros(size(Xcur))-origin(3,a)+pos(3,a))*1000;
    %Zcur=zeros(size(Xcur)); %z axis Center of Pressure
    %Calculation of the Torque Mz
    torque=cross((ones(size(Xcur,1),1)*pos(1:3,a)'-[Xcur Ycur Zcur]*.001),ForMom(:,1:3));
    tz=(ForMom(:,6)+torque(:,3)+toffset(:,3))*1000; %Mz output in Nmm
    %Gather all data and add to previous data for output
    Output=[Output ForMom(:,1:3) Xcur Ycur Zcur tz];
    %Add the Channel Names to a channel output
    ChanOut=[ChanOut 'FX' int2str(a) ' FY' int2str(a) ' FZ' int2str(a) ' X' int2str(a)	' Y' int2str(a) ' Z' int2str(a) ' MZ' int2str(a)];
%     end
    
%     %For treadmill forceplates
%     if (FPnumber >= 4)
%     Xcur=[];
%     Ycur=[];
%     a=FPnumber;
%     FPname=['F' int2str(a) 'X'];
%     FP_location=strmatch(FPname,ChannelsNew(:),'exact'); %Finds location of data
%     Trans=[R(:,:,a) zeros(size(R(:,:,a))); zeros(size(R(:,:,a))) R(:,:,a)]; %Builds transformation matrix
%     ForMom=-VDataNew(:,FP_location:(FP_location+5))*S(:,:,a)'*Trans; %Calulates forces and moments from the new data
%     % Allow for zero value of z force by setting threshold and setting COP
%     % to 0 if z force is below threshold
%     for jk=1:1:g
%     if  ForMom(jk,3)>=fythresh %|| ForMom(jk,3)<=-fythresh
%     x=((-ForMom(jk,5)+ForMom(jk,1)*.015)./ForMom(jk,3)+origin(1,a)+pos(1,a))*1000; %x axis Center of pressure
%     y=((ForMom(jk,4)+ForMom(jk,2)*.015)./ForMom(jk,3)+origin(2,a)+pos(2,a))*1000; %y axis Center of pressure
%     Xcur=[Xcur; x];
%     Ycur=[Ycur; y];
%     else
%     x=0;
%     y=0;
%     Xcur=[Xcur; x];
%     Ycur=[Ycur; y];
%     end
%     end
%     Zcur=(zeros(size(Xcur))-origin(3,a)+pos(3,a))*1000; %z axis Center of pressure Origin (-) to transform to lab reference frame
%     %Calculation of the Torque Mz
%     torque=ForMom(:,6)+Xcur.*ForMom(:,2)-Ycur.*ForMom(:,1);
%     %Gather all data and add to previous data for output
%     Output=[Output ForMom(:,1:3) Xcur Ycur Zcur torque];
%     %Add the Channel Names to a channel output
%     ChanOut=[ChanOut '  FX' int2str(a) '    FY' int2str(a) '    FZ' int2str(a) '    X' int2str(a)	'   Y' int2str(a) '     Z' int2str(a) '     MZ' int2str(a)];
%     end


  end


forces=zeros(size(Output));
% if length(FPused)==1 && FPused(1)==4
% fc_temp=find(Output(:,3)>=50);
% forces=zeros(size(Output));
% forces(fc_temp,1:7)=Output(fc_temp,1:7);
% %put 4 and 5 together for running
% %if gait=='y'
% elseif length(FPused)==2
% clear fc_comb;
% fc_temp1=find(Output(:,3)>=50);
% fc_temp2=find(Output(:,10)>=50);
% fc_comb=union(fc_temp1,fc_temp2);
% 
% forces(fc_temp1,4:7)=Output(fc_temp1,4:7);
% forces(fc_comb,1:3)=Output(fc_comb,1:3)+Output(fc_comb,8:10);
% forces(fc_temp2,4:7)=Output(fc_temp2,11:14);
% else
% %for walking use temp1 and temp2   

fc_temp1=find(Output(:,3)>=50);
fc_temp2=find(Output(:,10)>=50);
fc_temp3=find(Output(:,17)>=50);
fc_temp4=find(Output(:,24)>=50);

forces(fc_temp1,1:7)=Output(fc_temp1,1:7);
forces(fc_temp2,8:14)=Output(fc_temp2,8:14);
forces(fc_temp3,15:21)=Output(fc_temp3,15:21);
forces(fc_temp4,22:28)=Output(fc_temp4,22:28);

%JRF edit
% if length(FPused)==3
%     fc_temp3=find(Output(:,17)>=50);
%     forces(fc_temp3,15:21)=Output(fc_temp3,15:21);
% end




input_file = strrep(input_file, '.anc', '');

%Write Forces File
T=time(2)-time(1);
f=1/T;
[npts,nf]=size(Output);
nfp=nf/7;

 if isempty(directoryf)==1
    fid = fopen([input_file,'.forces'],'w');
 else
    fid = fopen([directoryf,input_file,'.forces'],'w');
 end

% Write the header
fprintf(fid,['[Force Data]']);
fprintf(fid,'\n');
fprintf(fid,['NumberOfForcePlates=',num2str(nfp)]);
fprintf(fid,'\n');
fprintf(fid,['SampleRate=%12.7f'],f);
fprintf(fid,'\n');
fprintf(fid,['NumberOfSamples=',num2str(npts)]);
fprintf(fid,'\n');
fprintf(fid,['#Sample ']);
fprintf(fid,'     %15s',ChanOut);
fprintf(fid,'\n');

% Write the data
 for i=1:npts
   fprintf(fid,[num2str(i)]);
   fprintf(fid,'\t%3.2f',forces(i,:));
   fprintf(fid,'\n');
 end
 
 
 mass=mean(forces(:,3)+forces(:,10))/9.81;

 
 
disp(['Wrote ',num2str(npts),' frames of force data to ',[input_file,'.forces']]);
fclose(fid);
done=1;
end
end









