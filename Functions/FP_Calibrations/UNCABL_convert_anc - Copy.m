% Analyze Forces from UNCABL

clear
clc
% 
files = [{'Trimmed_MUS_04_1.anc'},{'Trimmed_MUS_00_1.anc'},...
{'Trimmed_MUS_02_1.anc'}, {'Trimmed_MUS_20_1.anc'}, {'Trimmed_MUS_40_1.anc'}];

subjects = [{'Sub001'}];

FPcal_file = 'D:\ABL_Pipeline\MotionAnalysis_Pipeline\2_anc_to_forces\forcepla.cal'; % Lab Computer

for l=1:length(subjects)
    subject = char(subjects(l))
    odir = strcat(['E:\xxmisc\JumpData\', subjects{l},'\Generated_Files\']); % Update to your subject's Cortex Files
 
for t = 1:length(files)
    input_file = char(files(t))
    [mass(t)] = convertFPdata(input_file,FPcal_file,odir);

end
% submass(l) = mean(mass);
end
