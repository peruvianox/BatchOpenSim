% Analyze Forces from UNCABL

clear
clc
% 
% files = [{'A12_1.anc'}, {'I12_1.anc'}, {'W08_1.anc'}, {'W10_1.anc'}, {'W12_1.anc'}, {'W14_1.anc'}, {'W16_1.anc'}, {'W18_1.anc'}];
files = [{'A12_1.anc'}, {'I12_1.anc'}, {'W18_1.anc'}];
subjects = [{'Sub09'}];
% subjects = [{'Sub02'},{'Sub03'},{'Sub04'},{'Sub05'},{'Sub06'},{'Sub07'},{'Sub08'},{'Sub09'},{'Sub10'},{'Sub12'}];

FPcal_file = 'G:\Post\ABL_Pipeline\MotionAnalysis_Pipeline\2_anc_to_forces\forcepla_April2018.cal'; % Lab Computer

for l=1:length(subjects)
    subject = char(subjects(l))
    odir = strcat(['C:\Users\wihcl\Desktop\', subjects{l},'\Mocap\Generated_Files\']); % Update to your subject's Cortex Files

for t = 1:length(files)
    input_file = char(files(t))
    [mass(t)] = convertFPdata(input_file,FPcal_file,odir);
end

% submass(l) = mean(mass);
end
