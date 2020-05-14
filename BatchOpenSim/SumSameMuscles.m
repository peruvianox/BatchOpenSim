function [Hip, Knee, Ankle, AllMuscles] = SumSameMuscles(Hip, Knee, Ankle)

% Add Muscles together if they are just separate lines of action
% for example: glut_med1, glut_med2, glut_med3 = glut_med


% HIP
if isfield(Hip,'UMB_Muscles') && isfield(Hip,'UMB_Columns')
    for i = 1:length(Hip)
    Gmed_Ind = contains('glut_med', Hip(i).UMBColumns); 
    Gmed = sum(Hip(i).UMB_Muscles(:,Gmed_Ind)); 
    end
end
    
end