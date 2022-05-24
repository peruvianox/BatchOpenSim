function [scaleError] = ComputeScaleError(Subjects)

for subj = 1:length(Subjects)
    trial = contains({Subjects(subj).Trials(:).type}, 'static');
    scaleError(subj).Subj = Subjects(subj).name;
    scaleError(subj).TotalSqError = Subjects(subj).Trials(trial).ScaleErr.TotalSqErr;
    scaleError(subj).RMSError = Subjects(subj).Trials(trial).ScaleErr.RMSErr;
    scaleError(subj).MaxError = Subjects(subj).Trials(trial).ScaleErr.MaxErr;
    scaleError(subj).MaxMkr =  Subjects(subj).Trials(trial).ScaleErr.MaxMkr;
end
        
%% Plot Results
SubjColors = colormap(jet(length(Subjects))); 
close all; 
MkrSz = 20;
H = figure('Position',[100 100 800 600]); 

for subj = 1:length(Subjects)
    subplot(131); hold on;
    plot(1, scaleError(subj).TotalSqError,'.', 'Color', SubjColors(subj,:), 'MarkerSize',MkrSz);
    title('Total Square Error');
    ylabel('Error (m)'); 
    subplot(132); hold on;
    plot(1, scaleError(subj).RMSError,'.', 'Color', SubjColors(subj,:), 'MarkerSize',MkrSz);
    title('RMS Error');
    subplot(133); hold on;
    plot(1, scaleError(subj).MaxError,'.', 'Color', SubjColors(subj,:), 'MarkerSize',MkrSz);
    title('Max Error');
    text(1.2,scaleError(subj).MaxError, scaleError(subj).MaxMkr,  'Color', SubjColors(subj,:)); 
    text(0.2,scaleError(subj).MaxError, Subjects(subj).name,  'Color', SubjColors(subj,:)); 
end
saveas(H, 'ScaleError.png'); 


end
