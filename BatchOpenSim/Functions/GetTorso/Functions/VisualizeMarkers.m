function[] = VisualizeMarkers(Data, Markers2Plot, Video)

% display marker locations in 3D to visualize positions for quality
% assurance

%% extract markers of interest
if isfield('Data','MarkerData')
    MarkerData2Plot = Data.MarkerData; 
    for i = 1:length(MarkerData2Plot)
        if isempty(MarkerData2Plot(i).Trajectories) == 0
            Markers2Plot{i} = MarkerData2Plot(i).Name;
        end
    end
else
     MarkerData2Plot = GetMarkerTrajectories(Data.data, Data.colheaders, Markers2Plot);
end

%% plot first frame
MkrSz = 25;
i = 1;
figure('Position', [50 50 1200 800]);
NumMarkers = length(MarkerData2Plot);
hold on;
ToDel = zeros(NumMarkers,1); 
for j = 1:NumMarkers
    if isempty(MarkerData2Plot(j).Trajectories) == 0
        plot3(MarkerData2Plot(j).Trajectories(i,1), MarkerData2Plot(j).Trajectories(i,2), MarkerData2Plot(j).Trajectories(i,3), '.', 'MarkerSize',MkrSz);
    else
        ToDel(j) = 1;
    end
end
MarkerData2Plot(ToDel==1) = []; 
legend({MarkerData2Plot.Name});
NumMarkers = NumMarkers - sum(ToDel); 
% axis and title things
title({'Parent and Child positions at frame ' num2str(i)});
axis equal;
       
%% plot more frames?
if nargin < 3
    Video = questdlg('Would you like to see the rest of the frames as a video?');
end
if strcmp(Video, 'Yes')
    [m,~,~] = size(MarkerData2Plot(1).Trajectories);
    for i = 2:m
        clf; hold on;
        for j = 1:NumMarkers
            plot3(MarkerData2Plot(j).Trajectories(i,1), MarkerData2Plot(j).Trajectories(i,2), MarkerData2Plot(j).Trajectories(i,3), '.', 'MarkerSize',MkrSz);
        end
        legend({MarkerData2Plot.Name});
        % axis and title things
        title({'Parent and Child positions at frame ' num2str(i)});
        axis equal;
        hold off;
        pause(0.025);
    end
end
        
end


