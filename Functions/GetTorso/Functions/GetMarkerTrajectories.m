function[MarkerData] = GetMarkerTrajectories(data, colheaders, Markers)
%% GetMarkerTrajectories

% INPUTS
% data              data file loaded
% colheaders    name of the column headers for reference
% Markers        strings of marker names to load, input 'All' to load all

% OUTPUTS
% MarkerData     structure containing marker names, XYZ coorinate
%                       trajectories, and column of storage in raw data array

% Ricky Pimentel       November 2019
% Applied Biomechanics Lab      UNC - Chapel Hill

%% Define markers
if iscell(Markers)
    NumMarkers = length(Markers);
else
    if strcmp('All', Markers)
        Markers = colheaders;
        [m,~] = size(Markers);
        if m > 1
            Markers(2:m,:) = [];
        end
        for i = 1:length(Markers)
            if strcmp(Markers{i},' ')
                ToDel(i) = 1;
            elseif isempty(Markers{i})
                ToDel(i) = 1;
            else
                ToDel(i) = 0;
            end
        end
        Markers(logical(ToDel)) = [];
        NumMarkers = length(Markers);
    else
        NumMarkers = 1;
    end
end

%% loop through data array, extract data, and save in structure
for j = 1:NumMarkers
    for i = 1:length(colheaders)
        if strcmp(colheaders{1,i},Markers{j})
            MarkerData(j).Name = Markers{j};
            MarkerData(j).Trajectories = [data(:, i:i+2)];
            MarkerData(j).Col = i;
        end
    end
end

end
