function [AvgLegLength, Left, Right] = GetLegLength(MarkerData, Markers, LTimes, RTimes)

% this function will define the average leg length from marker trajectories
% for scaling factors. This leg length is not an anatomically accurate in
% terms of hip joint center to knee joint center to ankle joint center,
% just a quick check from ASIS to Knee to Ankle markers. 

if exist('LTimes') == 0
    LTimes = 1:length(MarkerData(1).Trajectories); 
end
if exist('RTimes') == 0
    RTimes = 1:length(MarkerData(1).Trajectories); 
end

%% Left Leg Length
% define markers
LASIS = mean(MarkerData(strcmp(Markers, 'L.ASIS')).Trajectories(LTimes,:), 1); 
LKnee = mean(MarkerData(strcmp(Markers, 'L.Knee')).Trajectories(LTimes,:),1); 
LAnkle = mean(MarkerData(strcmp(Markers, 'L.Ankle')).Trajectories(LTimes,:),1); 
% get lengths
Left.Thi = norm(LASIS - LKnee); 
Left.Shank = norm(LKnee - LAnkle); 
Left.Leg = Left.Thi + Left.Shank; 

%% Right Leg Length
% define markers
RASIS = mean(MarkerData(strcmp(Markers, 'R.ASIS')).Trajectories(RTimes,:), 1); 
RKnee = mean(MarkerData(strcmp(Markers, 'R.Knee')).Trajectories(RTimes,:),1); 
RAnkle = mean(MarkerData(strcmp(Markers, 'R.Ankle')).Trajectories(RTimes,:),1); 
% get lengths
Right.Thi = norm(RASIS - RKnee); 
Right.Shank = norm(RKnee - RAnkle); 
Right.Leg = Right.Thi + Right.Shank; 

%% Average Sides
AvgLegLength = mean([Left.Leg Right.Leg]);

end

