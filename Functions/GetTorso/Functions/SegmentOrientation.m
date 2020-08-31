function [Euler, Rot, Quat] = SegmentOrientation(Parent, Child, Viz)

%% Segment Orientation
% will generate euler angles and rotation matricies between two segments.
% The movement of the child segment will be referenced to the parent segment.

% INPUTS
% Parent    the 3D marker locations of at least 3 markers that define a segment
%              the first marker will define the primary axis (X), and so on (Y, Z)
% Child      the 3D marker locations of at least 3 markers that define a segment
%              the first marker will define the primary axis (X), and so on (Y, Z)
% Viz         whether or not to vizualize the orientations by plotting

% OUTPUTS
% Euler
% Rot
% Quat

% Ricky Pimentel       November 2019
% Applied Biomechanics Lab      UNC - Chapel Hill

%% create orientation of parent segment
% get mean location to treat as the center
ParentCenter = mean(Parent,3);

% use first point through center to define first axis
ParentAxis1 = ParentCenter - Parent(:,:,1);
ParentAxis1uv = ParentAxis1 ./ norm(ParentAxis1);

% use line between two other points to define second axis
ParentAxis2 = Parent(:,:,2) - Parent(:,:,3);
ParentAxis2uv = ParentAxis2 ./ norm(ParentAxis2);

% cross first two axes to get third axis
ParentAxis3 = cross(ParentAxis1uv, ParentAxis2uv);
ParentAxis3uv = ParentAxis3 ./ norm(ParentAxis3);

% Create unit vector
ParentVector = ParentAxis1 + ParentAxis2 + ParentAxis3;
ParentUV = ParentVector ./ norm(ParentVector);

% cross first and third axes to set axis 2 orthogonal
ParentAxis2 = cross(ParentAxis3, ParentAxis1);
ParentAxis2uv = ParentAxis2 ./ norm(ParentAxis2);


%% define orientation of child segment
% get mean location to treat as the center
ChildCenter = mean(Child,3);

% use first point through center to define first axis
ChildAxis1 = ChildCenter - Child(:,:,1);
ChildAxis1uv = ChildAxis1 ./ norm(ChildAxis1);

% use line between two other points to define second axis
ChildAxis2 = Child(:,:,2) - Child(:,:,3);
ChildAxis2uv = ChildAxis2 ./ norm(ChildAxis2);

% cross first two axes to get third axis
ChildAxis3 = cross(ChildAxis1uv, ChildAxis2uv);
ChildAxis3uv = ChildAxis3 ./ norm(ChildAxis3);

% Create unit vector
ChildVector = ChildAxis1 + ChildAxis2 + ChildAxis3;
ChildUV = ChildVector / norm(ChildVector);

% cross first and third axes to set axis 2 orthogonal
ChildAxis2 = cross(ChildAxis3, ChildAxis1);
ChildAxis2uv = ChildAxis2 ./ norm(ChildAxis2);

clearvars Dist Line

%% Plot parent orientations to check defined axes
if strcmp(Viz, 'Yes')
    figure;
    Tri = [1 2 3 1];
    VizFactor = 100;
    pLW = 3;
    cLW = 1;
    i = 1;
    clf;
    hold on;
    
    % parent
    plot3([0 ParentAxis1uv(i, 1)*VizFactor],...
        [0 ParentAxis1uv(i, 2)*VizFactor],...
        [0 ParentAxis1uv(i, 3)*VizFactor], 'r', 'LineWidth',pLW);
    plot3([0 ParentAxis2uv(i, 1)*VizFactor],...
        [0 ParentAxis2uv(i, 2)*VizFactor],...
        [0 ParentAxis2uv(i, 3)*VizFactor], 'g', 'LineWidth',pLW);
    plot3([0 ParentAxis3uv(i, 1)*VizFactor],...
        [0 ParentAxis3uv(i, 2)*VizFactor],...
        [0 ParentAxis3uv(i, 3)*VizFactor], 'b', 'LineWidth',pLW);
      
    % child
    plot3([0 ChildAxis1uv(i, 1)*VizFactor],...
        [0 ChildAxis1uv(i, 2)*VizFactor],...
        [0 ChildAxis1uv(i, 3)*VizFactor], 'r', 'LineWidth',cLW);
    plot3([0 ChildAxis2uv(i, 1)*VizFactor],...
        [0 ChildAxis2uv(i, 2)*VizFactor],...
        [0 ChildAxis2uv(i, 3)*VizFactor], 'g', 'LineWidth',cLW);
    plot3([0 ChildAxis3uv(i, 1)*VizFactor],...
        [0 ChildAxis3uv(i, 2)*VizFactor],...
        [0 ChildAxis3uv(i, 3)*VizFactor], 'b', 'LineWidth',cLW);
 
    axis equal;
    pause(0.1);
end

%% visualize
if strcmp(Viz, 'Yes')
    VizFactor = 100;
    Tri = [1 2 3 1];
    
    figure('Position', [50 50 1200 800]);
    [m,~,~] = size(Parent); 
    for i = 1:m
        clf;
        hold on;
        % plot parent
        plot3(ParentCenter(i,1), ParentCenter(i,2), ParentCenter(i,3), 'k*'); % centroid
        plot3(reshape([Parent(i,1, Tri)], [1 4]), reshape([Parent(i,2,Tri)], [1 4]), reshape([Parent(i,3, Tri)], [1 4]), '-k'); % points
        plot3([ParentCenter(i,1)  ParentCenter(i,1)+ParentAxis1uv(i,1).*VizFactor],...
            [ParentCenter(i,2) ParentCenter(i,2)+ParentAxis1uv(i,2).*VizFactor],...
            [ParentCenter(i,3) ParentCenter(i,3)+ParentAxis1uv(i,3).*VizFactor], '-r');
        plot3([ParentCenter(i,1) ParentCenter(i,1)+ParentAxis2uv(i,1).*VizFactor],...
            [ParentCenter(i,2) ParentCenter(i,2)+ParentAxis2uv(i,2).*VizFactor],...
            [ParentCenter(i,3) ParentCenter(i,3)+ParentAxis2uv(i,3).*VizFactor], '-g');
        plot3([ParentCenter(i,1) ParentCenter(i,1)+ParentAxis3uv(i,1).*VizFactor],...
            [ParentCenter(i,2) ParentCenter(i,2)+ParentAxis3uv(i,2).*VizFactor],...
            [ParentCenter(i,3) ParentCenter(i,3)+ParentAxis3uv(i,3).*VizFactor], '-b');
        
        % plot child
        plot3(ChildCenter(i,1), ChildCenter(i,2), ChildCenter(i,3), 'k*'); % centroid
        plot3(reshape([Child(i,1, Tri)], [1 4]), reshape([Child(i,2,Tri)], [1 4]), reshape([Child(i,3, Tri)], [1 4]), '-k'); % points
        plot3([ChildCenter(i,1)  ChildCenter(i,1)+ChildAxis1uv(i,1).*VizFactor],...
            [ChildCenter(i,2) ChildCenter(i,2)+ChildAxis1uv(i,2).*VizFactor],...
            [ChildCenter(i,3) ChildCenter(i,3)+ChildAxis1uv(i,3).*VizFactor], '-r');
        plot3([ChildCenter(i,1) ChildCenter(i,1)+ChildAxis2uv(i,1).*VizFactor],...
            [ChildCenter(i,2) ChildCenter(i,2)+ChildAxis2uv(i,2).*VizFactor],...
            [ChildCenter(i,3) ChildCenter(i,3)+ChildAxis2uv(i,3).*VizFactor], '-g');
        plot3([ChildCenter(i,1) ChildCenter(i,1)+ChildAxis3uv(i,1).*VizFactor],...
            [ChildCenter(i,2) ChildCenter(i,2)+ChildAxis3uv(i,2).*VizFactor],...
            [ChildCenter(i,3) ChildCenter(i,3)+ChildAxis3uv(i,3).*VizFactor], '-b');
        
        % plot std data as ellipses at each marker position?
        %     for k = 1:4
        %         % torso stds
        %         [x,y,z] = ellipsoid(data(i).Cycles.TorsoAvg(j,Xdata(k)), data(i).Cycles.TorsoAvg(j,Ydata(k)), data(i).Cycles.TorsoAvg(j,Zdata(k)),...
        %             data(i).Cycles.TorsoStd(j,Xdata(k)), data(i).Cycles.TorsoStd(j,Ydata(k)), data(i).Cycles.TorsoStd(j,Zdata(k)));
        %         surf(x,y,z);
        %         % pelvis stds
        %         [x,y,z] = ellipsoid(data(i).Cycles.PelvisAvg(j,Xdata(k)), data(i).Cycles.PelvisAvg(j,Ydata(k)), data(i).Cycles.PelvisAvg(j,Zdata(k)),...
        %             data(i).Cycles.PelvisStd(j,Xdata(k)), data(i).Cycles.PelvisStd(j,Ydata(k)), data(i).Cycles.PelvisStd(j,Zdata(k)));
        %         surf(x,y,z);
        %     end
        
        % axis and title things
        title({'Parent and Child positions at frame ' num2str(i)});
        axis equal;
        hold off;
        pause(0.1);
        
    end
end

%% Determine Rotation matricies, Euler angles, and Quaternions between parent and child
p0 = ParentUV;
p1 = ChildUV;

% calculate cross and dot products
C = cross(p0, p1) ;
D = dot(p0, p1) ;
NP0 = norm(p0) ; % used for scaling
if ~all(C==0) % check for colinearity
    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ;
    Rot = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
else
    Rot = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
end

% R is the rotation matrix from p0 to p1, so that (except for round-off errors) ...
% R * p0'      % ... equals p1
% inv(R) * p1' % ... equals p0

Quat = qGetQ(Rot);
Euler = Quat2Eul(Quat');

end