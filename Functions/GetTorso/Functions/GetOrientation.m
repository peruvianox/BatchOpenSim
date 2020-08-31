function[UnitVector, Vector, Coords] = GetOrientation(Matrix)

% INPUT
% Matrix   a nx3x3 array of marker positions

% OUTPUT 
% UnitVector    a unit vector of the orientation of the 3 markers
% Vector         the vector orientation
% Coords        axis coordinates for each of the orientation primary axes

%%
[m,~] = size(Matrix);
for i = 1:m
    % get mean location to treat as the center
    Center = mean(Matrix(i,:,:),3);
    
    % use first point through center to define first axis
    Axis1 = Center - Matrix(i,:,1);
    Axis1uv = Axis1 ./ norm(Axis1);
    
    % use line between two other points to define second axis
    Axis2 = Matrix(i,:,2) - Matrix(i,:,3);
    Axis2uv = Axis2 ./ norm(Axis2);
    
    % cross first two axes to get third axis
    Axis3 = cross(Axis1, Axis2);
    Axis3uv = Axis3 ./ norm(Axis3);
    
    % Create unit vector
    Vector(i,:,:) = Axis1 + Axis2 + Axis3;
    
    UnitVector(i,:,:) = Vector(i,:) ./ norm(Vector(i,:));
    
    Coords(i,:,1) = Axis1uv;
    Coords(i,:,2) = Axis2uv;
    Coords(i,:,3) = Axis3uv;
end

end