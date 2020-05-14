function [Eul] = Quat2Eul(q)
% Takes quaternions and calculates the equivalent Euler angles in degrees

% Inputs
% q = [q1,q2,q3,q4] = [q0,q1,q2,q3] quaternions with q0 being the "scalar" value

% Outputs
% bank,pitch,azimuth = Euler angles (rad) in 1,2,and 3 axis 
% bank,pitch,azimuth = roll,pitch,yaw

% Rename input data for matrix calculations
a = q(:,1);
b = q(:,2);
c = q(:,3);
d = q(:,4);

% Conversion!
m11 = 2 .* (b .* c + a .* d);
m12 = a .^ 2 + b .^ 2 - c .^ 2 - d .^ 2;
m21 = -2 .*( b .* d - a .*c);
m31 = 2 .* (c .* d + a .* b);
m32 = a .^ 2 - b .^ 2 - c .^ 2 + d .^ 2;

bank = atan2(m31,m32);
pitch1 = asin(m21);
azimuth = atan2(m11,m12);

% convert from radians to degrees
roll = bank*180/pi;
pitch = pitch1*180/pi;
yaw = azimuth*180/pi;

% Prepare for output
Eul = [roll,pitch,yaw];
return