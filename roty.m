function [Ry] = roty(angle)
% AUTHOR:   Adam Schonewille
% ABOUT:    Generates the 3x3 transformation matrix for a rotation of
%       `angle` about the y-axis
% INPUT: angle is an angle in radians

Ry = [cos(angle),  0, sin(angle);
      0,           1,          0;
      -sin(angle), 0, cos(angle)];
end