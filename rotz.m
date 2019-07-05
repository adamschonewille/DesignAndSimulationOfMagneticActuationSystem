function [Rz] = rotz(angle) 
% AUTHOR:   Adam Schonewille
% ABOUT:    Generates the 3x3 transformation matrix for a rotation of
%       `angle` about the z-axis

Rz = [cos(angle), -sin(angle),  0;
      sin(angle),  cos(angle),  0;
      0,                    0,  1,];
end