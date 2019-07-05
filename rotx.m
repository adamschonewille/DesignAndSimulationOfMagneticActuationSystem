function [Rx] = rotx(angle) 
% AUTHOR:   Adam Schonewille
% ABOUT:    Generates the 3x3 transformation matrix for a rotation of
%       `angle` about the x-axis

Rx = [1,           0,           0;
      0,  cos(angle), -sin(angle);
      0,  sin(angle),  cos(angle)];
end