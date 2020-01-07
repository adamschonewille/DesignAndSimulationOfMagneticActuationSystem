function [B] = dipoleField(position, mag_moment)
% AUTHOR: Adam Schonewille
% DATE: January 7th 2020
% ABOUT:    Finds the magnetic field at position for a magnet with magnetic
%           moment mag_moment
% 
%% INPUTS
% position   -  A vector of length 3
% mag_moment -  A vector of length 3
% 
%% RETURNED OUTPUTS
% B  - The 3 x 1 vector containing Bx, By, and Bz field components
%
%
% See also DIPOLEGRADIENT.

%% Function Code
% position and mag_moment should both be 3x1 vectors
% Define shorthand terms
p = position;
m = mag_moment;
p_hat = p/norm(p);

B = 1e-7/(norm(p)^3)*(3*p_hat*p_hat'-eye(3))*m;

end