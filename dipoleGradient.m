function [G] = dipoleGradient(position, mag_moment)
% AUTHOR: Adam Schonewille
% DATE: January 7th 2020
% ABOUT:    Finds the 5 unique gradient terms of a magnet with magnetic
% moment mag_moment at position assuming the dipole model is valid
% 
%% INPUTS
% position   -  A vector of length 3
% mag_moment -  A vector of length 3
% 
%% RETURNED OUTPUTS
% G  - The 5 x 1 vector of field gradient components including: 
%               dBxdx, dBxdy dBxdz dBydy dBydz
% 
% See also DIPOLEFIELD.

%% Function Code
%dipole pointing in z?

% position and mag_moment should both be 3x1 vectors

% constant out front, alpha, is mu_0 / (4 pi) , but mu_0 is 4pi e-7 so:
alpha = 1e-7;
coeff_x = [4 3 3]';
coeff_y = [3 4 3]';
p = position;
m = mag_moment;

B_submatrix = (3*p*p'-norm(p)^2*eye(3))*m;
% Derived Derivatives of B
dBxdx = alpha * ( ( dot(coeff_x.*p,m) )*( norm(p) )^(-5) + ...
    (-5)*p(1)*( norm(p) )^(-7)*( B_submatrix(1) ) );
dBxdy = alpha * ( ( 3*p(1)*m(2) - 2*m(1)*p(2) )*( norm(p) )^(-5) + ...
    (-5)*p(2)*( norm(p) )^(-7)*( B_submatrix(1) ) );
dBxdz = alpha * ( ( 3*p(1)*m(3) - 2*m(1)*p(3) )*( norm(p) )^(-5) + ...
    (-5)*p(3)*( norm(p) )^(-7)*( B_submatrix(1) ) );
dBydy = alpha * ( ( dot(coeff_y.*p,m) )*( norm(p) )^(-5) + ...
    (-5)*p(2)*( norm(p) )^(-7)*( B_submatrix(2) ) );
dBydz = alpha * ( ( 3*p(2)*m(3) - 2*m(2)*p(3) )*( norm(p) )^(-5) + ...
    (-5)*p(3)*( norm(p) )^(-7)*( B_submatrix(2) ) );

G = [dBxdx dBxdy dBxdz dBydy dBydz]';


end