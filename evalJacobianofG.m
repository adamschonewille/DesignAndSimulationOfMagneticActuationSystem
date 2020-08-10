function [J_G] = evalJacobianofG(z_0, EM_dim, mAct_sph, pAct_cartesion)
%evalJacobianofG  -  Numerically evaluates the jacobian of G at parameter
%                   inputs
% INPUTS
% z_0       - The nominal distance below the origin to the table surface
% EM_dim    - a vector containing the electromagnet properties in the form
%           of: EM_dim = [m, Diameter, Length]';
% mAct_sph  - A 2xn matrix containing the orientation information [rad] of the
%           EMs including beta (the azimuthal angle about the z-axis) and 
%           gamma (the polar angle from the z-axis).
% pAct_cartesion  - A 3xn matrix containing the cartesion coordinates [m] of
%                 each electromagnet
% OUTPUTS:
% J_G         - The NxN matrix of numerical gradient components of G 
%               evaluated at the input parameters.

% define parameters:
x = [ reshape(pAct_cartesion',1,[]), reshape(mAct_sph',1,[]) ]';
numActuators = size(pAct_cartesion,2);

% choose some really small number for numerically finding the derivative.
dx = 0.0001; 

G_nom = evalG(z_0, EM_dim, mAct_sph, pAct_cartesion);

N = length(G_nom);
numParams = length(x);

J_G = zeros(N,numParams);

for i = 1:numParams
    x_plus_dx = x;
    x_plus_dx(i,1) = x_plus_dx(i,1) + dx; 
    temp_params = reshape(x_plus_dx',numActuators,[])';
    pAct_cartesion_dx = temp_params(1:3,:);
    mAct_sph_dx = temp_params(4:5,:);
    G_dx = evalG(z_0, EM_dim, mAct_sph_dx, pAct_cartesion_dx);
    % Calc each colum of the jacobian using the definition of the
    % derivative
    J_G(:,i) = 1/dx * ( G_dx - G_nom );
    % this could be optimized as there are a lot of 0-valued entries in the
    % jacobian based on the theory.
end

% the jacobian should now be filled in.

end
