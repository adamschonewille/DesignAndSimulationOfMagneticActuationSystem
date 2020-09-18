function [U_3, S_3, V_3, U_8, S_8, V_8] = calcSystemSVD(mAct_sph, pAct_cartesion, EM_dim)
%calcSystemSVD
% INPUTS
% EM_dim    - a vector containing the electromagnet properties in the form
%           of: EM_dim = [m, Diameter, Length]';
% mAct_sph  - A 2xn matrix containing the orientation information [rad] of the
%           EMs including beta (the azimuthal angle about the z-axis) and 
%           gamma (the polar angle from the z-axis).
% pAct_cartesion  - A 3xn matrix containing the cartesion coordinates [m] of
%                 each electromagnet
% OUTPUTS:
% SVD_3     - The 3 matrices containing the singular value decomposition 
%           (SVD) of the control matrix whene only magnetic field
%           components are considered. This decomposition describes the
%           relations or transformations between input currents and out
%           magnetic field components and the sensitivity to certain field
%           components.
% SVD_8     - The 3 matrices containing the singular value decomposition 
%           (SVD) of the control matrix whene the magnetic field
%           components and gradient components are considered. This 
%           decomposition describes the relations or transformations 
%           between input currents and out magnetic field components and 
%           the sensitivity to certain field components.

numActuators = size(pAct_cartesion,2);
O = [0;0;0]; % Origin where the dipole field is evaluated at
z_hat = [0;0;1]; % the z-direction unit vector.

m_mag = EM_dim(1);
D = EM_dim(2);
L = EM_dim(3);
B_dip = zeros(3,numActuators);
B_dip_G = zeros(5,numActuators);
% Calc B and Z at the same time because both go through n loops.
for i = 1:numActuators
    % find EM magnetic moment axis using the orientation angles 
    Rzy = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat; 
    % Calc magnetic field and gradient control matrices
    B_dip(:,i) = dipoleField( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy);
    B_dip_G(:,i) = dipoleGradient( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy );
end
d = 0.012 + L/2;
Bmax = numActuators * norm( dipoleField([d;0;0],[m_mag;0;0]) );
Gmax = numActuators * norm( dipoleGradient([d;0;0],[m_mag;0;0]) );

[U_3,S_3,V_3] = svd(B_dip);
[U_8,S_8,V_8] = svd([B_dip/Bmax; B_dip_G/Gmax]);

end