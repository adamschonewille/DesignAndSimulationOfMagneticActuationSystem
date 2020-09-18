function [B, I] = calcMaxNonUniformBFields(mAct_sph, pAct_cartesion, EM_dim)
%calcMaxNonUniformBFields
% INPUTS
% EM_dim    - a vector containing the electromagnet properties in the form
%           of: EM_dim = [m, Diameter, Length]';
% mAct_sph  - A 2xn matrix containing the orientation information [rad] of the
%           EMs including beta (the azimuthal angle about the z-axis) and 
%           gamma (the polar angle from the z-axis).
% pAct_cartesion  - A 3xn matrix containing the cartesion coordinates [m] of
%                 each electromagnet
% OUTPUTS:
% B         - The 8x3 matrix with each column representing the maximum
%           magnetic  field that can be generated in that orientation
%           regardless of field uniformity

%% Setup Parameters
numActuators = size(pAct_cartesion,2);
% output is the max uniform B field in the 3 cardinal directions. Ideally:
%      | Bx 0  0  |
% B =  | 0  By 0  |
%      | 0  0  Bz | 
%      |          |
%      |    G     |
%      |          |
B = zeros(8,3);

O = [0;0;0]; % Origin where the dipole field is evaluated at
z_hat = [0;0;1]; % the z-direction unit vector.

%% Field Strength Objective Functions

m_mag = EM_dim(1);
D = EM_dim(2);
L = EM_dim(3);
B_dip = zeros(8,numActuators);
% Calc B and Z at the same time because both go through n loops.
for i = 1:numActuators
    % find EM magnetic moment axis using the orientation angles 
    Rzy = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat; 
    % Calc magnetic field and add to previous magnetic field contributions
    B_dip(1:3,i) = dipoleField( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy );
    B_dip(4:8,i) = dipoleGradient( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy );
end
% Should now have a 8xN matrix of field and gradient components for each
% actuator. In order to be able to control the gradient terms we will need
% at least 8 actuators. With exactly 8 actuators this is simple matrix
% inversion, but for these calculations we will use the psuedoinverse, so
% the control matrix does not need to be square.

B_desx = [1 0 0]'; % [T] this is large, but will be accounted for 
B_desy = [0 1 0]'; % [T] this is large, but will be accounted for
B_desz = [0 0 1]'; % [T] this is large, but will be accounted for
% Use the psuedo inverse of jsut the field components
B_dip_inv = pinv(B_dip(1:3,:));
I_x = B_dip_inv * B_desx;
I_y = B_dip_inv * B_desy;
I_z = B_dip_inv * B_desz;
% Find absolute maximum of I to find the component with limiting current
I_xMaxVal = max( abs(I_x) );
I_yMaxVal = max( abs(I_y) );
I_zMaxVal = max( abs(I_z) );
% these next steps are necessary as the Gadients will need to be calculated 
% as well
I_xMax = I_x/I_xMaxVal;
I_yMax = I_y/I_yMaxVal;
I_zMax = I_z/I_zMaxVal;
B = [B_dip*I_xMax B_dip*I_yMax B_dip*I_zMax];
I = [I_xMax I_yMax I_zMax];
% B now contains the max field able to be produced by these actuators. In
% order to maximize the field, G should be equal to 1/B since we are trying
% to set G=0 this will maximize B. Scaling is also needed to ensure the
% optimzation converges at a reasonable rate.

    
end