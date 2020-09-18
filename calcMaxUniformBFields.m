function [B, I] = calcMaxUniformBFields(mAct_sph, pAct_cartesion, EM_dim)
%calcMaxUniformBFields
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
%           uniform magnetic  field that can be generated in that orientation

%% Setup Parameters
numActuators = size(pAct_cartesion,2);
% output is the max uniform B field in the 3 cardinal directions. Ideally:
%      | Bx 0  0  |
% B =  | 0  By 0  |
%      | 0  0  Bz |  
%      |          |
%      |    G     |
%      |          |
% G is a 5x3 matrix of zeros (ideally)
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
% inversion.

if (numActuators == 8)
    B_desx = [1 0 0 0 0 0 0 0]'; % [T] this is large, but will be accounted for 
    B_desy = [0 1 0 0 0 0 0 0]'; % [T] this is large, but will be accounted for
    B_desz = [0 0 1 0 0 0 0 0]'; % [T] this is large, but will be accounted for
    B_dip_inv = inv(B_dip);
    I_x = B_dip_inv * B_desx;
    I_y = B_dip_inv * B_desy;
    I_z = B_dip_inv * B_desz;
    % Find absolute maximum of I to find the component with limiting current
    I_xMax = max( abs(I_x) );
    I_yMax = max( abs(I_y) );
    I_zMax = max( abs(I_z) );
    % these next steps are unecessary as the value of the field will just be 1
    % divided by the max current (since we chose 1 T as our desired field).
    % I = I_x/I_xMax;
    % B = B_dip*I;
    % Is the same as:
    % B = Bdes/I_xMax
    % since 
    % B = B_dip*I = B_dip*(I_x/I_xMax) = (B_dip*I_x)/I_xMax = Bdes/I_xMax
    % In other words, we get the samee results if with scale the input 
    % current by Imax than if we just scale the desired field by Imax.
    
    % Note: there are several scenarios here. The max current of the 8
    % actuators could be less than unity meaning that the actuators are not
    % saturated at the desired 1 T field. Dividing by the max current
    % will increase the acheivable output field beyond 1 T. This scenario
    % is less likely than the currents being much greater than unity and
    % the results being scaled back so that the actuators are not
    % oversaturated

    B = [B_desx/I_xMax B_desy/I_yMax B_desz/I_zMax];
    % B now contains the max field able to be produced by these actuators. In
    % order to maximize the field, G should be equal to 1/B since we are trying
    % to set G=0 this will maximize B. Scaling is also needed to ensure the
    % optimzation converges at a reasonable rate.
    I = [I_x/I_xMax I_y/I_yMax I_z/I_zMax];
else
    error('There are not 8 actuators. It is necessary to have at least 8 actuators to control uniform fields')
end
    
end