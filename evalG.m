function [G] = evalG(z_0, EM_dim, mAct_sph, pAct_cartesion)
%evalG
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
% G         - The Nx1 vector of G evaluated at the input parameters.

%% CODE

%% Setup Parameters
numActuators = size(pAct_cartesion,2);
B = zeros(3,1); 
Z = zeros(numActuators,1); % there are n actuators that must be below the table plane
V = zeros(numActuators*(numActuators-1)/2,1); % size is n(n-1)/2
DV = zeros(numActuators*(numActuators-1)/2,1); % size is n(n-1)/2
DA = zeros(numActuators*(numActuators-1)/2,1); % size is n(n-1)/2
O = [0;0;0]; % Origin where the dipole field is evaluated at
z_hat = [0;0;1]; % the z-direction unit vector.

%% Field Strength Objective Functions

m_mag = EM_dim(1);
D = EM_dim(2);
L = EM_dim(3);
theta = atan(D/L);
B_dip = zeros(3,numActuators);
B_dip_G = zeros(5,numActuators);
% Calc B and Z at the same time because both go through n loops.
for i = 1:numActuators
    % find EM magnetic moment axis using the orientation angles 
    Rzy = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat; 
    % Calc magnetic field and add to previous magnetic field contributions
    B(:,1) = B(:,1) + dipoleField( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy);
    B_dip(:,i) = dipoleField( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy);
    B_dip_G(:,i) = dipoleGradient( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy );
    
    % Calc the z height constraint.   
    phi = mod( abs(mAct_sph(2,i)), pi/2 );
    if ~(phi == mod( abs(mAct_sph(2,i)), pi ) )
        % if the component of the modulus is off the x axis then 
        % take the complement 
        phi = pi/2 - phi;
    end
    h = 1/2*sqrt(L^2+D^2)*sin(pi/2-(phi-theta));
    min_h = D/2; % Small radius of the cylinder
    % this check shoud never happen but okay:
    if ( abs(h) < min_h )
        h = min_h * sign(h);
    end
%     Z(i,1) = pAct_cartesion(3,i) + z_0 + h;
    Z(i,1) = pAct_cartesion(3,i) + h;
    Z(i,1) = -log( -Z(i,1) / z_0 );
end
% to maximize Bz set the optimization to 1/Bz -> 0
% B(3,1) = 1/B(3,1);
d = z_0 + L/2;
Bmax = numActuators * norm( dipoleField([d;0;0],[m_mag;0;0]) );
Gmax = numActuators * norm( dipoleGradient([d;0;0],[m_mag;0;0]) );
% B(3,1) = 1-B(3,1)/0.05;
B(3,1) = 0.2*Bmax/B(3,1);
B(1:2,1) = B(1:2,1)/Bmax;

% [U_3,S_3,V_3] = svd(B_dip)
[U_8,S_8,V_8] = svd([B_dip/Bmax; B_dip_G/Gmax])

isotropicFields = true;
if (isotropicFields)
    B_desx = [1 0 0 0 0 0 0 0]';
    B_desy = [0 1 0 0 0 0 0 0]';
    B_desz = [0 0 1 0 0 0 0 0]';
    B_dip_inv = inv( [B_dip; B_dip_G] );
    I_x = B_dip_inv * B_desx;
    I_y = B_dip_inv * B_desy;
    I_z = B_dip_inv * B_desz;
    % Find absolute maximum of I to find the limiting current
    I_xMax = max( abs(I_x) );
    I_yMax = max( abs(I_y) );
    I_zMax = max( abs(I_z) );
    % these next steps are unecessary as the value of the field will just be 1
    % divided by the max current.
    % I = I_x/I_xMax;
    % B = B_dip*I;

    B = [1/I_xMax; 1/I_yMax; 1/I_zMax]
    % B now contains the max field able to be produced by these actuators. In
    % order to maximize the field, G should be equal to 1/B since we are trying
    % to set G=0 this will maximize B. Scaling is also needed to ensure the
    % optimzation converges at a reasonable rate.
    B = 0.00001 * 1./B;
    B = 2*B;
    % Enforce isotropy (Bx = By and Bz = Bx)
    % This adds two new constraint/objective functions
%     B = [B; 0.025*(1-B(1)/B(3)); 0.025*(1-B(1)/B(2))];
    [MaxVal index] = max( abs(U_8(1:3,:)') );
    %Find singular values that are most sensitive to the b field components
    singular_vals = [S_8(index(1),index(1)), S_8(index(2),index(2)), S_8(index(3),index(3)) ]; 
    B = [B; 0.0001*(1/min(singular_vals)); 0.025*(1-max(singular_vals)/min(singular_vals))];
%     B = [B; 0.0001*(1/S_8(3,3)); 0.025*(1-S_8(1,1)/S_8(3,3))];
else
    B_desx = [1 0 0]';
    B_desy = [0 1 0]';
    B_desz = [0 0 1]';
    B_dip_inv = pinv(B_dip);
    I_x = B_dip_inv * B_desx;
    I_y = B_dip_inv * B_desy;
    I_z = B_dip_inv * B_desz;
    % Find absolute maximum of I to find the limiting current
    I_xMax = max( abs(I_x) );
    I_yMax = max( abs(I_y) );
    I_zMax = max( abs(I_z) );
    % these next steps are unecessary as the value of the field will just be 1
    % divided by the max current.
    % I = I_x/I_xMax;
    % B = B_dip*I;

    B = B_desx/I_xMax + B_desy/I_yMax + B_desz/I_zMax
    % B now contains the max field able to be produced by these actuators. In
    % order to maximize the field, G should be equal to 1/B since we are trying
    % to set G=0 this will maximize B. Scaling is also needed to ensure the
    % optimzation converges at a reasonable rate.
    B = 0.001 * 1./B;
    B = 2*B;
    % Enforce anisotropy (Bx = By and Bz = 2Bx)
    % This adds two new constraint/objective functions
    B = [B; 1-0.5*B(1)/B(3); 1-B(1)/B(2)];
end

% Z = 1./Z.^2;

% %% Isotropic or Uniformity Penalty functions
% % Add to this multiple Penalty Functions for the "uniformity" of the field
% a = 0.03; % bounding cube dist from origin
% % Test points are the vertices of a 0.06 m sidelength cube as well as where
% % the principle axis intersect the cube surface and origin (27 pts total)
% pTool_test = [a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a ;
%               0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a ;
%               0  0  0  0  0  0  0  0  a  a  a  a  a  a  a  a  a -a -a -a -a -a -a -a -a -a ];
% % B_test = zeros( size(pTool_test),3 );
% B_test = zeros( size(pTool_test,1),size(pTool_test,2),3 );
% 
% for j = 1:size(pTool_test,2)
%     B_dip = zeros(3,numActuators);
%     for i = 1:numActuators
%         % find EM magnetic moment axis using the orientation angles 
%         Rzy = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat; 
%         % Calc magnetic field and add to previous magnetic field contributions
%         B_dip(:,i) = dipoleField( ( pTool_test(:,j) - (pAct_cartesion(:,i)) ), m_mag*Rzy);
%     end 
%     %I = I_x/I_xMax;
%     % Using the currents that generate the max Bx By Bz fields at the
%     % origin, calculate what the fields are on the vertices of the bounding
%     % cube.
%     B_test(:,j,1) = B_dip*(I_x/I_xMax);
%     B_test(:,j,2) = B_dip*(I_y/I_yMax);
%     B_test(:,j,3) = B_dip*(I_z/I_zMax);
%     % The magnitude and direction of these values can be compared to the
%     % values at the origin to determine the uniformity of the max field in
%     % the space.
% end
% 
% % Comparing the magnitude and direction (2) of the origin field (3) and the 
% % field at the vertices (26) gives us 2 x 3 x 26 = 156 penalty functions
% % Magnitude Penalty is calculated as follows:
% % m_p_x = 1 - |B_testx|/|BxMax|
% % m_p_y = 1 - |B_testy|/|ByMax|
% % m_p_z = 1 - |B_testz|/|BzMax|
% % Direction Penalty is calculated as follows:
% % d_p_x = arccos( B_testx * BxMax / (|B_testx| |BxMax|) )
% % d_p_y = arccos( B_testy * ByMax / (|B_testy| |ByMax|) )
% % d_p_z = arccos( B_testz * BzMax / (|B_testz| |BzMax|) )
% % These penalty functions tend towards 0 as their conditions are met (they
% % become equal).
% mag_penalty = zeros(size(pTool_test,2)*3,1);
% dir_penalty = zeros(size(pTool_test,2)*3,1);
% for j = 1:size(pTool_test,2)
%     mag_penalty(3*j-2,1) = 1-norm(B_test(:,j,1))/(norm( B_desx/I_xMax ));
%     mag_penalty(3*j-1,1) = 1-norm(B_test(:,j,2))/(norm( B_desy/I_yMax ));
%     mag_penalty(3*j , 1) = 1-norm(B_test(:,j,3))/(norm( B_desz/I_zMax ));
%     
%     dir_penalty(3*j-2,1) = acos( B_test(:,j,1)'*B_desx/I_xMax / (norm(B_test(:,j,1))*norm(B_desx/I_xMax)) );
%     dir_penalty(3*j-1,1) = acos( B_test(:,j,2)'*B_desy/I_yMax / (norm(B_test(:,j,2))*norm(B_desy/I_yMax)) );
%     dir_penalty(3*j , 1) = acos( B_test(:,j,3)'*B_desz/I_zMax / (norm(B_test(:,j,3))*norm(B_desz/I_zMax)) );
% end
% 
% % mag_penalty = 1-norm(B_test(:,j,1))/(norm( B_desx/I_xMax ))
% % dir_penalty = acos( B_test(:,j,1)*B_desx/I_xMax / (norm(B_test(:,j,1))*norm(B_desx/I_xMax))
% 
% B_uniform = [mag_penalty; dir_penalty];
% % include these penalty functions to the B objective functions
% B = [B; 10*B_uniform];

%% Proximity Penalty Functions

% Calc V, the overlap in cylinder volumes approximated by multiple spheres
% Calc DV, the distance between in cylinder volumes approximated by multiple spheres
% Calc DA, the minimum distance between two cylinders when modeled as
% inifinite/finite axes.
index = 1;
for i = 1:numActuators-1
    Rzy0 = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat;
    for j = i+1:numActuators
        Rzy1 = rotz(mAct_sph(1,j))*roty(mAct_sph(2,j))*z_hat; 
%         V(index,1) = calcSphericalVolumeOverlap(pAct_cartesion(:,i), Rzy0, D/2, L, ...
%                                                 pAct_cartesion(:,j), Rzy1, D/2, L);
%         DV(index,1) = calcDistanceBetweenSpherical(pAct_cartesion(:,i), Rzy0, D/2, L, ...
%                                                   pAct_cartesion(:,j), Rzy1, D/2, L);
        DA(index,1) = calcMinDistBetweenCylinders(pAct_cartesion(:,i), Rzy0, D/2, L, ...
                                                  pAct_cartesion(:,j), Rzy1, D/2, L);                                      
        index = index + 1;
    end
    
end
% % Scale V by a constant
% Rzy = rotz(mAct_sph(1,1))*roty(mAct_sph(2,1))*z_hat;
% % sigma is the max overlapping volume possible
% sigma = calcSphericalVolumeOverlap(pAct_cartesion(:,1), Rzy, D/2, L, ...
%                                                 pAct_cartesion(:,1), Rzy, D/2, L);
% % V is the overlapping volume bounded from [0,1]
% V = V/sigma;
% G = [B; Z; V];

% % Take the inverse of the distance and scale appropriately
% sigma = 0.00001; 
% DV = sigma./DV;

% Take the inverse of the distance and scale appropriately
sigma = 0.001;
DA = sigma./DA;

%% All Objective functions together

G = [B; Z; DA];

end