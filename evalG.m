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

numActuators = size(pAct_cartesion,2);
B = zeros(3,1); 
Z = zeros(numActuators,1); % there are n actuators that must be below the table plane
V = zeros(numActuators*(numActuators-1)/2,1); % size is n(n-1)/2
DV = zeros(numActuators*(numActuators-1)/2,1); % size is n(n-1)/2
O = [0;0;0]; % Origin where the dipole field is evaluated at
z_hat = [0;0;1]; % the z-direction unit vector.

m_mag = EM_dim(1);
D = EM_dim(2);
L = EM_dim(3);
theta = atan(D/L);
% Calc B and Z at the same time because both go through n loops.
for i = 1:numActuators
    % find EM magnetic moment axis using the orientation angles 
    Rzy = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat; 
    % Calc magnetic field and add to previous magnetic field contributions
    B(:,1) = B(:,1) + dipoleField( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy);
    
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
% B(3,1) = 1-B(3,1)/0.05;
B(3,1) = 0.2*Bmax/B(3,1);
B(1:2,1) = B(1:2,1)/Bmax;


% Z = 1./Z.^2;

% Calc V, the overlap in cylinder volumes approximated by multiple spheres
% Calc DV, the distance between in cylinder volumes approximated by multiple spheres
index = 1;
for i = 1:numActuators-1
    Rzy0 = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat;
    for j = i+1:numActuators
        Rzy1 = rotz(mAct_sph(1,j))*roty(mAct_sph(2,j))*z_hat; 
%         V(index,1) = calcSphericalVolumeOverlap(pAct_cartesion(:,i), Rzy0, D/2, L, ...
%                                                 pAct_cartesion(:,j), Rzy1, D/2, L);
        DV(index,1) = calcDistanceBetweenSpherical(pAct_cartesion(:,i), Rzy0, D/2, L, ...
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

% Take the inverse of the distance and scale appropriately
sigma = 0.0001;
DV = sigma./DV;

G = [B; Z; DV];

end