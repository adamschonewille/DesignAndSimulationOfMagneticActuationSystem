function [mAct_sph pAct_cartesion] = optimizeEMPositions (numActuators, large_EM, small_EM, nomDist, figureNum)
% custum optimization algorithm to find the optimum conditions for 
% actuators below a nominal distance plane
% ___EM = [m_EM_large, d_EM_large, l_EM_large]';
% quadrants in 3D:
% #     x   y   z
% 1     +   +   +
% 2     -   +   +
% 3     -   -   +
% 4     +   -   +
% 5     +   +   -
% 6     -   +   -
% 7     -   -   -
% 8     +   -   -
%

%% CODE


%% Initialize the electromagnet configurations

% 6 actuator example
% mAct_sph = [ %actuator magnetization direction (axes) in spherical coordinates (beta and gamma)
%               0      0 1*2*pi/5 2*2*pi/5 3*2*pi/5 4*2*pi/5;    %Azimuth [deg]
%            pi/2 5*pi/6   5*pi/6   5*pi/6   5*pi/6   5*pi/6];    %Inclination [deg]
% pAct_sph = [ %actuator positions in spherical coordinates (alpha and phi)
%         r_mag r_mag r_mag r_mag r_mag r_mag;    %Radius [m]
%         0    0 1*2*pi/5 2*2*pi/5 3*2*pi/5 4*2*pi/5;    %Azimuth [deg]
%         0 pi/3     pi/3     pi/3     pi/3     pi/3];    %Inclination [deg]
%        

% Initialize to all zeros
% mAct_sph = zeros(2,numActuators);
% pAct_cartesion = zeros(3,numActuators);
% pAct_sph = zeros(3,numActuators);

% Randomize the Initial inputs:
% mAct_sph = rand(2,numActuators);
% pAct_cartesion = rand(3,numActuators);

% Randomize the Initial inputs (with scaling):
mAct_sph = [2*pi*rand(1,numActuators);
             pi/2*rand(1,numActuators)];
pAct_cartesion = [rand(1,numActuators)-0.5;
                  rand(1,numActuators)-0.5;
                  -(0.5-nomDist)*rand(1,numActuators)-0.15];

% Initial User Guess:
% mAct_sph = [  0   pi/2  pi  3*pi/2  0  0  0  0;    %Azimuth [rad]
%              pi/6 pi/6 pi/6   pi/6  0  0  0  0];    %Inclination [rad]
% pAct_cartesion = [ -0.15  0     0.15  0    -0.045   0.045   0.045  -0.045;  % X Position
%                     0    -0.15  0     0.15 -0.045  -0.045   0.045   0.045;	% Y Position
%                    -0.18 -0.18 -0.18 -0.18 -0.15   -0.15   -0.15   -0.15];% Z Position

% Initial User Guess (all vert):
% mAct_sph = [ 0  0  0  0  0  0  0  0;     %Azimuth [rad]
%              0  0  0  0  0  0  0  0];    %Inclination [rad]
% r_a = large_EM(2)/2 + sqrt( (large_EM(2)/2+small_EM(2)/2)^2 + (large_EM(2)/2)^2 );
% pAct_cartesion = [  0.0675  0.0675 -0.0675 -0.0675  r_a   0    -r_a   0   ; % X Position
%                     0.0675 -0.0675 -0.0675  0.0675  0     r_a   0    -r_a ;	% Y Position
%                    -0.15   -0.15   -0.15   -0.15   -0.15 -0.15 -0.15 -0.15];% Z Position

%% Setup Parameters

% Tethered Tool Parameters
pTool = [0 0 0]'; %tool position [m]
mTool_dir = [0 0 -1]'; %tool orientation
mTool_mag = 8.4e-3; %tool dipole moment magnetiude [Am^2]
mTool = mTool_dir/norm(mTool_dir)*mTool_mag; %tool dipole moment [Am^2]

% define x0
fittingParameters = [mAct_sph; pAct_cartesion];
Constants = [nomDist; large_EM; small_EM];

x_range = 1.0; % [m] from -0.5 to 0.5 
y_range = 1.0; % [m] from -0.5 to 0.5
z_range = -nomDist + 0.5; % [m] from -0.5 to -nomDist
beta_range = 2*pi;  % [rad] azimuthal angle about z from the x-axis
gamma_range = pi/2; % [rad] inclination angle about y from z-axis

%% Step Sizes and Opitimzation parameters
dp = 5* 0.001; % [m] about 1 mm.  p -> x, y, z
dtheta = 5* pi/360; % [rad] about 0.5 deg   theta -> beta, gamma
numSteps = 1e4; % number of steps to run without convergence.
convergence_tol = 1e-8; % convergence requirement 


%% Calculations 
% Calc initial config and enforce/check constraints 
% Calc initial B_vec

%     mAct_sph = [2*pi*rand(1,numActuators);
%                 pi/2*rand(1,numActuators)];
%     pAct_cartesion = [rand(1,numActuators)-0.5;
%                       rand(1,numActuators)-0.5;
%                     -(0.5-nomDist)*rand(1,numActuators)-0.15];

B = zeros(3,1);
x_hat = [1; 0; 0];
y_hat = [0; 1; 0];
z_hat = [0; 0; 1];
Rzy = rotz(mAct_sph(1,1))*roty(mAct_sph(2,1))*z_hat;
m_mag = large_EM(1);
O = [0;0;0];
if (numActuators == 1)
    B = dipoleField( (O-pAct_cartesion(:,1)), m_mag*Rzy);
end

%% Optimization Recursion
% Calc dB/dx, dB/dy, dB/dz, dB/dbeta, and dB/dgamma for each actuator
% For clarification:
%   dB/dx = B(x) - B(x+dX)   etc for the other variables also.
% This assumes a unit input of current, equivalent to a permanent magnet
% Determine the appropriate step to take for each actuator
% B_des = [0.25;
%          0.25;
%             1;];
% B_des = [0;
%          0;
%          1;];
B_des = [1;
         1;
         0;];
% numErrors = -numSteps;
fval = nan;
B_prev = [0;0;0];
az = 37.5; %for plotting
converge_bool = false;
for i = 1:numSteps
    % Find deltas
    dB_dx = B - dipoleField( ( O - (pAct_cartesion(:,1)+dp*x_hat) ), m_mag*Rzy);
    dB_dy = B - dipoleField( ( O - (pAct_cartesion(:,1)+dp*y_hat) ), m_mag*Rzy);
    dB_dz = B - dipoleField( ( O - (pAct_cartesion(:,1)+dp*z_hat) ), m_mag*Rzy);
    dRzy_dbeta = rotz(mAct_sph(1,1)+dtheta)*roty(mAct_sph(2,1))*z_hat;
    dRzy_dgamma = rotz(mAct_sph(1,1))*roty(mAct_sph(2,1)+dtheta)*z_hat;
    dB_dbeta = B - dipoleField( ( O - pAct_cartesion(:,1) ), m_mag*dRzy_dbeta);
    dB_dgamma = B - dipoleField( ( O - pAct_cartesion(:,1) ), m_mag*dRzy_dgamma);
    % we can create a 3 x 5 matrix relating 
    M_dels = [dB_dx dB_dy dB_dz dB_dbeta dB_dgamma];
    dels_dir = pinv(M_dels)*B_des;
    dels_dir = dels_dir/norm(dels_dir);
    B_prev = B;
%     if (i == numSteps-1)
%         disp("Last iter");
%     end
    % find new values for configuration:
    pAct_cartesion = pAct_cartesion + dp * dels_dir(1:3);
    mAct_sph = mAct_sph + dtheta * dels_dir(4:5);

    for n=1:numActuators
            theta = atan2(large_EM(2),large_EM(3));
            phi = mod( mAct_sph(2,n), pi/2 );
            if ~(phi == mod( mAct_sph(2,n), pi ) )
                % if the component of the modulus is off the x axis then 
                % take the complement 
                phi = pi/2 - phi;
            end
            h = sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2) * sin( pi/2 - (phi - theta) );
            % Ensure h is good for all angles
            h = abs(h);
            min_h = large_EM(2)/2; % Small radius of the cylinder
            if (h < min_h)
                h = min_h;
            end
            % Contrain to be only as close as the nomDist plane 
            if ( pAct_cartesion(3,n) > (0 - nomDist - h) )
                fprintf("Constraint broken. Z = %d Setting z valid value: %d \n", pAct_cartesion(3,n),(0 - nomDist - h) );
                pAct_cartesion(3,n) = 0 - nomDist - h;
                fprintf("h = %d \n", h);

            end

%             % Contrain to be as close to the nomDist plane as possible
%             pAct_cartesion(3,n) = 0 - nomDist - h;
            
%             % Constrain to a plane
%             pAct_cartesion(3,n) = 0 - nomDist - sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2);
    end
    
    
    % Recalc new B for next iteration
    Rzy = rotz(mAct_sph(1,1))*roty(mAct_sph(2,1))*z_hat;
    B = dipoleField( (O-pAct_cartesion(:,1)), m_mag*Rzy);
    if ( sum( abs(B-B_prev) ) < convergence_tol ) 
        converge_bool = true;
        disp("Simulation has reached convergance tolerance")
        break;
    end
    % Plot ?
    % duplicate for 1 EM case
    mAct_sph_plot = [mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph];
    pAct_plot = [pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion];
    
    plotEM(nomDist, large_EM, large_EM, rad2deg(mAct_sph_plot), pAct_plot, B(3), 25); % figure number arbitrarilly set to 25
    hold on
    [x y] = meshgrid(-0.75:0.1:0.75); % Generate x and y data
    z = -nomDist*ones(size(x, 1)); % Generate z data
    surf(x, y, z,'FaceAlpha',0.25) % Plot the surface
    zlim([-0.5,0.25])
    view([ az, 30 ])
    hold off
    az = az+0.5;
%     validConfig = false;
%     % keep finding configurations and checking collisions until the
%     % configuration is valid (no collisions)
%     while (~validConfig)
%         numErrors = numErrors + 1;
% %         disp("Generating new Configuration")
%         % Randomize Initial inputs (with scaling):
%         mAct_sph = [2*pi*rand(1,numActuators);  % [rad] Azimuth
%                      pi/2*rand(1,numActuators)];% [rad] Inclination
%         pAct_cartesion = [rand(1,numActuators)-0.5; % [m] x-axis from -0.5 to 0.5
%                           rand(1,numActuators)-0.5; % [m] y-axis from -0.5 to 0.5
%                           -(0.5-nomDist)*rand(1,numActuators)-0.15]; % [m] z-axis from -0.15 to -0.53 for r = 0.12
%         % Enforce constraint 1.
%         for n=1:numActuators
%             theta = atan2(large_EM(2),large_EM(3));
%             h = sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2) * sin( pi/2 - (mAct_sph(2,n) - theta) );
%            
% %             % Contrain to be only as close as the nomDist plane 
% %             if ( pAct_cartesion(3,n) > (0 - nomDist - h) )
% %                 pAct_cartesion(3,n) = 0 - nomDist - h;
% %             end
% 
%             % Contrain to be as close to the nomDist plane as possible
%             pAct_cartesion(3,n) = 0 - nomDist - h;
%             
% %             % Constrain to a plane
% %             pAct_cartesion(3,n) = 0 - nomDist - sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2);
%         end
%         % approximate the cylinders as spheres and check collision
%         % check large radius sphere first (sphere bounds cylinder):
%         radius_large = sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2) * ones(numActuators,1);
%         % if needed, check small radius sphere second (cylinder bounds sphere):
%         radius_small = large_EM(2)/2 * ones(numActuators,1);
%         if ( ~checkSphereCollision( pAct_cartesion, radius_large ) )
%             validConfig = true;
% %             disp("Valid Configuration")
%         elseif ( ~checkSphereCollision( pAct_cartesion, radius_small ) )
%                 % If large radii had collisions but small does not, then 
%                 % we are here.
%                 % Refine the search by seeing if the cylinders collide.
% %                 disp('Checking Cylinder Collision');
%                 if ( ~checkAllCylinderCollisions( mAct_sph, pAct_cartesion, radius_small, large_EM(3)*ones(numActuators,1) ) )
%                     % There is no collisions between cylinders
%                     validConfig = true;
% %                     disp("Valid Configuration")
%                 end 
%         end
%             
%         
%     end
%     % should have 1 new randomized collision-free configuration now.
%     % find the fitness of this random set-up


end

if (~converge_bool)
    disp("Simulation ended. Max number of steps reached.")
end






end