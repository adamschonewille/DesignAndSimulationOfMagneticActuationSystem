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

radiusVec = large_EM(2)/2*ones(1,numActuators);
lengthVec = large_EM(3)*ones(1,numActuators);


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


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% Randomize the Initial inputs:
% mAct_sph = rand(2,numActuators);
% pAct_cartesion = rand(3,numActuators);

% % Randomize the Initial inputs (with scaling):
mAct_sph = [2*pi*rand(1,numActuators);
             pi/2*rand(1,numActuators)];
pAct_cartesion = [rand(1,numActuators)-0.5;
                  rand(1,numActuators)-0.5;
                  -(0.5-nomDist)*rand(1,numActuators)-0.15];

% Ensure random inputs are valid
validConfig = false;
% keep finding configurations and checking collisions until the
% configuration is valid (no collisions)
while (~validConfig)

%         disp("Generating new Configuration")
    % Randomize Initial inputs (with scaling):
    mAct_sph = [2*pi*rand(1,numActuators);  % [rad] Azimuth
                 pi/2*rand(1,numActuators)];% [rad] Inclination
    pAct_cartesion = [rand(1,numActuators)-0.5; % [m] x-axis from -0.5 to 0.5
                      rand(1,numActuators)-0.5; % [m] y-axis from -0.5 to 0.5
                      -(0.5-nomDist)*rand(1,numActuators)-0.15]; % [m] z-axis from -0.15 to -0.53 for r = 0.12
    % Enforce constraint 1.
    for n=1:numActuators
        theta = atan2(large_EM(2),large_EM(3));
        h = sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2) * sin( pi/2 - (mAct_sph(2,n) - theta) );

%             % Contrain to be only as close as the nomDist plane 
%             if ( pAct_cartesion(3,n) > (0 - nomDist - h) )
%                 pAct_cartesion(3,n) = 0 - nomDist - h;
%             end

        % Contrain to be as close to the nomDist plane as possible
        pAct_cartesion(3,n) = 0 - nomDist - h;

%             % Constrain to a plane
%             pAct_cartesion(3,n) = 0 - nomDist - sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2);
    end
    % approximate the cylinders as spheres and check collision
    % check large radius sphere first (sphere bounds cylinder):
    radius_large = sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2) * ones(numActuators,1);
    % if needed, check small radius sphere second (cylinder bounds sphere):
    radius_small = large_EM(2)/2 * ones(numActuators,1);
    if ( ~checkSphereCollision( pAct_cartesion, radius_large ) )
        validConfig = true;
%             disp("Valid Configuration")
    elseif ( ~checkSphereCollision( pAct_cartesion, radius_small ) )
            % If large radii had collisions but small does not, then 
            % we are here.
            % Refine the search by seeing if the cylinders collide.
%                 disp('Checking Cylinder Collision');
            if ( ~checkAllCylinderCollisions( mAct_sph, pAct_cartesion, radius_small, large_EM(3)*ones(numActuators,1) ) )
                % There is no collisions between cylinders
                validConfig = true;
%                     disp("Valid Configuration")
            end 
    end


end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% Initial User Guess: Based off previous simulation
mAct_sph = [ 5.1509    2.6582    5.3473    2.5260    4.7252    2.3382    5.4930    4.7838;
             0.4719    0.1530    1.3648    1.5399    0.2293    1.1737    1.5134    1.3395];

pAct_cartesion = [ 0.3719   -0.1012    0.0897    0.4597    0.0040    0.2199   -0.0774   -0.5519;
                  -0.5164   -0.1186    0.3318    0.1936    0.2052   -0.0474    0.0523   -0.5468;
                  -0.3111   -0.3065   -0.2227   -0.1931   -0.3092   -0.2509   -0.1769   -0.2270];

% % Initial User Guess:
% mAct_sph = [  0   pi/2  pi  3*pi/2  0  0  0  0;    %Azimuth [rad]
%              pi/6 pi/6 pi/6   pi/6  0  0  0  0];    %Inclination [rad]
% pAct_cartesion = [ -0.25  0     0.25  0    -0.135   0.135   0.135  -0.135;  % X Position
%                     0    -0.25  0     0.25 -0.135  -0.135   0.135   0.135;	% Y Position
%                    -0.32 -0.32 -0.32 -0.32 -0.32   -0.32   -0.32   -0.32];% Z Position

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
dp = 1* 0.001; % [m] about 1 mm.  p -> x, y, z
dtheta = 1* pi/360; % [rad] about 0.5 deg   theta -> beta, gamma
% % irrational steps:
% dp = 4*pi* 0.001; % [m] about 1 mm.  p -> x, y, z
% dtheta = 4*pi* pi/360; % [rad] about 0.5 deg   theta -> beta, gamma

numSteps = 10000; % 1e6; % number of steps to run without convergence.
convergence_tol = 1e-9; % convergence requirement 
numCollisionChecks = numActuators-1;

%% Calculations 
% Calc initial config and enforce/check constraints 
% Calc initial B_vec

%     mAct_sph = [2*pi*rand(1,numActuators);
%                 pi/2*rand(1,numActuators)];
%     pAct_cartesion = [rand(1,numActuators)-0.5;
%                       rand(1,numActuators)-0.5;
%                     -(0.5-nomDist)*rand(1,numActuators)-0.15];
lengthB = 3+numCollisionChecks;
B_net = zeros(lengthB,1);
B_n = zeros(lengthB,numActuators);
x_hat = [1; 0; 0];
y_hat = [0; 1; 0];
z_hat = [0; 0; 1];
Rzy = rotz(mAct_sph(1,1))*roty(mAct_sph(2,1))*z_hat;
m_mag = large_EM(1);
O = [0;0;0];
if (numActuators == 1)
    Rzy = rotz(mAct_sph(1,1))*roty(mAct_sph(2,1))*z_hat;
    B_n(1:3) = dipoleField( (O-pAct_cartesion(:,1)), m_mag*Rzy);
%     B_n(4) = 0;
    B_net = B_n;
else
    for n = 1:numActuators
        Rzy = rotz(mAct_sph(1,n))*roty(mAct_sph(2,n))*z_hat;
        B_n(1:3,n) = dipoleField( (O-pAct_cartesion(:,n)), m_mag*Rzy);
    end
    B_net = sum(B_n,2);
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
B_des = [0;
         0;
         1];
% B_des = [1;
%          1;
%          0;];

B_des = [B_des; zeros(numCollisionChecks,1) ];
% numErrors = -numSteps;
fval = nan;
B_net_prev = [0;0;0;0];
az = 37.5; %for plotting
converge_bool = false;
% Setup derivative step matrices:
dB_dx = zeros(lengthB,numActuators);
dB_dy = zeros(lengthB,numActuators);
dB_dz = zeros(lengthB,numActuators);

dRzy_dbeta = zeros(3,numActuators);
dRzy_dgamma = zeros(3,numActuators);

dB_dbeta = zeros(lengthB,numActuators);
dB_dgamma = zeros(lengthB,numActuators);

M_dels = zeros(lengthB,5,numActuators);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
%%~~~~~~~~~~~~~~~~~~~~~ New Optimization Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~%%
x_params = [ reshape(pAct_cartesion',1,[]), reshape(mAct_sph',1,[]) ]';
% set arbitrary set size
mu = 0.001; 

% mu = [dp*ones(numActuators*3,1); dtheta*ones(numActuators*2,1)];

EM_dim = large_EM;
f_record = zeros(numSteps,0);
G = evalG(nomDist, EM_dim, mAct_sph, pAct_cartesion);
F = 1/2 * G'*G;

simulatedAnnealing = false;


for i = 1:numSteps
    % update F
    f_record(i) = F;
    
    temp_params = reshape(x_params',numActuators,[])';
    pAct_cartesion = temp_params(1:3,:);
    mAct_sph = temp_params(4:5,:);
    J_G = evalJacobianofG(nomDist, EM_dim, mAct_sph, pAct_cartesion);
    G = evalG(nomDist, EM_dim, mAct_sph, pAct_cartesion);
%     G';
    x_prev = x_params;
    % Update F for the next loop
    F_prev = F;
    F = 1/2 * G'*G;
    G_prev = G;
    % Find new parameters, ensure that F is decreasing.
    improvement_bool = false;
%     mu = 0.01 * rand(1,1);
    
    x_params = x_prev - mu .* J_G' * G_prev;
    
%     while ~( improvement_bool )
%         % find new parameters
%         x_params = x_prev - mu .* J_G' * G_prev;
%         temp_params = reshape(x_params',numActuators,[])';
%         pAct_cartesion = temp_params(1:3,:);
%         mAct_sph = temp_params(4:5,:);
%         % Calc G to calc F and ensure that it is decreasing.
%         G = evalG(nomDist, EM_dim, mAct_sph, pAct_cartesion);
%         F = 1/2 * G'*G;
%         diff = F_prev - F
%         % If the function F doesn't decrease, then decrease the step size
%         % and repeat this loop.
%         if ( (F_prev - F) < 0 )
%             mu = mu/10;
%         elseif ( ( (F_prev - F) < convergence_tol ) && (mu > 1) )
%             % if the difference is less than the convergence tolerance and
%             % the step size is very small then confirm convergence
%             converge_bool = true;
%             disp("Simulation has reached convergance tolerance")
%             improvement_bool = true;
%             break;            
%         else
%             % If there is a positive difference (improvement) to the cost 
%             % function then end this loop and continue optimization
%             improvement_bool = true;
%             break;
%         end
%     end


%     x_params = x_prev - mu .* J_G' * G;
   
    % Check convergence:
%     F = 1/2 * G'*G;
%     if ( abs(F_prev - F) < convergence_tol )
%         converge_bool = true;
%         disp("Simulation has reached convergance tolerance")
%         break;
%     end
        
    figure(5)
    semilogy(f_record, 'ok')
    hold on
    set(gcf, 'Position',  [0, 50, 450, 300])
    title('Convergence')
    hold off

    figure(10)
    plot(G, 'o-k')
    hold on
    set(gcf, 'Position',  [0, 450, 450, 300])
    title('Evaluation of G')
    hold off
    
    
    
    % Add simulated annealing noise to the configuration every x steps.
    if (simulatedAnnealing)
        if ( mod(i, 100) == 0 )
                mAct_sph_noise = dtheta * (rand(2*numActuators,1)-0.5);
                pAct_cartesion_noise = dp * (rand(3*numActuators,1)-0.5);
                % Noise should decrease over the iterations
                mAct_sph_noise = mAct_sph_noise * numSteps/100 / i;
                pAct_cartesion_noise = pAct_cartesion_noise * numSteps/100 / i;
                % Add noise
                x_params = x_params + [pAct_cartesion_noise; mAct_sph_noise];
        end
    end
    
    
        % Plot ?
    if (numActuators == 1)
        % duplicate for 1 EM case
        mAct_sph_plot = [mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph];
        pAct_plot = [pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion];
    elseif (numActuators == 2)
        % duplicate for 2 EM case
        mAct_sph_plot = [mAct_sph mAct_sph mAct_sph mAct_sph];
        pAct_plot = [pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion];
    else
        % same same for 8 EM case
        mAct_sph_plot = mAct_sph;
        pAct_plot = pAct_cartesion;
    end
    
    plotEM(nomDist, large_EM, large_EM, rad2deg(mAct_sph_plot), pAct_plot, F, 25); % figure number arbitrarilly set to 25
    hold on
    [x y] = meshgrid(-0.75:0.1:0.75); % Generate x and y data
    z = -nomDist*ones(size(x, 1)); % Generate z data
    surf(x, y, z,'FaceAlpha',0.25) % Plot the surface
    xlim([-0.75,0.75])
    ylim([-0.75,0.75])
    zlim([-0.5,0.25])
    view([ az, 30 ])
    hold off
%     az = az+0.5;

    


    
    
end

if (~converge_bool)
    disp("Simulation ended. Max number of steps reached.")
end

return
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

for i = 1:numSteps
       
    % Find deltas for each actuator
    for n = 1:numActuators
        Rzy = rotz(mAct_sph(1,n))*roty(mAct_sph(2,n))*z_hat;
        dB_dx(1:3,n) = B_n(1:3,n) - dipoleField( ( O - (pAct_cartesion(:,n)+dp*x_hat) ), m_mag*Rzy);
        dB_dy(1:3,n) = B_n(1:3,n) - dipoleField( ( O - (pAct_cartesion(:,n)+dp*y_hat) ), m_mag*Rzy);
        dB_dz(1:3,n) = B_n(1:3,n) - dipoleField( ( O - (pAct_cartesion(:,n)+dp*z_hat) ), m_mag*Rzy);
        dRzy_dbeta(:,n) = rotz(mAct_sph(1,n)+dtheta)*roty(mAct_sph(2,n))*z_hat;
        dRzy_dgamma(:,n) = rotz(mAct_sph(1,n))*roty(mAct_sph(2,n)+dtheta)*z_hat;
        dB_dbeta(1:3,n) = B_n(1:3,n) - dipoleField( ( O - pAct_cartesion(:,n) ), m_mag*dRzy_dbeta(:,n));
        dB_dgamma(1:3,n) = B_n(1:3,n) - dipoleField( ( O - pAct_cartesion(:,n) ), m_mag*dRzy_dgamma(:,n));
        % find how the changes in position and orientation affect
        % collisions by calculating the distance between actuators
        if (numActuators == 1)
%             B_n(4) = 0;
            % dB/d__ terms should all stay 0 no matter what. 
            % no need to set them to 0 again.
        else
            % calc gradient change terms for collisions
            
        end
    
        
        % we can create a 3 x 5 matrix relating 
        M_dels(:,:,n) = [dB_dx(:,n) dB_dy(:,n) dB_dz(:,n) dB_dbeta(:,n) dB_dgamma(:,n)];
        dels_dir = pinv(M_dels(:,:,n))*B_des;
        dels_dir = dels_dir/norm(dels_dir);
        
        %is this the right place for this call?
        B_net_prev = B_net;

%         if (i == numSteps-1)
%             disp("Last iter");
%         end
        
        % find new values for configuration:
        pAct_cartesion(:,n) = pAct_cartesion(:,n) + dp * dels_dir(1:3);
        mAct_sph(:,n) = mAct_sph(:,n) + dtheta * dels_dir(4:5);
        
    end
    
    % Add simulated annealing noise to the configuration every x steps.
    if ( mod(i, 250) == 0 )
            mAct_sph_noise = [dtheta*rand(1,numActuators);
                              dtheta*rand(1,numActuators)];
            pAct_cartesion_noise = [dp*rand(1,numActuators);
                                     dp*rand(1,numActuators);
                                     dp*rand(1,numActuators)];
            % Noise should decrease over the iterations
            mAct_sph_noise = mAct_sph_noise * numSteps/i;
            pAct_cartesion_noise = pAct_cartesion_noise * numSteps/i;
            % Add noise
            mAct_sph = mAct_sph + mAct_sph_noise;
            pAct_cartesion = pAct_cartesion + pAct_cartesion_noise;
    end
    
    
    
    % enforce constraints: Table plane limit on z
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

%         % Contrain to be as close to the nomDist plane as possible
%         pAct_cartesion(3,n) = 0 - nomDist - h;
% 
%         % Constrain to a plane
%         pAct_cartesion(3,n) = 0 - nomDist - sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2);
    end
    
    % enforce constraints: collisions between EMs
    % if two actuators collide find the vector between them and seperate by
    % a proportional amount in the x-y projection of that vector.
    if (numActuators > 1)
        W = zeros(size(pAct_cartesion));
        % set up the unit vectors of the axis for each cylinder.
        for n = 1:numActuators
            Rzy = rotz(mAct_sph(1,n)) * roty(mAct_sph(2,n));
            W(:,n) =  Rzy * z_hat;
        end
        
        validConfig = false;
        
        while ~(validConfig)
            validConfig = true; 
            % set to true. if a collision is found, the config is
            % invalidated and will have to recheck collision on another
            % cycle to ensure that the change was sufficient to avoid
            % collision.
            % if no collisions were found the configuration remains valid.
            for i = 1:numActuators-1
                for j = i+1:numActuators
                    if ( checkCylinderCollision(pAct_cartesion(:,i), W(:,i), radiusVec(i), lengthVec(i), ...
                                                pAct_cartesion(:,j), W(:,j), radiusVec(j), lengthVec(j) ) )
    %                     cond = true;
    %                     break;
                        validConfig = false;
                        p_diff = pAct_cartesion(:,j) - pAct_cartesion(:,i);
                        p_diff_hat = p_diff/norm(p_diff);
                        p_diff_hat_noZ = p_diff_hat;
                        p_diff_hat_noZ(3) = 0;
                        pAct_cartesion(:,j) = pAct_cartesion(:,j) + dp * p_diff_hat_noZ;
                        pAct_cartesion(:,i) = pAct_cartesion(:,i) - dp * p_diff_hat_noZ;
                    end

                end

            end
        end
    end
    
    
    % Recalc new B for next iteration
    for n = 1:numActuators
        Rzy = rotz(mAct_sph(1,n))*roty(mAct_sph(2,n))*z_hat;
        B_n(1:3,n) = dipoleField( (O-pAct_cartesion(:,n)), m_mag*Rzy);
    end
    B_net = sum(B_n,2);
%     Rzy = rotz(mAct_sph(1,1))*roty(mAct_sph(2,1))*z_hat;
%     B = dipoleField( (O-pAct_cartesion(:,1)), m_mag*Rzy);
    if ( sum( abs(B_net-B_net_prev) ) < convergence_tol ) 
        converge_bool = true;
        disp("Simulation has reached convergance tolerance")
        break;
    end
    
    % Plot ?
    if (numActuators == 1)
        % duplicate for 1 EM case
        mAct_sph_plot = [mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph mAct_sph];
        pAct_plot = [pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion];
    elseif (numActuators == 2)
        % duplicate for 2 EM case
        mAct_sph_plot = [mAct_sph mAct_sph mAct_sph mAct_sph];
        pAct_plot = [pAct_cartesion pAct_cartesion pAct_cartesion pAct_cartesion];
    else
        % same same for 8 EM case
        mAct_sph_plot = mAct_sph;
        pAct_plot = pAct_cartesion;
    end
    
    plotEM(nomDist, large_EM, large_EM, rad2deg(mAct_sph_plot), pAct_plot, B_net(3), 25); % figure number arbitrarilly set to 25
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