function [mAct_sph pAct_cartesion] = MonteCarloEMPositions(numActuators, large_EM, small_EM, nomDist, figureNum, quadrants)
% finds optimum conditions for actuators below nomDist plane

% small_EM = [m_EM_small, d_EM_small, l_EM_small]';
% large_EM = [m_EM_large, d_EM_large, l_EM_large]';

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% 6 actuator example
% mAct_sph = [ %actuator magnetization direction (axes) in spherical coordinates (beta and gamma)
%               0      0 1*2*pi/5 2*2*pi/5 3*2*pi/5 4*2*pi/5;    %Azimuth [deg]
%            pi/2 5*pi/6   5*pi/6   5*pi/6   5*pi/6   5*pi/6];    %Inclination [deg]
% pAct_sph = [ %actuator positions in spherical coordinates (alpha and phi)
%         r_mag r_mag r_mag r_mag r_mag r_mag;    %Radius [m]
%         0    0 1*2*pi/5 2*2*pi/5 3*2*pi/5 4*2*pi/5;    %Azimuth [deg]
%         0 pi/3     pi/3     pi/3     pi/3     pi/3];    %Inclination [deg]
%        

%% ~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTRAINTS   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% 1. The actuators should be as close to the operating plane as possible.
%    This means: The cylinders should be touching the operating plane.


%% ~~~~~~~~~~~~~~~~~~~ INITIALIZING PARAMETERS  ~~~~~~~~~~~~~~~~~~~~~~~~ %%
% Initialize all to zero:
% mAct_sph = zeros(2,numActuators);
% pAct_cartesion = zeros(3,numActuators);
%%% pAct_sph = zeros(3,numActuators);
% Randomize the Initial inputs:
% mAct_sph = rand(2,numActuators);
% pAct_cartesion = rand(3,numActuators);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Randomize Initial inputs (with scaling):
% mAct_sph = [2*pi*rand(1,numActuators);  % [rad] Azimuth
%              pi/2*rand(1,numActuators)];% [rad] Inclination
% pAct_cartesion = [rand(1,numActuators)-0.5; % [m] x-axis from -0.5 to 0.5
%                   rand(1,numActuators)-0.5; % [m] y-axis from -0.5 to 0.5
%                   -(0.5-nomDist)*rand(1,numActuators)-0.15]; % [m] z-axis from -0.15 to -0.53 for r = 0.12

% All actuators vertical:
% mAct_sph = [ 0  0  0  0  0  0  0  0;     % Azimuth [rad]
%              0  0  0  0  0  0  0  0];    % Inclination [rad]
% r_a = large_EM(2)/2 + sqrt( (large_EM(2)/2+small_EM(2)/2)^2 + (large_EM(2)/2)^2 ) +0.25;
% h_o = large_EM(3)/2 + nomDist;
% pAct_cartesion = [  0.0675  0.0675 -0.0675 -0.0675  r_a   0    -r_a   0   ;  % X Position
%                     0.0675 -0.0675 -0.0675  0.0675  0     r_a   0    -r_a ;  % Y Position
%                    -h_o    -h_o    -h_o    -h_o    -h_o  -h_o  -h_o  -h_o];  % Z Position
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%             
% Initialize to a User Guess:
mAct_sph = [  0   pi/2  pi  3*pi/2  0  0  0  0;    % Azimuth [rad]
             pi/6 pi/6 pi/6   pi/6  0  0  0  0];   % Inclination [rad]
pAct_cartesion = [ -0.25  0     0.25  0    -0.0675   0.0675   0.0675  -0.0675; % X Position
                    0    -0.25  0     0.25 -0.0675  -0.0675   0.0675   0.0675; % Y Position
                   -0.33 -0.33 -0.33 -0.33 -0.3   -0.3   -0.3   -0.3]; % Z Position
% pAct_cartesion = [ -0.15  0     0.15  0    -0.045   0.045   0.045  -0.045; % X Position
%                     0    -0.15  0     0.15 -0.045  -0.045   0.045   0.045; % Y Position
%                    -0.18 -0.18 -0.18 -0.18 -0.15   -0.15   -0.15   -0.15]; % Z Position

% All actuators vertical:
% mAct_sph = [ 0  0  0  0  0  0  0  0;     % Azimuth [rad]
%              0  0  0  0  0  0  0  0];    % Inclination [rad]
% r_a = large_EM(2)/2 + sqrt( (large_EM(2)/2+small_EM(2)/2)^2 + (large_EM(2)/2)^2 );
% pAct_cartesion = [  0.0675  0.0675 -0.0675 -0.0675  r_a   0    -r_a   0   ;  % X Position
%                     0.0675 -0.0675 -0.0675  0.0675  0     r_a   0    -r_a ;  % Y Position
%                    -0.15   -0.15   -0.15   -0.15   -0.15 -0.15 -0.15 -0.15]; % Z Position

% define x0
fittingParameters = [mAct_sph; pAct_cartesion];
Constants = [nomDist; large_EM; small_EM];

% The best result from the Monte Carlo sim will be recorded and returned at
% the end of the function call:
% Initialize a basis fitness consisting of all vertical Electromagnets
mAct_sph_best = mAct_sph;
pAct_cartesion_best = pAct_cartesion;
best_fitness = calcFitness(fittingParameters, Constants);
fprintf("Fitness to beat is: %d \n", best_fitness);

fitness_record = [ best_fitness ];

% % Plotting the baseline configuration:
% figure(figureNum)
% hold on
% plotEM(nomDist, large_EM, small_EM, rad2deg(mAct_sph), pAct_cartesion, best_fitness, figureNum);
% hold on
% [x y] = meshgrid(-0.75:0.1:0.75); % Generate x and y data
% z = -nomDist*ones(size(x, 1)); % Generate z data
% surf(x, y, z,'FaceAlpha',0.25) % Plot the surface
% hold off

%% ~~~~~~~~~~~~~~~~~~~ MONTE CARLO SIMULATION  ~~~~~~~~~~~~~~~~~~~~~~~~~ %%
numLoops = 10000;
numErrors = -numLoops;
fval = nan;
for i = 1:numLoops
    validConfig = false;
    % keep finding configurations and checking collisions until the
    % configuration is valid (no collisions)
    while (~validConfig)
        numErrors = numErrors + 1;
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
    % should have 1 new randomized collision-free configuration now.
    % find the fitness of this random set-up
    fittingParameters = [mAct_sph; pAct_cartesion];
    Constants = [nomDist; large_EM; small_EM];
    fval = calcFitness(fittingParameters,Constants);
    if ( fval < best_fitness )
        mAct_sph_best = mAct_sph;
        pAct_cartesion_best = pAct_cartesion;
        best_fitness = fval;
        fitness_record = [fitness_record, best_fitness];
    end
    if ( mod(i,100) == 0)
        fprintf("Loop %d complete \n", i)
    end
end 
fprintf("Number of Regenerated Configurations: %d \n", numErrors);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
mAct_sph = mAct_sph_best;
pAct_cartesion = pAct_cartesion_best;

figure(figureNum)
hold on
plot(fitness_record);
hold off
% plotEM(nomDist, large_EM, small_EM, rad2deg(mAct_sph), pAct_cartesion, best_fitness, figureNum+1);
% hold on
% [x y] = meshgrid(-0.75:0.1:0.75); % Generate x and y data
% z = -nomDist*ones(size(x, 1)); % Generate z data
% surf(x, y, z,'FaceAlpha',0.25) % Plot the surface
% hold off
end

% % Tethered Tool Parameters
% pTool = [0 0 0]'; %tool position [m]
% mTool_dir = [0 0 -1]'; %tool orientation in negative z axis
% mTool_mag = 8.4e-3; %tool dipole moment magnetiude [Am^2]
% mTool = mTool_dir/norm(mTool_dir)*mTool_mag; %tool dipole moment [Am^2]
% 
% 
% 
% fun_opt = @(x) calcFitness(x,Constants);
% fun_con = @(x) calcConstraint(x,Constants);
% % Maximize for Gradient and Field strength
% options = optimset('PlotFcns','optimplotfval','TolX',1e-9,'TolCon',1e-3,'MaxIter',10e4,'MaxFunEvals', 200e4) %'Algorithm','active-set'
% % [x,fval] = fminsearch(fun,fittingParameters,options)
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = [   0    0    0    0    0    0    0    0;
%          0    0    0    0    0    0    0    0;
%       -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5;
%       -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5;
%       -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5];
% ub = [ 2*pi 2*pi 2*pi 2*pi 2*pi 2*pi 2*pi 2*pi;  % 360  360  360  360  360  360  360  360;
%        pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2 pi/2;  % 90   90   90   90   90   90   90   90;
%        0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5;
%        0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5;
%       -nomDist -nomDist -nomDist -nomDist -nomDist -nomDist -nomDist -nomDist];
% % nonlcon = @fun_con;
% [x,fval] = fmincon(fun_opt,fittingParameters,A,b,Aeq,beq,lb,ub,fun_con,options)
% %https://www.mathworks.com/help/optim/ug/fmincon.html
% figureNum = figureNum+1;
% 
% mAct_sph = x(1:2,:);
% pAct = x(3:5,:);
% 
% % ~~~~ These are found already in the plotEM function ~~~~ %
% % % Define large cylinder for plotting
% % [X_l,Y_l,Z_l] = cylinder();
% % Z_l = -Z_l.*large_EM(3);
% % X_l =  X_l.*large_EM(2)/2;
% % Y_l =  Y_l.*large_EM(2)/2;
% % % Define small cylinder for plotting
% % [X_s,Y_s,Z_s] = cylinder();
% % Z_s = -Z_s.*small_EM(3);
% % X_s =  X_s.*small_EM(2)/2;
% % Y_s =  Y_s.*small_EM(2)/2;
% 
% plotEM(nomDist, large_EM, small_EM, rad2deg(mAct_sph), pAct, fval, figureNum+1);
% 
% end