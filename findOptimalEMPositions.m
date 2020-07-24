function [mAct_sph pAct] = findOptimalEMPositions(numActuators, large_EM, small_EM, nomDist, figureNum, quadrants)
% finds optimum conditions for actuators below nomDist plane

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
%%% pAct_sph = zeros(3,numActuators);
% Randomize the Initial inputs:
% mAct_sph = rand(2,numActuators);
% pAct_cartesion = rand(3,numActuators);
% Randomize the Initial inputs (with scaling):

% mAct_sph = [2*pi*rand(1,numActuators);
%              pi/2*rand(1,numActuators)];
% pAct_cartesion = [rand(1,numActuators)-0.5;
%                   rand(1,numActuators)-0.5;
%                   -(0.5-nomDist)*rand(1,numActuators)-0.15];

% Initial User Guess:
% mAct_sph = [  0   pi/2  pi  3*pi/2  0  0  0  0;    %Azimuth [rad]
%              pi/6 pi/6 pi/6   pi/6  0  0  0  0];    %Inclination [rad]
% pAct_cartesion = [ -0.15  0     0.15  0    -0.045   0.045   0.045  -0.045;  % X Position
%                     0    -0.15  0     0.15 -0.045  -0.045   0.045   0.045;	% Y Position
%                    -0.18 -0.18 -0.18 -0.18 -0.15   -0.15   -0.15   -0.15];% Z Position
% Initial User Guess (all vert):
mAct_sph = [ 0  0  0  0  0  0  0  0;     %Azimuth [rad]
             0  0  0  0  0  0  0  0];    %Inclination [rad]
r_a = large_EM(2)/2 + sqrt( (large_EM(2)/2+small_EM(2)/2)^2 + (large_EM(2)/2)^2 );
pAct_cartesion = [  0.0675  0.0675 -0.0675 -0.0675  r_a   0    -r_a   0   ; % X Position
                    0.0675 -0.0675 -0.0675  0.0675  0     r_a   0    -r_a ;	% Y Position
                   -0.15   -0.15   -0.15   -0.15   -0.15 -0.15 -0.15 -0.15];% Z Position

% Tethered Tool Parameters
pTool = [0 0 0]'; %tool position [m]
mTool_dir = [0 0 -1]'; %tool orientation
mTool_mag = 8.4e-3; %tool dipole moment magnetiude [Am^2]
mTool = mTool_dir/norm(mTool_dir)*mTool_mag; %tool dipole moment [Am^2]

% define x0
fittingParameters = [mAct_sph; pAct_cartesion];
Constants = [nomDist; large_EM; small_EM];

fun_opt = @(x) calcFitness(x,Constants);
fun_con = @(x) calcConstraint(x,Constants);
% Maximize for Gradient and Field strength
options = optimset('PlotFcns','optimplotfval','TolX',1e-9,'TolCon',1e-3,'MaxIter',10e4,'MaxFunEvals', 200e4) %'Algorithm','active-set'
% [x,fval] = fminsearch(fun,fittingParameters,options)
A = [];
b = [];
Aeq = [];
beq = [];
lb = [   0    0    0    0    0    0    0    0;
         0    0    0    0    0    0    0    0;
      -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5;
      -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5;
      -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5];
% Create the upper bound as close to the patient plane as possible for any
% EM orientation:
h_l = sqrt((large_EM(2)/2)^2+(large_EM(3)/2)^2);
h_s = sqrt((small_EM(2)/2)^2+(small_EM(3)/2)^2);
ub_z_l = -nomDist - h_l;
ub_z_s = -nomDist - h_s;
ub = [ 2*pi   2*pi   2*pi   2*pi   2*pi   2*pi   2*pi   2*pi;  % 360  360  360  360  360  360  360  360;
       pi/2   pi/2   pi/2   pi/2   pi/2   pi/2   pi/2   pi/2;  % 90   90   90   90   90   90   90   90;
       0.5    0.5    0.5    0.5    0.5    0.5    0.5    0.5;
       0.5    0.5    0.5    0.5    0.5    0.5    0.5    0.5;
       ub_z_l ub_z_l ub_z_l ub_z_l ub_z_s ub_z_s ub_z_s ub_z_s];
% nonlcon = @fun_con;
% fmincon with matrix inputs
% [x,fval] = fmincon(fun_opt,fittingParameters,A,b,Aeq,beq,lb,ub,fun_con,options);
%https://www.mathworks.com/help/optim/ug/fmincon.html
% fmincon with vectorized inputs
[x,fval] = fmincon(fun_opt,reshape(fittingParameters',1,[]),A,b,Aeq,beq,reshape(lb',1,[]),reshape(ub',1,[]),fun_con,options);
isVec = true;

% [x,fval] = ga(fun_opt,numel(fittingParameters),A,b,Aeq,beq,reshape(lb',1,[]),reshape(ub',1,[]),fun_con); %,options

figureNum = figureNum+1;

if isVec
    x = reshape(x,8,[])' ;
end

mAct_sph = x(1:2,:);
pAct = x(3:5,:);


% Define large cylinder for plotting
[X_l,Y_l,Z_l] = cylinder();
Z_l = -Z_l.*large_EM(3);
X_l = X_l.*large_EM(2)/2;
Y_l = Y_l.*large_EM(2)/2;
% Define small cylinder for plotting
[X_s,Y_s,Z_s] = cylinder();
Z_s = -Z_s.*small_EM(3);
X_s = X_s.*small_EM(2)/2;
Y_s = Y_s.*small_EM(2)/2;

plotEM(nomDist, large_EM, small_EM, rad2deg(mAct_sph), pAct, fval, figureNum+1);

end
% 
% %% Plotting 
% % Plot cylinders after rotating and translating
% figure(figureNum+1)
% hold on
% for l=1:4
%     % 3 x 3 Transformation (Rotation Only) Matrix 
%     RotM = rotz(deg2rad(mAct_sph(1,l)))*roty(deg2rad(mAct_sph(2,l)));
%     % Rotated cylinder points:
%     % | RotM   0  |     | X1      | X2
%     % |           |  x  | Y1  ... | Y2 
%     % |  0   RotM |     | Z1      | Z2
%     temp_cylinder = [RotM zeros(size(RotM));zeros(size(RotM)) RotM] * [X_l(1,:);Y_l(1,:);Z_l(1,:);X_l(2,:);Y_l(2,:);Z_l(2,:)];
%     X_mod = [temp_cylinder(1,:);temp_cylinder(4,:)] + pAct(1,l);
%     Y_mod = [temp_cylinder(2,:);temp_cylinder(5,:)] + pAct(2,l);
%     Z_mod = [temp_cylinder(3,:);temp_cylinder(6,:)] + pAct(3,l);
%     hold on
%     surf(X_mod,Y_mod,Z_mod) 
% end
% for l=5:8
%         % 3 x 3 Transformation (Rotation Only) Matrix 
%     RotM = rotz(deg2rad(mAct_sph(1,l)))*roty(deg2rad(mAct_sph(2,l)));
%     % Rotated cylinder points:
%     % | RotM   0  |     | X1      | X2
%     % |           |  x  | Y1  ... | Y2 
%     % |  0   RotM |     | Z1      | Z2
%     temp_cylinder = [RotM zeros(size(RotM));zeros(size(RotM)) RotM] * [X_s(1,:);Y_s(1,:);Z_s(1,:);X_s(2,:);Y_s(2,:);Z_s(2,:)];
%     X_mod = [temp_cylinder(1,:);temp_cylinder(4,:)] + pAct(1,l);
%     Y_mod = [temp_cylinder(2,:);temp_cylinder(5,:)] + pAct(2,l);
%     Z_mod = [temp_cylinder(3,:);temp_cylinder(6,:)] + pAct(3,l);
%     hold on
%     surf(X_mod,Y_mod,Z_mod) 
% end
% axis equal
% [X,Y,Z] = sphere();
% X = X.*nomDist;
% Y = Y.*nomDist;
% Z = Z.*nomDist;
% surf(X,Y,Z)
% xlabel("X")
% ylabel("Y")
% zlabel("Z")
% title("Plot actuator positions with a fitness of: " + num2str(fval))
% hold off
% figureNum = figureNum + 1;
% 
% end