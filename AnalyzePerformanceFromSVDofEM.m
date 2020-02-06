%% Evaluate different EM configurations for best performance using 
%               Singular Value Decomposition (SVD)
%
%  Author: Adam Schonewille
%  MASc Student, University of Toronto
%
%  Start Date: Feb 6th 2020
%  Last Editted: Feb 6th 2020
%

clear all; close all; clc; % ensure there are no variables or open plots

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');

%% ~~~~~~~~~~ ~~~~~~~~~~ GLOBAL CONSTANTS ~~~~~~~~~~ ~~~~~~~~~~ %%
% Essential:
mu_0 = pi*4e-7;         % [H/m or N/A^2 (amp-turns)]

% For Surgical applications:
% Head Dimensions taken from the 99th Percentile as an upper limit (max)
% https://en.wikipedia.org/wiki/Human_head
max_Head_Width  = 0.165; % [m] Ear to Ear
max_Head_Length = 0.217; % [m] Forehead (brows) to Back of Head
max_Head_Height = 0.255; % [m] Bottom of Chin to Top of Head
% Because actuators are from the bottom of the patient (in the supine
% postion) we will use the patient head length
r = 0.12; % [m] Distance from posterior of head to the center of the brain

%% ~~~~~~~~~~ ~~~~~~~~~~ OPERATING PARAMETERS ~~~~~~~~~~ ~~~~~~~~~~ %%
% Approximate the EM as a dipole source at its face
m_EM_large = 2*pi*(0.12)^3*(0.11407)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 63.8 mT
m_EM_small = 2*pi*(0.12)^3*(0.04846)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 18.9 mT
d_EM_large = 2*(45 + 22.5) / 1000;  % [m]
d_EM_small = 2*(30 + 15) / 1000;    % [m]
l_EM_large = 360 / 1000;            % [m]
l_EM_small = 240 / 1000;            % [m]
% Accumulate EM data
large_EM = [m_EM_large, d_EM_large, l_EM_large]';
small_EM = [m_EM_small, d_EM_small, l_EM_small]';

% Tethered Tool Parameters
pTool = [0 0 0]'; %tool position [m]
mTool_dir = [0 0 -1]'; %tool orientation
mTool_mag = 8.4e-3; %tool dipole moment magnetiude [Am^2] = [Nm/T]
mTool = mTool_dir/norm(mTool_dir)*mTool_mag; %tool dipole moment [Am^2]

%% ~~~~~~~~~~ ~~~~~~~~~~ SIMULATION TARGETS ~~~~~~~~~~ ~~~~~~~~~~ %%
targetMaxB_Field = 0.100; %[T] or higher
targetMaxForce   = 0.100; %[N] or higher
targetMaxField_G = targetMaxForce/mTool_mag; % [T/m]
% Create normalization matrix D0. This matrix with map unity inputs to unity
% desired outputs. Outputs greater than unity mean that fields greater than
% desired are possible. Inputs greater than unity higher currents are
% needed to produce a given desired output. Singular values less than unity
% should correspond to saturated input actuators.
DB = 1/targetMaxB_Field*eye(3);
DG = 1/targetMaxField_G*eye(5);
D0 = [DB, zeros(3,5); zeros(5,3), DG];

% for plotting
figureNum = 1;

%% ~~~~~~~~~~ ~~~~~~~~~~ SIMULATION INPUTS ~~~~~~~~~~ ~~~~~~~~~~ %%
mAct = [m_EM_large m_EM_large m_EM_large m_EM_large m_EM_small m_EM_small m_EM_small m_EM_small];
r_l = d_EM_large/2;
r_s = d_EM_small/2;
alpha = asin(sqrt(2)*r_s*sin(pi/4)/(r_l+r_s));
beta = pi-pi/4-alpha;
x_sqr = (r_l+r_s)*sin(beta)/sin(pi/4);
s = x_sqr/sqrt(2);
mAct_sph = [  pi/4 3*pi/4 5*pi/4 7*pi/4 0  0  0  0;    %Azimuth [rad]
             -pi/6  -pi/6  -pi/6  -pi/6 0  0  0  0];    %Inclination [rad]
pAct_cartesion = [    s     -s     -s      s  sqrt(2)*r_s         0     -sqrt(2)*r_s          0    ;    % X Position
                      s      s     -s     -s         0     sqrt(2)*r_s          0     -sqrt(2)*r_s ;	% Y Position
                   -0.12-r_l*sin(pi/6)  -0.12-r_l*sin(pi/6)  -0.12-r_l*sin(pi/6)  -0.12-r_l*sin(pi/6)     -0.12        -0.12         -0.12         -0.12 ];% Z Position

%% ~~~~~~~~~~ SIMULATION: EVALUATE SINGULAR VALUES ~~~~~~~~~~ %%              
nAct = size(pAct_cartesion,2);
RzyAct = zeros(3,1,nAct); 
% This is the unit vector of the magnetic moment for Euler angle rotations
% beta and gamma about Z and Y axes, respectively
% Assumes that the actuators dipole moment are alligned with Z axis
% originally
for i=1:nAct 
    RzyAct(:,:,i) = [ cosd(mAct_sph(1,i))*sind(mAct_sph(2,i));
                      sind(mAct_sph(1,i))*sind(mAct_sph(2,i));
                     -cosd(mAct_sph(2,i)) ];        
end
% Now we should have variables:
% pAct: xyz coordinates of actuator magnet centers (m)
% RzyAct: ZY Euler angle rotation matrix for each actuator magnet
% mAct: magnetic dipole moment of actuators [Am^2]
pTool = [0,0,0]';
% We want to test a variety of tool positions and assume that if it is
% controllable at extrema, it will also be controllable within extrema
a = 0.03; % bounding cube dist from origin
% Test points are the vertices of a 0.06 m sidelength cube as well as where
% the principle axis intersect the cube surface and origin (27 pts total)
pTool = [0  a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a ;
         0  0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a ;
         0  0  0  0  0  0  0  0  0  a  a  a  a  a  a  a  a  a -a -a -a -a -a -a -a -a -a ];
% The test points along the z axis corresponding to the mean, worst, and 
% best situations are at indices 1, 10, and 19, respectively

nOutputs = 8; % 3 for field and 5 for gradient
singularValues = zeros(1,size(pTool,2)); % For fitness evaluation (smaller sum is best)
K = zeros(nOutputs,nAct,size(pTool,2));    % Acuation matrix. Relates normalized currents to output field
U = zeros(nOutputs,nOutputs,size(pTool,2)); % m x m  Orthonormal basis of Km 
S = zeros(nOutputs,nAct,size(pTool,2));     % m x n  Singular Values (Scaling)
V = zeros(nAct,nAct,size(pTool,2));         % n x n  Orthonormal basis of Kn
% Iterate over all the different tested tool positions in the workspace
for j=1:size(pTool,2)
    Bpr = zeros(3,nAct);
    Gpr = zeros(5,nAct);
    % Iterate over all actuators and calculate their individual
    % contribution to the field at the workspace point of interest
    for i=1:size(pAct_cartesion,2)
        p = -pAct_cartesion(:,i) + pTool(:,j);    
        p_hat = p/norm(p);

        KB =   mAct(i)*1e-7 / (norm(p)^3); %note mu/(4*pi) = 1*e-7
        KG = 3*mAct(i)*1e-7 / (norm(p)^4);

    %     Bpr(:,2*i-1:2*i) = KB*( 3* (p_hat*p_hat') - eye(3) ) * RzyAct(:,:,i);
        Bpr(:,i) = KB*( 3* (p_hat*p_hat') - eye(3) ) * RzyAct(:,:,i);

        Gp = [
            3*p_hat(1)-5*p_hat(1)^3         p_hat(2)-5*p_hat(1)^2*p_hat(2)  p_hat(3)-5*p_hat(1)^2*p_hat(3);
            p_hat(2)-5*p_hat(1)^2*p_hat(2)  p_hat(1)-5*p_hat(1)*p_hat(2)^2  -5*p_hat(1)*p_hat(2)*p_hat(3);
            p_hat(3)-5*p_hat(1)^2*p_hat(3)  -5*p_hat(1)*p_hat(2)*p_hat(3)   p_hat(1)-5*p_hat(1)*p_hat(3)^2;
            p_hat(1)-5*p_hat(1)*p_hat(2)^2  3*p_hat(2)-5*p_hat(2)^3         p_hat(3)-5*p_hat(2)^2*p_hat(3);
            -5*p_hat(1)*p_hat(2)*p_hat(3)   p_hat(3)-5*p_hat(2)^2*p_hat(3)  p_hat(2)-5*p_hat(2)*p_hat(3)^2;
            ];

    %     Gpr(:,2*i-1:2*i) = KG * Gp * RzyAct(:,:,i);
        Gpr(:,i) = KG * Gp * RzyAct(:,:,i);
    end

    K(:,:,j) = [Bpr;Gpr];
    K0 = D0*K(:,:,j);
    [U(:,:,j),S(:,:,j),V(:,:,j)] = svd( K0 );
%     [U,S,V] = svd(Bpr);
%     B = sum(Bpr,2)
%     G = sum(Gpr,2);
    max3values = zeros(2,3);
    for m=1:size(U,2)
        magnitude = norm(U(1:3,m,j));
        if (magnitude > max3values(1,1) )
            if (magnitude > max3values(1,2) )
                if (magnitude > max3values(1,3) )                    
                    % shift all over
                    max3values(1,1) = max3values(1,2);
                    max3values(1,2) = max3values(1,3);
                    max3values(1,3) = magnitude;
                    % shift indices too
                    max3values(2,1) = max3values(2,2);
                    max3values(2,2) = max3values(2,3);
                    max3values(2,3) = m;
                else
                    % shift all but 3
                    max3values(1,1) = max3values(1,2);
                    max3values(1,2) = magnitude;
                    % shift indices too
                    max3values(2,1) = max3values(2,2);
                    max3values(2,2) = m;
                end
            else
                max3values(1,1) = magnitude;
                max3values(2,1) = m;
            end
        end
                    
    end
%     U(:,:,j)
    maxSV = S(max3values(2,3),max3values(2,3));
    minSV = S(max3values(2,1),max3values(2,1));
    singularValues(j) = 1/minSV;% + maxSV/minSV;
    printout = sprintf('Point %d complete',j);
    disp(printout);
end

plotEM(r, large_EM, small_EM, rad2deg(mAct_sph), pAct_cartesion, sum(singularValues), figureNum);
figureNum = figureNum + 1;