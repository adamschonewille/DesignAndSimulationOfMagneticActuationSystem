%% Single dipole vs Multiple dipoles to Approximate EM Fields
%
%  Author: Adam Schonewille
%  MASc Student, University of Toronto
%
%  Start Date: Feb 11th 2020
%  Last Editted: Feb 11th 2020
%

clear all; close all; clc; % ensure there are no variables or open plots

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');

%% ~~~~~~~~~~ ~~~~~~~~~~ GLOBAL CONSTANTS ~~~~~~~~~~ ~~~~~~~~~~ %%
% Essential:
mu_0 = pi*4e-7;         % [H/m or N/A^2 (amp-turns)]

% Distance used in ANSYS simulations
r = 0.12; % [m] Distance from posterior of head to the center of the brain

%% ~~~~~~~~~~ ~~~~~~~~~~ OPERATING PARAMETERS ~~~~~~~~~~ ~~~~~~~~~~ %%
% Approximate the EM as a single dipole source at its face
B_target_l = 0.11407; % [T] for J = 24 A/mm^2. At J = 6 A/mm^2 -> 63.8 mT
B_target_s = 0.04846; % [T] for J = 24 A/mm^2. At J = 6 A/mm^2 -> 18.9 mT
m_EM_large = 2*pi*(0.12)^3*(0.11407)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 63.8 mT
m_EM_small = 2*pi*(0.12)^3*(0.04846)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 18.9 mT

d_EM_large = 2*(45 + 22.5) / 1000;  % [m]
d_EM_small = 2*(30 + 15) / 1000;    % [m]
r_l = d_EM_large/2;
r_s = d_EM_small/2;
l_EM_large = 360 / 1000;            % [m]
l_EM_small = 240 / 1000;            % [m]
% Accumulate EM data
large_EM = [m_EM_large, d_EM_large, l_EM_large]';
small_EM = [m_EM_small, d_EM_small, l_EM_small]';

% for plotting
figureNum = 1;

%% ~~~~~~~~~~ ~~~~~~~~~~ SIMULATION INPUTS ~~~~~~~~~~ ~~~~~~~~~~ %%
nAzimuthal = 4:4:120;
nAxial = 4:2:62;
nRadial = 2:1:31;
nDipoles = nAzimuthal.*nRadial.*nAxial;
nSimulations = length(nDipoles);
m = zeros(size(nDipoles));

%% Large EM Simulation
dx = 0.001;
x = dx:dx:0.2; % Points to evaluate field at
% Create a 3 x 1 x m x n matrix 
% Where:  3 x 1 is the Bx, By, Bz at each mth point for the nth simulation
B = zeros(3,1,length(x),nSimulations);

val_index = 0;
tol = 1e-9;
for i=1:length(x)
    if ( abs(x(i)-r) < tol )
        % Found index with critical value in it
        val_index = i;
    end
    if ( (i == length(x)) && (val_index == 0) )
        % Could not find value
        disp('Error: the vector "x" does not contain the critical search position');
        return;
    end
end

% Find the magnetic moment value such that all dipoles have the same m 
% and the field strength at 'r' matches the ANSYS Simulations
for n = 1:nSimulations
    % Initialize the positions of the dipoles at the start of each
    % simulation
    pDipoles = zeros(3,nDipoles(n));    
    % Find slice locations for the dipoles to be located at along the EM
    % axis
    slicesL = linspace(0,-l_EM_large,nAxial(n));
    slicesR = linspace(0,r_l,nRadial(n));
    dphi = 2*pi / nAzimuthal(n);
    index = 1;
    % Discretize dipole positions and save into vector:
    for i = 1:nAxial(n)
%         % initialize first point in center
%         pDipoles(:,index) = [slicesL(i) 0 0]';
%         index = index + 1;
        for j = 1:nAzimuthal(n)
            for k = 1:nRadial(n)
                pDipoles(:,index) = [slicesL(i) slicesR(k)*sin(j*dphi) slicesR(k)*cos(j*dphi)]';
                index = index + 1;
            end
        end    
    end
    % All positions should now be encoded/calculated and stored in pDipoles
    % For all dipoles we should calculate the B-field at r with |m| = 1
    % Assume all small dipoles, dm, have the same magnetic moment, m, and
    % that they are all parallel and in the axial direction of the EM
    mdir = [1 0 0]';
    pTool = [r 0 0]';
    % B_norm is the normalized B field from the dipole contributions with
    % |m| = 0
    B_norm = [0 0 0]';
    for i = 1:size(pDipoles,2)
        B_norm = B_norm + dipoleField([x(val_index) 0 0]'-pDipoles(:,i), mdir);
    end
    % find magnetic moment of this simulation:
    % Currently |m| = 1  --> B_norm
    % so B_norm * m = B_target
    m(n) = B_target_l/B_norm(1);
    % print out B_norm, the y and z values should be 0
    B_norm
    % While we already have the m and pDipoles calculated, we should
    % calculate the fields at some points along the axis:
    for j = 1:length(x)
        % at each j points calculate through all dipoles
        for i = 1:size(pDipoles,2)
            % zeros(3,1,length(x),nSimulations);
            B(:,:,j,n) = B(:,:,j,n) + dipoleField([x(j) 0 0]'-pDipoles(:,i), m(n)*mdir);
        end
    end
    printout = sprintf('Simulation %d complete',n);
    disp(printout);
end

%% Time to plot yo
% Plot yo numbas

figure(1)
hold on
for i=1:nSimulations
    Bx = zeros(size(x));
    for j = 1:length(x)
        Bx(j) = B(1,1,j,i);
    end
    plot(x,Bx,'Color', [1 1-nDipoles(i)/nDipoles(end) 0])
end
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',...
        '16','17','18','19','20','21','22','23','24','25','26','27','28',...
        '29','30','NorthWest')
title('On-Axis Magnetic field for different numbers of dipoles in an EM')
xlabel('x [m]')
ylabel('Bx [T]')
hold off
% Plot difference between simulations
return;
figure(2)
hold on
for i=2:nSimulations
    Bx = zeros(size(x));
    for j = 1:length(x)
        Bx(j) = B(1,1,j,i);
    end
    Bx_prev = zeros(size(x));
    for j = 1:length(x)
        Bx_prev(j) = B(1,1,j,i-1);
    end
    plot(x,Bx,'Color', [1 1-nDipoles(i)/nDipoles(end) 0])
end
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',...
        '16','17','18','19','20','21','22','23','24','25','26','27','28',...
        '29','30','NorthWest')
title('On-Axis Magnetic field for different numbers of dipoles in an EM')
xlabel('x [m]')
ylabel('Bx [T]')
hold off


return;

mAct = [m_EM_large m_EM_large m_EM_large m_EM_large 0 0 0 0];
r_l = d_EM_large/2;
r_s = d_EM_small/2;
alpha = asin(sqrt(2)*r_l*sin(pi/4)/(r_l+r_s));
beta = pi-pi/4-alpha;
x_sqr = (r_l+r_s)*sin(beta)/sin(pi/4);
s = x_sqr/sqrt(2);
mAct_sph = [  0     pi/2  pi   3*pi/2 0  0  0  0;    %Azimuth [rad]
             -pi/6 -pi/6 -pi/6 -pi/6  0  0  0  0];    %Inclination [rad]
pAct_cartesion = [ sqrt(2)*r_l        0    -sqrt(2)*r_l          0     s    -s    -s     s  ;  % X Position
                          0    sqrt(2)*r_l         0     -sqrt(2)*r_l  s     s    -s    -s  ;	% Y Position
                       -0.15       -0.15        -0.15         -0.15  -0.12 -0.12 -0.12 -0.12 ];% Z Position

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



