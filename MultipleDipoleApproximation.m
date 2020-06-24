%% Single dipole vs Multiple dipoles to Approximate EM Fields
%
%  Author: Adam Schonewille
%  MASc Student, University of Toronto
%
%  Start Date: Feb 11th 2020
%  Last Editted: Mar 6th 2020
%
% This method is shown to be inaccurate compared to ANSYS simulations due
% to the magnetization of the core being nonuniform. This script assumes
% that the magnetized core acts like a permanent magnet with uniform
% magnetization.

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

t_l =  22.5 / 1000;         % [m] the thickness of the wire layers wrapping the coil
r_core_l =   45 / 1000;     % [m] Inner Coil radius
d_l = 2*(r_core_l + t_l);   % [m] Outer Coil diameter
% L_l = 360 / 1000;           % [m] Length of the large coil
L_l = r_core_l*2;
J = 24;                     % [A/mm^2]
packing_factor = pi/4;
A_l = t_l * L_l * packing_factor;   % [m^2]
i_l = J * A_l * 1e6;                % [A]
% Using the equation from Erics paper for setting the theoretical magnetization
% m_l = pi/3 * (3*r_core_l^2+3*r_core_l*t_l+t_l^2) * i_l;
% dipoleField([r-(-L_l/2) 0 0]', [m_l 0 0]')
% dipoleField([r 0 0]', [m_l 0 0]')
% m_l = pi*r_core_l^2 * i_l
% dipoleField([r-(-L_l/2) 0 0]', [m_l 0 0]')    
% dipoleField([r 0 0]', [m_l 0 0]')

% for near fields to match the ANSYS model, use this:
% m_l = 0.1750 * 2*pi*(0.12+L_l/2)^3*(0.11407)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 63.8 mT
% for far fields to match the ANSYS model, use this:
% m_l = 0.4555 * 2*pi*(0.12+L_l/2)^3*(0.11407)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 63.8 mT

% do m_l for dipole at center EM with B = 0.1 T at r=0.12

% from ANSYS:
% m_l = 2*pi*(0.12+L_l/2)^3*(0.0299)/mu_0;

% from COMSOL:
m_l = 2*pi*(0.12+L_l/2)^3*(0.0227)/mu_0;

% V_l = pi*r_core_l^2*L_l; % [m^3] volume of the magnetized core 
V_l = pi*(r_core_l)^2*L_l; % [m^3] volume of the whole magnetized EM
M_l = m_l/V_l; % [A/m]

%% NOTE THIS m_l is way too small. Field value at ANSYS is roughly 10% 

% d_EM_small = 2*(30 + 15) / 1000;    % [m]
% r_l = d_l/2;
% r_s = d_EM_small/2;
% l_EM_small = 240 / 1000;            % [m]
% Accumulate EM data
large_EM = [m_EM_large, d_l, L_l]';
% small_EM = [m_EM_small, d_EM_small, l_EM_small]';



%% ~~~~~~~~~~ ~~~~~~~~~~ SIMULATION INPUTS ~~~~~~~~~~ ~~~~~~~~~~ %%
% nAzimuthal = 12:12:360;
nAzimuthal = 10:6:184;
nAxial = 4:2:62;
nRadial = 2:1:31;
nDipoles = nAzimuthal.*nRadial.*nAxial;
nSimulations = length(nDipoles);
m = zeros(size(nDipoles)); % not used anymore?

largest_s   = 1/nAzimuthal(end) * 2 * pi * r_core_l;
largest_l = 1/nAxial(end) * L_l;
largest_r = 1/nRadial(end) * r_core_l;

%% Large EM Simulation
dx = 0.001;
x = dx:dx:0.2; % Points to evaluate field at
% Create a 3 x 1 x m x n matrix 
% Where:  3 x 1 is the Bx, By, Bz at each mth point for the nth simulation
B = zeros(3,1,length(x),nSimulations);


%% Check that the point exists in the array for calibrating the overall m
val_index = 0;
tol = 1e-9;
% This section is no longer relevant
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

%% Simulation for dividing up the volume into multiple dipoles further and further until convergence
choice = menu("Run Simulation: Division of volume into multiple dipoles until convergence?",'Yes','No');
 % 0 is Exit
 % 1 is Yes (option 1)
 % 2 is No (option 2)
if (choice == 1)

    % Find the magnetic moment value such that all dipoles have the same m 
    % and the field strength at 'r' matches the ANSYS Simulations
    tic
    for n = 1:nSimulations
        % Initialize the positions of the dipoles at the start of each
        % simulation
        pDipoles = zeros(3,nDipoles(n)); % [m] Positions     
        mDipoles = zeros(1,nDipoles(n)); % [Am^2] Magnetic moments  
        % Find slice locations for the dipoles to be located at along the EM
        % axis
        slicesL = linspace(0,-L_l,nAxial(n)+1);
        dl = abs(slicesL(2)-slicesL(1));
    %     slicesR = linspace(0,r_l,nRadial(n)+1);
%         slicesR = linspace(0,r_core_l,nRadial(n)+1); % To core radius since the core acts as a permanent magnet
        slicesR = linspace(0,r_core_l,nRadial(n)+1); % To EM radius
        dr = abs(slicesR(2)-slicesR(1));
        dphi = 2*pi / nAzimuthal(n); % nAzimuthal(n) does not need a + 1 (circular divisions)
        index = 1;
        % Discretize dipole positions and save into vector:
        for i = 1:nAxial(n)
    %         % initialize first point in center
    %         pDipoles(:,index) = [slicesL(i) 0 0]';
    %         index = index + 1;
            for j = 1:nAzimuthal(n)
                for k = 1:nRadial(n)
                    % Set pDipoles as the centroids of the volumes they
                    % represent:
                    Cr = 2/dphi*sin(dphi/2)*(2*slicesR(k)^2 + 2*slicesR(k)*dr + 2/3*dr^2)/(2*slicesR(k)+dr);
                    p_x = slicesL(i) - dl/2;
                    p_y = Cr * sin(j*dphi); % Using dphi straight as it is just a symmetric rotation about the EM axis to add dphi/2 to it
                    p_z = Cr * cos(j*dphi);
                    pDipoles(:,index) = [p_x p_y p_z]';
                    % Calc volume
                    A_tmp = ( 2*slicesR(k)+dr ) * dr*dphi/2;
                    V_tmp = A_tmp * dl;
                    % Calc m from M and V
                    mDipoles(:,index) = V_tmp * M_l;
                    index = index + 1;
                end
            end         
        end
        % All positions should now be encoded/calculated and stored in pDipoles
        % For all dipoles we should calculate the B-field at r with |m| = 1
        % Assume all small dipoles, dm, have the same magnetic moment, m, and
        % that they are all parallel and in the axial direction of the EM
        mdir = [1 0 0]';
        pTool = [x(val_index) 0 0]';

    %     % B_norm is the normalized B field from the dipole contributions with
    %     % |m| = 0
    %     B_norm = [0 0 0]';
    %     for i = 1:size(pDipoles,2)
    %         B_norm = B_norm + dipoleField(pTool-pDipoles(:,i), mdir.*mDipoles(1,i));
    %     end
    %     % find magnetic moment of this simulation:
    %     % Currently |m| = 1  --> B_norm
    %     % so B_norm * m = B_target
    %     m(n) = B_target_l/B_norm(1);
    %     % print out B_norm, the y and z values should be 0
    %     B_norm


        % While we already have the m and pDipoles calculated, we should
        % calculate the fields at some points along the axis:
        for j = 1:length(x)
            % at each j points calculate through all dipoles
            for i = 1:size(pDipoles,2)
                % zeros(3,1,length(x),nSimulations);
    %             B(:,:,j,n) = B(:,:,j,n) + dipoleField([x(j) 0 0]'-pDipoles(:,i), m(n)*mdir);
                pTest = [x(j) 0 0]';
                B(:,:,j,n) = B(:,:,j,n) + dipoleField(pTest-pDipoles(:,i), mDipoles(i)*mdir);
            end
        end
        printout = sprintf('Simulation %d complete',n);
        disp(printout);
    end
    toc

    %% Time to plot yo
    % Plot yo numbas

    % for plotting
    figureNum = 1;
    lw = 1.5;
    ms = 1.5;
    figure(figureNum)
    hold on
    Bx = zeros(size(x));
    for i=1:nSimulations
        Bx = zeros(size(x));
        for j = 1:length(x)
            Bx(j) = B(1,1,j,i);
        end
        plot(x,Bx,'Color', [1 1-nDipoles(i)/nDipoles(end) 0])
    end

    title('On-Axis Magnetic field for different numbers of dipoles in an EM')
    xlabel('x [m]')
    ylabel('Bx [T]')
    B_init = zeros(3,length(x));
    for i=1:length(x)
        B_init(:,i) = dipoleField([x(i)-(-L_l/2) 0 0]', [m_l 0 0]');
    end
    plot(x,B_init(1,:),'k','LineWidth',lw,'MarkerSize',ms)
    B_fit = zeros(3,length(x));
%     for i=1:length(x)
%         mm = 2*pi*(0.12)^3*(0.11407)/mu_0;
%         B_fit(:,i) = dipoleField([x(i) 0 0]', [mm 0 0]');
%     end
%     plot(x,B_fit(1,:),'b','LineWidth',lw,'MarkerSize',ms)
    legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',...
            '16','17','18','19','20','21','22','23','24','25','26','27','28',...
            '29','30','Dipole EM center, m = 7014','NorthWest')
    hold off
    figureNum = figureNum + 1;
    % Plot difference between simulations

end
%%
lw = 1.5;
ms = 1.5;
    
x_sims = [0.001; 0.002; 0.003; 0.004; 0.005; 0.006; 0.007; 0.008; ...
          0.009; 0.010; 0.011; 0.012; 0.013; 0.014; 0.015; 0.016; ...
          0.017; 0.018; 0.019; 0.020; 0.021; 0.022; 0.023; 0.024; ...
          0.025; 0.026; 0.027; 0.028; 0.029; 0.030; 0.040; 0.050; ...
          0.060; 0.070; 0.080; 0.090; 0.100; 0.110; 0.120; 0.130; ...
          0.140; 0.150; 0.160; 0.170; 0.180; 0.190; 0.200; 0.210; ...
          0.220; 0.230; 0.240; 0.250];

B_coils_ANSYS = [0.2363; 0.2337; 0.2310; 0.2284; 0.2257; 0.2231; ...
                 0.2204; 0.2177; 0.2151; 0.2124; 0.2094; 0.2062; ...
                 0.2030; 0.1998; 0.1963; 0.1926; 0.1881; 0.1836; ...
                 0.1792; 0.1748; 0.1706; 0.1664; 0.1622; 0.1580; ...
                 0.1538; 0.1497; 0.1456; 0.1414; 0.1373; 0.1331; ...
                 0.1013; 0.0809; 0.0627; 0.0505; 0.0409; 0.0336; ...
                 0.0272; 0.0224; 0.0189; 0.0162; 0.0137; 0.0118; ...
                 0.0105; 0.0092; 0.0082; 0.0075; 0.0068; 0.0063; ...
                 0.0059; 0.0057; 0.0055; 0.0055]; 

B_CoredCoils_ANSYS = [0.4463; 0.4447; 0.4430; 0.4414; 0.4398; 0.4382; ...
                      0.4365; 0.4349; 0.4333; 0.4316; 0.4293; 0.4264; ...
                      0.4236; 0.4207; 0.4169; 0.4120; 0.4054; 0.3987; ...
                      0.3921; 0.3857; 0.3789; 0.3715; 0.3641; 0.3567; ...
                      0.3494; 0.3423; 0.3351; 0.3279; 0.3207; 0.3134; ...
                      0.2480; 0.2014; 0.1582; 0.1282; 0.1044; 0.0863; ...
                      0.0700; 0.0577; 0.0488; 0.0418; 0.0354; 0.0307; ...
                      0.0272; 0.0238; 0.0214; 0.0194; 0.0177; 0.0163; ...
                      0.0154; 0.0148; 0.0144; 0.0142];

                  
                  
% i_l = i_l * pi/4;
% i_l = 33929.20066; % was 38170

B_funct = @(z) mu_0*i_l/(L_l*2)*(cos(atan2(r_core_l,z+L_l)) - cos(atan2(r_core_l,z)));

figure(5)
hold on
% Look at coil contribution
fplot(B_funct,[-L_l 0.25],'LineWidth',lw,'MarkerSize',ms,'Color',[235/255, 113/255, 52/255])
plot(x_sims,B_coils_ANSYS,'LineWidth',lw,'MarkerSize',ms,'Color',[186/255, 7/255, 7/255])
% Look at core contribution
plot(x,Bx,'LineWidth',lw,'MarkerSize',ms,'Color',[6/255, 127/255, 184/255])
plot(x_sims,B_CoredCoils_ANSYS-B_coils_ANSYS,'LineWidth',lw,'MarkerSize',ms,'Color',[6/255, 18/255, 184/255])
% Plot original dipole approx
plot(x,B_init(1,:),'LineWidth',lw,'MarkerSize',ms,'Color',[0/255, 125/255, 23/255])

title("On-axis Magnetic Field Contributions from Coil and Core")
xlabel("x Distance [m]")
ylabel("Bx [T]")
xlim([0 0.2])
ylim([0 0.6])
legend("Modeled Analytical Coils","ANSYS Coils","Multiple Dipoles Core","ANSYS Core","Single dipole")
hold off



%%
return;

choice = menu("Run Simulation: Determine the field for the converged configuration?",'Yes','No');
 % 0 is Exit
 % 1 is Yes (option 1)
 % 2 is No (option 2)
if (choice == 1)
    % Define the B to H relationship within the core material:
    % Obtained from Supermendure on ANSYS workbench
    % Field Intensity  Flux Density
    %       [A/m]   [T]
    B_H = [ 0       0
            22.2	1.8
            24.5	2
            31.7	2.055
            47.6	2.085
            63.8	2.11
            79.4	2.13
            158.7	2.19
            317.5	2.236
            476.2	2.27
            634.9	2.28
            793.7	2.3
            1587	2.32
            3174	2.35
            4762	2.375
            6348	2.39
            7937	2.395
            8730	2.4
            17460	2.4274
            34920	2.4538
            69840	2.5459
            139680	2.69
            279360	2.91
            558720	3.29]
    
    
end

return;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(figureNum)
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
figureNum = figureNum + 1;

return;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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



