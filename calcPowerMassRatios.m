function [] = calcPowerMassRatios(figureNum)
% AUTHOR:   Adam Schonewille
% CREATED:     July 5th 2019
% LAST EDITED: July 15th 2019
% ABOUT:    Calculates and plots the actuator output power vs mass
%
%% INPUTS
% 
%
%% RETURNED OUTPUTS
% 

%% Function Code

mu_0 = pi*4e-07;
Z = 0.120; % [m] Distance from closest actuator face to the center of the workspace

% For Electromagmets:
% maxCurrent = 20; % [A]
currentDensity = 3; % [A/mm^2]
Power = 1800; % [W]
rho_Cu = 8960; %[kg/m^3] Density of Copper
rho_Supermendur = 8120; %[kg/m^3] Density of Core material
packing_factor = pi/4;

% For Permanent Magnets 
avgRemanence = 1.45; % [T]
avgCoercivity = 1.25e6; % [A/m]
rho_NdFeB = 7500; % [kg/m^3] Density of Neodymium Magnet (7300-7700)

% outputArg1 = inputArg1;
% outputArg2 = inputArg2;

numPoints = 10001;
mass = linspace(0.01, 100, numPoints); % From 10 grams to 10 kg
B_PermanentMagnet = zeros(length(mass),1);
m = zeros(length(mass),1); % Debugging purposes

%% Permanent Magnet
for i=1:numPoints
    volume = mass(i)/rho_NdFeB; % [m^3]
    sideLength = volume^(1/3);  % [m]
    m(i) = volume*avgRemanence/mu_0;
    B_PermanentMagnet(i) = mu_0*m(i)/(2*pi*(Z+sqrt(2)*sideLength/2)^3); % [T]
end

mass_discrete = 25.4/1000*(0.5:0.5:9); % [m]
mass_discrete = mass_discrete.^3 * rho_NdFeB; % [kg]

B_Magnet_Sim = [0.000221760;
                0.001444600;
                0.004118200;
                0.008095800;
                0.013009000;
                0.019897000;
                0.026141000;
                0.032443000;
                0.039003000;
                0.049828000;
                0.056226000;
                0.063676000;
                0.072515000;
                0.083478000;
                0.085689000;
                0.096825000;
                0.102670000;
                0.111460000];


markerSize = 5;
lineWidth = 2;
            
figure(figureNum)
plot(mass,B_PermanentMagnet*1000,"MarkerSize", markerSize,"LineWidth",lineWidth)
hold on 

plot(mass_discrete,B_Magnet_Sim*1000,'xb',"MarkerSize", markerSize,"LineWidth",lineWidth)

xlabel("Mass [kg]")
ylabel("On-axis Magnetic Field [mT] at a distance of " + num2str(Z) + " [m]" )
title("On-axis Magnetic Field at "+ num2str(Z) + " [m] vs Mass of Actuators (Q = 8 for Electromagnets)")
% legend("Permanent Magnet", "Location","SouthEast")
% hold off
% figureNum = figureNum +1;

%% Electromagnet from Experiments
% All from a specific experiment with Q = 8 for all cases
% Oct 24th 2019
L =      1/1000 * (80:40:480); % [m] The radius of the electromagnet core
R_core = 1/1000 * (10:5:60);   % [m] The length of the electromagnet
r =      1/1000 * (5:2.5:30);  % [m] The thickness of the coil surrounding the core
B_Electromagnet_J_3 = [0.000318219725;
                     0.001085626073;
                     0.006689077855;
                     0.011647687264;
                     0.017939724572;
                     0.023171295615;
                     0.022473120811;
                     0.030596100000;
                     0.049745500000;
                     0.050088900000;
                     0.060794600000];
B_Electromagnet_J_6 = [ 0.000526720000;   %0.000667750000;
                        0.002149450000;   %0.002030000000;
                        0.004811550000;   %0.005000200000;
                        0.009445750000;   %0.012819000000;
                        0.018949000000;   %0.016811000000;
                        0.032585750000;
                        0.039335500000;
                        0.063814000000;
                        0.089053000000;
                        0.096940666667;
                        0.118056666667 ];
B_Electromagnet_J_24 = [ 0.001918833;
                         0.008242367;
                         0.014959333;
                         0.022812667;
                         0.048461000;
                         0.062727667;
                         0.075294333;
                         0.114074667;
                         0.140411667;
                         0.158966667;
                         0.203443333 ];
                 
B_Electromagnet_P = [0.004302180428;
                     0.009229306070;
                     0.018202944721;
                     0.025955937284;
                     0.031707106895;
                     0.043617000000;  % 0.071116433753; After changing mesh
                     0.058249931377;
                     0.068153800000;
                     0.072124000000;  % 0.095997800000; After changing mesh
                     0.084183600000;
                     0.089674700000];

mass_electromagnet = ones(length(L),1);
                 
for i = 1:length(L)
    mass_electromagnet(i) = L(i)*pi*R_core(i)^2*rho_Supermendur + ...
        packing_factor*L(i)*pi*((r(i)+R_core(i))^2-R_core(i)^2)*rho_Cu;
end

% figure(figureNum)
hold on
plot(mass_electromagnet,B_Electromagnet_J_3*1000,'ok',"MarkerSize", markerSize,"LineWidth",lineWidth);
plot(mass_electromagnet,B_Electromagnet_J_6*1000,'xk',"MarkerSize", markerSize,"LineWidth",lineWidth); 
plot(mass_electromagnet,B_Electromagnet_J_24*1000,'^k',"MarkerSize", markerSize,"LineWidth",lineWidth);
plot(mass_electromagnet,B_Electromagnet_P*1000,'or',"MarkerSize", markerSize,"LineWidth",lineWidth);
% xlabel("Mass [kg]")
% ylabel("On-axis Magnetic Field [mT] at a distance of " + num2str(Z) + " [m]" )
% title("Magnetic Field vs Mass of Actuator")
% legend("Const Current Density 3 A/mm$^2$","Const Power of 1.8 kW", "Location","SouthEast")

% Plot dipole approximation breakdown:
thresholdMass = (Z/2)^3*rho_NdFeB
plot([thresholdMass, thresholdMass],[0 205],"-.","MarkerSize", markerSize,"LineWidth",lineWidth)

legend("Permanent Magnet Dipole (Theory)","Permanent Magnet Simulations","Const. Current Density 3 A/mm$^2$","Const. Current Density 6 A/mm$^2$", "Const. Current Density 24 A/mm$^2$","Const. Power of 1.8 kW","Threshold for dipole model approx.", "Location","SouthEast")
% xlim([0 12.5])
ylim([0 205])



hold off

end

