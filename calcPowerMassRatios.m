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
maxCurrent = 20; % [A]

% For Permanent Magnets 
avgRemanence = 1.45; % [T]
rho_NdFeB = 7500; % [kg/m^3] Density of Neodymium Magnet

% outputArg1 = inputArg1;
% outputArg2 = inputArg2;

numPoints = 10001;
mass = linspace(0.01, 10, numPoints); % From 10 grams to 10 kg
B_PermanentMagnet = zeros(length(mass),1);
m = zeros(length(mass),1); % Debugging purposes
B_Electromagnet = zeros(length(mass),1);

%% Permanent Magnet
for i=1:numPoints
    volume = mass(i)/rho_NdFeB; % [m^3]
    sideLength = volume^(1/3);  % [m]
    m(i) = volume*avgRemanence/mu_0;
    B_PermanentMagnet(i) = mu_0*m(i)/(2*pi*(Z+sqrt(2)*sideLength/2)^3); % [T]
end

figure(figureNum)
plot(mass,B_PermanentMagnet*1000)
hold on 

xlabel("Mass [kg]")
ylabel("On-axis Magnetic Field [mT] at a distance of " + num2str(Z) + " [m]" )
title("Magnetic Field vs Mass of Actuator")
legend("Permanent Magnet", "Location","SouthEast")
hold off


end

