%% Simulations for designing the next novel actuation system.
%
%  Author: Adam Schonewille
%  MASc Student, University of Toronto
%
%  Start Date: July 5th 2019
%  Last Editted: Feb 6th 2020
%

clear all; close all; clc; % ensure there are no variables or open plots

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');


%% ~~~~~~~~~~ ~~~~~~~~~~ GLOBAL CONSTANTS ~~~~~~~~~~ ~~~~~~~~~~ %%
% Essential:
mu_0 = pi*4e-7;         % [H/m or N/A^2 (amp-turns)]
epsilon_0 = 8.854e-12;  % [F/m or N/V^2]
rel_Temp = 293;         % [K] Standard temperature

% Non-essential:
% maxTemperature_Cu = 198.889; % [Degrees C], MP of Copper
% mu_rel_Cu = 0.999994;   % [ ] 
alpha_Cu = 0.00404;     % [1/K] Temperature coefficient
k_Cu = 401;             % [W/(m*K)] Thermal conductivity @ 293 K
rho_Cu = 1.68e-8;       % [Ohm*m] Resistivity @ 293 K
N2_temp = 77;           % [K] (Liquid Nitrogen)

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
magnetSize = 0.0254; % [m]
Remenance = 1.45; % [T]

m_EM_large = 2*pi*(0.12)^3*(0.11407)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 63.8 mT
m_EM_small = 2*pi*(0.12)^3*(0.04846)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 18.9 mT
d_EM_large = 2*(45 + 22.5) / 1000;  % [m]
l_EM_large = 360 / 1000;            % [m]
d_EM_small = 2*(30 + 15) / 1000;    % [m]
l_EM_small = 240 / 1000;            % [m]

large_EM = [m_EM_large, d_EM_large, l_EM_large]';
small_EM = [m_EM_small, d_EM_small, l_EM_small]';

%% ~~~~~~~~~~ ~~~~~~~~~~ SIMULATION TARGETS ~~~~~~~~~~ ~~~~~~~~~~ %%
isotropicFields = true;
MaxB_Field = 0.040; %[T] or higher

% for plotting
figureNum = 1;

%% ~~~~~~~~~~ ~~~~~~~~~~ POWER/WEIGHT SIMULATION ~~~~~~~~~~ ~~~~~~~~~~ %%
choice = menu("Run Simulation for: Finding the Power/Mass ratio of Actuators?",'Yes','No');
 % 0 is Exit
 % 1 is Yes (option 1)
 % 2 is No (option 2)
if (choice == 1)
    calcPowerMassRatios(figureNum);   
    figureNum = figureNum +1;
end % Any other option just skips/does nothing


%% ~~~~~~~~~~ ~~~~~~~~~~ ANSYS RESULTS ANALYSIS SIMULATION ~~~~~~~~~~ ~~~~~~~~~~ %% 
choice = menu("Run Simulation for: Finding spacial magnetic field of an Electromagnet?",'Yes','No');
 % 0 is Exit
 % 1 is Yes (option 1)
 % 2 is No (option 2)
if (choice == 1)
    calcElectromagnetField(figureNum);   
    figureNum = figureNum +1;
end % Any other option just skips/does nothing

%% ~~~~~~~~~~ ~~~~~~~~~~ POSITION OF ACTUATORS SIMULATION ~~~~~~~~~~ ~~~~~~~~~~ %% 
choice = menu("Run Simulation for: Finding Position of Actuators?",'Yes','No');
 % 0 is Exit
 % 1 is Yes (option 1)
 % 2 is No (option 2)
if (choice == 1)
    numActuators = 8;
    quadrants = [5 6 7 8];
    findOptimalEMPositions(numActuators, large_EM, small_EM, r, figureNum); 
    figureNum = figureNum + 1;
    findOptimalEMPositions(numActuators, large_EM, small_EM, r, figureNum);
    figureNum = figureNum + 1;
    findOptimalEMPositions(numActuators, large_EM, small_EM, r, figureNum);  
    figureNum = figureNum + 1;
end % Any other option just skips/does nothing
% checkCylinderCollision(centerPoint0, W0, r0, h0, centerPoint1, W1, r1, h1)

%checkCylinderCollision([1 1 1], [0 0 1], 0.07, 0.5, [-1 -1 -1], [0 0 1], 0.07, 0.5)

