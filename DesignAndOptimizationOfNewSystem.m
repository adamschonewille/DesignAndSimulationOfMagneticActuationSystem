%% Simulations for designing the next novel actuation system.
%
%  Author: Adam Schonewille
%  MASc Student, University of Toronto
%
%  Start Date: July 5th 2019
%  Last Editted: July 5th 2019
%

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
rho_Cu = 1.68e-8;       % [Ohm*m] Resisticity @ 293 K
N2_temp = 77;           % [K] (Liquid Nitrogen)

% For Surgical applications:
% Head Dimensions taken from the 99th Percentile as an upper limit (max)
% https://en.wikipedia.org/wiki/Human_head
max_Head_Width  = 0.165; % [m] Ear to Ear
max_Head_Length = 0.217; % [m] Forehead (brows) to Back of Head
max_Head_Height = 0.255; % [m] Bottom of Chin to Top of Head

%% ~~~~~~~~~~ ~~~~~~~~~~ OPERATING PARAMETERS ~~~~~~~~~~ ~~~~~~~~~~ %%
magnetSize = 0.0254; % [m]
Remenance = 1.45; % [T]

%% ~~~~~~~~~~ ~~~~~~~~~~ SIMULATION TARGETS ~~~~~~~~~~ ~~~~~~~~~~ %%
isotropicFields = true;
MaxB_Field = 0.040; %[T]


% for plotting
figureNum = 1;




