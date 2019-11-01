%% ConstantCurrentDensityAndGeometry
% A script for finding the input parameters for the ANSYS simulation based
% off of a constant coil geometry and constant current density J.

clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');

% Coil Geometry for Testing
L = 0.210;      % [m]
r_i = 0.021;    % [m]
r_o = 0.0315;   % [m]
r = r_o - r_i;  % [m]

max_Area = L*r*pi/4; % pi/4 is the max packing factor possible

% Set Current Density J
J = 3;          % [A/mm^2] recommended to use 3 A/mm^2 to avoid excessive heating
                % with a cooling system this number can be 6 A/mm^2 (Octomag)
J = J * 1e6;    % [A/m^2]
                
                
test_gauges = 6:1:28; % Gauges 6 through 28 are the arbitrary gauges tested
numPoints = length(test_gauges);

% Initialize vectors:
number_Of_Turns = zeros(numPoints,1);       % [] also known as N
effective_Total_Area = zeros(numPoints,1);  % [m^2]
total_Wire_Length = zeros(numPoints,1);     % [m]
resistance = zeros(numPoints,1);            % [Ohm]
applied_Current = zeros(numPoints,1);       % [A] I = J * A
power_Loss = zeros(numPoints,1);            % [W] P = I^2 R
power_Required = zeros(numPoints,1);        % [W] P = I^2 R + LI^2

wireDiam = zeros(numPoints,1);              % [m]
resistivity = zeros(numPoints,1);           % [Ohm*m]

for i = 1 : numPoints
    [wireDiam(i) resistanceLength] = gaugeLookup( test_gauges(i) );
    wireDiam(i) = wireDiam(i)/1000; % mm to m
    resistanceLength = resistanceLength/1000; % mOhm to Ohm
    
    resistivity(i) = resistanceLength * ( pi/4*wireDiam(i)^2 ); % rho = R / l * A
    n_a = floor( L/wireDiam(i) );
    n_r = floor( r/wireDiam(i) );
    number_Of_Turns(i) = n_a * n_r; % cubic patterning of wires when stacked
    % Cross-sectional Area of one wire times by N
    effective_Total_Area(i) = pi/4 * wireDiam(i)^2 * number_Of_Turns(i);
    % find length of wire needed and resulting resistance
    total_Wire_Length(i) = 2*pi*n_r*n_a*(r_i+n_r*wireDiam(i)/2); % see notes
    resistance(i) = resistivity(i)*total_Wire_Length(i) / ( pi/4*wireDiam(i)^2 ); % R = rho*l/A
    % Current required based on current density and Power necessary
    applied_Current(i) = J * pi/4*wireDiam(i)^2;
    power_Loss(i) = applied_Current(i)^2 * resistance(i); % P = I^2 R
    
end

% Needed for ANSYS Simulation:
number_Of_Turns
effective_Total_Area
applied_Current

figureNum = 1;
figure(figureNum)
hold on
plot(test_gauges, power_Loss,'.-',"MarkerSize",15, "Color", [119/255, 41/255, 83/255]); %[119/255, 41/255, 83/255]
title("Power consumption vs Wire Gauge [AWG]");
xlabel("Wire gauge [AWG]");
ylabel("Heat power loss [W]");
hold off
figureNum = figureNum + 1;
% Ubuntu Purple (CANONICAL AUBERGINE):  119, 41, 83
% Ubuntu Orange:                        221, 72, 20     
% Ubuntu Terminal Purple (Darker):       48, 10, 36

%% Post-ANSYS Experiment Data Analysis:

% [file_inductance,path_inductance] = uigetfile('*.txt');
% if isequal(file_inductance,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(path_inductance,file_inductance)]);
% end


uiwait(msgbox('Select the file containing experimental data for: Inductance','Constant current dennsity and Geometry Experiments','modal'));
inductance = uiimport();
inductance = struct2array(inductance);


uiwait(msgbox('Select the file containing experimental data for: X','Constant current dennsity and Geometry Experiments','modal'));
X = uiimport();
X = struct2array(X);


uiwait(msgbox('Select the file containing experimental data for: Y','Constant current dennsity and Geometry Experiments','modal'));
Y = uiimport();
Y = struct2array(Y);


uiwait(msgbox('Select the file containing experimental data for: B-Field','Constant current dennsity and Geometry Experiments','modal'));
B_mag = uiimport();
B_mag = struct2array(B_mag);


%% Analysis and Plotting
bandwidth = resistance./(2*pi*inductance);

figure(figureNum)
hold on
plot(wireDiam, (resistance.*applied_Current.^2 + 1/2*inductance.*applied_Current.^2),'.',"MarkerSize",15, "Color", [209/255, 19/255, 51/255]); %[119/255, 41/255, 83/255]
title("Total Power vs Wire Diameter [m]");
xlabel("Wire Diameter [m]");
ylabel("Power [W]");
ylim([0 50])
hold off
figureNum = figureNum + 1;

figure(figureNum)
hold on
plot(wireDiam, 100*1/2*(inductance.*applied_Current.^2)./(resistance.*applied_Current.^2 + 1/2*inductance.*applied_Current.^2),'.',"MarkerSize",15, "Color", [209/255, 19/255, 51/255]); %[119/255, 41/255, 83/255]
title("Efficiency vs Wire Diameter [m]");
xlabel("Wire Diameter [m]");
ylabel("Efficiency [\%]");
hold off
figureNum = figureNum + 1;

figure(figureNum)
hold on
plot(wireDiam, resistance,'.',"MarkerSize",15, "Color", [209/255, 19/255, 51/255]); %[119/255, 41/255, 83/255]
title("Resistance vs Wire Diameter [m]");
xlabel("Wire Diameter [m]");
ylabel("Resistance [Ohm]");
hold off
figureNum = figureNum + 1;



figure(figureNum)
hold on
plot(wireDiam, inductance,'.',"MarkerSize",15, "Color", [66/255, 135/255, 245/255]); %[119/255, 41/255, 83/255]
title("Inductance vs Wire Diameter [m]");
xlabel("Wire Diameter [m]");
ylabel("Inductance [H]");
hold off
figureNum = figureNum + 1;


figure(figureNum)
hold on
plot(wireDiam, bandwidth,'.',"MarkerSize",15, "Color", [66/255, 135/255, 245/255]); %[119/255, 41/255, 83/255]
title("Bandwidth vs Wire Diameter [m]");
xlabel("Wire Diameter [m]");
ylabel("Cutoff Frequency [Hz]");
hold off
figureNum = figureNum + 1;

figure(figureNum)
hold on
plot(wireDiam, 100*effective_Total_Area./max_Area,'.',"MarkerSize",15, "Color", [66/255, 135/255, 245/255]); %[119/255, 41/255, 83/255]
title("Effective coil winding area used vs Wire Diameter [m]");
xlabel("Wire Diameter [m]");
ylabel("Efficient area used [\%]");
hold off
figureNum = figureNum + 1;

figure(figureNum)
hold on
for i = 1:numPoints
    plot(Y(:,i)-L, B_mag(:,i),'.',"MarkerSize",15 );
    legend
end
title("Magnetic Field vs Axial Distance [m]");
xlabel("Distance along Y [m]");
ylabel("By Magnetic Field Component [T]");
ylim([0 0.15])
legend show
hold off
figureNum = figureNum + 1;

%
figure(figureNum)
% Normalized to account for wire wrappings that are not ideal for the
% constant geometry of the experiment
hold on
plot(wireDiam, bandwidth./(effective_Total_Area./max_Area),'.',"MarkerSize",15, "Color", [66/255, 135/255, 245/255]); %[119/255, 41/255, 83/255]
title("Normalized Bandwidth vs Wire Diameter [m]");
xlabel("Wire Diameter [m]");
ylabel("Cutoff Frequency (Normalized) [Hz]");
hold off
figureNum = figureNum + 1;

return;

%% Determine the reason for the curved results from the Total Resistive Power

wireDiameter_continuous = 0.050:0.00050:10; %[mm] from 50 microns to 10 mm
avg_resistivity = 1.7248e-08; %[Ohm m] from the table previously, also 1.72e-8 online
numPoints_continuous = length(wireDiameter_continuous);

% Initialize vectors:
number_Of_Turns_continuous = zeros(numPoints_continuous,1);       % [] also known as N
total_Wire_Length_continuous = zeros(numPoints_continuous,1);     % [m]
resistance_continuous = zeros(numPoints_continuous,1);            % [Ohm]
applied_Current_continuous = zeros(numPoints_continuous,1);       % [A] I = J * A
power_Loss_continuous = zeros(numPoints_continuous,1);            % [W] P = I^2 R


for i = 1 : numPoints_continuous
    wireDiameter_continuous(i) = wireDiameter_continuous(i)/1000; % mm to m
    n_a = floor( L/wireDiameter_continuous(i) );
    n_r = floor( r/wireDiameter_continuous(i) );
    number_Of_Turns_continuous(i) = n_a * n_r; % cubic patterning of wires when stacked
    % find length of wire needed and resulting resistance
    total_Wire_Length_continuous(i) = 2*pi*n_r*n_a*(r_i+n_r*wireDiameter_continuous(i)/2); % see notes
    resistance_continuous(i) = avg_resistivity*total_Wire_Length_continuous(i) / ( pi/4*wireDiameter_continuous(i)^2 ); % R = rho*l/A
    % Current required based on current density and Power necessary
    applied_Current_continuous(i) = J * pi/4*wireDiameter_continuous(i)^2;
    power_Loss_continuous(i) = applied_Current_continuous(i)^2 * resistance_continuous(i); % P = I^2 R
end

figure(figureNum)
% Shows the reason for the bump-like trends in the plotted graphs. Aliasing
% occurs for a set geometry and changing wire diameter. The Wire sizes that
% were determined by AWG are coincidental.
hold on
plot(wireDiameter_continuous,power_Loss_continuous,'-',"MarkerSize",15, "Color", [119/255, 41/255, 83/255]); %[119/255, 41/255, 83/255])
plot(wireDiam, power_Loss,'.',"MarkerSize",15, "Color", [66/255, 135/255, 245/255]); %[119/255, 41/255, 83/255]
title("Comparison of Power loss vs Wire Diameter [m] for AWG and continuous wire");
xlabel("Wire Diameter [m]");
ylabel("Power Consumption due to Joule Heating [W]");
ylim([0 45]);
hold off
figureNum = figureNum + 1;

% Double check that this approximation below is correct
power_Loss_max = (max_Area*J)^2 * (avg_resistivity * (2*pi*(r_i+r/2)) / ( max_Area ) ); % P = I^2 R