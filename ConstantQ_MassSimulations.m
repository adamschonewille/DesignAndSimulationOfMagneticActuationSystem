%% ConstantQ Mass Simulations
% A script for analysing the results of ANSYS simulations for a constant Q
% ratio (L/r_core) of 8 with r_core varying from 10 mm to 60 mm by 5 mm
% increments. 
% Two different experiments were run simultaneously:
%   Constant Current density: J = 3 A / mm^3 (limited Joule Heating)
%       As the size of the electromagnets grow throughout the experiments
%       the cross-sectional area gets larger and the overall current 
%       increases. This is not limited by power in any regard.
%   Constant Power input: 1800 W (max power output of an AC socket 120V 15 A)
%       Here the current density is only limited by the amount of power put
%       into the system. The current density values tested typically are
%       not feasible due to the extreme temperatures that would be produced
%       even with an advanced cooling system attempting to mitigate them.
%       NOTE: this power is assuming all power purely goes into resistive
%       heat losses and does not account for the power going into the
%       magnetic field ( 1/2 L I^2 ). From previous experiments it was
%       found that the magnetic field power was only about 5% efficient so
%       there is some validity to this assumption.

% Due to way the coils were modeled, there was only one turn and large
% currents. This was done to circumvent the need to chose wire diameters
% and potentially skewing results based on geometric fitting. As a result,
% no meaningful inductance values were produced or recorded.

clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
             
figureNum = 1;

%% Post-ANSYS Experiment Data Analysis:
% Import synthesized data

% [file_inductance,path_inductance] = uigetfile('*.txt');
% if isequal(file_inductance,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(path_inductance,file_inductance)]);
% end

uiwait(msgbox('Select the file containing experimental data for: X','Constant Current Density','modal'));
X = uiimport();
X = struct2array(X);


uiwait(msgbox('Select the file containing experimental data for: Y','Constant Current Density','modal'));
Y = uiimport();
Y = struct2array(Y);


uiwait(msgbox('Select the file containing experimental data for: B-Field','Constant Current Density','modal'));
B_mag = uiimport();
B_mag = struct2array(B_mag);

% return;
%% Analysis and Plotting
numPoints = size(X,2);
L = 0;  % starting length of the electromagnet from the origin

colors = 1/255*[235,  64,  52;
                235, 131,  52;
                235, 232,  52;
                162, 235,  52;
                 52, 235, 159;
                 52, 229, 235;
                 52, 140, 235;
                 70,  52, 235;
                131,  52, 235;
                223,  52, 235;
                235,  52, 110];

figure(figureNum)
hold on
for i = 1:numPoints
    % Sort points to be ascending in Y
    [Y(:,i),indices] = sort( Y(:,i) );
    B_temp = B_mag(:,i);
    for j = 1 : size(B_mag,1)
        B_mag(j,i) = B_temp(indices(j));
    end
    plot(Y(:,i)-L, B_mag(:,i),'-o','MarkerSize',2.5,'LineWidth',1.5,'Color',colors(i,:));
%     legend
end
title("Magnetic Field vs Axial Distance [m] for a constant current density of 3 A/mm$^2$");
xlabel("Distance along Y [m]");
ylabel("By Magnetic Field Component [T]");
% ylim([0 0.15])
legend show
legend('R=10mm','R=15mm','R=20mm','R=25mm','R=30mm','R=35mm','R=40mm','R=45mm','R=50mm','R=55mm','R=60mm')
hold off
figureNum = figureNum + 1;

return;

%% Next set
clear X;
clear Y;
clear B_mag;

uiwait(msgbox('Select the file containing experimental data for: X','Constant Power','modal'));
X = uiimport();
X = struct2array(X);


uiwait(msgbox('Select the file containing experimental data for: Y','Constant Power','modal'));
Y = uiimport();
Y = struct2array(Y);


uiwait(msgbox('Select the file containing experimental data for: B-Field','Constant Power','modal'));
B_mag = uiimport();
B_mag = struct2array(B_mag);


numPoints = size(X,2);
L = 0;  % starting length of the electromagnet from the origin

figure(figureNum)
hold on
for i = 1:numPoints
    % Sort points to be ascending in Y
    [Y(:,i),indices] = sort( Y(:,i) );
    B_temp = B_mag(:,i);
    for j = 1 : size(B_mag,1)
        B_mag(j,i) = B_temp(indices(j));
    end
    plot(Y(:,i)-L, B_mag(:,i),'-o','MarkerSize',2.5,'LineWidth',1.5 );
    legend
end
title("Magnetic Field vs Axial Distance [m] For a constant power of 1.8 kW");
xlabel("Distance along Y [m]");
ylabel("By Magnetic Field Component [T]");
% ylim([0 0.15])
legend show
legend('R=10mm','R=15mm','R=20mm','R=25mm','R=30mm','R=35mm','R=40mm','R=45mm','R=50mm','R=55mm','R=60mm')
hold off
figureNum = figureNum + 1;