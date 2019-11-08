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


%% Analysis and Plotting
numExperiments = size(X,2);
numPoints = size(X,1);
L = 0;  % starting length of the electromagnet from the origin
% arbitrary rainbow colours:
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
% Sort points to be ordered and then plot
for i = 1:numExperiments
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


%% Average Values in the plot for same points

figure(figureNum)
hold on
Y_avg = [];
B_mag_avg = [];
% Take sorted values and average the values that are on the same point
prev_Val = NaN;
initialize_flag = 1; 

% cycle through columns where each column is an experiment
for n = 1:numExperiments
    index = 1;
    % For each experiment cycle through the datapoints
    for i = 1:numPoints
        
        count = 0;
        B_acc = 0;
        % If the data point is the same as the previous or it is NaN 
        % then skip it         
        if ( ( Y(i,n) ~= prev_Val ) && (isnan(Y(i,n)) == 0) )
            % look ahead through the data and check to see if a value
            % repeats. Accumulate these values.
            for j = i:numPoints
                if ( abs( Y(j,n) - Y(i,n) ) <= 1e-6 )
                    B_acc = B_acc + B_mag(j,n);
                    count = count + 1;
                end
            end
            % Combine identical points to be one point with an average
            % value
            Y_avg(index,n) = Y(i,n);
            B_mag_avg(index,n) = B_acc/count;
            index = index + 1;
            % Set prevVal as this number to skip it in the future. The data
            % is ordered in ascending order so this should suffice
            prev_Val = Y(i,n);
        end
    end
    plot(Y_avg(:,n)-L, B_mag_avg(:,n),'-o','MarkerSize',2.5,'LineWidth',1.5,'Color',colors(n,:));
    if (initialize_flag)
        % set only once, ensures that the matrix ends in NaN for blank
        % results and not in zeros that will be plotted.
        Y_avg = [Y_avg NaN*ones(length(Y_avg),numExperiments-1)];
        B_mag_avg = [B_mag_avg NaN*ones(length(B_mag_avg),numExperiments-1)];
        initialize_flag = 0;
    end
end
title("Averaged Magnetic Field vs Axial Distance [m] for a constant current density of 3 A/mm$^2$");
xlabel("Distance along Y [m]");
ylabel("Average By Magnetic Field Component [T]");
% ylim([0 0.15])
legend show
legend('R=10mm','R=15mm','R=20mm','R=25mm','R=30mm','R=35mm','R=40mm','R=45mm','R=50mm','R=55mm','R=60mm')
hold off
figureNum = figureNum + 1;


%% Next set

uiwait(msgbox('Select the file containing experimental data for: X','Constant Power','modal'));
X_P = uiimport();
X_P = struct2array(X_P);


uiwait(msgbox('Select the file containing experimental data for: Y','Constant Power','modal'));
Y_P = uiimport();
Y_P = struct2array(Y_P);


uiwait(msgbox('Select the file containing experimental data for: B-Field','Constant Power','modal'));
B_mag_P = uiimport();
B_mag_P = struct2array(B_mag_P);


numExperiments_P = size(X_P,2);
numPoints_P = size(X_P,1);
L = 0;  % starting length of the electromagnet from the origin

figure(figureNum)
hold on
for i = 1:numExperiments_P
    % Sort points to be ascending in Y
    [Y_P(:,i),indices] = sort( Y_P(:,i) );
    B_temp = B_mag_P(:,i);
    for j = 1 : size(B_mag_P,1)
        B_mag_P(j,i) = B_temp(indices(j));
    end
    plot(Y_P(:,i)-L, B_mag_P(:,i),'-o','MarkerSize',2.5,'LineWidth',1.5 );
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

%% Average the Power results
Y_P_avg = [];
B_mag_P_avg = [];
% Take sorted values and average the values that are on the same point
prev_Val = NaN;
initialize_flag = 1; 

figure(figureNum)
hold on
% cycle through columns where each column is an experiment
for n = 1:numExperiments_P
    index = 1;
    % For each experiment cycle through the datapoints
    for i = 1:numPoints_P
        
        count = 0;
        B_acc = 0;
        % If the data point is the same as the previous or it is NaN 
        % then skip it         
        if ( ( Y_P(i,n) ~= prev_Val ) && (isnan(Y_P(i,n)) == 0) )
            % look ahead through the data and check to see if a value
            % repeats. Accumulate these values.
            for j = i:numPoints
                if ( abs( Y_P(j,n) - Y_P(i,n) ) <= 1e-6 )
                    B_acc = B_acc + B_mag_P(j,n);
                    count = count + 1;
                end
            end
            % Combine identical points to be one point with an average
            % value
            Y_P_avg(index,n) = Y_P(i,n);
            B_mag_P_avg(index,n) = B_acc/count;
            index = index + 1;
            % Set prevVal as this number to skip it in the future. The data
            % is ordered in ascending order so this should suffice
            prev_Val = Y_P(i,n);
        end
    end
    plot(Y_avg(:,n)-L, B_mag_P_avg(:,n),'-o','MarkerSize',2.5,'LineWidth',1.5,'Color',colors(n,:));
    if (initialize_flag)
        % set only once
        Y_P_avg = [Y_P_avg NaN*ones(length(Y_P_avg),numExperiments_P-1)];
        B_mag_P_avg = [B_mag_P_avg NaN*ones(length(B_mag_P_avg),numExperiments_P-1)];
        initialize_flag = 0;
    end
end
title("Averaged Magnetic Field vs Axial Distance [m] for a constant Power of 1.8kW");
xlabel("Distance along Y [m]");
ylabel("Average By Magnetic Field Component [T]");
% ylim([0 0.15])
legend show
legend('R=10mm','R=15mm','R=20mm','R=25mm','R=30mm','R=35mm','R=40mm','R=45mm','R=50mm','R=55mm','R=60mm')
hold off
figureNum = figureNum + 1;