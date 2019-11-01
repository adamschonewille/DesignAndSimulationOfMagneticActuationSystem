%% Process / Synthesis ANSYS Results
% Requires data points across the line X = 0 and takes a 2D slice in the XY
% plane as a input. If the data is in a different plane it can be
% explicitly changed in this code.

% Saves 3 data files each in a matrix format where the columns are
% individual experiments and the rows are the results from each experiment
%   X: The X value of the values along X = 0. As such this should be an NxM
%   matrix of zeros. Measured in meters
%   Y: The Y values along the axis where X = 0 (y-axis). These points are 
%   the from the face of the electromagnet outwards. The starting point can
%   be explicitly set below. Measured in meters
%   B_mag: The component of the magnetic field in the y-direction along the
%   y-axis. Measured in Teslas

clear all; clc; close all;

%% Initialize these variables based on the ANSYS simulation 
start_length = 0.0; %[m] The end of the electromagnet. Values greater than 
                    %this along the axis will be recorded
electromagnet_Axis = 1; % X = 1, Y = 2, Z = 3
perpendicular_Axis = 3; % axis that forms the plane

%% Set up parameters and start synthesizing data

X = [];
Y = [];
MagneticField_y_Dir = [];
format_str = "";    % for printing data in it's respective colums

% Prompt for dataset
choice = menu("Choose a Dataset?",'Yes','No');
 % 0 is Exit
 % 1 is Yes (option 1)
 % 2 is No (option 2)
dataset_Num = 1; % initalize the dataset number to record experiments
while (choice == 1) 
    format_str = format_str + "%.8f";
    uiimport(); % opens the import data menu
    disp("Choose ANSYS data file for analysis. Press any key to continue ... ");
    pause(); % press any button to continue after the data is imported
    
    X_temp = [];
    Y_temp = [];
    MagneticField_y_Dir_temp = [];
    
    % EXTRACTION OF AXIAL MAGNETIC FIELD VALUES FROM LARGER DATASET
    % Looks for only data points where x = 0 and makes a temporary vector
    
    % Different cases depending on the setup of the ANSYS simulation:
    if ( electromagnet_Axis == 1 && perpendicular_Axis == 3 )
        for n=1:length(XLocationm)
            if (ZLocationm(n) == 0) && (XLocationm(n) >= start_length)
                X_temp = [X_temp;  ZLocationm(n)];
                Y_temp = [Y_temp;  XLocationm(n)];
                MagneticField_y_Dir_temp = [MagneticField_y_Dir_temp; DirectionalMagneticFluxDensityT(n) ];
            end
        end        
        
    elseif ( electromagnet_Axis == 2 && perpendicular_Axis == 1 )
        for n=1:length(XLocationm)
            if (XLocationm(n) == 0) && (YLocationm(n) >= start_length)
                X_temp = [X_temp;  XLocationm(n)];
                Y_temp = [Y_temp;  YLocationm(n)];
                MagneticField_y_Dir_temp = [MagneticField_y_Dir_temp; DirectionalMagneticFluxDensityT(n) ];
            end
        end
    end
    
    elements = length(X_temp);
    
    % Add the entire temporary vector to the total vector
    % check that the vector is a compactible size to the matrix and zeropad
    % if it is not.
    if ( length(X_temp) < size(X,1) )
        X(:,dataset_Num) = [X_temp; NaN*zeros(size(X,1)-length(X_temp),1)];
        Y(:,dataset_Num) = [Y_temp; NaN*zeros(size(Y,1)-length(Y_temp),1)];
        MagneticField_y_Dir(:,dataset_Num) = [MagneticField_y_Dir_temp; NaN*zeros(size(Y,1)-length(Y_temp),1) ];
    else 
        X(:,dataset_Num) = X_temp;
        Y(:,dataset_Num) = Y_temp;
        MagneticField_y_Dir(:,dataset_Num) = MagneticField_y_Dir_temp;
    end 
    disp("Dateset " + num2str(dataset_Num) + " analyzed.");
    disp("Number of elements: " + num2str(elements));
    dataset_Num = dataset_Num+1;
    clear MagneticField_y_Dir_temp
    clear XLocationm;
    clear YLocationm;
    clear DirectionalMagneticFluxDensityT;
    % Not necessary but do anyways to avoid build up of variables:
    clear NodeNumber;
    clear ZLocationm;
    
    choice = menu("Choose another Dataset?",'Yes','No');
    if (choice == 1)
        format_str = format_str + "     "; % Add a tab for spacing, else add a new line
    end
end % Any other option just skips/does nothing
format_str = format_str + "\n";     % move to the next line

disp("All datesets collected.")

return;
%% Write to file 

NumStr = 1:1:size(X,2);
X_header = "X_exp_"+NumStr;
Y_header = "Y_exp_"+NumStr;
B_header = "B_exp_"+NumStr;
% Colums are datasets for each experiment
outputFileName_X = strcat(date, '_outputData_X.txt');
outputFileName_Y = strcat(date, '_outputData_Y.txt');
outputFileName_B_Field = strcat(date, '_outputData_B_Field.txt');

disp("Writing to file now ...");
fid = fopen(outputFileName_X,'wt');
% fprintf(fid,'%.8f\n',X_header); % I don't think that the header is necessary
fprintf(fid,format_str,X'); 
% writematrix(X,outputFileName_X,'Delimiter','tab')
fclose(fid);
clear fid


fid = fopen(outputFileName_Y,'wt');
% fprintf(fid,'%.8f\n',X_header); % I don't think that the header is necessary
fprintf(fid,format_str,Y');
fclose(fid);
clear fid

fid = fopen(outputFileName_B_Field,'wt');
% fprintf(fid,'%.8f\n',X_header); % I don't think that the header is necessary
fprintf(fid,format_str,MagneticField_y_Dir');
fclose(fid);
clear fid

disp("File written.");