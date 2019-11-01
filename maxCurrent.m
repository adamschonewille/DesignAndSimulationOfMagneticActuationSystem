function [Imax, Bmax, L_coil] = optimizeCoil(availablePower, ohms_per_Length, Ri, Ro, wireDiam, Z)
% Author: Adam Schonewille
% Code Adapted from: Sajad Salmanipour
% Date Created: September 21 2018
% Last Edited:  July 15 2019
%% INPUTS
% availablePower    -- [W] The power of the electrical system in Watts
% ohms_per_Length   -- [ohm/m] The wire resistivity in Ohm per meter
% Ri                -- [m] The inner radius of the entire coil
% Ro                -- [m] The outer radius of the entire coil
% wireDiam          -- [m] The wire diameter
% Z                 -- [m] The distance to center of workspace from the
%                      coil face
%% RETURNED OUTPUTS
% Imax              -- [A] The current at the optimal power and largest
%                   produced B-field
% Bmax              -- [T] The largest magnetic field that can be generated
%                   when operating at max power
% L_coil            -- [m] The length of the coil at these conditions

%% Function Code        
% magnetic permeability in free space
u0 = 4*pi*10^(-7);
numPoints = 100; % number of test points (N*dI) for current optimization
maximumCurrent = 20; % [A]
% Look at all current scenarios from 1 A to 20 A
I_vector = linspace(0.5, maximumCurrent, numPoints);
% Corresponding B-fields achieved for each Current
B_vector = zeros(1,numPoints);
% Nturns_radial - is the number of turns radially in the coil (indep. of current)
Nturns_radial = ceil( (Ro-Ri)/wireDiam );
% Sum up the circumference of each radial loop to get the length of wire
% for one axial slice
length_of_axial_slice = 0;
for n = 0:Nturns_radial-1
    % go through each loop that make up one slice and add up their
    % circumferences to get the total length in a radial slice
    length_of_axial_slice = length_of_axial_slice + 2*pi*(Ri + wireDiam/2 + n*wireDiam);
end

% Avg radius for quicker calculations
rAvg = (Ro+Ri)/2;
% loop
L_vector = zeros(1,numPoints);  % recorded length vectors for each optimal electromagnet
for j=1:numPoints
        % select current from span of values
        I_test = I_vector(j);   
        % Calculate the total Resistance from the power available - Ohms
        coil_resistance = availablePower/I_test^2;
        % Length of wire - m
        L_wire = coil_resistance/ohms_per_Length;
        % find # of axial turns
        Nturns_axial = L_wire/length_of_axial_slice;
%         % N - Total number of turns (Averaging method)
%         N = coil_resistance/(2*pi*((Ri+Ro)/2)*ohms_per_Length);

        % L_coil - resulting length of the coil
        L_coil = ceil(Nturns_axial)*wireDiam;
        L_vector(j) = L_coil;
        if( Nturns_axial < 1 ), break,  end
        % Add up all the contributions to the on-axis Magnetic Field
        Bsol = 0;
        for m=0:ceil(Nturns_axial)-1
            % for each axial layer of the coil compute the distance to the
            % workspace (working backwards from the face of the coil)
            d = Z + wireDiam/2 + m*wireDiam;
            % Assume all the coils are located on top of each other at the
            % average coil radius
            Bsol  = Bsol + 0.5*u0*Nturns_radial*I_test*rAvg^2*(rAvg^2+d^2)^(-1.5);
        end
        % save the resulting coil Bfield in the vector
        B_vector(j)= Bsol;
        % why is this increased 1000 fold
end
% find the maximum Bfield configuration.
[Bmax,Index] = max(B_vector);
%     Bmax
Imax = I_vector(Index); 
L_coil = L_vector(Index);

% figure(1)
% hold on
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
% set(0,'defaulttextInterpreter','latex');
% %ylim([0 1.05])
% %xlim([1-d_lambda 1+d_lambda])
% %legend('Calculated Reflectance','Location','SouthEast')
% plot(I_vector,B_vector/B_vector(Index))
% plot(I_vector,L_vector/max(L_vector))
% xlabel("Test Current (A)")
% ylabel("Resulting B-Field $$\left| B \right|$$ (T)")
% title("Normalized B-field strength and Coil Length")
% hold off

disp("Ri = " + num2str(Ri) + " and Ro = " + num2str(Ro) + " Complete")

end