function [wireDiam resistanceLength] = gaugeLookup( gauge )
% Author: Adam Schonewille
% Date: September 21 2018

%% Description of Inputs
% gauge  - The wire gauge that is being designed for

%% Returned Values
% wireDiam         - The wire diameter corresponding to the gauge in mm
% resistanceLength - The resistance in mOhm/m of the wire gauge

%% Function Code
% if the wire gauge is not an integer then an error is thrown
if ( (gauge <= 0) | gauge > 40 | (mod(gauge,1) ~= 0) )
    error("Wire Gauge must be a positive Integer value between 1-40.")
end

% Table created from wikipedia wore gauge chart:
% https://en.wikipedia.org/wiki/American_wire_gauge
gaugeTable =   [1, 7.348, 0.4066;
                2, 6.544, 0.5127;
                3, 5.827, 0.6465;
                4, 5.189, 0.8152;
                5, 4.621, 1.028;
                6, 4.115, 1.296;
                7, 3.665, 1.634;
                8, 3.264, 2.061;
                9, 2.906, 2.599;
                10,2.588, 3.277;
                11,2.305, 4.132;
                12,2.053, 5.211;
                13,1.828, 6.571;
                14,1.628, 8.286;
                15,1.450, 10.45;
                16,1.291, 13.17;
                17,1.150, 16.61;
                18,1.024, 20.95;
                19,0.912, 26.42;
                20,0.812, 33.31;
                21,0.723, 42.00;
                22,0.644, 52.96;
                23,0.573, 66.79;
                24,0.511, 84.22;
                25,0.455, 106.2;
                26,0.405, 133.9;
                27,0.361, 168.9;
                28,0.321, 212.9;
                29,0.286, 268.5;
                30,0.255, 338.6;
                31,0.227, 426.9;
                32,0.202, 538.3;
                33,0.180, 678.8;
                34,0.160, 856.0;
                35,0.143, 1079;
                36,0.127, 1361;
                37,0.113, 1716;
                38,0.101, 2164;
                39,0.0897, 2729;
                40,0.0799, 3441];

wireDiam = gaugeTable(gauge,2);
resistanceLength = gaugeTable(gauge,3);
end