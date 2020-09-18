function [T] = torqueMatrix(mag_moment)
% AUTHOR: Adam Schonewille
% DATE: September 14th 2020
% ABOUT:    Finds the equivalent matrix of a 3 x 1 or 1 x 3 vector
%           for finding the torque. This is the same as the skew symmetric 
%           matrix 
% 
%% INPUTS
% mag_moment -  A vector of size 3 x 1
%                        or
%               A vector of size 1 x 3
% 
%% RETURNED OUTPUTS
% F  - The 3 x 3 Matrix for Torque Calculations (similar to Skew-Symmetric
% Matrix)
%
%
% See also FORCEMATRIX.

%% Function Code
% Checks that the size indexes are 1 and 3 or that they are 3 and 1.
% Each index length much match summing the logical values to 2 and checking
% against 2 for both.
m = mag_moment;
    if( (sum(size(m) == [1,3]) == 2) || (sum(size(m) == [3,1]) == 2) )
        % If vector is the right length then create the Matrix
        T = [   0   -m(3)   m(2);
              m(3)   0     -m(1);
             -m(2)   m(1)   0   ];
    else 
        error("The input vector is not a 3 x 1 or 1 x 3 vector.")
    end

end