function [F] = forceMatrix(mag_moment)
% AUTHOR: Adam Schonewille
% DATE: January 6th 2020
% ABOUT:    Finds the equivalent matrix of a 3 x 1 or 1 x 3 vector
%           for finding the force 
% 
%% INPUTS
% mag_moment -  A vector of size 3 x 1
%                        or
%               A vector of size 1 x 3
% 
%% RETURNED OUTPUTS
% F  - The 3 x 5 Matrix for Force Calculations (similar to Skew-Symmetric
% Matrix)

%% Function Code
% Checks that the size indexes are 1 and 3 or that they are 3 and 1.
% Each index length much match summing the logical values to 2 and checking
% against 2 for both.
m = mag_moment;
    if( (sum(size(m) == [1,3]) == 2) || (sum(size(m) == [3,1]) == 2) )
        % If vector is the right length then create the Matrix
        F = [ m(1)  m(2)  m(3)    0     0 ;
                0   m(1)    0   m(2)  m(3);
             -m(3)    0   m(1) -m(3)  m(2)];
    else 
        error("The input vector is not a 3 x 1 or 1 x 3 vector.")
    end

end