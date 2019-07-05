function [skewSymmetricMatrix] = skew(X)
% AUTHOR: Adam Schonewille
% DATE: November 1st 2018
% ABOUT:    Finds the skew-symmetric matric of a 3 x 1 or 1 x 3 vector
% 
%% INPUTS
% X           - A vector of size 3 x 1
%                        or
%               A vector of size 1 x 3
% 
%% RETURNED OUTPUTS
% skewSymmetricMatrix  - The 3 x 3 Skew-Symmetric Matrix

%% Function Code
% Checks that the size indexes are 1 and 3 or that they are 3 and 1.
% Each index length much match summing the logical values to 2 and checking
% against 2 for both.
    if( (sum(size(X) == [1,3]) == 2) || (sum(size(X) == [3,1]) == 2) )
        % If vector is the right length then skew the vector
        skewSymmetricMatrix = [   0  -X(3)  X(2); 
                                X(3)    0  -X(1);
                               -X(2)  X(1)    0 ];
    else 
        error("The input vector is not a 3 x 1 or 1 x 3 vector.")
    end

end

