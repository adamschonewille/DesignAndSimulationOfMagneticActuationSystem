function [dist] = calcMinDistBetweenCylinders(p0, m0, r0, L0, p1, m1, r1, L1)
%calcDistanceBetweenSpherical
% INPUTS
% p0 - center position of the first EM
% m0 - orientation unit vector of the first EM
% r0 - radius of the first EM
% L0 - length of the first EM
% p1 - center position of the second EM
% m1 - orientation unit vector of the second EM
% r1 - radius of the second EM
% L1 - length of the second EM
% OUTPUTS
% dist - the sum of each distance between cylinders when modeled as 
% finite axes.

lambda = cross(m0, m1);

% check first if the two lines are parallel since the algorithms only work
% for skewed line pairs. The two lines will be parallel if their cross
% product is 0
if ( norm(lambda) == 0 )
    diff = p1-p0;
    dist = norm( cross(m0,diff)/norm(m0) );
    % If the minimum distance is = 0 this means that the two lines are
    % colinear and the cylinders are either pointing towards one another
    % from a distance, pointing the same direction, or completely
    % overlapping.
    
    % This current representation of the min distance for parallel is not
    % ideal by any means, but we assume that in a numerical simulation this
    % loop will rarely be called due to the precise requirements needed.
else
    diff = p1-p0;
    dist = norm( dot( lambda/norm(lambda), diff) );
    
    % ensure that this distance is between two pointst that fall on each
    % cylinder or else we must take the distance from the end of one or 
    % more of the cylinder ends.
    lambda0 = cross(m0,lambda);
    lambda1 = cross(m1,lambda);
    c0 = p0 + dot(diff, lambda1)/dot(m0, lambda1) * m0;
    c1 = p1 + dot(-diff,lambda0)/dot(m1, lambda0) * m1;
    
    cond0 = norm(c0-p0);
    cond1 = norm(c1-p1);
    
    if (cond0 > L0/2)
        if (cond1 > L1/2)
            % both points lay outside the cylinders
            % disregard the distance value found ad test all permutations
            % of distances from cylinder 0 end to cylinder 1 end:
            distances = [norm( (p0+L0/2*m0) - (p1+L1/2*m1) );
                         norm( (p0+L0/2*m0) - (p1-L1/2*m1) );
                         norm( (p0-L0/2*m0) - (p1-L1/2*m1) );
                         norm( (p0-L0/2*m0) - (p1+L1/2*m1) )];
            dist = min(distances);
        else
            % only point 0 is outside of cylinder 0
            % Find the new minimum between two ends and point 1
            distances = [norm( (c1) - (p0+L0/2*m0) );
                         norm( (c1) - (p0-L0/2*m0) )];
            dist = min(distances);
        end
    elseif (cond1 > L1/2)
        % no point checking if cond0 is satisfied again    
        % only point 1 is outside of cylinder 1
        % Find the new minimum between two ends and point 0
        distances = [norm( (c0) - (p1+L1/2*m1) );
                     norm( (c0) - (p1-L1/2*m1) )];
        dist = min(distances);
    end
    % if neither of the conditions were not met then the minimum distance
    % will not have changed.
    
end

% currently the minimum distance is the mimimum distance between two axes,
% to get the distance between cylinder surfaces we assume that the mim
% distance vector will be perpendicular to m and therefore be also
% perpendicular to the cylinder surface. This distance must then have the
% radii of the two cylinders subtracted from it.

dist = dist - r0 - r1;


end