function [dist] = calcDistanceBetweenSpherical( p0, m0, r0, L0, p1, m1, r1, L1 )
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
% dist - the sum of each distance between cylinders when approximated as 
% several (5) stacked spheres. 
% NOTE: the spheres that approximate the cylinders make the cylinder seem
% 'pill-shaped' as there is a sphere located at each end of the EM.


numSpheres = 7; % this is fixed? THIS MUST BE >= 2 !!!

spheres_pos0 = zeros(3,numSpheres);
spheres_pos1 = zeros(3,numSpheres);

% Initialize the positions of all the spheres that make up EM 1 and EM 2
% move the starting point from the center of each hEM to their ends
p_end0 = p0 - L0/2*m0; 
p_end1 = p1 - L1/2*m1; 
% Find separation distances between evenly spaced spheres on each EM
dL0 =  L0/(numSpheres-1);
dL1 =  L1/(numSpheres-1);

% assign sphere positions by adding dL to the starting point in the
% unit vector direction of the EM (m0 and m1)
for n = 1:numSpheres
    spheres_pos0(:,n) = p_end0 + (n-1)*dL0*m0;
    spheres_pos1(:,n) = p_end1 + (n-1)*dL1*m1;
end

% Time to find the total distances across each EM
dist = 0;
for i = 1:numSpheres
    for j = 1:numSpheres
        % distance between spheres is the difference in center position
        % minus their combined radii.
        temp = norm(spheres_pos0(:,i)-spheres_pos1(:,j))-(r0+r1);
        % If any of the distances are close to 0 we want them to blowup so
        % we add the distances in parallel to exagerate the smaller values
        dist = dist + 1/temp;
        
    end
    
end

% Distances should be calculated now :) 
dist = 1/dist;

end