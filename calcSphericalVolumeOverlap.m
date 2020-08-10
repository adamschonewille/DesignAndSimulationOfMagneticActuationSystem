function [volume] = calcSphericalVolumeOverlap( p0, m0, r0, L0, p1, m1, r1, L1 )
%calcSphericalVolumeOverlap
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
% volume - the total volume of the overlapping EM when approximated as 
% several (5) stacked spheres. If the EMs do not overlap, then the function returns 0.
% NOTE: the spheres that approximate the cylinders make the cylinder seem
% 'pill-shaped' as there is a sphere located at each end of the EM.


% check if actuators are really far apart. If their bounding spheres are
% not intersecting then they do not overlap
R0 = 1/2*sqrt(L0^2+(r0*2)^2);
R1 = 1/2*sqrt(L1^2+(r1*2)^2); 
if ( 0 == findSphereOverlap( p0, R0, p1, R1 ) )
    volume = 0;
    return;
end

% if initial quickcheck fails then calc everything

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

% Time to find the total overlapping volume across each EM
volume = 0;
for i = 1:numSpheres
    for j = 1:numSpheres
        volume = volume + findSphereOverlap( spheres_pos0(:,i), r0, ...
                                             spheres_pos1(:,j), r1 );
    end
    
end

% volume should be calculated now :) 
% volume;

end