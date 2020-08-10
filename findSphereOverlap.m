function [volume] = findSphereOverlap( p0, r0, p1, r1 )
%findSphereOverlap
% INPUTS
% p0 - position of first sphere
% r0 - radius of first sphere
% p1 - position of second sphere
% r1 - radius of second sphere
% OUTPUTS
% volume - the volume of the overlapping spheres. If the spheres do not
% overlap, then the function returns 0.


% set the same as wolfram: https://mathworld.wolfram.com/Sphere-SphereIntersection.html
R = r0;
r = r1;

% Find the distance between them d
d = norm(p0-p1);
if ( d > (r0+r1) )
    % if the seperating distance is larger than the sum of the radii then
    % there is no way the spheres overlap.
    volume = 0;
elseif ( r == R )  
    % use the simpler formula if r = R (and isn't prone to dividing by zero
    % errors when d = 0)
    volume = 1/12 * pi * (4*R + d) * (4*R - d)^2;
else
    volume = pi * (R+r-d)^2 * (d^2 + 2*d*r - 3*r^2 + 2*d*R + 6*r*R - 3*R^2) / (12 * d);
end


end