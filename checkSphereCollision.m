function [cond] = checkSphereCollision(pAct_cartesion, radiusVec)
%checkSphereCollision
% Returns true if there is a collision between any two spheres. A collision
% is determined to exist if the distance between two electromagnets centers
% are less than the sum of the radii of their bounding spheres.
%% Parameters:
% pAct_cartesion - 3 x n matrix containing the cartesion coords of n
%               actuators
% radiusVec - 1 x n vector containing the radii of the spheres that bound
%               the actuators
%% Output:
% cond - true if a collision exists and false otherwise

%% Code:
% initialize:
cond = false;
numActuators = size(pAct_cartesion,2);

% Check that there is more than 1 Actuator:
if (numActuators == 1)
    return
end

% search through all actuator pairs and measure the distance between     
for i = 1:numActuators-1 
    for j = i+1:numActuators
        % determine if the distance between sphere origins is less than the
        % sum of the radius distances of each sphere:
        if ( norm(pAct_cartesion(:,i)-pAct_cartesion(:,j)) < (radiusVec(i)+radiusVec(j)) )
            cond = true;
            break;
        end
    end
end

end