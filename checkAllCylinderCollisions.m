function [cond] = checkAllCylinderCollisions( mAct_sph, pAct_cartesion, radiusVec, lengthVec )
%checkAllCylinderCollisions
% Checks all combinations of cylinders until a collision is found. As soon
% as one collision is found the function returns true. If no collisions are
% found then the function returns false.
%% Parameters:
% mAct_sph      - 2 x n matrix containing the azimuthal and polar angles to 
%               define the unit vector of the axis of each n cylinders.
% pAct_cartesion- 3 x n matrix containing the [x y z] cartesion coords of
%               the centerpoint of each n cylinder.
% radiusVec     - n x 1 vector containing the radii of each n cylinder.
% lengthVec     - n x 1 vector containing the overall length of each n cylinder.
%% Output:
% cond - true if there is a collision

%% Code:
% Initialize result as no collisions 
cond = false;
z_norm = [0 0 1]';
W = zeros(size(pAct_cartesion));
numActuators = size(pAct_cartesion, 2);
% set up the unit vectors of the axis for each cylinder.
for n = 1:numActuators
    Rzy = rotz(mAct_sph(1,n)) * roty(mAct_sph(2,n));
    W(:,n) =  Rzy * z_norm;
end
% run the check cylinder collision code for all cylinders
for i = 1:numActuators-1
    for j = i+1:numActuators
        if ( ~~checkCylinderCollision(pAct_cartesion(:,i), W(:,i), radiusVec(i), lengthVec(i), ...
                                     pAct_cartesion(:,j), W(:,j), radiusVec(j), lengthVec(j) ) )
            cond = true;
            break;
        end
       
    end
    
end
                                     
                                     

end