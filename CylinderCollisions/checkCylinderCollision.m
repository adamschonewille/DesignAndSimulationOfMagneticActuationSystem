function [cond] = checkCylinderCollision(centerPoint0, W0, r0, h0, ...
                                         centerPoint1, W1, r1, h1)
%checkCylinderCollision
%% Parameters:
% centerPoint0 - [x y z] center of mass of a uniformly distributed cylinder
% W0 - unit vector direction of the cylinder axis of rotational symmetry
% r0 - cylinder radius
% h0 - overall length of the cylinder
%% Output:
% cond - ? bool

%% Code:

delta = centerPoint1 - centerPoint0;
W0xW1 = cross(W0,W1);
lenW0xW1 = norm(W0xW1);
h0Div2 = h0 / 2;
h1Div2 = h1 / 2;
rSum = r0 + r1;
true = 1;
false = 0;
cond = true;

if ( lenW0xW1 > 0 ) 
    % Test for separation by W0.
    if ( (r1 * lenW0xW1 + h0Div2 + h1Div2 * norm(dot(W0,W1)) - norm(dot(W0, delta)) ) < 0 ) 
        cond = false;
        return;
    end
    % Test for separration by W1.
    if ( (r0 * lenW0xW1 + h0Div2 * norm(dot(W0,W1)) + h1Div2 - norm(dot(W1, delta)) ) < 0 ) 
        cond = false;
        return;
    end
    % Test for separation by W0xW1.
    if ( (rSum*lenW0xW1 - norm(dot(W0xW1, delta)) ) < 0 )
        cond = false;
        return;
    end
    % Test for separation by directions perpendicular to W0.
    if ( separatedByCylinderPerpendiculars(centerPoint0 ,W0, r0 , h0 , centerPoint1 ,W1, r1 , h1) ) 
        cond = false;
        return;
    end
    % Test for separation by directions perpendicular to W1.
    if ( separatedByCylinderPerpendiculars(centerPoint1 ,W1, r1 , h1 , centerPoint0 ,W0, r0 , h0 ) )
        cond = false;
        return;
    end
    % Test for separation by other directions.
    if ( separatedByOtherDirections(W0, r0 , h0 ,W1, r1 , h1 , delta ) )
        cond = false;
        return;
    end
    
else
    % Test for separation by height.
    if ( ( h0Div2 + h1Div2 - norm(dot(W0, delta)) ) < 0 ) 
        cond = false;
        return;
    end
    % Test for separation radially.
    if ( ( rSum - norm( delta - dot(W0, delta)*W0 ) ) < 0 ) 
        cond = false;
        return;
    end
    % If parallel cylinders are not separated by height or radial distance,
    % then the cylinders must overlap.
end

% If none of the functions return false, then cond will still be true.

end
