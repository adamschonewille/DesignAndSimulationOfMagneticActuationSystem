function [cond] = separatedByCylinderPerpendiculars(centerPoint0, W0, r0, h0, ...
                                                    centerPoint1, W1, r1, h1)
delta = centerPoint1 - centerPoint0;
c1 = dot(W0,W1) ;
b1 = sqrt( 1 - c1*c1 );
V0 = (W1 - c1*W0) / b1;
U0 = cross(V0 , W0);
h1b1Div2 = h1*b1 / 2;
c1sqr = c1*c1;
a2 = dot( delta, U0 );
b2 = dot( delta, V0 );

cond = false;

% Test directions ( 1 - t )*U0 + t * V0
if ( F(0, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
    % U0 is a separating direction
    cond = true;
    return;
end

if ( F(1, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
    % V0 is a separating direction
    cond = true;
    return;
end

if( FDer(0, r0, r1, h1b1Div2, c1sqr, a2, b2) >= 0 )
    % no separation by perpendicular directions
    cond = false;
    return;
end

if ( FDer(1, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
    % no separation by perpendicular directions
    cond = false;
    return;
end                                

% Next checks

% Use bisection to locatie t-bar for which F(t-bar) is a minimum. The upper
% bound maxIterations may be chosen to guarantee a specified number of digits
% of precision in the t?variable.
maxIterations = 1000;
fd0 = 0; fd1 = 0; tmid = 0; fdmid = 0;
t0 = 0;
t1 = 1;

for i=0:maxIterations  
    tmid = 0.5 * ( t0 + t1 );
    if ( F(tmid, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
        % ( 1 - t )*U0 + t * V0 is a se parating direction
        cond = true;
        return;
    end
    fdmid = FDer(tmid, r0, r1, h1b1Div2, c1sqr, a2, b2);
    if ( fdmid > 0 )
        t1 = tmid ;
    elseif ( fdmid < 0 )
        t0 = tmid ;
    else
        break ;
    end
end

% Test directions ( 1-t )*(-U0) + t * V0.
a2 = -a2;
if ( F(0, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
    % U0 is a separating direction
    cond = true;
    return;
end

if ( F(1, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
    %V0 is a separating direction
    cond = true; 
    return;
end

if ( FDer(0, r0, r1, h1b1Div2, c1sqr, a2, b2) >= 0 )
    % no separation by perpendicular directions
    cond = false; 
    return;
end

if ( FDer(1, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
    % no separation by perpendicular directions
    cond = false;
    return;
end

% Use bisection to locatie t-bar for which F(t-bar) is a minimum. The upper
% bound maxIterations may be chosen to guarantee a specified number of digits
% of precision in the t?variable.
t0 = 0 ;
t1 = 1 ;

for i = 0:maxIterations
    tmid = 0.5 * ( t0 + t1 );
    if ( F(tmid, r0, r1, h1b1Div2, c1sqr, a2, b2) <= 0 )
        % ( 1-t )*U0 + t * V0 is a separating direction
        cond = true;
        return;
    end
    fdmid = FDer(tmid, r0, r1, h1b1Div2, c1sqr, a2, b2);
    if ( fdmid > 0 )
        t1 = tmid;
    elseif ( fdmid < 0 )
        t0 = tmid;
    else
        break;
    end
end


end