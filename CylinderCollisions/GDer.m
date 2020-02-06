function [gradient] = GDer(s, t, r0, h0Div2, r1 , h1Div2, a0, b0, c0, a1, b1, c1, lenDelta)
omsmt = 1 - s - t;
ssqr = s*s;
tsqr = t*t;
omsmtsqr = omsmt*omsmt;
temp = ssqr + tsqr + omsmtsqr;
L0 = a0*s + b0*t + c0*omsmt;
L1 = a1*s + b1*t + c1*omsmt;
Q0 = temp - L0*L0;
Q1 = temp - L1*L1;

diffS = s - omsmt;
diffT = t - omsmt;
diffa0c0 = a0 - c0;
diffa1c1 = a1 - c1;
diffb0c0 = b0 - c0;
diffb1c1 = b1 - c1;
halfQ0s = diffS - diffa0c0 * L0;
halfQ1s = diffS - diffa1c1 * L1 ;
halfQ0t = diffT - diffb0c0 * L0;
halfQ1t = diffT - diffb1c1 * L1 ;
factor0 = r0 / sqrt(Q0);
factor1 = r1 / sqrt(Q1);
signL0 = sign(L0);
signL1 = sign(L1);
gradient = [0,0];
gradient(1) = gradient(1) + halfQ0s * factor0;
gradient(1) = gradient(1) + halfQ1s * factor1;
gradient(1) = gradient(1) + h0Div2 * diffa0c0 * signL0;
gradient(1) = gradient(1) + h1Div2 * diffa1c1 * signL1;
gradient(1) = gradient(1) + lenDelta;
gradient(2) = gradient(2) + halfQ0t * factor0;
gradient(2) = gradient(2) + halfQ1t * factor1;
gradient(2) = gradient(2) + h0Div2 * diffb0c0 * signL0;
gradient(2) = gradient(2) + h1Div2 * diffb1c1 * signL1;
gradient(2) = gradient(2) + lenDelta;

% val0 = gradient(1);
% val1 = gradient(2);
gradient;
end