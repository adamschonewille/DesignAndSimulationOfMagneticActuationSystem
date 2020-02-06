function [val] = F(t, r0, r1, h1b1Div2, c1sqr, a2, b2)
% Based on the paper obtained from:
% https://www.geometrictools.com/Documentation/IntersectionOfCylinders.pdf
omt = 1 - t ;
tsqr = t * t ;
omtsqr = omt*omt ;
term0 = r0 * sqrt( omtsqr + tsqr ) ;
term1 = r1 * sqrt( omtsqr + c1sqr * tsqr ) ;
term2 = h1b1Div2 * t ;
term3 = norm( omt*a2 + t * b2 );
val = term0 + term1 + term2 - term3 ;
end