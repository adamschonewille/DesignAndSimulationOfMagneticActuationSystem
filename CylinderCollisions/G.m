function [fval, fder] = G(s, t, r0, h0Div2, r1, h1Div2, a0, b0, c0, a1, b1, c1, lenDelta)
omsmt = 1 - s - t;
ssqr = s*s;
tsqr = t*t;
omsmtsqr = omsmt*omsmt;
temp = ssqr + tsqr + omsmtsqr;
L0 = a0*s + b0*t + c0*omsmt;
L1 = a1*s + b1*t + c1*omsmt;
Q0 = temp - L0*L0;
Q1 = temp - L1*L1;
fval = (r0 * sqrt(Q0) + r1 * sqrt(Q1) + h0Div2 * norm(L0) + h1Div2 * norm(L1) - omsmt * lenDelta);

if nargout > 1 % gradient required
    fder = GDer(s, t, r0, h0Div2, r1 , h1Div2, a0, b0, c0, a1, b1, c1, lenDelta);
end


end