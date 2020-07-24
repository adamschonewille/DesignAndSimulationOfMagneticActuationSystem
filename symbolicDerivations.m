% Derivations
clear all; close all; clc;

syms x y z mx my mz;
p = [x; y; z];
p_hat = 1/norm(p)*p;
m = [mx; my; mz];

alpha = 1e-7;

B = alpha / norm(p)^3 * ( 3*p_hat*( dot(p_hat,m) ) - m )

Bxx = diff( B(1,1), x)
Bxy = diff( B(1,1), y)
Bxz = diff( B(1,1), z)

Byy = diff( B(2,1), y)
Byz = diff( B(2,1), z)

-(mx - (3*x*((mx*conj(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(1/2) + (my*conj(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(1/2) + (mz*conj(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(1/2)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(1/2))/(10000000*(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2))

