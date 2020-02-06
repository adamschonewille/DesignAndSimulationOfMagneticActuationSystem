function [cond] = separatedByOtherDirections(W0, r0, h0, ...
                                             W1, r1, h1, delta)
% Minimize G(s,t) subject to s >= 0, t >= 0, and s + t <= 1 . If at any
% iterate you find a value for which G <= 0, return ’true’.
% If no separating directions have been found at the end of the 
% minimization, return ’false’.

% s0 = rand(1);
% t0 = rand(1);
s0 = 0.5;
t0 = 0.5;


% define x0
% fittingParameters = [mAct_sph; pAct_cartesion];
% Constants = [nomDist; large_EM; small_EM];

h0Div2 = h0 / 2;
h1Div2 = h1 / 2;
% delta = centerPoint1 - centerPoint0;
lenDelta = norm(delta);
W0_hat = W0./norm(W0);
W1_hat = W1./norm(W1);
a0 = W0_hat(1);
b0 = W0_hat(2);
c0 = W0_hat(3);
a1 = W1_hat(1);
b1 = W1_hat(2);
c1 = W1_hat(3);

%% Optimize via fmincon
% 
% fun_opt = @(x) G(x(1), x(2), r0, h0Div2, r1, h1Div2, a0, b0, c0, a1, b1, c1, lenDelta);
% fun_con = @(x) cylinderConstraints(x); % Just return 0 always for c, and ceq
% %  s >= 0,  t >= 0, and s + t <= 1
% % -s <= 0, -t <= 0, and 
% % Minimize function G(s,t)
% options = optimset('Display','off','TolX',1e-4,'MaxIter',10e4); % display = none or off
% % options = optimset('PlotFcns','optimplotfval','TolX',1e-3,'MaxIter',10e4); % display = none or off
% A = []; %A = [1,1];
% b = []; %b = [1];
% Aeq = [];
% beq = [];
% lb = [0 0];
% ub = [1,1];
% % nonlcon = @fun_con;
% [inputs,G_val] = fmincon(fun_opt,[s0,t0],A,b,Aeq,beq,lb,ub,fun_con,options);%,nonlcon,options);


%% Optimize via own algorithm:
G_val = findGMinimum(s0, t0, r0, h0Div2, r1, h1Div2, a0, b0, c0, a1, b1, c1, lenDelta);



%% Check condition and return appropriate value
% If resulting global minima (G is convex apparently) is less than or = 0
% then collision is "true"
if (G_val <= 0)
    cond = true;
else
    cond = false;
end
% s, t, r0, h0Div2, r1 , h1Div2, a0, b0, c0, a1, b1, c1, lenDelta

end