function [Gmin] = findGMinimum(s, t, r0, h0Div2, r1, h1Div2, a0, b0, c0, a1, b1, c1, lenDelta)

[Gmin, Ggrad] = G(s, t, r0, h0Div2, r1, h1Div2, a0, b0, c0, a1, b1, c1, lenDelta);

gamma = 0.001; % Step size (slow convergence)

a_0 = [s;t];
a_n = a_0 + gamma*Ggrad;
maxiter = 1000;
Gprev = Gmin;

tolerance = 1e-3;

for n=1:maxiter
    % Ensure that the new values are within bounds
    if ( a_n(1) < 0 )
        a_n(1) = 0;
    end
    
    if ( a_n(2) < 0 )
        a_n(2) = 0;
    end
    
    if ( (a_n(1) + a_n(2)) > 1 )
        a_n(1) = a_n(1) / (a_n(1) + a_n(2));
        a_n(2) = a_n(2) / (a_n(1) + a_n(2));
    end
    
    % Update the gradient and current G value
    [Gmin, Ggrad] = G(a_n(1), a_n(2), r0, h0Div2, r1, h1Div2, a0, b0, c0, a1, b1, c1, lenDelta);
    % Update next parameters
    a_n = a_0 + gamma*Ggrad;    
    % Check if Gmin falls below zero for the criterion of cylinder
    % collision
    if (Gmin <=0)
        return;
    end
    
    % Check convergence to a tolerance
    if ( abs(Gprev-Gmin) < tolerance )
        return;
    end
    Gprev = Gmin;
    
end

disp("Max Iterations Reached")

end