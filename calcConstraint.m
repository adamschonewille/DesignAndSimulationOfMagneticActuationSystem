function [c, ceq] = calcConstraint(fittingParameters, Constants)
% [c,ceq]
% c = ...     % Compute nonlinear inequalities at x.
% c(x) <= 0 for all entries of c.
%
% ceq = ...   % Compute nonlinear equalities at x.


% For Genetic Algorithm:
% input matrix is changed into a vector mx1

isVec = false;
if ( (size(fittingParameters,1) == 1) || (size(fittingParameters,2) == 1) )
    % if the input is a vector then reorganize into a matrix.
    fittingParameters = reshape(fittingParameters,8,[])' ;
    isVec = true;
end

minDist = Constants(1);
d_large = Constants(3);
d_small = Constants(6);
l_large = Constants(4);
l_small = Constants(7);

nAct = size(fittingParameters,2);

mAct_sph = rad2deg(fittingParameters(1:2,:));
RzyAct = zeros(3,nAct);
for i=1:nAct 
    RzyAct(:,i) = [ cosd(mAct_sph(1,i))*sind(mAct_sph(2,i));
                    sind(mAct_sph(1,i))*sind(mAct_sph(2,i));
                    cosd(mAct_sph(2,i)) ];  %-cosd(mAct_sph(2,i)) ];       
end
W = RzyAct;
EM_size = [l_large/2 l_large/2 l_large/2 l_large/2 l_small/2 l_small/2 l_small/2 l_small/2];
EM_size = [EM_size;EM_size;EM_size];
% This is to shift the location of the dipole for a dipole located at the 
% face of the EM instead of the center
% Cn = fittingParameters(3:5,:) - EM_size.*W; 

Cn = fittingParameters(3:5,:);

c = zeros(1,nAct);
ceq = c;

extra_c = [];
count = 1;
for i = 1:length(c)
    % Ensure that all the EM are below the minDist plane regardless of
    % angle:
    if (i<=4)
        c(i) = minDist + fittingParameters(5,i) + d_large/2 * sin(fittingParameters(2,i));
    else
        c(i) = minDist + fittingParameters(5,i) + d_small/2 * sin(fittingParameters(2,i));
    end
    
    if (i <= length(c)-1)
        for j = i+1:length(c)
            % Only check if the faces of EM are spaced far enough apart:
%             if (i <= 4) && (j <= 4)
%                 extra_c(count) = d_large - norm(fittingParameters(3:5,i)-fittingParameters(3:5,j) );
%                 count = count+1;
%             elseif (i <= 4) && (j > 4)
%                 extra_c(count) = d_large/2 + d_small/2 - norm(fittingParameters(3:5,i)-fittingParameters(3:5,j) );
%                 count = count+1;
%             elseif (i > 4) && (j > 4)
%                 extra_c(count) = d_small - norm(fittingParameters(3:5,i)-fittingParameters(3:5,j) );
%                 count = count+1;
%             end
            % Check if EMs collide with one another
            if (i <= 4) && (j <= 4)
                extra_c(count) = checkCylinderCollision(Cn(:,i), W(:,i), d_large/2, l_large, ...
                                                        Cn(:,j), W(:,j), d_large/2, l_large);
                count = count+1;
            elseif (i <= 4) && (j > 4)
                extra_c(count) = checkCylinderCollision(Cn(:,i), W(:,i), d_large/2, l_large, ...
                                                        Cn(:,j), W(:,j), d_small/2, l_small);
                count = count+1;
            elseif (i > 4) && (j > 4)
                extra_c(count) = checkCylinderCollision(Cn(:,i), W(:,i), d_small/2, l_small, ...
                                                        Cn(:,j), W(:,j), d_small/2, l_small);
                count = count+1;
            end
        end
    end
end

% original constraint function
% c = [c, 10*extra_c];
% ceq = zeros(size(c));

% new constraint function, the collisions must be = 0 (false)

ceq = 10*extra_c;
% zeros(size(c));
c = [c, zeros(1,size(ceq,2)-size(c,2))];

if (isVec)
    c = reshape(c',1,[]);
    ceq = reshape(ceq',1,[]);    
end

% For ceq, ensure that the values are the same magnitude

% minDist = Constants(1);
% mAct = [Constants(2) Constants(2) Constants(2) Constants(2) ...
%         Constants(5) Constants(5) Constants(5) Constants(5)];
% mAct_sph = rad2deg(fittingParameters(1:2,:));
% pAct = fittingParameters(3:5,:);
% 
% nAct = size(pAct,2);
% %% Axes & constants
% RzyAct = zeros(3,1,nAct); 
% % This is the unit vector of the magnetic moment for Euler angle rotations
% % beta and gamma about Z and Y axes, respectively
% % Assumes that the actuators dipole moment are alligned with Z axis
% % originally
% for i=1:nAct 
%     RzyAct(:,:,i) = [ cosd(mAct_sph(1,i))*sind(mAct_sph(2,i));
%                       sind(mAct_sph(1,i))*sind(mAct_sph(2,i));
%                      -cosd(mAct_sph(2,i)) ];        
% end
% 
% % pAct: xyz coordinates of actuator magnet centers (m)
% % RzyAct: ZY Euler angle rotation matrix for each actuator magnet
% % mAct: magnetic dipole moment of actuators [Am^2]
% pTool = [0,0,0]';
% 
% Bpr = zeros(3,nAct);
% Gpr = zeros(5,nAct);
% 
% for i=1:size(pAct,2)
%     p=-pAct(:,i)+pTool;    
%     p_hat=p/norm(p);
%     
%     KB =   mAct(i)*1e-7 / (norm(p)^3); %note mu/(4*pi) = 1*e-7
%     KG = 3*mAct(i)*1e-7 / (norm(p)^4);
% 
%     Bpr(:,i) = KB*( 3* (p_hat*p_hat') - eye(3) ) * RzyAct(:,:,i);
%     
%     Gp = [
%         3*p_hat(1)-5*p_hat(1)^3         p_hat(2)-5*p_hat(1)^2*p_hat(2)  p_hat(3)-5*p_hat(1)^2*p_hat(3);
%         p_hat(2)-5*p_hat(1)^2*p_hat(2)  p_hat(1)-5*p_hat(1)*p_hat(2)^2  -5*p_hat(1)*p_hat(2)*p_hat(3);
%         p_hat(3)-5*p_hat(1)^2*p_hat(3)  -5*p_hat(1)*p_hat(2)*p_hat(3)   p_hat(1)-5*p_hat(1)*p_hat(3)^2;
%         p_hat(1)-5*p_hat(1)*p_hat(2)^2  3*p_hat(2)-5*p_hat(2)^3         p_hat(3)-5*p_hat(2)^2*p_hat(3);
%         -5*p_hat(1)*p_hat(2)*p_hat(3)   p_hat(3)-5*p_hat(2)^2*p_hat(3)  p_hat(2)-5*p_hat(2)*p_hat(3)^2;
%         ];
% 
%     Gpr(:,i) = KG * Gp * RzyAct(:,:,i);
% end

% [U,S,V] = svd([Bpr;Gpr]);

% Equality Function       
% ceq(1) = abs(1-S(1,1)/S(3,3));

% ceq(2) = 10^6*(8-rank([Bpr;Gpr]));

