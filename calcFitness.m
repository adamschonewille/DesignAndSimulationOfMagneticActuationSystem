function f = calcFitness(fittingParameters, Constants)
%minimize this function

%% Rewritten by Adam from field_force_gradient_J.m
%fittingParameters = [mAct_sph; pAct_cartesion];
minDist = Constants(1);
m_large = Constants(2);
d_large = Constants(3);
m_small = Constants(5);
d_small = Constants(6);
mAct = [Constants(2) Constants(2) Constants(2) Constants(2) ...
        Constants(5) Constants(5) Constants(5) Constants(5)];
mAct_sph = rad2deg(fittingParameters(1:2,:));
pAct = fittingParameters(3:5,:);

mTool_mag = 8.4e-3; %tool dipole moment magnetiude [Am^2] = [Nm/T]
targetMaxB_Field = 0.100; %[T] or higher
targetMaxForce   = 0.100; %[N] or higher
targetMaxField_G = targetMaxForce/mTool_mag; % [T/m]
DB = 1/targetMaxB_Field*eye(3);
DG = 1/targetMaxField_G*eye(5);
D0 = [DB, zeros(3,5); zeros(5,3), DG];

nAct = size(pAct,2);
% Check that the z position of the actuators is in the correct bound:
% This has been accounted fo rin the constraints

RzyAct = zeros(3,1,nAct); 
% This is the unit vector of the magnetic moment for Euler angle rotations
% beta and gamma about Z and Y axes, respectively
% Assumes that the actuators dipole moment are alligned with Z axis
% originally
for i=1:nAct 
    RzyAct(:,:,i) = [ cosd(mAct_sph(1,i))*sind(mAct_sph(2,i));
                      sind(mAct_sph(1,i))*sind(mAct_sph(2,i));
                     -cosd(mAct_sph(2,i)) ];        
end
% Now we should have variables:
% pAct: xyz coordinates of actuator magnet centers (m)
% RzyAct: ZY Euler angle rotation matrix for each actuator magnet
% mAct: magnetic dipole moment of actuators [Am^2]
pTool = [0,0,0]';
% We want to test a variety of tool positions and assume that if it is
% controllable at extrema, it will also be controllable within extrema
a = 0.03; % bounding cube dist from origin
% Test points are the vertices of a 0.06 m sidelength cube as well as where
% the principle axis intersect the cube surface and origin (27 pts total)
pTool = [0  a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a  ;
         0  0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a ;
         0  0  0  0  0  0  0  0  0  a  a  a  a  a  a  a  a  a -a -a -a -a -a -a -a -a -a ];

% 
% Bpr = zeros(3,nAct);
% Gpr = zeros(5,nAct);
singularValues = zeros(1,size(pTool,2));
for j=1:size(pTool,2)
    Bpr = zeros(3,nAct);
    Gpr = zeros(5,nAct);
    for i=1:size(pAct,2)
        p=-pAct(:,i)+pTool(:,j);    
        p_hat=p/norm(p);

        KB =   mAct(i)*1e-7 / (norm(p)^3); %note mu/(4*pi) = 1*e-7
        KG = 3*mAct(i)*1e-7 / (norm(p)^4);

    %     Bpr(:,2*i-1:2*i) = KB*( 3* (p_hat*p_hat') - eye(3) ) * RzyAct(:,:,i);
        Bpr(:,i) = KB*( 3* (p_hat*p_hat') - eye(3) ) * RzyAct(:,:,i);

        Gp = [
            3*p_hat(1)-5*p_hat(1)^3         p_hat(2)-5*p_hat(1)^2*p_hat(2)  p_hat(3)-5*p_hat(1)^2*p_hat(3);
            p_hat(2)-5*p_hat(1)^2*p_hat(2)  p_hat(1)-5*p_hat(1)*p_hat(2)^2  -5*p_hat(1)*p_hat(2)*p_hat(3);
            p_hat(3)-5*p_hat(1)^2*p_hat(3)  -5*p_hat(1)*p_hat(2)*p_hat(3)   p_hat(1)-5*p_hat(1)*p_hat(3)^2;
            p_hat(1)-5*p_hat(1)*p_hat(2)^2  3*p_hat(2)-5*p_hat(2)^3         p_hat(3)-5*p_hat(2)^2*p_hat(3);
            -5*p_hat(1)*p_hat(2)*p_hat(3)   p_hat(3)-5*p_hat(2)^2*p_hat(3)  p_hat(2)-5*p_hat(2)*p_hat(3)^2;
            ];

    %     Gpr(:,2*i-1:2*i) = KG * Gp * RzyAct(:,:,i);
        Gpr(:,i) = KG * Gp * RzyAct(:,:,i);
    end

    K = [Bpr;Gpr];
    K0 = D0*K;
    [U,S,V] = svd(K0);
%     [U,S,V] = svd(Bpr);
%     B = sum(Bpr,2)
%     G = sum(Gpr,2);
    max3values = zeros(2,3);
    for m=1:size(U,2)
        magnitude = norm(U(1:3,m));
        if (magnitude > max3values(1,1) )
            if (magnitude > max3values(1,2) )
                if (magnitude > max3values(1,3) )                    
                    % shift all over
                    max3values(1,1) = max3values(1,2);
                    max3values(1,2) = max3values(1,3);
                    max3values(1,3) = magnitude;
                    % shift indices too
                    max3values(2,1) = max3values(2,2);
                    max3values(2,2) = max3values(2,3);
                    max3values(2,3) = m;
                else
                    % shift all but 3
                    max3values(1,1) = max3values(1,2);
                    max3values(1,2) = magnitude;
                    % shift indices too
                    max3values(2,1) = max3values(2,2);
                    max3values(2,2) = m;
                end
            else
                max3values(1,1) = magnitude;
                max3values(2,1) = m;
            end
        end
                    
    end
    U
    maxSV = S(max3values(2,3),max3values(2,3));
    minSV = S(max3values(2,1),max3values(2,1));
    singularValues(j) = 1/minSV;% + maxSV/minSV;
end
    
% G_reshape = [G(1:3)';
%              G(2)  G(4)  G(5);
%              G(3)  G(5) -G(1)-G(4) ];

% Fitness / Cost Function         
% f =  K*norm(B-[0.1 0.1 0.1]')^2+(1-K)*norm(F)^2 + 1/abs(S(3,3)) + S(1,1)/S(3,3);
% f =  1/abs(S(3,3)) + abs(1-S(1,1)/S(3,3)); %norm(B-[0.25 0.25 0.25]')^2 + 

f = sum(singularValues);
%Constants = [nomDist; large_EM; small_EM];
% plotEM(minDist, Constants(2:4), Constants(5:7), mAct_sph, pAct, f, 10);

