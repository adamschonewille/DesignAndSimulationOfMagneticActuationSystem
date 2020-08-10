function [] = plotEM(nomDist, large_EM, small_EM, mAct_sph, pAct, fval, figureNum)

% Define large cylinder for plotting
[X_l,Y_l,Z_l] = cylinder();
Z_l = -Z_l.*large_EM(3) + large_EM(3)/2;
X_l = X_l.*large_EM(2)/2;
Y_l = Y_l.*large_EM(2)/2;
% Define small cylinder for plotting
[X_s,Y_s,Z_s] = cylinder();
Z_s = -Z_s.*small_EM(3) + small_EM(3)/2;
X_s = X_s.*small_EM(2)/2;
Y_s = Y_s.*small_EM(2)/2;

%% Plotting 
% Plot cylinders after rotating and translating
figure(figureNum)
for l=1:4
    % 3 x 3 Transformation (Rotation Only) Matrix 
    RotM = rotz(deg2rad(mAct_sph(1,l)))*roty(deg2rad(mAct_sph(2,l)));
    % Rotated cylinder points:
    % | RotM   0  |     | X1      | X2
    % |           |  x  | Y1  ... | Y2 
    % |  0   RotM |     | Z1      | Z2
    temp_cylinder = [RotM zeros(size(RotM));zeros(size(RotM)) RotM] * [X_l(1,:);Y_l(1,:);Z_l(1,:);X_l(2,:);Y_l(2,:);Z_l(2,:)];
    X_mod = [temp_cylinder(1,:);temp_cylinder(4,:)] + pAct(1,l);
    Y_mod = [temp_cylinder(2,:);temp_cylinder(5,:)] + pAct(2,l);
    Z_mod = [temp_cylinder(3,:);temp_cylinder(6,:)] + pAct(3,l);
    
    X_mod = real(X_mod);
    Y_mod = real(Y_mod);
    Z_mod = real(Z_mod);
    
    surf(X_mod,Y_mod,Z_mod) 
    hold on
end
for l=5:8
        % 3 x 3 Transformation (Rotation Only) Matrix 
    RotM = rotz(deg2rad(mAct_sph(1,l)))*roty(deg2rad(mAct_sph(2,l)));
    % Rotated cylinder points:
    % | RotM   0  |     | X1      | X2
    % |           |  x  | Y1  ... | Y2 
    % |  0   RotM |     | Z1      | Z2
    temp_cylinder = [RotM zeros(size(RotM));zeros(size(RotM)) RotM] * [X_s(1,:);Y_s(1,:);Z_s(1,:);X_s(2,:);Y_s(2,:);Z_s(2,:)];
    X_mod = [temp_cylinder(1,:);temp_cylinder(4,:)] + pAct(1,l);
    Y_mod = [temp_cylinder(2,:);temp_cylinder(5,:)] + pAct(2,l);
    Z_mod = [temp_cylinder(3,:);temp_cylinder(6,:)] + pAct(3,l);
    
    X_mod = real(X_mod);
    Y_mod = real(Y_mod);
    Z_mod = real(Z_mod);
    
    hold on
    surf(X_mod,Y_mod,Z_mod) 
end
[X,Y,Z] = sphere();
X = X.*nomDist;
Y = Y.*nomDist;
Z = Z.*nomDist;
surf(X,Y,Z)
xlabel("X")
ylabel("Y")
zlabel("Z")
title("Plot actuator positions with a fitness of: " + num2str(fval))
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5,0.12])
axis equal
hold off

end