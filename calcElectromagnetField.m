function [] = calcElectromagnetField(figureNum)

uiimport();
disp("Choose ANSYS data file for analysis (XLocationm, YLocationm, and DirectionalMagneticFluxDensityT). Press any key to continue ... ");
pause();

figure(figureNum)
hold on
maxB = max(DirectionalMagneticFluxDensityT);
minB = min(DirectionalMagneticFluxDensityT);
% c = zeros(length(DirectionalMagneticFluxDensityT),1);
c = DirectionalMagneticFluxDensityT; % color is directly proportional to height
sz = 75;

% for n = 1:length(DirectionalMagneticFluxDensityT)
%     norm_val = (DirectionalMagneticFluxDensityT(n) - minB)/(maxB - minB); % normalizes values from 0 to 1
%     c(n) = norm_val*9+1;
%     if (norm_val < 0.5)
% %         plot3(XLocationm(n),YLocationm(n),DirectionalMagneticFluxDensityT(n),'.', "Color", [1 norm_val*2 0], "MarkerSize", 50 );
%          
%     else
% %         plot3(XLocationm(n),YLocationm(n),DirectionalMagneticFluxDensityT(n),'.', "Color", [(2-2*norm_val) 1 0], "MarkerSize", 50 ); 
%     end
% end

% (DirectionalMagneticFluxDensityT - minB)/(maxB - minB); % normalizes values from 0 to 1
scatter3(XLocationm, YLocationm, DirectionalMagneticFluxDensityT, sz, c, 'filled');
colorbar();
% surf(XLocationm,YLocationm,DirectionalMagneticFluxDensityT)
plot([0.03 0.03 -0.03 -0.03 0.03],[0 0.21 0.21 0 0],"Color", [0 0 0])
pbaspect([1 1 1])
title("Y component of the B-field produced by an Electromagnet");
xlabel("X [m]");
ylabel("Y [m]");
zlabel("By [T]");
hold off
figureNum = figureNum +1;

% INTERPOLATION CODE %
numPoints = 501;
inter_Points = linspace(0,0.45,numPoints);

%                    z1 z2 z3 z4 Interpolated z
inter_BField = zeros(numPoints, 5);
% Find closest 3 to 5 points

%                    x1 y1 dist1 x2 y2 dist2 x3 y3 dist3 x4 y4 dist4
nearest_Points = ones(numPoints, 12);
for i = 1:numPoints
    % Initialize first 4 points
    tempDist1 = sqrt( ( 0-XLocationm(1) )^2 + ( inter_Points(i)-YLocationm(1) )^2 );
    tempDist2 = sqrt( ( 0-XLocationm(2) )^2 + ( inter_Points(i)-YLocationm(2) )^2 );
    tempDist3 = sqrt( ( 0-XLocationm(3) )^2 + ( inter_Points(i)-YLocationm(3) )^2 );
    tempDist4 = sqrt( ( 0-XLocationm(4) )^2 + ( inter_Points(i)-YLocationm(4) )^2 );
    nearest_Points(i,1) = XLocationm(1); % x1
    nearest_Points(i,2) = YLocationm(1); % y1
    nearest_Points(i,3) = tempDist1;     % dist1
    nearest_Points(i,4) = XLocationm(2); % x2
    nearest_Points(i,5) = YLocationm(2); % y2
    nearest_Points(i,6) = tempDist2;     % dist2
    nearest_Points(i,7) = XLocationm(3); % x3
    nearest_Points(i,8) = YLocationm(3); % y3
    nearest_Points(i,9) = tempDist3;     % dist3
    nearest_Points(i,10) = XLocationm(4);% x4 
    nearest_Points(i,11) = YLocationm(4);% y4
    nearest_Points(i,12) = tempDist4;    % dist4
    inter_BField(i,:) = [DirectionalMagneticFluxDensityT(1), DirectionalMagneticFluxDensityT(2), DirectionalMagneticFluxDensityT(3), DirectionalMagneticFluxDensityT(4), 0];
    for j = 5:length(XLocationm)  
        % Calc distance from test point to measurement point
        tempDist = sqrt( ( 0-XLocationm(j) )^2 + ( inter_Points(i)-YLocationm(j) )^2 );
        % Check if distance is less than other calcs and store result
        if (tempDist < nearest_Points(i,12)) || (tempDist < nearest_Points(i,9)) || (tempDist < nearest_Points(i,6)) || (tempDist < nearest_Points(i,3))
            % If the test distance is smaller than any previous test
            % distance than the largest test distance should be replaced.
            if (nearest_Points(i,12) > nearest_Points(i,9)) && (nearest_Points(i,12) > nearest_Points(i,6)) && (nearest_Points(i,12) > nearest_Points(i,3))
                nearest_Points(i,10) = XLocationm(j);
                nearest_Points(i,11) = YLocationm(j);
                nearest_Points(i,12) = tempDist; 
                inter_BField(i,4) = DirectionalMagneticFluxDensityT(j);
            elseif (nearest_Points(i,9) > nearest_Points(i,12)) && (nearest_Points(i,9) > nearest_Points(i,6)) && (nearest_Points(i,9) > nearest_Points(i,3))
                nearest_Points(i,7) = XLocationm(j);
                nearest_Points(i,8) = YLocationm(j);
                nearest_Points(i,9) = tempDist; 
                inter_BField(i,3) = DirectionalMagneticFluxDensityT(j);
            elseif (nearest_Points(i,6) > nearest_Points(i,12)) && (nearest_Points(i,6) > nearest_Points(i,9)) && (nearest_Points(i,6) > nearest_Points(i,3))
                nearest_Points(i,4) = XLocationm(j);
                nearest_Points(i,5) = YLocationm(j);
                nearest_Points(i,6) = tempDist;
                inter_BField(i,2) = DirectionalMagneticFluxDensityT(j);                
            elseif (nearest_Points(i,3) > nearest_Points(i,12)) && (nearest_Points(i,3) > nearest_Points(i,6)) && (nearest_Points(i,3) > nearest_Points(i,9))
                nearest_Points(i,1) = XLocationm(j);
                nearest_Points(i,2) = YLocationm(j);
                nearest_Points(i,3) = tempDist;
                inter_BField(i,1) = DirectionalMagneticFluxDensityT(j);                
            end
        end
    end
    disp("Point " + i + " done.");
    % Check that the 3 points surround the point of interest
    % https://www.geeksforgeeks.orgcheck-whether-a-given-point-lies-inside-a-triangle-or-not/
    if (nearest_Points(i,12) >= nearest_Points(i,9)) && (nearest_Points(i,12) >= nearest_Points(i,6)) && (nearest_Points(i,12) >= nearest_Points(i,3))
        % 4th point is furthest away
        % A = [x1(y2 – y3) + x2(y3 – y1) + x3(y1-y2)]/2
        A = abs( nearest_Points(i,1)*(nearest_Points(i,5)-nearest_Points(i,8)) + ... 
                 nearest_Points(i,4)*(nearest_Points(i,8)-nearest_Points(i,2)) + ...
                 nearest_Points(i,7)*(nearest_Points(i,2)-nearest_Points(i,5)) ) / 2;
        % Calculate individual triangles     
        A1 = abs( 0*(nearest_Points(i,5)-nearest_Points(i,8)) + ... 
                  nearest_Points(i,4)*(nearest_Points(i,8)-inter_Points(i)) + ...
                  nearest_Points(i,7)*(inter_Points(i)-nearest_Points(i,5)) ) / 2; 
        A2 = abs( nearest_Points(i,1)*(inter_Points(i)-nearest_Points(i,8)) + ... 
                  0*(nearest_Points(i,8)-nearest_Points(i,2)) + ...
                  nearest_Points(i,7)*(nearest_Points(i,2)-inter_Points(i)) ) / 2;
        A3 = abs( nearest_Points(i,1)*(nearest_Points(i,5)-inter_Points(i)) + ... 
                  nearest_Points(i,4)*(inter_Points(i)-nearest_Points(i,2)) + ...
                  0*(nearest_Points(i,2)-nearest_Points(i,5)) ) / 2;
        if (A1+A2+A3 == A)
            disp("Point is inside");
        else
            disp("Point is outside");
        end
            
    elseif (nearest_Points(i,9) > nearest_Points(i,12)) && (nearest_Points(i,9) > nearest_Points(i,6)) && (nearest_Points(i,9) > nearest_Points(i,3))
        % 3rd point is furthest away 

    elseif (nearest_Points(i,6) > nearest_Points(i,12)) && (nearest_Points(i,6) > nearest_Points(i,9)) && (nearest_Points(i,6) > nearest_Points(i,3))
        nearest_Points(i,4) = XLocationm(j);
        nearest_Points(i,5) = YLocationm(j);
        nearest_Points(i,6) = tempDist; 
    elseif (nearest_Points(i,3) > nearest_Points(i,12)) && (nearest_Points(i,3) > nearest_Points(i,6)) && (nearest_Points(i,3) > nearest_Points(i,9))
        nearest_Points(i,1) = XLocationm(j);
        nearest_Points(i,2) = YLocationm(j);
        nearest_Points(i,3) = tempDist;         
    end
end

% EXTRACTION OF AXIAL MAGNETIC FIELD VALUES
X = [];
Y = [];
MagneticField_y_Dir = [];
for n=1:length(XLocationm)
    if (XLocationm(n) == 0) && (YLocationm(n) >= 0.21)
       X = [X  XLocationm(n)];
       Y = [Y  YLocationm(n)];
       MagneticField_y_Dir = [MagneticField_y_Dir DirectionalMagneticFluxDensityT(n) ];
    end
end

figure(figureNum)
hold on
plot(Y,MagneticField_y_Dir,'o')
title("Y component of the B-field produced by an Electromagnet");
xlabel("Y [m]");
ylabel("By [mT]");
hold off


end