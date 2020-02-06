function [pAct,RzyAct,axesAct] = setup_parameters(pAct_sph, mAct_sph)
% Inputs:
% pAct_sph = actuator magnet positions in spherical coordinates [radius(m); azimuth(deg); inclination (deg)];
% rAct_sph = actuator magnet rotational axes in spherical coordinates [azimuth(deg); inclination (deg)];

% Outputs:
% pAct: xyz coordinates of actuator magnet centers (m)
% RzyAct: ZY Euler angle rotation matrix for each actuator magnet
% axesAct: direction of axis of rotation in xyz 

%[pAct,R,axesAct] = setup_constants(dist, az, incl, choose, rot_ax, actM_mag) 

min_r = 0.12;
d_large = 2*(45 + 22.5) / 1000;
d_small = 2*(30 + 15) / 1000;

nAct=size(pAct_sph,2); %number of magnets

%% Magnet locations
pAct=[
    pAct_sph(1,:).*sind(pAct_sph(3,:)).*cosd(pAct_sph(2,:));
    pAct_sph(1,:).*sind(pAct_sph(3,:)).*sind(pAct_sph(2,:));
    pAct_sph(1,:).*cosd(pAct_sph(3,:))
    ];

for i=1:nAct
    if (pAct(3,i) > -min_r-d_large/2*sind())
              
    end
    
end

%% Axes & constants

axesAct=[
    sind(rAct_sph(2,:)).*cosd(rAct_sph(1,:));
    sind(rAct_sph(2,:)).*sind(rAct_sph(1,:));
    cosd(rAct_sph(2,:))
    ];    

RzyAct = zeros(3,2,nAct);
for i=1:nAct 
    RzyAct(:,:,i) = [
        sind(rAct_sph(1,i)) cosd(rAct_sph(2,i))*cosd(rAct_sph(1,i));
       -cosd(rAct_sph(1,i)) cosd(rAct_sph(2,i))*sind(rAct_sph(1,i));
        0                  -sind(rAct_sph(2,i));
        ];        
end

