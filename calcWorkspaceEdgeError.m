function [Error_field, Error_grad] = calcWorkspaceEdgeError(mAct_sph, pAct_cartesion, EM_dim, a, I)
%calcWorkspaceEdgeError
% INPUTS
% I         - The nx1 vector of input currents for the electromagnets that  
%           are present.
% a         - The distance from the origin to the face of a bounding cube
% EM_dim    - a vector containing the electromagnet properties in the form
%           of: EM_dim = [m, Diameter, Length]';
% mAct_sph  - A 2xn matrix containing the orientation information [rad] of the
%           EMs including beta (the azimuthal angle about the z-axis) and 
%           gamma (the polar angle from the z-axis).
% pAct_cartesion  - A 3xn matrix containing the cartesion coordinates [m] of
%                 each electromagnet
% OUTPUTS:
% Error_field   - The 2 x 26 (9+9+8 points) matrix containing the magnitude
%               error in the first row and the directional error in the
%               second row for the magnetic field components at each test
%               point relative to the field at the center of the workspace
%               for the same input current. These values are given as
%               percent errors of the desired value.
% Error_grad    - The 5 x 26 matrix containing the difference in magnetic 
%               gradient field components at each test point relative to 
%               the gradient components at the center of the workspace
%               for the same input current. ( B_probe - B_bench )


%% Isotropic or Uniformity Penalty functions



%% Field at the center of the workspace
% First find the field a the center of the workspace that will act as the
% benchmark for the error.
B_bench = zeros( 8,1 );
numActuators = size(pAct_cartesion,2);
O = [0;0;0]; % Origin where the dipole field is evaluated at
z_hat = [0;0;1]; % the z-direction unit vector.
m_mag = EM_dim(1);
B_dip = zeros(3,numActuators);
B_dip_G = zeros(5,numActuators);
for i = 1:numActuators
    % find EM magnetic moment axis using the orientation angles 
    Rzy = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat; 
    % Calc magnetic field and gradient control matrices
    B_dip(:,i) = dipoleField( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy);
    B_dip_G(:,i) = dipoleGradient( ( O - (pAct_cartesion(:,i)) ), m_mag*Rzy );
end
B_bench = [B_dip;B_dip_G]*I;

%% Find the magnetic field at all test points on the edges of the workspace
% a = 0.03; % bounding cube dist from origin
% Test points are the vertices of a 0.06 m sidelength cube as well as where
% the principle axis intersect the cube surface and origin (27 pts total)
pTool_test = [a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a  0  a  a  0 -a -a -a  0  a ;
              0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a  0  0  a  a  a  0 -a -a -a ;
              0  0  0  0  0  0  0  0  a  a  a  a  a  a  a  a  a -a -a -a -a -a -a -a -a -a ];
pTool = [0 0 0]';
% The magnetic field will be probed at all 'tool test locations', located
% at the edges of the workspace described by a bounding cube.
B_probe = zeros( 8,size(pTool_test,2) );
for j = 1:size(pTool_test,2)
    B_dip = zeros(3,numActuators);
    B_dip_G = zeros(5,numActuators);
    for i = 1:numActuators
        % find EM magnetic moment axis using the orientation angles 
        Rzy = rotz(mAct_sph(1,i))*roty(mAct_sph(2,i))*z_hat; 
        % Calc magnetic field and add to previous magnetic field contributions
        B_dip(:,i) = dipoleField( ( pTool_test(:,j) - (pAct_cartesion(:,i)) ), m_mag*Rzy);
        B_dip_G(:,i) = dipoleGradient( ( pTool_test(:,j) - (pAct_cartesion(:,i)) ), m_mag*Rzy);
    end 
    %I = I_x/I_xMax;
    % Using the currents that generate the max Bx By Bz fields at the
    % origin, calculate what the fields are on the vertices of the bounding
    % cube.
    B_probe(:,j) = [B_dip;B_dip_G]*I;
    % The magnitude and direction of these values can be compared to the
    % values at the origin to determine the uniformity of the max field in
    % the space.
end

%% Calculating error

mag_error = zeros(1,size(pTool_test,2));
dir_error = zeros(1,size(pTool_test,2));
Error_grad = zeros(5,size(pTool_test,2));
for j = 1:size(pTool_test,2)
    mag_error(1,j) = ( norm(B_probe(1:3,j)) - norm(B_bench(1:3)) )/(norm( B_bench(1:3) ));
    dir_error(1,j) = acos( B_probe(1:3,j)'*B_bench(1:3) / ( norm(B_probe(1:3,j))*norm(B_bench(1:3)) ) );
    
    Error_grad(1:5,j) = B_probe(4:8,j) - B_bench(4:8);
end

Error_field = [mag_error; dir_error];


end