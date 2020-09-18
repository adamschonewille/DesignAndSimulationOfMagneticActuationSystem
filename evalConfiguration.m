% EVALUATE CONFIGURATION
clear all; clc; close all;

% mAct_sph =    [-0.0000   -0.0000   -0.0000    0.0000    0.0004   -0.0003    0.0003   -0.0004;
%                -0.0010   -0.0000    0.0040    0.0004    0.0284   -0.0272   -0.0266    0.0280];
% 
% 
% pAct_cartesion = [-0.5004   -0.0002    0.5007   -0.0002   -0.1699    0.1670    0.1682   -0.1711;
%                   -0.0000   -0.5027   -0.0001    0.5034   -0.1659   -0.1670    0.1637    0.1625;
%                   -0.3001   -0.3000   -0.3003   -0.3000   -0.3018   -0.3018   -0.3017   -0.3018];
              
              
mAct_sph =    [0.0000   0.0000   0.0000    0.0000    0.0000    0.0000    0.0000    0.0000;
               0.0000   0.0000   0.0000    0.0000    0.0000    0.0000    0.0000    0.0000];

pAct_cartesion = [-0.5000    0.0000    0.5000    0.0000   -0.1670    0.1670    0.1670   -0.1670;
                   0.0000   -0.5000    0.0000    0.5000   -0.1670   -0.1670    0.1670    0.1670;
                  -0.3000   -0.3000   -0.3000   -0.3000   -0.3000   -0.3000   -0.3000   -0.3000];

mAct_sph =    [3.8337    1.7582    1.7906    2.2138    4.8549    4.9023   -0.0025    2.2042;
    0.1038    1.2821    0.2663    0.0407    0.8914    0.0550    1.0818    0.5916];


pAct_cartesion = [ 0.0827   -0.0567   -0.0759    0.3317   -0.2452    0.0802   -0.4175   -0.1105;
    0.0864   -0.2169    0.0624   -0.1040   -0.0851   -0.0762    0.0967    0.4022;
   -0.3033   -0.2351   -0.3086   -0.3025   -0.2851   -0.3031   -0.2640   -0.3070];
              
              
mu_0 = pi*4e-7;         % [H/m or N/A^2 (amp-turns)]
r = 0.12; % [m] Distance from posterior of head to the center of the brain

%% ~~~~~~~~~~ ~~~~~~~~~~ OPERATING PARAMETERS ~~~~~~~~~~ ~~~~~~~~~~ %%
magnetSize = 0.0254; % [m]
Remenance = 1.45; % [T]

m_EM_large = 2*pi*(0.12)^3*(0.11407)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 63.8 mT
% m_EM_small = 2*pi*(0.12)^3*(0.04846)/mu_0; % [A m^2] (for 24 A/mm^2) at 6 A/mm^2 -> 18.9 mT
d_EM_large = 2*(45 + 22.5) / 1000;  % [m]
% d_EM_small = 2*(30 + 15) / 1000;    % [m]
l_EM_large = 360 / 1000;            % [m]
% l_EM_small = 240 / 1000;            % [m]
large_EM = [m_EM_large, d_EM_large, l_EM_large]';

%% Evaluations

G = evalG(r, large_EM, mAct_sph, pAct_cartesion)
F = 1/2 * G'*G

a = 0.030; % [m] 1/2 sidelength of bounding cube for system workspace
[B_maxNonUniform, I_maxNonUniform] = calcMaxNonUniformBFields(mAct_sph, pAct_cartesion, large_EM);
B_maxNonUniform
[Error_Bx_nonUniform, Error_GBx_nonUniform] = calcWorkspaceEdgeError(mAct_sph, pAct_cartesion, large_EM, a, I_maxNonUniform(:,1))
[Error_By_nonUniform, Error_GBy_nonUniform] = calcWorkspaceEdgeError(mAct_sph, pAct_cartesion, large_EM, a, I_maxNonUniform(:,2))
[Error_Bz_nonUniform, Error_GBz_nonUniform] = calcWorkspaceEdgeError(mAct_sph, pAct_cartesion, large_EM, a, I_maxNonUniform(:,3))

[B_maxUniform, I_maxUniform] = calcMaxUniformBFields(mAct_sph, pAct_cartesion, large_EM);
B_maxUniform
[Error_Bx_Uniform, Error_GBx_Uniform] = calcWorkspaceEdgeError(mAct_sph, pAct_cartesion, large_EM, a, I_maxNonUniform(:,1))
[Error_By_Uniform, Error_GBy_Uniform] = calcWorkspaceEdgeError(mAct_sph, pAct_cartesion, large_EM, a, I_maxNonUniform(:,2))
[Error_Bz_Uniform, Error_GBz_Uniform] = calcWorkspaceEdgeError(mAct_sph, pAct_cartesion, large_EM, a, I_maxNonUniform(:,3))

[U_3, S_3, V_3, U_8, S_8, V_8] = calcSystemSVD(mAct_sph, pAct_cartesion, large_EM)

