function [ssModel] = stateSpace(params,fltcond)
%% DESCRIPTION
% State Space model for the aeroelastic system of a typical section of a wing
% Author: Filipe Franca, Universidade Federal de Uberlandia

% The system can have 2 or 3 DOF, being the first two the plunge (heave)
% and pitch movements, and the later the deflection of a control surface 
% The unsteady aerodynamics are approximated using the Jones exponential
% approx. for the Wagner function, which includes 2 lag states in the model

%% 
% Input Parameters
nDOF = params.nDOF;
nStates = 2*nDOF + 2;

% Call Roger RFA for the specific aerodynamic modelling or data input
[aerodynamics] = aerodynamics3DOF(params,fltcond);

% State-space model matrices
Mt   = params.M - aerodynamics.A1;              % Mass matrix (Structural Mass - Aerodynamics) - terms that multiply x_dot_dot
Ct   = params.B - aerodynamics.A2;              % Damping matrix - terms that multiply x_dot
Kt   = params.K - aerodynamics.A3;              % Stiffness matrix  - terms that multiply x

A4   = aerodynamics.A4;
B1   = aerodynamics.B1;
B2   = aerodynamics.B2;
B3   = aerodynamics.B3;
B4   = aerodynamics.B4;

% state (or system) matrix
state = [zeros(nDOF)        eye(nDOF)           zeros(nDOF,2); 
         -inv(Mt)*Kt        -inv(Mt)*Ct         inv(Mt)*A4;
        B3-B1*inv(Mt)*Kt    B2-B1*inv(Mt)*Ct    B4+B1*inv(Mt)*A4];


% input (or control) matrix
input = zeros(nStates,1);

% output matrix
output = zeros(2*nDOF,nStates);

% feedforward matrix
feedfoward = zeros(2*nDOF,1);

% Output Parameters
ssModel = ss(state,input,output,feedfoward);

end