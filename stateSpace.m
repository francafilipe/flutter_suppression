function [ssModel] = stateSpace(params,fltcond,simcond)
%% DESCRIPTION
% State Space model for the aeroelastic system of a typical section of a wing
% Author: Filipe Franca, Universidade Federal de Uberlandia

% The system can have 2 or 3 DOF, being the first two the plunge (heave)
% and pitch movements, and the later the deflection of a control surface 
% The unsteady aerodynamics are approximated via the Roger's RFA approx.,
% using nLAG number of lag parameters, which introduces nLAG*nDOF more
% states to the system

%% 

% Input Parameters
nDOF = params.nDOF;
b    = params.b;
Voo  = fltcond.V;
nLAG = simcond.nLAG;
gamma = simcond.gamma;

% Call Roger RFA for the specific aerodynamic modelling or data input
[RFAs] = rogerRFA(params,fltcond,simcond);

% State-space model matrices
M   = params.M - RFAs.A2*(params.b/fltcond.V)^2;     % Mass matrix (Structural Mass - Aerodynamics) - terms that multiply x_dot_dot
B   = params.B - RFAs.A1*(params.b/fltcond.V);       % Damping matrix - terms that multiply x_dot
K   = params.K - RFAs.A0;                            % Stiffness matrix  - terms that multiply x

% state (or system) matrix
nStates = 2*nDOF + nDOF*nLAG;
state = zeros(nStates);
state(1:nDOF,(nDOF+1):2*nDOF)   = eye(nDOF);                                                            % first three rows
state(nDOF+1:2*nDOF,:)          = cat(2, -M\K, -M\B, M\RFAs.A3, M\RFAs.A4, M\RFAs.A5, M\RFAs.A6);       % dynamic equations (4th - 6th rows)
for nlag=1:nLAG
state(((1+nlag)*nDOF+1):(2+nlag)*nDOF,(nDOF+1):2*nDOF)                  = eye(nDOF);                    % first temr of the lag parameter derivative equations (4th -6th columns)
state(((1+nlag)*nDOF+1):(2+nlag)*nDOF,((1+nlag)*nDOF+1):(2+nlag)*nDOF)  = -(Voo/b)*gamma(nlag)*eye(3);  % lag terms (second term of the derivative equation for the lags)
end

% input (or control) matrix
input = zeros(nStates,1);
Ba = [0 0 params.K(nDOF,nDOF)]';
input(1:nDOF) = -inv(M)*Ba;

% output matrix
output = zeros(2*nDOF,nStates);
for k=1:2*nDOF
    output(k,k) = 1;
end

% feedforward matrix
feedfoward = zeros(2*nDOF,1);


%% 
% Insert actuator's dynamics
% tau = 5;
% state(3,:) = zeros(1,length(state(3,:)));
% state(3,3) = -1/tau;
% 
% input = zeros(nStates,1);
% input(3,1) = 1/tau;

%% 
% Output Parameters
ssModel = ss(state,input,output,feedfoward);

end