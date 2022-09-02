function [output] = controlDesignLQR(params,fltcond,samplingTime,weights)
%% DESCRIPTION
% Control law design using LQR for a space state model 
% Author: Filipe Franca, Universidade Federal de Uberlandia

%% INPUTS & DEFINITIONS
% Number of states
nDOF = params.nDOF;
nLAG = params.RFAs.nLAG;
nControl = 3;
nStates  = nDOF*(2+nLAG)+nControl;

% Aeroservoelastic State-Space Model Representation
sys = ssAeroservoelastic(params,fltcond);

% Discretizando o sitema
discreteSys = c2d(sys, samplingTime, 'zoh');

% Matriz de peso dos estados Qstate
Qstate   = weights.State.*eye(nStates);
Qcontrol = weights.Control;

% LQR Control Design
[K,S,lamb] = dlqr(discreteSys.A, discreteSys.B, Qstate, Qcontrol);

% Evaluate Closed Loop Stability

% Output Parameters
output.Gain = K;            % Control Gain
output.S    = S;            % Infinite horizon solution of the associated discrete-time Riccati equation
output.CLPoles = lamb;      % Closed loop poles for the system
output.discreteSystem = discreteSys; % Discrete System
end