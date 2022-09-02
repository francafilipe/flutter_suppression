function [ssModel] = ssAeroservoelastic(params,fltcond)
%% DESCRIPTION
% State Space model for the aeroelastic system
% Author: Filipe Franca, Universidade Federal de Uberlandia

% The unsteady aerodynamics are approximated via the Roger's RFA approx.,
% using nLAG number of lag parameters, which introduces nLAG*nDOF more
% states to the system

% Input variables:
% params        - systems parameters, such as the modal matrices, AICs etc
% flcond        - flight conditions for the modal determination

%% INPUT PARAMETERS
% Actuator's dynamics parameters
freqActuator   = 10;            % [Hz]  Actuator's frequency
coeffsActuator = [27 13.5 4.5];

% Turbulence definition
Lg    = 0;                      % [m]   scale of turbulence
Voo   = fltcond.Voo;            % [m/s] undisturbed flow velocity
sigma = 0;                      % [m/s]   gust velocity RMS value
a     = 0;                      % [rad/s]  break frequency of the low-pass filter

% No of states definition
nDOF = params.nDOF;
nLAG = params.RFAs.nLAG;
nControl = 3;

% Compute state space models for individual systems
ssAe = ssAeroelastic(params,fltcond);
ssAc = ssActuator(freqActuator, coeffsActuator);
ssGt = ssTurbulence(Lg,Voo,sigma,a);

% Complete state space representation
A  = [          ssAe.A             ssAe.B;
        zeros(3,nDOF*(2+nLAG))     ssAc.A   ];

B  = [  zeros(12,1); ssAc.B];

C  = eye(nDOF*(2+nLAG)+3);

D  = zeros(nDOF*(2+nLAG)+3,1);

ssModel = ss(A, B, C, D);

end

