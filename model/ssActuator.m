function [ssModel] = ssActuator(frequency,coeffs)
%% DESCRIPTION
% State Space model for the control surface actuator
% Author: Filipe Franca, Universidade Federal de Uberlandia

% 

% Input variables:
% frequency     - control surface actuator's frequency
% coeffs        - coefficients for the actuators dynamics

% ------------------------------- IMPORTANT -------------------------------
% OBS: The algorithm is implemented to use only 4 lag parameters and is yet
% to be updated to take as an input the number of lag terms
% -------------------------------------------------------------------------

%% CONTROL SURFACE ACTUATOR MODEL
omega = 2*pi*frequency;
A0 = coeffs(1)*omega^3;
A1 = coeffs(2)*omega^2;
A2 = coeffs(3)*omega;

% State Space model matrices for the actuators dynamics
A = [0   1   0;
     0   0   1;
   -A0 -A1 -A2];
 
B = [0; 0; A0];
C = eye(3);
D = zeros(3,1);

% State Space model
ssModel = ss(A,B,C,D);

end