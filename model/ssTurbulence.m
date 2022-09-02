function [ssModel] = ssTurbulence(Lg,Voo,sigma,a)
%% DESCRIPTION
% State Space model for the control surface actuator
% Author: Filipe Franca, Universidade Federal de Uberlandia

% 

% Input variables:
% L             - scale of turbulence
% Voo           - coefficients for the actuators dynamics
% sigma         - RMS value of the gust velocity
% a             - break frequency for the low-pass filter

%% 
tau = Lg/Voo;   % time constant for the turbulence model

% State Space representation matrices
Ag = zeros(3,3);
Bg = zeros(3,1);
Cg = zeros(2,3);

ssModel = ss(Ag, Bg, Cg, zeros(2,1));
end