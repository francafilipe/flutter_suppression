%% DESCRIPTION
% 

clear; close all; clc

%% INPUTS & MODEL PARAMETERS
% Inputs
Voo = 564;                  % [in/s]    undisturbed flow speed
rho = 0.000044256;          % [lb/in3]  air density
x0  = [5.0 0.0 0.0];        % Initial Condition {h theta beta}

dt  = 1e-2;                 % [sec] Sampling time
T   = 3;                    % [sec] Simulation Time
k   = floor(T/dt);          % [-]   Sampling indice (last value)

% Load model parameters & flight condition
params  = modelParameters();
fltcond = struct('V',Voo,'rho',rho);

% Load state-space model
sys = stateSpace(params,fltcond);

%% SIMULATION
options_ode45 = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);

x = zeros(8,k+1);    % State Vector (for all sampling time steps)
input = zeros(1,k+1); % Control action value (for all sampling time steps)
x(1:3,1) = x0;

for i=1:k
    % Showing simulation iteration
    fprintf('Starting %d° iteration\n\n',i);
    
    % Set control action value
    input(i) = 0;
    % Evolving the system dynamics
    y0 = x(:,i);
    [t, y] = ode45(@(t,y) system(t,y,sys,input(i)), [0.0 dt], y0);
    x(:,i+1) = y(length(t),:);
end

% Plot Results
plotResults(T,k,x,input);

%% Dynamic model of the system
function xdot = system(~,y,sys,input)
    % Input variables
    A    = sys.A;
    B    = sys.B;
    u    = input;
    % State-Space systems equations
    xdot = A*y + B.*u;
end
