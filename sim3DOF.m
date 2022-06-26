%% DESCRIPTION
% 

clear; close all; clc;

%% INPUTS & MODEL PARAMETERS
% Inputs
Voo = 1300;                 % [in/s]    undisturbed flow speed (Vflutter3DOF = 460 in/s)
rho = 0.000044256;          % [lb/in3]  air density
x0  = [5.0 0.0 0.0];        % Initial Condition {h theta beta}

dt  = 1e-2;                 % [sec] Sampling time
T   = 2;                    % [sec] Simulation Time

k   = floor(T/dt);          % [-]   Sampling indice (last value)
tControl = 1;               % [sec] Time at which the control is activated

reducedFreq = [0.05:0.05:0.5, 0.6:0.1:1.0, 1.2:0.2:2.0];    % set of reduced frequencies
gamma       = [0.2, 0.4, 0.6, 0.8];                         % lag parameters

% Load model parameters & flight condition
params  = modelParameters();
fltcond = struct('V',Voo,'rho',rho);
simcond = struct('k',reducedFreq,'gamma',gamma,'nLAG',length(gamma),'samplingTime',dt,'SimulationTime',T,'noIterations',k);

% Load state-space model
[sys,~,~] = stateSpaceControl(params,fltcond,simcond);

% Design control law 
K = controlDesign(params,fltcond,simcond);

%% SIMULATION
options_ode45 = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);

x = zeros(15,k+1);      % State Vector (for all sampling time steps)
input = zeros(1,k+1);   % Control action value (for all sampling time steps)
x(1:3,1) = x0;

for i=1:k
    % Showing simulation iteration
    fprintf('Starting %d° iteration\n\n',i);
    
    % Set control action value
    if i >= tControl/dt
        input(i) = -15;%K*x(:,i);
        % input(i) = min(input(i), 15*pi/180);
        % input(i) = max(input(i), -15*pi/180);
    else
        input(i) = 0;
    end

    % Evolving the system dynamics
    y0 = x(:,i);
    [t, y] = ode45(@(t,y) system(t,y,sys,input(i)), [0.0 dt], y0);
    x(:,i+1) = y(length(t),:);
end

% Plot Results
plotResults(T,k,x,input);

%% 
% Dynamic model of the system
function xdot = system(~,y,sys,input)
    % Input variables
    A    = sys.A;
    B    = sys.B;
    u    = input;

    % State-Space systems equations
    xdot = A*y + B.*u;
end
