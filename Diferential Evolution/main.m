%% DESCRIPTION 
% 

clear; close all; clc
format long

%% MODEL & SIMULATION PARAMETERS

% 3DOF Simulation params
Voo = 1300;                 % [in/s]    undisturbed flow speed (Vflutter3DOF = 460 in/s)
rho = 0.000044256;          % [lb/in3]  air density
x0  = [5.0 0.0];            % Initial Condition {h theta beta}

dt  = 1e-2;                 % [sec] Sampling time
T   = 0.5;                  % [sec] Simulation Time

k   = floor(T/dt);          % [-]   Sampling indice (last value)
tControl = 0;               % [sec] Time at which the control is activated

reducedFreq = [0.05:0.05:0.5, 0.6:0.1:1.0, 1.2:0.2:2.0];    % set of reduced frequencies
gamma       = [0.2, 0.4, 0.6, 0.8];                         % lag parameters

% Load model parameters & flight condition
params  = modelParameters();
fltcond = struct('V',Voo,'rho',rho);
simcond = struct('k',reducedFreq,'gamma',gamma,'nLAG',length(gamma),'samplingTime',dt,'SimulationTime',T,'noIterations',k,'initialCondition',x0,'controlActivationTime',tControl);

% Create state-space model
[sys,~,~] = stateSpaceControl(params,fltcond,simcond);

% Create discrete model
discreteSys = c2d(sys, simcond.samplingTime, 'zoh');

%%

       %  x1  x2  x3  x4
XVmin = [  0   0   0   0];
XVmax = [1e4 1e4 1e4 1e4];     

D   = length(XVmin);        % No of optimization variables
NP  = 120*D;                % Population size
itermax = 100;              % Max number of iterations (generations)
VTR = 1.e-20;               % Stoping Criteria - error for the last population

% Input parameters variable (model & simulation params)
y = struct('modelParams', params, 'flightConditions', fltcond, 'simulationConditions', simcond, 'continuousSys', sys, 'discreteSys', discreteSys); 

F = 0.8;    % taxa de perturbacao (0-2)
CR = 0.8;   % probabilidade de cruzamento (0-1)
strategy=7; % estrategia - veja a minha tese para mais detalhes ...
refresh=10; % nivel de impressao na tela - nao precisa alterar esse valor ...

[X,FO,NF] = differential_evolution('eval_objective',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);

% X e o projeto final obtido
% FO e o valor da funcao objetivo
% NF e o numero de avaliacoes da funcao objetivo

%% 
