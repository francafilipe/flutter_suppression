function f = eval_objective(x,optParams)
%% DESCRIPTION

%% INPUTS & MODEL PARAMETERS
tic
%  Optimization Parameters
alpha = x(1);
beta  = x(2);
gamma = x(3);
delta = x(4);

% Load model parameters & flight condition
params  = optParams.modelParams;
fltcond = optParams.flightConditions;
simcond = optParams.simulationConditions;

% Load state-space model
sys = optParams.continuousSys;
discreteSys = optParams.discreteSys;


%% CONTROL DESIGN 
% Define control design conditon
fltcond.V = params.Vflutter;

% Matriz de peso dos estados Qstate
Qstate = eye(15).*[alpha alpha beta beta gamma gamma gamma gamma gamma gamma gamma gamma gamma gamma gamma];

% Projetando o controlador de espaco de estados
[K,~,~] = dlqr(discreteSys.A, discreteSys.B, Qstate, delta);

%% NUMERIC SIMULATION (Closed-loop)
dt  = simcond.samplingTime;                 % [sec] Sampling time
T   = simcond.SimulationTime;               % [sec] Simulation Time

k   = simcond.noIterations;                 % [-]   Sampling indice (last value)
tControl = simcond.controlActivationTime;   % [sec] Time at which the control is activated

x = zeros(15,k+1);                          % State Vector (for all sampling time steps)
input = zeros(1,k+1);                       % Control action value (for all sampling time steps)
x(1:2,1) = simcond.initialCondition;        % Initial Condition {h theta beta}

for i=1:k
    % Showing simulation iteration
     fprintf('Starting %d° iteration\n\n',i);
    
    % Set control action value
    if i >= tControl/dt
        input(i) = -K*x(:,i);
        input(i) = min(input(i), 15*pi/180);
        input(i) = max(input(i), -15*pi/180);
    else
        input(i) = 0;
    end

    % Evolving the system dynamics
    y0 = x(:,i);
    [t, y] = ode45(@(t,y) (sys.A*y+sys.B.*input(i)), [0.0 dt], y0);
    x(:,i+1) = y(length(t),:);
end

%%
% Função Objetivo
f = norm(x(1,:)) + norm(x(2,:)) + norm(x(3,:)) + norm(x(4,:)) + ...
     norm(x(end-2,:)) + norm(x(end-1,:)) + norm(x(end,:));

%% 
% Display
fprintf('--------------------------------------------------------------\n')
fprintf('Vetor x: %d, %d, %d, %d \n',alpha,beta,gamma,delta)
fprintf('FO: %d\n', f)
fprintf('Tempo execucao: %d\n', toc)
fprintf('--------------------------------------------------------------\n\n\n')
end
