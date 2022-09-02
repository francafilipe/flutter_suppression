%% DESCRIPTION & INTRODUCTION
% 

clear; close all; clc;
SetPaths;

%% INPUTS & MODEL PARAMETERS
% Load model parameters
params = Airfoil2DOF();

% Load | Define flight conditions
fltcond.Voo = 1165;
fltcond.rho = 0.000044256;
fltcond.qoo = (fltcond.rho*fltcond.Voo^2)/2;

% Load aerodynamics and allocate
load('RFA_Airfoil2DOF.mat');
params.RFAs = RFAs;

% Calculate state space model
sys = ssAeroelastic(params,fltcond);

%% SIMULATION
% Simulation inputs
x0  = [5.0 0.0 0.0 0.0];    % Initial Condition {h theta hDot thetaDot}
dt  = 1e-3;                 % [sec] Sampling time
T   = 1;                    % [sec] Simulation Time
k   = floor(T/dt);          % [-]   Sampling indice (last value)

% Variable allocation
x = zeros(length(sys.A),k+1);       % State Vector (for all sampling time steps)
input = zeros(3,k+1);               % Control action value (for all sampling time steps)
x(1:4,1) = x0;

for i=1:k
    % Showing simulation iteration
    fprintf('Starting %d° iteration\n\n',i);
    
    % Set control action value
    input(:,i) = [0.0 0.0 0.0];

    % Evolving the system dynamics
    y0 = x(:,i);
    [t, y] = ode45(@(t,y) (sys.A*y+sys.B*input(:,i)), [0.0 dt], y0);
    x(:,i+1) = y(length(t),:);
end

%% POST-PROCESSING
t = linspace(0,T,k+1);

% Plot Results
figure(1)
tiledlayout(4, 1, 'TileSpacing', 'compact')

nexttile; hold on; box on;
plot(t,x(1,:)/39.37,'-k','linewidth',1.5)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
% xlabel('$t \quad [sec]$','interpreter','latex','FontSize',12);
ylabel('$h$ [m]','interpreter','latex','FontSize',12);
% title('Modo de Translacao','interpreter','latex','FontSize',12);

nexttile; hold on; box on;
plot(t,x(3,:)/39.37,'-k','linewidth',1.5)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
% xlabel('$t \quad [sec]$','interpreter','latex','FontSize',12); 
ylabel('$\dot{h}$ [m/s]','interpreter','latex','FontSize',12); 
% title('Modo de Translacao - Derivada','interpreter','latex','FontSize',12);


nexttile; hold on; box on;
plot(t,x(2,:)*180/pi,'-k','linewidth',1.5)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
% xlabel('$t \quad [sec]$','interpreter','latex','FontSize'r,12); 
ylabel('$\theta$ [deg]','interpreter','latex','FontSize',12); 
% title('Modo de Rotacao','interpreter','latex','FontSize',12);

nexttile; hold on; box on;
plot(t,x(4,:)*180/pi,'-k','linewidth',1.5)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
xlabel('$t \quad [sec]$','interpreter','latex','FontSize',12); 
ylabel('$\dot{\theta}$ [deg/s]','interpreter','latex','FontSize',12); 
% title('Modo de Rotacao - Derivada','interpreter','latex','FontSize',12);
