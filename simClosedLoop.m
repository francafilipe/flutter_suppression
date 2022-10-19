%% DESCRIPTION & INTRODUCTION
% 

clear; close all; clc;
SetPaths;

%% INPUTS & MODEL PARAMETERS
% Load model parameters
params = Airfoil2DOF();

% Define flight conditions for simulation
fltcond.Voo = 1165;
fltcond.rho = 0.000044256;
fltcond.qoo = (fltcond.rho*fltcond.Voo^2)/2;

% Load aerodynamics and allocate
load('RFA_Airfoil2DOF.mat');
params.RFAs = RFAs;

% Aeroservoelastic State-Space Model Representation
sys = ssAeroservoelastic(params,fltcond);


% ---------------------- Compute | Load Control Law -----------------------
% Control Design Parameters
samplingTime = 1e-2;
controlWeights.State    = 1*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
controlWeights.Control  = 1;

% Flight Condition at which the control law will be designed
fltcondControl.Voo = 1165;
fltcondControl.rho = 0.000044256;
fltcondControl.qoo = (fltcondControl.rho*fltcondControl.Voo^2)/2;

control = controlDesignLQR(params,fltcondControl,samplingTime,controlWeights);
K = control.Gain;

%% SIMULATION
% Simulation inputs
x0  = [5.0 0.0 0.0 0.0];   % Initial Condition {h theta hDot thetaDot}
dt  = 1e-2;                 % [sec] Sampling time
T   = 1.5;                    % [sec] Simulation Time
k   = floor(T/dt);          % [-]   Sampling indice (last value)
tControl = 0;             % Control activation timestamp
controlON = 0;

% Variable allocation
x = zeros(length(sys.A),k+1);       % State Vector (for all sampling time steps)
input = zeros(1,k+1);               % Control action value (for all sampling time steps)
x(1:4,1) = x0;
xOL = x;

for i=1:k
    % Showing simulation iteration
    fprintf('Starting %d° iteration\n\n',i);
    
    % Set control action value
    if (i>=tControl/dt && x(1,i)>=0)
        controlON = 1;
    end
    
    input(i) = -K*x(:,i)*controlON;

    % Evolving the system dynamics
    y0 = x(:,i);
    [t, y] = ode45(@(t,y) (sys.A*y+sys.B*input(i)), [0.0 dt], y0);
    x(:,i+1) = y(length(t),:);
    
    % Evolving the system dynamics
    y0OL = xOL(:,i);
    [t, yOL] = ode45(@(t,y) (sys.A*y+sys.B*0), [0.0 dt], y0OL);
    xOL(:,i+1) = yOL(length(t),:);
end

%% POST-PROCESSING
t = linspace(0,T,k+1);
input(k+1) = input(k);

% ----------------------------- Plot Results ------------------------------
% -------------------------------------------------------------------------
figure(1); % Modes | DOFs Dynamic Response
tiledlayout(4, 1, 'TileSpacing', 'compact')

nexttile; hold on; box on;
plot(t,x(1,:)/39.37,'-k','linewidth',2)
plot(t,xOL(1,:)/39.37,'-ok','linewidth',1)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
% xlabel('$t \quad [sec]$','interpreter','latex','FontSize',12);
ylabel('$h$ [m]','interpreter','latex','FontSize',12);

lg = legend('$Control \ On$','$Control \ Off$','interpreter','latex','FontSize',13);
lg.Orientation = 'horizontal';
lg.Location = 'Northoutside';

nexttile; hold on; box on;
plot(t,x(3,:)/39.37,'-k','linewidth',2)
plot(t,xOL(3,:)/39.37,'-ok','linewidth',1)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
% xlabel('$t \quad [sec]$','interpreter','latex','FontSize',12); 
ylabel('$\dot{h}$ [m/sec]','interpreter','latex','FontSize',12); 


nexttile; hold on; box on;
plot(t,x(2,:)*180/pi,'-k','linewidth',2)
plot(t,xOL(2,:)*180/pi,'-ok','linewidth',1)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
% xlabel('$t \quad [sec]$','interpreter','latex','FontSize'r,12); 
ylabel('$\theta$ [deg]','interpreter','latex','FontSize',12); 

nexttile; hold on; box on;
plot(t,x(4,:)*180/pi,'-k','linewidth',2)
plot(t,xOL(4,:)*180/pi,'-ok','linewidth',1)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
xlabel('$t$ [s]','interpreter','latex','FontSize',12); 
ylabel('$\dot{\theta}$  [deg/s]','interpreter','latex','FontSize',12); 

% -------------------------------------------------------------------------
figure(2); % Control Input & Control Surface States
tiledlayout(3, 1, 'TileSpacing', 'compact')

nexttile; hold on; box on;
plot(t,x(end-2,:)*180/pi,'-k','linewidth',2)
plot(t,input*180/pi,'--ok','linewidth',0.5)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
ylabel('$\delta$ [deg]','interpreter','latex','FontSize',12);
lg = legend('$\delta$','$u_{c}$','interpreter','latex','FontSize',13);
lg.Orientation = 'horizontal';
lg.Location = 'Northoutside';

nexttile; hold on; box on;
plot(t,x(end-1,:)*180/pi,'-k','linewidth',1.5)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
ylabel('$\dot{\delta}$ [deg/s]','interpreter','latex','FontSize',12); 

nexttile; hold on; box on;
plot(t,x(end,:)*180/pi,'-k','linewidth',1.5)
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
xlabel('$t \ $[s]','interpreter','latex','FontSize',12); 
ylabel('$\ddot{\delta}$  [deg/s2]','interpreter','latex','FontSize',12); 

