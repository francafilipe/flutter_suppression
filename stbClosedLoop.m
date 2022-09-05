%% DESCRIPTION
% 

clear; close all; clc;
SetPaths; 

%% ANALYSIS CNODITIONS
% Define flight conditions | analysis points
fltcond.rho = 0.000044256;          % [lb/in3]  air density
Voo = [1165];       % [in/s]    undisturbed flow speed
% Voo = [100 1165 1235 1300];
nVoo = length(Voo);

% Load model parameters
params = Airfoil2DOF();

% Load aerodynamics and allocate
load('RFA_Airfoil2DOF.mat');
params.RFAs = RFAs;

% Define poles array
nPoles = params.nDOF*(2+params.RFAs.nLAG)+3;
poles = zeros(nVoo,nPoles);

% ---------------------- Compute | Load Control Law -----------------------
% Control Design Parameters
samplingTime = 1e-2;
controlWeights.State    = 1e3*[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
controlWeights.Control  = 1;

% Flight Condition at which the control law will be designed
fltcondControl.Voo = 1165;
fltcondControl.rho = 0.000044256;
fltcondControl.qoo = (fltcond.rho*fltcondControl.Voo^2)/2;

control = controlDesignLQR(params,fltcondControl,samplingTime,controlWeights);
K = control.Gain;

for i=1:nVoo
    % Define flight condition
    fltcond.Voo = Voo(i);
    fltcond.qoo = (fltcond.rho*fltcond.Voo^2)/2;

    % Calculate Aeroelastic State Space model
    sys = ssAeroservoelastic(params,fltcond);
    dSys = c2d(sys, samplingTime, 'zoh');
    currPoles = eig(dSys.A-dSys.B*K);
    
    if i==1
        sortedPoles = currPoles;
    else
        lastPoles = poles(i-1,:);
        for r = 1:nPoles
            [minValue,closestIndex] = min(abs(currPoles-lastPoles(r)));
            sortedPoles(r) = currPoles(closestIndex);
            currPoles(closestIndex) = [];
        end
    end
    
    if Voo(i)==1165
        VOLF_poles = sortedPoles;
    elseif Voo(i)==1235
        VCLF_poles = sortedPoles;
    end

    poles(i,:) = sortedPoles;
end

poles'
%% POST-PROCESSING
xCircle = 1.005*cos(0:pi/180:2*pi)/(2*pi); 
yCircle = 1.005*sin(0:pi/180:2*pi)/(2*pi);

% Plot Results
plt = figure(1); hold on; box on;
plot(real(poles)/(2*pi),imag(poles)/(2*pi),'-ok','linewidth',1.5)
plot(xCircle,yCircle,'--k'); plot([-1 1]*1.2/(2*pi),[0 0], '--k'); plot([0 0],[-1 1]*1.2/(2*pi), '--k')
% plot(real(VOLF_poles)/(2*pi),imag(VOLF_poles)/(2*pi),'or','linewidth',1.5,'linestyle','none')
% plot(real(VCLF_poles)/(2*pi),imag(VCLF_poles)/(2*pi),'ob','linewidth',1.5,'linestyle','none')
xlim([-0.05 0.2]); ylim([-0.2 0.2]);
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
xlabel('$Real \ [Hz]$','interpreter','latex','FontSize',12); ylabel('$Imag \ [Hz]$','interpreter','latex','FontSize',12); 

% Indicate V_OLF condition 
% VOLFannotation = annotation('textarrow', [0.63 0.687], [0.78 0.695], 'String','$V_{OLF}$','interpreter','latex','FontSize',12);
% VCLFannotation = annotation('textarrow', [0.76 0.71], [0.76 0.697], 'String','$V_{CLF}$','interpreter','latex','FontSize',12);
% Mod1annotation = annotation('textarrow', [0.66 0.685], [0.44 0.375], 'String','$2nd \ Mode$','interpreter','latex','FontSize',12);
% Mod2annotation = annotation('textarrow', [0.55 0.62], [0.365 0.365], 'String','$1st \ Mode \ $','interpreter','latex','FontSize',12);
