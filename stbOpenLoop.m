%% DESCRIPTION
% 

clear; close all; clc;
SetPaths; 

%% ANALYSIS CNODITIONS
% Define flight conditions | analysis points
fltcond.rho = 0.000044256;          % [lb/in3]  air density
Voo = [700];       % [in/s]    undisturbed flow speed
% Voo = [0.1 10:20:1200];
nVoo = length(Voo);

% Load model parameters
params = load('modal-Goland-Wing-6modes.mat');

% Load aerodynamics and allocate
load('RFA-Goland-Wing-Ma05-6Modes.mat');
params.RFAs = RFAs;

% Define poles array
nPoles = params.nDOF*(2+params.RFAs.nLAG);
poles = zeros(nVoo,nPoles);

for i=1:nVoo
    % Define flight condition
    fltcond.Voo = Voo(i);
    fltcond.qoo = (fltcond.rho*fltcond.Voo^2)/2;

    % Calculate Aeroelastic State Space model
    sys = ssAeroelastic(params,fltcond);
    currPoles = eig(sys.A);
    
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
    
    poles(i,:) = sortedPoles;
end

%% POST-PROCESSING
limyMax = max(imag(poles),[],'all');
limyMin = min(imag(poles),[],'all');

% Plot Results
plt = figure(1); hold on; box on;
plot(real(poles)/(2*pi),imag(poles)/(2*pi),'-ok','linewidth',1.5)
plot([-2.5 0.5],[0 0], '--k')
plot([0 0],[limyMin limyMax]*1.2/(2*pi), '--k')
set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
xlabel('$Real \ [Hz]$','interpreter','latex','FontSize',12); ylabel('$Imag \ [Hz]$','interpreter','latex','FontSize',12); 

% Indicate V_OLF condition 
VOLFannotation = annotation('textarrow', [0.82 0.785], [0.65 0.75], 'String','$V_{OLF}$','interpreter','latex','FontSize',12);
Mod1annotation = annotation('textarrow', [0.38 0.36], [0.68 0.76], 'String','$Primeiro \ Modo \ Aeroelastico$', 'interpreter','latex','FontSize',12);
Mod2annotation = annotation('textarrow', [0.60 0.64], [0.40 0.33], 'String','$Segundo \ Modo \ Aeroelastico$','interpreter','latex','FontSize',12);