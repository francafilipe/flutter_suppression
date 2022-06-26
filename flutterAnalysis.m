%% FLUTTER SPEED ANALYSIS 
%  Script that analyze the freq and damping of the modes for different speed
%  Author: Filipe Franca, Universidade Federal de Uberlandia


clear; close all; clc

%% INPUTS & MODEL PARAMETERS
% Inputs
Voo = 0:100:1500;                   % [in/s]    undisturbed flow speed
rho = 0.000044256;              % [lb/in3]  air density

% Load model parameters & flight condition
params  = modelParameters();
fltcond = struct('V',Voo,'rho',rho);
simcond = struct('k',{},'gamma',{},'nLAG',{});

% Load state-space model
% sys = stateSpace(params,fltcond,simcond);


%% Free Vibration Analasis
[PHI,lamb]=eig(params.K,params.M);	% PHI - matriz modal, lamb - matriz cuja diagonal s�o os autovalores

% Modal matrices
Mm = PHI'*params.M*PHI;
Km = PHI'*params.K*PHI;
Cm = PHI'*params.B*PHI;

% Natural frequencies
omegaN    = sqrt(diag(lamb));                   % [rad/s] Natural frequencies
fN        = omegaN/(2*pi);                      % [Hz]
Qci       = diag(Cm)./(2*diag(Mm).*omegaN);     % Damping factor
omegaD    = sqrt(1-Qci.^2).*omegaN;             % Damped natural frequencies


%% FREQUENCY MATCHING (p-k method)
wMatch      = zeros(params.nDOF,length(Voo));
wMatch(:,1) = omegaN; % Wind-Off frequency for all modes (freq. for V = 0)
g           = zeros(params.nDOF,length(Voo));
g(:,1)      = Qci;    % Wind-Off damping ratio

for i=2:length(Voo)
    
    fprintf('-----------------------------------------------------------\n');
    fprintf('Calculating for Speed = %d m/s\n\n',Voo(i))
    fltcond.V = Voo(i);
    
    for n=1:params.nDOF

        % Initial Guess for the frequency of the mode n for i-th Velocity
        wMatch(n,i) = wMatch(n,i-1);
        
        while true
            k     = (wMatch(n,i)*params.b/fltcond.V);

            % Aerodynamics loads
            aero = aeroMatrices(params,fltcond,k);
            Q    = aero.AIC;      % AIC matrix
            Qm   = PHI'*Q*PHI;    % AIC matrix modal form

            % Eigenvalues determination
            [V,A] = eig(Km+Qm,Mm);
            A     = sort(diag(A));
            omega = sqrt(A(n));

            % Check for convergence of omega
            if abs(omega-wMatch(n,i))>1e-6
                wMatch(n,i) = omega;
            else
                break
            end
        end
        
        % Note corresponding damping ratio for the converged freqs.
        g(n,i) = Cm(n,n)/(2*Mm(n,n)*wMatch(n,i));     % Damping factor

    end
end

%%
% Plot Resultados 
figure(1);

subplot(2,1,1); hold on; grid minor;
plot(Voo,wMatch(1,:),'-o','linewidth',1.5)
plot(Voo,wMatch(2,:),'-o','linewidth',1.5)
plot(Voo,wMatch(3,:),'-o','linewidth',1.5)
ylabel('Frequencia ~ rad/s ~')
ylim([0 150])

subplot(2,1,2); hold on; grid minor;
plot(Voo,g(1,:),'-o','linewidth',1.5)
plot(Voo,g(2,:),'-o','linewidth',1.5)
plot(Voo,g(3,:),'-o','linewidth',1.5)
ylabel('Amortecimento ~ g ~')
xlabel('Velocidade ~ m/s ~')
