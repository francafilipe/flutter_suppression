%% DESCRIPTION
% 

clear; close all; clc

%% INPUTS & MODEL PARAMETERS
% Inputs
Voo = 540;                  % [in/s]    undisturbed flow speed
rho = 0.000044256;          % [lb/in3]  air density
x0  = [5.0 0.0 0.0];        % Initial Condition {h theta beta}

reducedFreq = [0.05:0.05:0.5, 0.6:0.1:1.0, 1.2:0.2:2.0];    % set of reduced frequencies
gamma       = [0.2, 0.4, 0.6, 0.8];                         % lag parameters

% Load model parameters & flight condition
params  = modelParameters();
fltcond = struct('V',Voo,'rho',rho);
simcond = struct('k',reducedFreq,'gamma',gamma,'nLAG',length(gamma));

%% COMPUTING
% Aerodynamics matrices (reduced frequency dependent AICs)
[aero] = aeroMatrices(params,fltcond,simcond.k);

% Roger's RFAs
[RFAs] = rogerRFA(params,fltcond,simcond);
for k=1:length(reducedFreq)
Q(:,:,k) =  RFAs.A0 + RFAs.A1*(1i*reducedFreq(k)) + RFAs.A2*(1i*reducedFreq(k))^2 + ...
            RFAs.A3*(1i*reducedFreq(k))/(1i*reducedFreq(k)+gamma(1)) + ...
            RFAs.A4*(1i*reducedFreq(k))/(1i*reducedFreq(k)+gamma(2)) + ...
            RFAs.A5*(1i*reducedFreq(k))/(1i*reducedFreq(k)+gamma(3)) + ...
            RFAs.A6*(1i*reducedFreq(k))/(1i*reducedFreq(k)+gamma(4));

    % Reduced-frequency dependent indices of rfa
    Q11(k) = Q(1,1,k);
    Q12(k) = Q(1,2,k);
    Q13(k) = Q(1,3,k);
    Q21(k) = Q(2,1,k);
    Q22(k) = Q(2,2,k);
    Q23(k) = Q(2,3,k);
    Q31(k) = Q(3,1,k);
    Q32(k) = Q(3,2,k);
    Q33(k) = Q(3,3,k);

    % Reduced-frequency dependent indices of AIC
    AIC11(k) = aero(k).AIC(1,1);
    AIC12(k) = aero(k).AIC(1,2);
    AIC13(k) = aero(k).AIC(1,3);
    AIC21(k) = aero(k).AIC(2,1);
    AIC22(k) = aero(k).AIC(2,2);
    AIC23(k) = aero(k).AIC(2,3);
    AIC31(k) = aero(k).AIC(3,1);
    AIC32(k) = aero(k).AIC(3,2);
    AIC33(k) = aero(k).AIC(3,3);

end

%% EVALUATE RESULTS
figure(1)
tiledlayout(3,3)

nexttile; hold on; grid minor;
plot(real(Q11),imag(Q11),'-k','linewidth',1.5)
plot(real(AIC11),imag(AIC11),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{11}');
legend('RFA','AIC','location','southeast')

nexttile; hold on; grid minor;
plot(real(Q12),imag(Q12),'-k','linewidth',1.5)
plot(real(AIC12),imag(AIC12),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{12}');
legend('RFA','AIC','location','southwest')

nexttile; hold on; grid minor;
plot(real(Q13),imag(Q13),'-k','linewidth',1.5)
plot(real(AIC13),imag(AIC13),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{13}');
legend('RFA','AIC','location','southwest')


nexttile; hold on; grid minor;
plot(real(Q21),imag(Q21),'-k','linewidth',1.5)
plot(real(AIC21),imag(AIC21),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{21}');
legend('RFA','AIC','location','southeast')

nexttile; hold on; grid minor;
plot(real(Q22),imag(Q22),'-k','linewidth',1.5)
plot(real(AIC22),imag(AIC22),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{22}');
legend('RFA','AIC','location','southeast')

nexttile; hold on; grid minor;
plot(real(Q23),imag(Q23),'-k','linewidth',1.5)
plot(real(AIC23),imag(AIC23),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{23}');
legend('RFA','AIC','location','southeast')


nexttile; hold on; grid minor;
plot(real(Q31),imag(Q31),'-k','linewidth',1.5)
plot(real(AIC31),imag(AIC31),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{31}');
legend('RFA','AIC','location','southeast')

nexttile; hold on; grid minor;
plot(real(Q32),imag(Q32),'-k','linewidth',1.5)
plot(real(AIC32),imag(AIC32),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{32}');
legend('RFA','AIC','location','southwest')

nexttile; hold on; grid minor;
plot(real(Q33),imag(Q33),'-k','linewidth',1.5)
plot(real(AIC33),imag(AIC33),'ok','markersize',5)
xlabel('Real'); ylabel('Imag'); title('Q_{33}');
legend('RFA','AIC','location','southwest')


