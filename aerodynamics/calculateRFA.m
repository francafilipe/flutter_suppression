%% DESCRIPTION & INTRODUCTION
% 

clear; close all; clc
SetPaths;

%% 
% Analysis Configuration
plotar = 1;
salvar = 1;
fileNameSaving  = '/home/francafilipe/Desktop/Research/scripts/flutter_suppression/data-results/RFA_Airfoil2DOF.mat';
% fileNameLoading = '/home/francafilipe/Desktop/Research/scripts/flutter_suppression/data-results/AIC.mat';

% Inputs & Model Parameters
k = [0.05:0.05:0.5, 0.6:0.1:1.0, 1.2:0.2:2.0];      % set of reduced frequencies
lag = [0.2 0.4 0.6 0.8];                            % lag parameters

params  = Airfoil2DOF();

Voo = params.V_OLF;                        % [in/s]    undisturbed flow speed (Vflutter3DOF = 460 in/s)
rho = 0.000044256;                         % [lb/in3]  air density
qoo = rho*(Voo^2)/2;

fltcond = struct('Voo',Voo,'rho',rho,'qoo',qoo);
nDOF = params.nDOF;

% Definition of the tabular data for reduced-freq dependent aerodynamics
aerodynamic = aeroMatrices(params,fltcond,k);
% AIC = zeros(length(k),nDOF,nDOF+1);
for i=1:length(k)
    AIC(i,:,:) = aerodynamic(i).AIC(1:nDOF,1:nDOF+1);
end

% RFAs Determination by Roger's method
RFAs = rogerRFA(k,lag,nDOF,AIC);

%% SAVING & EVALUATING RESULTS
if salvar
    save(fileNameSaving,'RFAs');
end

if plotar
    % Calculate the Q matrix given the RFA in each reduced frequency
    Qqq = zeros(length(k),nDOF,nDOF);
    Qqc = zeros(length(k),nDOF,1);

    for index=1:length(k)
       Qqq(index,:,:) = RFAs.modal.A0 + ...
                        RFAs.modal.A1*(1i*k(index)) + ...
                        RFAs.modal.A2*(1i*k(index))^2 + ...
                        RFAs.modal.A3*(1i*k(index))/(1i*k(index)+lag(1)) + ...
                        RFAs.modal.A4*(1i*k(index))/(1i*k(index)+lag(2)) + ...
                        RFAs.modal.A5*(1i*k(index))/(1i*k(index)+lag(3)) + ...
                        RFAs.modal.A6*(1i*k(index))/(1i*k(index)+lag(4));

       Qqc(index,:,:) = RFAs.control.A0 + ...
                        RFAs.control.A1*(1i*k(index)) + ...
                        RFAs.control.A2*(1i*k(index))^2 + ...
                        RFAs.control.A3*(1i*k(index))/(1i*k(index)+lag(1)) + ...
                        RFAs.control.A4*(1i*k(index))/(1i*k(index)+lag(2)) + ...
                        RFAs.control.A5*(1i*k(index))/(1i*k(index)+lag(3)) + ...
                        RFAs.control.A6*(1i*k(index))/(1i*k(index)+lag(4));
    end
    
    % Plot modes (or DOF) RFAs and AICs for each reduced frequency
    figure(1)
    tiledlayout(nDOF,nDOF, 'TileSpacing', 'compact')

    for i=1:nDOF
        for j=1:nDOF
        nexttile; hold on; box on;
        plot(real(Qqq(:,i,j)),imag(Qqq(:,i,j)),'-k','linewidth',1.5)
        plot(real(AIC(:,i,j)),imag(AIC(:,i,j)),'ok','markersize',5.5)
        set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
        xlabel('Real','interpreter','latex','FontSize',12); ylabel('Imag','interpreter','latex','FontSize',12); 
        tit = ['$Q_{',num2str(i),num2str(j),'}$'];
        title(tit,'interpreter','latex','FontSize',12);
        end
    end
    lg = legend('RFA','AIC','interpreter','latex','FontSize',13);
    lg.Orientation = 'horizontal';
    lg.Layout.Tile = 'North';
    
    % Plot control RFAs and AICs for each reduced frequency
    figure(2)
    tiledlayout(1,nDOF, 'TileSpacing', 'compact')
    for j=1:nDOF
        nexttile; hold on; box on;
        plot(real(Qqc(:,j)),imag(Qqc(:,j)),'-k','linewidth',1.5)
        plot(real(AIC(:,j,nDOF+1)),imag(AIC(:,j,nDOF+1)),'ok','markersize',5.5)
        set(gca,'FontName','cmr12'); set(gca,'FontSize',10)
        xlabel('Real','interpreter','latex','FontSize',12); ylabel('Imag','interpreter','latex','FontSize',12); 
        tit = ['$Q_{3',num2str(j),'}$'];
        title(tit,'interpreter','latex','FontSize',12);
    end
    lg = legend('RFA','AIC','interpreter','latex','FontSize',13);
    lg.Orientation = 'horizontal';
    lg.Layout.Tile = 'North';
end