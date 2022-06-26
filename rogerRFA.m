function [RFAs] = rogerRFA(params,fltcond,simcond)
%% Description
% RFA for the unsteady aerodynamics of a 3DOF aeroelastic system
% Author: Filipe Franca, Universidade Federal de Uberlandia

% This function uses the Roger's method to approximate the reduced frequency
% dependent aerodynamic matrices using a rational function approximation (RFA)
% thus allowing a time domain modelling of the aeroelastic problem

% Input parameters:
% params  - struct contain the parameters of the system (geometric & mass carachteristics)
% fltcond - flight condition for which the simulation is performed (velocity and air density)
% simcond - simulation conditions, such as the set of reduced frequency (k) for which the 
            % aerodynamics are evaluated and the lag parameters (gamma) used for the approximation of the unsteady aerodynamic

%%
% Attribute variables from input parameters
k       = simcond.k;
gamma   = simcond.gamma;

% Aerodynamic matrices dependent on the reduced frequency
Qaero = aeroMatrices(params,fltcond,k);

% Auxiliar matrix used in the rational function aproximation
B = [ones(1,length(k));     1i*k;     -k.^2;     1i*k./(1i*k+gamma(1));     1i*k./(1i*k+gamma(2));     1i*k./(1i*k+gamma(3));     1i*k./(1i*k+gamma(4))].';
sumB = real(B)'*real(B) + imag(B).'*imag(B);

for r=1:params.nDOF
    for s=1:params.nDOF

        % Get the rs terms of all k-th AIC matrices
        Qrs = zeros(length(k),1);
        for i=1:length(k)
            Qrs(i,1) = Qaero(i).AIC(r,s);
        end

        % Solve linear system for the terms of the RFA approximation
        sumQB = real(B)'*real(Qrs) + imag(B).'*imag(Qrs);
        ars = linsolve(sumB,sumQB);

        % Alocating the calculated terms given approximation in the RFA matrices
        RFAs.A0(r,s) = ars(1);
        RFAs.A1(r,s) = ars(2);
        RFAs.A2(r,s) = ars(3);
        RFAs.A3(r,s) = ars(4);
        RFAs.A4(r,s) = ars(5);
        RFAs.A5(r,s) = ars(6);
        RFAs.A6(r,s) = ars(7);

    end
end


end