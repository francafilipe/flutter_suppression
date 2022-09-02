function [RFAs] = rogerRFA(k,lags,nDOF,AIC)
%% DESCRIPTION & INTRODUCTION
% RFA for the unsteady aerodynamics of a 3DOF aeroelastic system
% Author: Filipe Franca, Universidade Federal de Uberlandia

% This function uses the Roger's method to approximate the reduced frequency
% dependent aerodynamic matrices using a rational function approximation (RFA)
% thus allowing a time domain modelling of the aeroelastic problem

% Input parameters:
% k         - reduced frequencies used to define the unsteady aerodynamics
% gamma     - array of lag parameters used in the RFA solution ( no_lags = length(gamma) )
% nDOF      - number of DOF or modes that are being used in the solution
% AIC       - reduced frequency dependent aerodynamics matrices
              % the matrix dimension is (i = reduced frequency, r = no of rows, s = no of columns)

% ------------------------------- IMPORTANT -------------------------------
% OBS: The algorithm is implemented to use only 4 lag parameters and is yet
% to be updated to take as an input the number of lag terms
% -------------------------------------------------------------------------

%% CALCULATE RFAs
nk     = length(k);                    % number of reduced frequencies used
nLAG   = length(lags);                % number of lag parameters used

% Auxiliar matrix used in the rational function aproximation
B = [ones(1,length(k));     1i*k;     -k.^2;     1i*k./(1i*k+lags(1));     1i*k./(1i*k+lags(2));     1i*k./(1i*k+lags(3));     1i*k./(1i*k+lags(4))].';
sumB = real(B)'*real(B) + imag(B).'*imag(B);

for r=1:nDOF
    for s=1:nDOF+1

        % Get the rs terms of all k-th AIC matrices
        Qrs = zeros(length(k),1);
        for i=1:length(k)
            Qrs(i,1) = AIC(i,r,s);
        end

        % Solve linear system for the terms of the RFA approximation
        sumQB = real(B)'*real(Qrs) + imag(B).'*imag(Qrs);
        ars = linsolve(sumB,sumQB);

        % Alocating the calculated terms given approximation in the RFA matrices
        A0(r,s) = ars(1);
        A1(r,s) = ars(2);
        A2(r,s) = ars(3);
        A3(r,s) = ars(4);
        A4(r,s) = ars(5);
        A5(r,s) = ars(6);
        A6(r,s) = ars(7);
    end
end

%% OUTPUT VARIABLE
% Elastic modes RFAs
RFAs.modal.A0   = A0(1:nDOF,1:nDOF);
RFAs.modal.A1   = A1(1:nDOF,1:nDOF);
RFAs.modal.A2   = A2(1:nDOF,1:nDOF);
RFAs.modal.A3   = A3(1:nDOF,1:nDOF);
RFAs.modal.A4   = A4(1:nDOF,1:nDOF);
RFAs.modal.A5   = A5(1:nDOF,1:nDOF);
RFAs.modal.A6   = A6(1:nDOF,1:nDOF);

% Control RFAs
RFAs.control.A0 = A0(1:nDOF,nDOF+1);
RFAs.control.A1 = A1(1:nDOF,nDOF+1);
RFAs.control.A2 = A2(1:nDOF,nDOF+1);
RFAs.control.A3 = A3(1:nDOF,nDOF+1);
RFAs.control.A4 = A4(1:nDOF,nDOF+1);
RFAs.control.A5 = A5(1:nDOF,nDOF+1);
RFAs.control.A6 = A6(1:nDOF,nDOF+1);

% General attribute for the RFA
RFAs.lags = lags;
RFAs.nLAG = nLAG;
RFAs.k    = k;
RFAs.nk   = nk;
end