function [ssModel] = ssAeroelastic(params,fltcond)
%% DESCRIPTION
% State Space model for the aeroelastic system
% Author: Filipe Franca, Universidade Federal de Uberlandia

% The unsteady aerodynamics are approximated via the Roger's RFA approx.,
% using nLAG number of lag parameters, which introduces nLAG*nDOF more
% states to the system

% Input variables:
% params        - systems parameters, such as the modal matrices, AICs etc
% flcond        - flight conditions for the modal determination

% ------------------------------- IMPORTANT -------------------------------
% OBS: The algorithm is implemented to use only 4 lag parameters and is yet
% to be updated to take as an input the number of lag terms
% -------------------------------------------------------------------------

%% INPUT PARAMETERS

% General variables
nDOF = params.nDOF;         % number of modes or DOFs of the system
b    = params.b;            % characteristic lenght
lags = params.RFAs.lags;    % lag parameters used for the RFA
nLAG = params.RFAs.nLAG;    % number of lag terms used for the RFA
Voo  = fltcond.Voo;         % flight speed
rho  = fltcond.rho;         % air density
qoo  = fltcond.qoo;         % dynamic pressure

% Modal matrices
Mqq = params.Mqq;           % Modal mass matrix (or system mass)
Dqq = zeros(size(Mqq));           % Modal damping matrix
Kqq = params.Kqq;           % Modal stiffness
% Mqc = params.Mqc;           % inertial coupling matrix btw the elastic modes and control surface
Mqc = zeros(nDOF,1);

% Aerodynamic approximation matrices
RFAs = params.RFAs;

%% AEROELASTIC MODEL
% Dynamic equation matrices
M = Mqq - qoo*RFAs.modal.A2*(b/Voo)^2;
D = Dqq - qoo*RFAs.modal.A1*(b/Voo);
K = Kqq - qoo*RFAs.modal.A0;

% Mc = Mqc - qoo*RFAs.control.A2*(b/Voo)^2;
Mc = zeros(nDOF,1);

% dynamic matrix
A = [   zeros(nDOF,nDOF)     eye(nDOF)       zeros(nDOF,nDOF)          zeros(nDOF,nDOF)          zeros(nDOF,nDOF)          zeros(nDOF,nDOF);
            -M\K               -M\D            qoo*inv(M)                  qoo*inv(M)               qoo*inv(M)                qoo*inv(M);
        zeros(nDOF,nDOF)   RFAs.modal.A3  -(Voo/b)*lags(1)*eye(nDOF)   zeros(nDOF,nDOF)          zeros(nDOF,nDOF)          zeros(nDOF,nDOF);
        zeros(nDOF,nDOF)   RFAs.modal.A4     zeros(nDOF,nDOF)       -(Voo/b)*lags(2)*eye(nDOF)   zeros(nDOF,nDOF)          zeros(nDOF,nDOF);
        zeros(nDOF,nDOF)   RFAs.modal.A5     zeros(nDOF,nDOF)          zeros(nDOF,nDOF)       -(Voo/b)*lags(3)*eye(nDOF)   zeros(nDOF,nDOF);
        zeros(nDOF,nDOF)   RFAs.modal.A6     zeros(nDOF,nDOF)          zeros(nDOF,nDOF)          zeros(nDOF,nDOF)     -(Voo/b)*lags(4)*eye(nDOF);
        ];

% input matrix
% B = [   zeros(nDOF,1)           zeros(nDOF,1)               zeros(nDOF,1);
%    qoo*M\RFAs.control.A0  qoo*(b/Voo)*M\RFAs.control.A1        -M\Mc;
%         zeros(nDOF,1)           RFAs.control.A3             zeros(nDOF,1);
%         zeros(nDOF,1)           RFAs.control.A4             zeros(nDOF,1);
%         zeros(nDOF,1)           RFAs.control.A5             zeros(nDOF,1);
%         zeros(nDOF,1)           RFAs.control.A6             zeros(nDOF,1)       ];

B = zeros(length(A),1);

% output matrix
% C = [ones(1,2*nDOF) zeros(1,nDOF*nLAG)].*eye(nDOF*(nDOF+nLAG));
C = eye(length(A));

% feedforward matrix
D = zeros(size(B));

% State Space model
ssModel = ss(A, B, C, D);
 
end