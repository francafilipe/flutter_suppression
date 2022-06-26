function [ssComplete,ssAeroelastic,ssActuator] = stateSpaceControl(params,fltcond,simcond)
%% DESCRIPTION
% State Space model for the aeroelastic system of a typical section of a
% wing for a 2 modes (plunge and pitch) and a control surface action
% Author: Filipe Franca, Universidade Federal de Uberlandia

% The unsteady aerodynamics are approximated via the Roger's RFA approx.,
% using nLAG number of lag parameters, which introduces nLAG*nDOF more
% states to the system

%% 

% Input Parameters
nDOF = params.nDOF;
b    = params.b;
Voo  = fltcond.V;
nLAG = simcond.nLAG;
gamma = simcond.gamma;

% Call Roger RFA for the specific aerodynamic modelling or data input
[RFAs] = rogerRFA(params,fltcond,simcond);


%% 
%  Transforming 3DOF model and aerodynamic matrices to its 2DOF and control
%  input equivalents

% 2DOF Model matrices
Mmode = params.M(1:2,1:2);
Bmode = params.B(1:2,1:2);
Kmode = params.K(1:2,1:2);

% 2DOF Aerodynamics matrices
As0 = RFAs.A0(1:2,1:2);
As1 = RFAs.A1(1:2,1:2);
As2 = RFAs.A2(1:2,1:2);
As3 = RFAs.A3(1:2,1:2);
As4 = RFAs.A4(1:2,1:2);
As5 = RFAs.A5(1:2,1:2);
As6 = RFAs.A6(1:2,1:2);

% Control surface Model matrices
Mc = params.M(1:2,3);

% Control surface Aerodynamics matrices
Ac0 = RFAs.A0(1:2,3);
Ac1 = RFAs.A1(1:2,3);
Ac2 = RFAs.A2(1:2,3);
Ac3 = RFAs.A3(1:2,3);
Ac4 = RFAs.A4(1:2,3);
Ac5 = RFAs.A5(1:2,3);
Ac6 = RFAs.A6(1:2,3);

%%
%  Equations of Motions for a 2DOF sys with a control surface input

% Sys matrices corrected for aerodynamics approximation terms
M   = Mmode - As2*(params.b/fltcond.V)^2;     % Mass matrix (Structural Mass - Aerodynamics) - terms that multiply x_dot_dot
B   = Bmode - As1*(params.b/fltcond.V);       % Damping matrix - terms that multiply x_dot
K   = Kmode - As0;                            % Stiffness matrix  - terms that multiply x

% State Space model matrices for the dynamic system
Aae = [zeros(2,2)       eye(2)          zeros(2,2)          zeros(2,2)       zeros(2,2)          zeros(2,2);
       -M\K             -M\B            M\As3               M\As4            M\As5               M\As6;
       zeros(2,2)       eye(2)  -(Voo/b)*gamma(1)*eye(2)    zeros(2,2)       zeros(2,2)          zeros(2,2);
       zeros(2,2)       eye(2)          zeros(2,2)  -(Voo/b)*gamma(2)*eye(2) zeros(2,2)          zeros(2,2);
       zeros(2,2)       eye(2)          zeros(2,2)          zeros(2,2)  -(Voo/b)*gamma(3)*eye(2) zeros(2,2);
       zeros(2,2)       eye(2)          zeros(2,2)          zeros(2,2)       zeros(2,2)  -(Voo/b)*gamma(4)*eye(2)   ];
   
Bae = [zeros(2,1)        zeros(2,1)       zeros(2,1);
      -M\Ac0   -M\Ac1*(params.b/fltcond.V) -M\(Mc-Ac2*(params.b/fltcond.V)^2);
       zeros(2,1)        ones(2,1)        zeros(2,1);
       zeros(2,1)        ones(2,1)        zeros(2,1);
       zeros(2,1)        ones(2,1)        zeros(2,1);
       zeros(2,1)        ones(2,1)        zeros(2,1)       ];

% Output (ss model)
ssAeroelastic = ss(Aae,Bae,eye(12),zeros(12,3));

%%
%  Actuators Dynamics

w = params.omegaBeta;

% State Space model matrices for the actuators dynamics
Aac = [0         1         0;
       0         0         1;
      -27*w^3   -13.5*w^2  -4.5*w^2 ];

Bac = [0; 0; 27*w^3];

% Output (ss model)
ssActuator = ss(Aac,Bac,eye(3),zeros(3,1));

%% 
%  Complete State Space model 

%  State (or system) matrix
state  = [Aae Bae; zeros(3,12) Aac];
input  = [zeros(12,1); Bac];
output = eye(15);
feedfwd = zeros(15,1);

% Output (ss model)
ssComplete = ss(state,input,output,feedfwd);

end
