function [output] = aeroMatrices(params,fltcond,k)
%% DESCRIPTION 
% Generalized aerodynamic forces for aeroelastic analysis
% Author: Filipe Franca, Universidade Federal de Uberlandia

% This function calculates the aerodynamics matrices of mass, damping and
% stiffness for the generalized dynamic aeroelastic system for a set of
% reduced frequencies k (input)

% The output variables are:
% A1        - Aerodynamic-mass matrices for each reduced frequency k
% A2        - Aerodynamic-damping matrices for each reduced frequency k
% A3        - Aerodynamic-stiffness matrices for each reduced frequency k

% The input variables are:
% k         - set of reduced frequencies at which the matrices are calculated, 
% params    - systems paramters 
% fltcond   - flight condition given for the simulation

%% 

% Input parameters
Voo = fltcond.V;        % Voo - Undisturbed flow speed, m/s
rho = fltcond.rho;      % rho - air density, kg/m3
b   = params.b;         % b   - airfoil half chord, m
a   = params.a;         % a   - location of the elastic axis from half chord positon relative to the half chord value , none
e   = params.e;         % e   - location of the control surface hinge axis from half chord positon relative to the half chord value , none

% Theodorsen function for the reduced-frequencies
Ck = theodorsen(k);

%% Generalized Forces for each reduced matrices

output = struct('k',{},'A1',{},'A2',{},'A3',{},'AIC',{});

for i=1:length(k)

% Alocating the reduced frequency in the output variable
output(i).k = k(i);

% ----------------- Oscilatory Aerodynamic Derivatives --------------------
% Lift
derivative.lift.h           = 0;
derivative.lift.hDot        = 2*pi*rho*Voo*b*Ck(i);
derivative.lift.hDot2       = pi*rho*(b^2);

derivative.lift.alpha       = 2*pi*rho*(Voo^2)*b*Ck(i);
derivative.lift.alphaDot    = 2*pi*rho*Voo*(b^2)*(0.5-a)*Ck(i) + pi*rho*Voo*(b^2);
derivative.lift.alphaDot2   = -pi*rho*a*(b^3);

derivative.lift.beta        = 2*rho*(Voo^2)*b*(sqrt(1-e^2) + acos(e))*Ck(i);
derivative.lift.betaDot     = rho*Voo*(b^2)*( ((1-2*e)*acos(e) + (2-e)*sqrt(1-e^2))*Ck(i) - (e*sqrt(1-e^2) - acos(e)) );
derivative.lift.betaDot2    = -rho*(b^3)*( e*acos(e) - (1/3)*(2 + e^2)*sqrt(1-e^2) );


% Pitching moment
derivative.pitchMoment.h            = 0;
derivative.pitchMoment.hDot         = pi*rho*Voo*(b^2) - 2*pi*rho*Voo*(b^2)*(0.5-(a+0.5)*Ck(i));
derivative.pitchMoment.hDot2        = pi*rho*a*(b^3);

derivative.pitchMoment.alpha        = pi*rho*(Voo^2)*(b^2) - 2*pi*rho*(Voo^2)*(b^2)*(0.5-(a+0.5)*Ck(i));
derivative.pitchMoment.alphaDot     = -2*pi*rho*Voo*(b^3)*(0.5-a)*(0.5-(a+0.5)*Ck(i));
derivative.pitchMoment.alphaDot2    = -pi*rho*(b^4)*(a^2 + 1/8);

derivative.pitchMoment.beta         = -rho*(Voo^2)*(b^2)*(e*sqrt(1-e^2)-acos(e)) - 2*rho*(Voo^2)*(b^2)*(sqrt(1-e^2)+acos(e))*(0.5-(a+0.5)*Ck(i));
derivative.pitchMoment.betaDot      = -rho*Voo*(b^3)*( (1/3)*(e^2-1)*sqrt(1-e^2) - (e-a)*(e*sqrt(1-e^2)-acos(e)) )  -  rho*Voo*(b^3)*( (1-2*e)*acos(e) + (2-e)*sqrt(1-e^2) )*( 0.5 - (a+0.5)*Ck(i) );
derivative.pitchMoment.betaDot2     = -rho*(b^4)*(  ((1/8)+e^2)*acos(e)  -  (1/8)*e*(7+2*e^2)*sqrt(1-e^2)  +  (e-a)*( (1/3)*(2+e^2)*sqrt(1-e^2) - e*acos(e) )   );


% Hinge moment
derivative.hingeMoment.h            = 0;
derivative.hingeMoment.hDot         = -rho*Voo*(b^2)*(e*sqrt(1-e^2)-acos(e)) - 2*rho*Voo*(b^2)*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e))*Ck(i) );
derivative.hingeMoment.hDot2        =  rho*(b^3)*(e*acos(e)-(1/3)*(2+e^2)*sqrt(1-e^2));

derivative.hingeMoment.alpha        = -rho*(Voo^2)*(b^2)*(e*sqrt(1-e^2)-acos(e)) - 2*rho*(Voo^2)*(b^2)*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e))*Ck(i) );
derivative.hingeMoment.alphaDot     =  rho*Voo*(b^3)*( a*(e*sqrt(1-e^2)-acos(e)) + (1/3)*(sqrt(1-e^2))^3 - (1/3)*(2+e^2)*sqrt(1-e^2) + e*acos(e) ) - 2*rho*Voo*(b^3)*(0.5-a)*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e))*Ck(i) );
derivative.hingeMoment.alphaDot2    = -rho*(b^4)*( (1/8+e^2)*acos(e) - (1/8)*e*(7+2*e^2)*sqrt(1-e^2) + (e-a)*( (1/3)*(2+e^2)*sqrt(1-e^2) - e*acos(e) )  );

derivative.hingeMoment.beta         = -(rho*(Voo^2)*(b^2)/pi)*( 2*e*sqrt(1-e^2)*acos(e) - (1-e^2) - (acos(e))^2 )  - (2/pi)*rho*(Voo^2)*(b^2)*(sqrt(1-e^2)+acos(e))*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e))*Ck(i) );
derivative.hingeMoment.betaDot      = -(1/pi)*rho*Voo*(b^3)*((1-2*e)*acos(e)+(2-e)*sqrt(1-e^2))*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e))*Ck(i) );
derivative.hingeMoment.betaDot2     =  (rho*(b^4)/pi)*( (1/4)*e*(7+2*e^2)*sqrt(1-e^2)*acos(e) - (1/8+e^2)*(acos(e))^2 - (1/8)*(1-e^2)*(5*e^2+4) );


% -------------------------- Aerodynamic matrices -------------------------

% Aerodynamic-mass matrix
A1 = [-derivative.lift.hDot2        -derivative.lift.alphaDot2         -derivative.lift.betaDot2;
       derivative.pitchMoment.hDot2  derivative.pitchMoment.alphaDot2   derivative.pitchMoment.betaDot2;
       derivative.hingeMoment.hDot2  derivative.hingeMoment.alphaDot2   derivative.hingeMoment.betaDot2];

% Aerodynamic-damping matrix
A2 = [-derivative.lift.hDot         -derivative.lift.alphaDot         -derivative.lift.betaDot;
       derivative.pitchMoment.hDot   derivative.pitchMoment.alphaDot   derivative.pitchMoment.betaDot;
       derivative.hingeMoment.hDot   derivative.hingeMoment.alphaDot   derivative.hingeMoment.betaDot];
 
% Aerodynamic-stiffness matrix
A3 = [-derivative.lift.h         -derivative.lift.alpha         -derivative.lift.beta;
       derivative.pitchMoment.h   derivative.pitchMoment.alpha   derivative.pitchMoment.beta;
       derivative.hingeMoment.h   derivative.hingeMoment.alpha   derivative.hingeMoment.beta];


% Alocating matrices in struct variable
output(i).A1 = A1;
output(i).A2 = A2;
output(i).A3 = A3;

% Calculating & alocating AIC matrix given the aerodynamics matrices
AIC = ((Voo/b)^2)*A1*(1i*k(i))^2 + (Voo/b)*A2*1i*k(i) + A3;
output(i).AIC = AIC;

end
end