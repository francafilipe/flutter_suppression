function [output] = aerodynamics3DOF(params,fltcond)
%% DESCRIPTION 
% Generalized aerodynamic forces for aeroelastic analysis
% Author: Filipe Franca, Universidade Federal de Uberlandia

% This function calculates the aerodynamics and lag matrices
% for the generalized dynamic aeroelastic system using Jones's exponential
% approximation for the Wagner function

% The output variables are:
% Ai        - Aerodynamic loads relative to DOFs and its derivatives 
% Bi        - Lag parameters derivatives relative to DOFs and its derivatives

% The input variables are:
% params    - systems paramters
% fltcond   - flight condition given for the simulation

%% 

% Input parameters
Voo = fltcond.V;        % Voo - Undisturbed flow speed, m/s
rho = fltcond.rho;      % rho - air density, kg/m3
b   = params.b;         % b   - airfoil half chord, m
a   = params.a;         % a   - location of the elastic axis from half chord positon relative to the half chord value , none
e   = params.e;         % e   - location of the control surface hinge axis from half chord positon relative to the half chord value , none

%% Generalized Forces for each reduced matrices

% ----------------- Oscilatory Aerodynamic Derivatives --------------------
% Lift
derivative.lift.h           = 0;
derivative.lift.hDot        = 2*pi*rho*Voo*b;
derivative.lift.hDot2       = pi*rho*(b^2);

derivative.lift.alpha       = 2*pi*rho*(Voo^2)*b;
derivative.lift.alphaDot    = 2*pi*rho*Voo*(b^2)*(0.5-a) + pi*rho*Voo*(b^2);
derivative.lift.alphaDot2   = -pi*rho*a*(b^3);

derivative.lift.beta        = 2*rho*(Voo^2)*b*(sqrt(1-e^2) + acos(e));
derivative.lift.betaDot     = rho*Voo*(b^2)*( ((1-2*e)*acos(e) + (2-e)*sqrt(1-e^2)) - (e*sqrt(1-e^2) - acos(e)) );
derivative.lift.betaDot2    = -rho*(b^3)*( e*acos(e) - (1/3)*(2 + e^2)*sqrt(1-e^2) );


% Pitching moment
derivative.pitchMoment.h            = 0;
derivative.pitchMoment.hDot         = pi*rho*Voo*(b^2) - 2*pi*rho*Voo*(b^2)*(0.5-(a+0.5));
derivative.pitchMoment.hDot2        = pi*rho*a*(b^3);

derivative.pitchMoment.alpha        = pi*rho*(Voo^2)*(b^2) - 2*pi*rho*(Voo^2)*(b^2)*(0.5-(a+0.5));
derivative.pitchMoment.alphaDot     = -2*pi*rho*Voo*(b^3)*(0.5-a)*(0.5-(a+0.5));
derivative.pitchMoment.alphaDot2    = -pi*rho*(b^4)*(a^2 + 1/8);

derivative.pitchMoment.beta         = -rho*(Voo^2)*(b^2)*(e*sqrt(1-e^2)-acos(e)) - 2*rho*(Voo^2)*(b^2)*(sqrt(1-e^2)+acos(e))*(0.5-(a+0.5));
derivative.pitchMoment.betaDot      = -rho*Voo*(b^3)*( (1/3)*(e^2-1)*sqrt(1-e^2) - (e-a)*(e*sqrt(1-e^2)-acos(e)) )  -  rho*Voo*(b^3)*( (1-2*e)*acos(e) + (2-e)*sqrt(1-e^2) )*( 0.5 - (a+0.5) );
derivative.pitchMoment.betaDot2     = -rho*(b^4)*(  ((1/8)+e^2)*acos(e)  -  (1/8)*e*(7+2*e^2)*sqrt(1-e^2)  +  (e-a)*( (1/3)*(2+e^2)*sqrt(1-e^2) - e*acos(e) )   );


% Hinge moment
derivative.hingeMoment.h            = 0;
derivative.hingeMoment.hDot         = -rho*Voo*(b^2)*(e*sqrt(1-e^2)-acos(e)) - 2*rho*Voo*(b^2)*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e)) );
derivative.hingeMoment.hDot2        =  rho*(b^3)*(e*acos(e)-(1/3)*(2+e^2)*sqrt(1-e^2));

derivative.hingeMoment.alpha        = -rho*(Voo^2)*(b^2)*(e*sqrt(1-e^2)-acos(e)) - 2*rho*(Voo^2)*(b^2)*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e)) );
derivative.hingeMoment.alphaDot     =  rho*Voo*(b^3)*( a*(e*sqrt(1-e^2)-acos(e)) + (1/3)*(sqrt(1-e^2))^3 - (1/3)*(2+e^2)*sqrt(1-e^2) + e*acos(e) ) - 2*rho*Voo*(b^3)*(0.5-a)*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e)) );
derivative.hingeMoment.alphaDot2    = -rho*(b^4)*( (1/8+e^2)*acos(e) - (1/8)*e*(7+2*e^2)*sqrt(1-e^2) + (e-a)*( (1/3)*(2+e^2)*sqrt(1-e^2) - e*acos(e) )  );

derivative.hingeMoment.beta         = -(rho*(Voo^2)*(b^2)/pi)*( 2*e*sqrt(1-e^2)*acos(e) - (1-e^2) - (acos(e))^2 )  - (2/pi)*rho*(Voo^2)*(b^2)*(sqrt(1-e^2)+acos(e))*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e)) );
derivative.hingeMoment.betaDot      = -(1/pi)*rho*Voo*(b^3)*((1-2*e)*acos(e)+(2-e)*sqrt(1-e^2))*( 0.5*(acos(e)-e*sqrt(1-e^2)) + ((1+e/2)*sqrt(1-e^2)-(e+0.5)*acos(e)) );
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


% ------------------------- Lag Parameters Matrices -----------------------
A4 = [-2*pi*Voo*b*rho             -2*pi*Voo*b*rho;
      2*pi*Voo*b^2*rho*(a + 1/2)  2*pi*Voo*b^2*rho*(a + 1/2);
      (-2*Voo*b^2*rho*((e/2 + 1)*(1 - e^2)^(1/2) - acos(e)*(e + 1/2))) (-2*Voo*b^2*rho*((e/2 + 1)*(1 - e^2)^(1/2) - acos(e)*(e + 1/2)))];

B1 = [-0.1650   (0.1650*b*(a - 1/2))    ((0.0825*b*((1 - e^2)^(1/2)*(e - 2) + acos(e)*(2*e - 1)))/(pi));
      -0.335    (0.335*b*(a - 1/2))     ((0.1675*b*((1 - e^2)^(1/2)*(e - 2) + acos(e)*(2*e - 1)))/(pi))];
    
B2 = [0     -(0.1650*Voo)     (-(0.1650*Voo*(acos(e) + (1 - e^2)^(1/2)))/(pi));
      0     -(0.335*Voo)      (-(0.335*Voo*(acos(e) + (1 - e^2)^(1/2)))/(pi))];
    
B3 = [0 0 0;
      0 0 0];
    
B4 = [-(0.041*Voo)/(b)          0;
            0           -(0.32*Voo)/(b)];  

% Alocating matrices in struct variable
output = struct('A1',A1,'A2',A2,'A3',A3,'A4',A4,'B1',B1,'B2',B2,'B3',B3,'B4',B4);

end