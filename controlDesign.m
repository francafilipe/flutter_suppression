function [gain] = controlDesign(params,fltcond,simcond)
%% DESCRIPTION
% Control law design function given a system in space state format
% Author: Filipe Franca, Universidade Federal de Uberlandia

%%

% Define control design conditon
fltcond.V = params.Vflutter;

% Load state-space model
[sys,~,~] = stateSpaceControl(params,fltcond,simcond);

sysControlability = ctrb(sys.A,sys.B);
sysRank    = rank(sysControlability);

% Discretizando o sitema
discreteSys = c2d(sys, simcond.samplingTime, 'zoh');

% Matriz de peso dos estados Qstate
Qstate = [100 100 100 100 1 1 1 1 1 1 1 1 1 1 1].*eye(15);

% Projetando o controlador de espaco de estados
[K,s,e] = dlqr(discreteSys.A, discreteSys.B, Qstate, 0);

% Output
gain = K;

end