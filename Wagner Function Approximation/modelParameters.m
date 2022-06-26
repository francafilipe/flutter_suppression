function params = modelParameters()

    % Geometric parameters
    b =  5;              % [m] half-chord value
    a = -0.15;           % [-] Position of the elastic axis in relation to the half-chord
    e =  0.6;            % [-] Position of the hinge axis in relation to the half-chord

    % Mass properties
    m      = 0.2642;
    Stheta = 0.3302;
    Itheta = 2.5627;
    Sbeta  = 0.0396;
    Ibeta  = 0.0528;

    % Uncopled natural frequencies and damping
    fHeave =  4.78;		% [Hz] Heave uncopled natural frequency
    fTheta = 10.20;     % [Hz] Pitch uncopled natural frequency
    fBeta  =  6.37;		% [Hz] Flap deflection uncopled natural frequency

    cHeave = 0.05;	 	% [-] Heave structural damping
    cTheta = 0.05;	 	% [-] Pitch structural damping
    cBeta  = 0.05;	 	% [-] Flap deflection structural damping


    omegaHeave = fHeave*2*pi;
    omegaTheta = fTheta*2*pi;
    omegaBeta  = fBeta*2*pi;

    % Mass, damping & stiffness matrices
    M = [	  m   		Stheta                  Sbeta
            Stheta      Itheta          Ibeta+b*(e-a)*Sbeta
            Sbeta   Ibeta+b*(e-a)*Sbeta         Ibeta 	     ]; 

    K = [ m*omegaHeave^2        0                   0
            0           Itheta*omegaTheta^2         0
            0                   0           Ibeta*omegaBeta^2  ];

    C = [ cHeave*(2*m*omegaHeave)           0                       0
            0               cTheta*(2*Itheta*omegaTheta)            0
            0                               0           cBeta*(2*Ibeta*omegaBeta)    ];

    % Output parameters
    params.b = b;
    params.a = a;
    params.e = e;
    params.M = M;
    params.B = C;
    params.K = K;
    params.nDOF = 3;

end