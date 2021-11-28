# Importing dependencies
import matplotlib.pyplot as plt
from model import *
from scipy.integrate import odeint 
from numpy import *
from math import pi

# Simulation and Modelling parameters and conditions
Voo = 8.3                           # [m/s] Flight speed

Ts = 1e-3                           # [sec] Sampling Time
T  = 5                              # [sec] Simulation total time
k  = int(T/Ts)                      # [-] Sampling indice (last value)

x = zeros((18,k))                   # State Vector
input = zeros(k)                    # Input (Control Action) Vector

x[0,0] = 0.261799                   # [rad] Initial Condition in pitch angle state (15°)

# Get state space LTI system matrices
matrices = ss_matrices(Voo)

for i in range(k-1):
    # Define the control action value
    input[i] = 0

    # Evolving the system dynamics
    x0 = x[:,i]                     # Initial condition
    y  = odeint(system,x0,[0.0, Ts],args=(matrices,input,)) # Run integral solver for the dyanmic solution
    x[:,i+1] = y[-1]


# Show Results
figure1, plots = plt.subplots(3,2)

plots[0,0].plot(arange(0, T, Ts), x[0,:])
plots[0,1].plot(arange(0, T, Ts), x[1,:])
plots[1,0].plot(arange(0, T, Ts), x[2,:])
plots[1,1].plot(arange(0, T, Ts), x[3,:])
plots[2,0].plot(arange(0, T, Ts), x[4,:])
plots[2,1].plot(arange(0, T, Ts), x[5,:])

plt.suptitle('State Variables Response\n (V_init = 7 m/s, alpha_init = 0°)', fontweight='bold')
plots[0,0].set_ylabel('h [m]'); plots[1,0].set_ylabel('\N{GREEK SMALL LETTER ALPHA} [rad]'); plots[2,0].set_ylabel('\N{GREEK SMALL LETTER BETA} [rad]')
plots[0,1].set_ylabel('h [m]'); plots[1,1].set_ylabel('alpha [rad/s]'); plots[2,1].set_ylabel('delta [rad/s]')
plots[2,0].set_xlabel('t [sec]'); plots[2,1].set_xlabel('t [sec]')

plt.show()
