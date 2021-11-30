# Importing dependencies
import matplotlib.pyplot as plt
from model import *
from scipy.signal import cont2discrete
from scipy.integrate import odeint 
from numpy import *
from math import pi

# Simulation and Modelling parameters and conditions
Voo = 8.3                           # [m/s] Flight speed

dt = 1e-3                           # [sec] Sampling Time
T  = 5                              # [sec] Simulation total time
k  = int(T/dt)                      # [-] Sampling indice (last value)

x = zeros((18,k))                   # State Vector
input = zeros(k)                    # Input (Control Action) Vector

x[2,0] = 0.261799                   # [rad] Initial Condition in pitch angle state (15°)

# Define State Space model
cont_sys = ss_matrices(Voo)         # Get continuous state space matrices
disc_sys = cont2discrete((cont_sys.A,cont_sys.B,cont_sys.C,cont_sys.D), dt, method='zoh', alpha=None)   
                                    # Sampling of state space matrices for discrete system 

#for i in range(k-1):
# Define the control action value
input[0] = 0

# Evolving the system dynamics
x0 = x[:,0]                     # Initial condition
y  = odeint(system,x0,arange(0,T,dt),args=(disc_sys,input[0],)) # Run integral solver for the dyanmic solution
x = y[-1]


# Show Results
figure1, plots = plt.subplots(3,2)

plots[0,0].plot(arange(0, T, dt), x[0,:])
plots[0,1].plot(arange(0, T, dt), x[1,:])
plots[1,0].plot(arange(0, T, dt), x[2,:])
plots[1,1].plot(arange(0, T, dt), x[3,:])
plots[2,0].plot(arange(0, T, dt), x[4,:])
plots[2,1].plot(arange(0, T, dt), x[5,:])

plt.suptitle('State Variables Response\n (V_init = 7 m/s, alpha_init = 0°)', fontweight='bold')
plots[0,0].set_ylabel('h [m]'); plots[1,0].set_ylabel('\N{GREEK SMALL LETTER ALPHA} [rad]'); plots[2,0].set_ylabel('\N{GREEK SMALL LETTER BETA} [rad]')
plots[0,1].set_ylabel('h [m]'); plots[1,1].set_ylabel('alpha [rad/s]'); plots[2,1].set_ylabel('delta [rad/s]')
plots[2,0].set_xlabel('t [sec]'); plots[2,1].set_xlabel('t [sec]')

plt.show()
