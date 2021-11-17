# Importing dependencies
import matplotlib.pyplot as plt
from model import *
from scipy.integrate import odeint 
from numpy import *
from math import pi

# Simulation and Modelling parameters and conditions
Voo = 8.0                           # [m/s] Flight speed

t = arange(0,5,0.0001)              # [sec] Time range
x = zeros(18)                       # State Vector Initial Condition
x[4] = -0.261799

# Get state space LTI system matrices
matrices = ss_matrices(Voo)

# Define the control action value
input = 0

# Run integral solver
y = odeint(system,x,t,args=(matrices,input,))

# Show Results
figure1, plots = plt.subplots(3,2)

plots[0,0].plot(t, y[:,0])
plots[0,1].plot(t, y[:,1])
plots[1,0].plot(t, y[:,2])
plots[1,1].plot(t, y[:,3])
plots[2,0].plot(t, y[:,4])
plots[2,1].plot(t, y[:,5])

plt.suptitle('State Variables Response\n (V_init = 7 m/s, alpha_init = 0°)', fontweight='bold')
plots[0,0].set_ylabel('h [m]'); plots[1,0].set_ylabel('\N{GREEK SMALL LETTER ALPHA} [rad]'); plots[2,0].set_ylabel('\N{GREEK SMALL LETTER BETA} [rad]')
plots[0,1].set_ylabel('h [m]'); plots[1,1].set_ylabel('alpha [rad/s]'); plots[2,1].set_ylabel('delta [rad/s]')
plots[2,0].set_xlabel('t [sec]'); plots[2,1].set_xlabel('t [sec]')

plt.show()
