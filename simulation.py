# Importing dependencies
import matplotlib.pyplot as plt
from numpy.linalg import eig, matrix_rank
from model import *
from controller import *
from control import lqr, ctrb, ss
from control.matlab import initial
from scipy.integrate import odeint, solve_ivp
from numpy import *
from math import pi

# Simulation and Modelling parameters and conditions
Voo = 8                           # [m/s] Flight speed
V_flutter = 8.3                     # [m/s] Flutter speed

dt = 1e-3                           # [sec] Sampling Time
T  = 5                              # [sec] Simulation total time
k  = int(T/dt)                      # [-] Sampling indice (last value)

x = zeros((18,k))                   # State Vector
input = zeros(k)                    # Input (Control Action) Vector

x[2,0] = 0.261799                   # [rad] Initial Condition in pitch angle state (15°)

# Define State Space model
cont_sys = ss_matrices(Voo)         # Get continuous state space matrices
disc_sys = cont_sys.sample(dt, method='zoh', alpha=None)   # Sampling of state space matrices for discrete system 

# Calculate Control Law
K, S, E = LQR_(V_flutter,dt)

# Evaluate MF dynamics
y, t = MF_analysis(Voo,cont_sys,K,x[:,0],plot=False)

for i in range(k-1):
    # Define the control action value
    input[i] = 0.1#-dot(K,x[:,i])
    # Evolving the system dynamics
    x0 = x[:,i]                     # Initial condition
    sol  = solve_ivp(system,[0.0, dt],x0,args=(cont_sys,input[i],),method='RK45') # Run integral solver for the dyanmic solution
    x[:,i+1] = sol.y[:,-1]


# Show Results
figure1, plots = plt.subplots(3,2)

plots[0,0].plot(arange(0, T, dt), x[0,:])
plots[0,1].plot(arange(0, T, dt), x[1,:])
plots[1,0].plot(arange(0, T, dt), x[2,:])
plots[1,1].plot(arange(0, T, dt), x[3,:])
plots[2,0].plot(arange(0, T, dt), x[4,:])
plots[2,1].plot(arange(0, T, dt), x[5,:])

plt.suptitle('State Variables Response', fontweight='bold')
plots[0,0].set_ylabel('h [m]'); plots[1,0].set_ylabel('\N{GREEK SMALL LETTER ALPHA} [rad]'); plots[2,0].set_ylabel('\N{GREEK SMALL LETTER BETA} [rad]')
plots[0,1].set_ylabel('h [m/s]'); plots[1,1].set_ylabel('\N{GREEK SMALL LETTER ALPHA} [rad/s]'); plots[2,1].set_ylabel('\N{GREEK SMALL LETTER BETA} [rad/s]')
plots[2,0].set_xlabel('t [sec]'); plots[2,1].set_xlabel('t [sec]')

plt.show()
