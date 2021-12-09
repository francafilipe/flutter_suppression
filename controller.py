# Importing dependencies
import matplotlib.pyplot as plt
from numpy.linalg import eig, matrix_rank
from model import *
from control import lqr, ctrb, ss
from control.matlab import initial
from scipy.integrate import odeint, solve_ivp
from numpy import *
from math import pi


def LQR_(Voo,dt):
    # Define State Space model for flight condition
    cont_sys = ss_matrices(Voo)         # Get continuous state space matrices
    disc_sys = cont_sys.sample(dt, method='zoh', alpha=None)   # Sampling of state space matrices for discrete system 

    # Calculate Control Law
    Q = identity((18))
    R = 1
    K, S, E = lqr(cont_sys.A,cont_sys.B, Q, R)

    return K, S, E


def MF_analysis(Voo,sys,K,x0,plot=False):
    # Define Closed Loop State Space Model
    A = sys.A - dot(sys.B,K)        # system dynamics matrix
    B = zeros(sys.B.shape)          # input matrix
    C = sys.C                       # output matrix
    D = sys.D                       # feedfoward matrix

    MFsys = ss(A,B,C,D)

    y, t = initial(MFsys,T=5,X0=x0)

    if plot:
        # Show Results
        figure1, plots = plt.subplots(3,2)

        plots[0,0].plot(t, y[:,0])
        plots[0,1].plot(t, y[:,1])
        plots[1,0].plot(t, y[:,2])
        plots[1,1].plot(t, y[:,3])
        plots[2,0].plot(t, y[:,4])
        plots[2,1].plot(t, y[:,5])

        plt.suptitle('Closed Loop Dynamics', fontweight='bold')
        plots[0,0].set_ylabel('h [m]'); plots[1,0].set_ylabel('\N{GREEK SMALL LETTER ALPHA} [rad]'); plots[2,0].set_ylabel('\N{GREEK SMALL LETTER BETA} [rad]')
        plots[0,1].set_ylabel('h [m]'); plots[1,1].set_ylabel('alpha [rad/s]'); plots[2,1].set_ylabel('delta [rad/s]')
        plots[2,0].set_xlabel('t [sec]'); plots[2,1].set_xlabel('t [sec]')

        plt.show()

    return y, t

