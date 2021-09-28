# Importing dependencies
from matplotlib.pyplot import *
from numpy.linalg.linalg import inv
from scipy.integrate import odeint 
from numpy import *

def model_nocontrol(x,t,k):
    # Importing the model's parameters
    import parameters as params
    import aerodynamics

    # Call Roger RFA for the specific aerodynamic modelling or data input
    RFA = aerodynamics.roger_RFA(aerodynamics.theodorsen,k,params.gamma)

    # State-space model matrices
    q = 0.5*params.rho*(params.Voo**2)
    M = params.M - 0.5*params.rho*(params.b**2)*RFA[2,:,:]
    B = -0.5*params.rho*params.b*params.Voo*RFA[1,:,:]
    K = params.K - q*RFA[0,:,:]
    A = q*RFA[3:params.nLAG+3,:,:]
    print(params.K,'\n',K)
    print(RFA[0,:,:])
    # state matrix
    state = zeros((2*params.nDOF+params.nDOF*params.nLAG,2*params.nDOF+params.nDOF*params.nLAG))                                                            # first column
    state[0:params.nDOF,params.nDOF:2*params.nDOF] = identity(params.nDOF)                                                                                  # first row (identity in the nDOF-th position)
    state[params.nDOF:2*params.nDOF,:] = concatenate((-dot(inv(M),K), -dot(inv(M),B), dot(inv(M),A[0,:,:]), dot(inv(M),A[1,:,:]), dot(inv(M),A[2,:,:]), dot(inv(M),A[3,:,:])),axis=1)       # second row (model equation of motion)
    for nlag in range(params.nLAG):
        n = params.nDOF
        state[(2+nlag)*n:(3+nlag)*n,params.nDOF:2*params.nDOF] = identity(params.nDOF)                                  # second column (identities starting after second row)
        state[(2+nlag)*n:(3+nlag)*n,(2+nlag)*n:(3+nlag)*n] = -(params.Voo/params.b)*params.gamma[nlag]*identity(3)      # lag parameters
    savetxt('state_matrix.txt', state, fmt='%10.4f')
    return dot(state,x)


from math import pi

# Define time range and initial conditions
t = arange(0,10,0.001)
x = zeros(18)
x[1] = 0.1

k = array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])

# Run integral solver
y = odeint(model_nocontrol,x,t,args=(k,))
