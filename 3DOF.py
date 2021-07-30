# Importing dependencies
from numpy.linalg.linalg import inv
from scipy import *
from numpy import *

def model_nocontrol(t,x,k):
    # Importing the model's parameters
    import parameters as params
    import aerodynamics

    # Call Roger RFA for the specific aerodynamic modelling or data input
    RFA = aerodynamics.roger_RFA(aerodynamics.theodorsen,k,params.gamma)

    # State-space model matrices
    mu = params.rho*(params.b**2)*(params.Voo**2)
    M = params.M - mu*RFA[2,:,:]*(params.b/params.Voo)**2
    B = -mu*RFA[1,:,:]*(params.b/params.Voo)
    K = params.K - mu*RFA[0,:,:]
    A = mu*RFA[3:params.nLAG+3,:,:]

    # state matrix
    state = zeros((2*params.nDOF+params.nDOF*params.nLAG,2*params.nDOF+params.nDOF*params.nLAG))
    state[0:params.nDOF,params.nDOF:2*params.nDOF] = identity(params.nDOF)                                                                                  # first row (identity in the nDOF-th position)                                                                                         # first column
    state[params.nDOF:2*params.nDOF,:] = concatenate((inv(M)*B, inv(M)*K, inv(M)*A[0,:,:], inv(M)*A[1,:,:], inv(M)*A[2,:,:], inv(M)*A[3,:,:]),axis=1)       # second row (model equation of motion)
    for nlag in range(params.nLAG):
        n = params.nDOF
        state[(2+nlag)*n:(3+nlag)*n,params.nDOF:2*params.nDOF] = identity(params.nDOF)                                  # second column (identities starting after second row)
        state[(2+nlag)*n:(3+nlag)*n,(2+nlag)*n:(3+nlag)*n] = -(params.Voo/params.b)*params.gamma[nlag]*identity(3)      # lag parameters

    return None

from math import pi

t = array([0, 10])
x = array([0, 0, 0, 0, 0, 0])

k = array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])

test = model_nocontrol(t,x,k)
