# Importing Dependencies
from math import pi
from control import ss
import numpy
import parameters as params
from numpy import cos, sqrt, array, complex, zeros, ones, linalg, conjugate, dot, real, imag, identity, concatenate
from numpy.linalg.linalg import inv


def system(t,x,matrices,input):
    # Description of the function

    # Input variables
    xp = zeros((18,))               # state derivatives vector initialize
    A = matrices.A                  # state (or system) matrix
    B = matrices.B.flatten()        # input (or control) matrix
    u = input                       # input (or control) value

    # State Space equation for the system
    xp[:] = dot(A,x)+dot(B,u)
    
    return xp


def ss_matrices(Voo):
    # Description of the function 

    # Call Roger RFA for the specific aerodynamic modelling or data input
    RFA = roger_RFA(theodorsen,params.k,params.gamma)

    # State-space model matrices
    q = 0.5*params.rho*(Voo**2)
    M = params.M - 0.5*params.rho*(params.b**2)*RFA[2,:,:]
    B = -0.5*params.rho*params.b*Voo*RFA[1,:,:]
    K = params.K - q*RFA[0,:,:]
    A = q*RFA[3:params.nLAG+3,:,:]

    # state (or system) matrix
    state = zeros((2*params.nDOF+params.nDOF*params.nLAG,2*params.nDOF+params.nDOF*params.nLAG))                                                            # first column
    state[0:params.nDOF,params.nDOF:2*params.nDOF] = identity(params.nDOF)                                                                                  # first row (identity in the nDOF-th position)
    state[params.nDOF:2*params.nDOF,:] = concatenate((-dot(inv(M),K), -dot(inv(M),B), dot(inv(M),A[0,:,:]), dot(inv(M),A[1,:,:]), dot(inv(M),A[2,:,:]), dot(inv(M),A[3,:,:])),axis=1)       # second row (model equation of motion)
    for nlag in range(params.nLAG):
        n = params.nDOF
        state[(2+nlag)*n:(3+nlag)*n,params.nDOF:2*params.nDOF] = identity(params.nDOF)                                  # second column (identities starting after second row)
        state[(2+nlag)*n:(3+nlag)*n,(2+nlag)*n:(3+nlag)*n] = -(Voo/params.b)*params.gamma[nlag]*identity(3)      # lag parameters

    # input (or control) matrix
    input = zeros((2*params.nDOF+params.nDOF*params.nLAG,1))
    Ba = numpy.array([0, 0, (params.rb**2)*(params.omega_b**2)])
    input[params.nDOF:2*params.nDOF,0] = -dot(inv(M),Ba)

    # output matrix
    output = zeros((2*params.nDOF,2*params.nDOF+params.nDOF*params.nLAG))
    for k in range(2*params.nDOF):
        output[k,k] = 1

    # feedforward matrix
    feedforward = zeros((2*params.nDOF,1))

    matrices = ss(state,input,output,feedforward)

    return matrices


def roger_RFA(aeroFunction,k,gamma):
    j = complex(0,1)
    z = aeroFunction(k)
    print(z)
    B = -array([ones(len(k)), j*k, -k**2, j*k/(j*k+gamma[0]), j*k/(j*k+gamma[1]), j*k/(j*k+gamma[2]), j*k/(j*k+gamma[3])]).T
    Q = zeros(len(k))
    RFA = zeros((3+params.nLAG,params.nDOF,params.nDOF), dtype='complex') # this is the 3D matrix containing all matrices of the RallFA approximation

    sumB = dot(real(B).T,real(B)) + dot(imag(B).T,imag(B))

    for r in range(params.nDOF):
        for s in range(params.nDOF):
            Q = z[:,r,s][:,None]
            sumQB = dot(real(B).T,real(Q)) + dot(imag(B).T,imag(Q))

            aij = - dot(linalg.inv(sumB),sumQB)
            RFA[:,r,s] = aij[:,0]

    return RFA


def theodorsen(k):
    # Values of Theodordsen function parameters
    T1 = (params.c*cos(params.c)**(-1)) - (2+params.c**2)*sqrt(1-params.c**2)/3
    T2 = 0
    T3 = 0.25*params.c*(7+2*params.c**2)*sqrt(1-params.c**2)*cos(params.c)**(-1)  - (1/8+params.c**2)*(cos(params.c)**2)**(-1) - (1-params.c**2)*(5*params.c**2+4)/8
    T4 = params.c*sqrt(1-params.c**2) - cos(params.c)**(-1) 
    T5 = 2*params.c*sqrt(1-params.c**2)*cos(params.c)**(-1) - (1-params.c**2) - (cos(params.c)**2)**(-1)
    T6 = 0
    T7 = params.c*(7+2*params.c**2)*sqrt(1-params.c**2)/8 - (1/8+params.c**2)*cos(params.c)**(-1)
    T8 = params.c*cos(params.c)**(-1) - 1/3*(1+2*params.c**2)*sqrt(1-params.c**2)
    T9 = 0.5*(params.a*T4 + sqrt(1-params.c**2)*(1-params.c**2)/3)
    T10 = sqrt(1-params.c**2) + cos(params.c)**(-1)
    T11 = (2-params.c)*sqrt(1-params.c**2) + (1-2*params.c)*cos(params.c)**(-1)
    T12 = (2+params.c)*sqrt(1-params.c**2) - (1+2*params.c)*cos(params.c)**(-1)
    T13 = -0.5*(T7 + (params.c-params.a)*T1)
    T14 = 0
    T15 = T4 + T10
    T16 = T1 - (params.c-params.a)*T4 - T8 + 0.5*T11
    T17 = -T1 + (params.a-0.5)*T4 - 2*T9
    T18 = T5 - T4*T10
    T19 = -0.5*T4*T11

    # Matrices
    Mnc = array([[-pi, pi*params.a, T1], [pi*params.a, -pi*(1/8+params.a**2), -2*T13], [T1, -2*T13, T3/pi]])        # Non-Circulatory Aerodynamic mass
    Knc = array([[0, 0, 0], [0, 0, -T15], [0, 0, -T18/pi]])                                                         # Non-Circulatory Aerodynamic stiffness
    Bnc = array([[0, -pi, T4], [0, pi*(params.a-1/2), -T16], [0, -T17, -T19/pi]])                                   # Non-Circularoty Aerodynamic damping (T4 value should be negative, review later)
    S1  = array([0, 1, T10/pi])
    S2  = array([1, (0.5-params.a), T11/(2*pi)])                # (the first value should be 0 instead of 1 - review later)
    R   = array([[-2*pi], [2*pi*(params.a+0.5)], [-T12]])

    # Theodorsen function for the k set of reduced frequency
    j = complex(0,1)
    Ck  = 0.5 + 0.0075/(j*k + 0.0455) + 0.10055/(j*k+0.3)

    z = zeros((len(k),params.nDOF,params.nDOF), dtype='complex')
    for i in range(len(k)):
        z[i,:,:]  = 2*(Mnc*(j*k[i])**2 + (Bnc+Ck[i]*R*S2)*(j*k[i]) + Knc+Ck[i]*R*S1)

    return z
