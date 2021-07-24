# Dependencies
from math import pi
import parameters as params
from numpy import arccos, sqrt, array, complex, zeros, ones, linalg, conjugate, dot

def theodorsen(k):
    # Values of Theodordsen function parameters
    T1 = (params.c*arccos(params.c)) - (2+params.c**2)*sqrt(1-params.c**2)
    T2 = 0
    T3 = 0.25*params.c*(7+2*params.c**2)*sqrt(1-params.c**2)*arccos(params.c) - (1/8+params.c**2)*arccos(params.c)**2 - (1-params.c**2)*(5*params.c**2+4)/8
    T4 = params.c*sqrt(1-params.c**2) - arccos(params.c)
    T5 = 2*params.c*sqrt(1-params.c**2)*arccos(params.c) - (1-params.c**2) - arccos(params.c)**2
    T6 = 0
    T7 = params.c*(7+2*params.c**2)*sqrt(1-params.c**2)/8 - (1/8+params.c**2)*arccos(params.c)
    T8 = params.c*arccos(params.c) - 1/3*(1+2*params.c)*sqrt(1-params.c**2)
    T9 = 0.5*(params.a*T4 + sqrt(1-params.c**2)*(1-params.c**2)/3)
    T10 = sqrt(1-params.c**2) + arccos(params.c)
    T11 = (2-params.c)*sqrt(1-params.c**2) + (1-2*params.c)*arccos(params.c)
    T12 = (2+params.c)*sqrt(1-params.c**2) - (1+2*params.c)*arccos(params.c)
    T13 = 0.5*(T7 + (params.c-params.a)*T1)
    T14 = 0
    T15 = T4 + T10
    T16 = T1 - (params.c-params.a)*T4 - T8 + 0.5*T11
    T17 = -T1 + (params.a-0.5)*T4 - 2*T9
    T18 = T5 - T4*T10
    T19 = 0.5*T4*T11

    # Matrices
    Mnc = array([[-pi, pi*params.a, T1], [pi*params.a, -pi*(1/8+params.a**2), -2*T13], [T1, -2*T13, T3/pi]])        # Non-Circulatory Aerodynamic mass
    Knc = array([[0, 0, 0], [0, 0, -T15], [0, 0, -T18/pi]])                                                           # Non-Circulatory Aerodynamic stiffness
    Bnc = array([[0, -pi, -T4], [0, pi*(params.a-1/2), -T16], [0, -T17, -T19/pi]])                                  # Non-Circularoty Aerodynamic damping
    S1  = array([0, 1, T10/pi])
    S2  = array([0, (0.5-params.a), T11/(2*pi)])
    R   = array([[-2*pi], [2*pi*(1+0.5)], [-10]])

    # Theodorsen function for the k set of reduced frequency
    j = complex(0,1)
    Ck  = 0.5 + 0.0075/(j*k + 0.0455) + 0.10055/(j*k+0.3)

    z = zeros((len(k),params.nDOF,params.nDOF), dtype='complex')
    for i in range(len(k)):
        z[i,:,:]  = 2*(Mnc*(j*k[i])**2 + (Bnc+Ck[i]*R*S2)*(j*k[i]) + Knc+Ck[i]*R*S1)

    return z


def roger_RFA(aeroFunction,k,gamma):
    Aap = zeros((3+params.nLAG,params.nDOF,params.nDOF), dtype='complex') # this is the 3D matrix containing all matrices of the RFA approximation
    z = aeroFunction(k)
    j = complex(0,1)

    AIC = zeros(len(k))
    for r in range(params.nDOF):
        for s in range(params.nDOF):
            AIC = z[:,r,s][:,None]
            F = array([ones(len(k)), j*k, -k**2, 1/(j*k+gamma[0]), 1/(j*k+gamma[1]), 1/(j*k+gamma[2]), 1/(j*k+gamma[3])]).T

            aij = dot(linalg.inv(dot(F.T,conjugate(F))+dot(conjugate(F.T),F)),(dot(conjugate(F.T),conjugate(AIC))+dot(F.T,conjugate(AIC))))
            Aap[:,r,s] = aij[:,0]

    return Aap
