"""
This function defines all geometric and aerodynamic parameters for the 3DOF typical wing section
The section's 3DOFs are: heave (h), pitch angle (a) and flap deflection angle (d)

Author: Filipe França
"""

# Dependencies
from math import pi
import numpy

# Simulation and modelling parameters and conditions
nDOF = 3                    # [-] Number of DOFs
nLAG = 4                    # [-] Number of lag aerodynamic terms
gamma = numpy.array([0.2, 0.4, 0.6, 0.8]) # [-] Aerodynamic lag term poles (RFA adjustment parameters)

k = numpy.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
                            # [-] Set of reduced frequencies


# Geometric and Inertial parameters
freq_h = 1.54               # [Hz] Natural Frequency of the Heave Motion
freq_a = 3.26               # [Hz] Natural Frequency of the Pitch Motion
freq_b = 8.01               # [Hz] Natural Frequency of the Flap Deflection

omega_h = freq_h*2*pi       # [rad/s] Natural Frequency of the Heave Motion
omega_a = freq_a*2*pi       # [rad/s] Natural Frequency of the Pitch Motion
omega_b = freq_b*2*pi       # [rad/s] Natural Frequency of the Flap Deflection

a = -0.5                    # [-] Position of the elastic axis in relation to the mid-chord value
c = 0.5                     # [-] Position of the flap hinge point in relation to the mid-chord value
b = 0.125                   # [-] Wing section mid-chord

xa = 0.66                   
xb = 0.0028
ra = 0.7303
rb = 0.0742

mu = 28.3467
rho = 0.002378              # [kg/m3] Complex Density
m  = pi*rho*mu*b**2         # [kg] Wing section mass

# System matrices
M = m*numpy.array([[1, xa, xb], [xa, ra**2, rb**2+xb*(c-a)], [xb, rb**2+xb*(c-a), rb**2]])          # Mass
K = m*numpy.array([[omega_h**2, 0, 0], [0, (ra**2)*(omega_a**2), 0], [0, 0, rb**2*omega_b**2]])     # Stiffness
# B = numpy.array([[25.0423, 6.2304, 0.0264], [6.2304, 4.7756, 0.0784], [0.0264, 0.0784, 0.0059]])    # Structural damping
B = numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])                                                  # Structural damping
