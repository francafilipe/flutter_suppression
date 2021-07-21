"""
This function defines all geometric and aerodynamic parameters for the 3DOF typical wing section
The section's 3DOFs are: heave (h), pitch angle (a) and flap deflection angle (d)

Author: Filipe França
"""

# Dependencies
from math import pi

freq_h = 1.54               # [Hz] Natural Frequency of the Heave Motion
freq_a = 1.54               # [Hz] Natural Frequency of the Pitch Motion
freq_b = 1.54               # [Hz] Natural Frequency of the Flap Deflection

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


