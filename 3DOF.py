# Importing dependencies
from scipy import *
from numpy import *

def model_nocontrol(t,x,k):
    # Importing the model's parameters
    import parameters
    import aerodynamics

    # Call Roger RFA for the specific aerodynamic modelling or data input
    Aap = aerodynamics.roger_RFA(aerodynamics.theodorsen,k,parameters.gamma)
    #print(Aap)

    return None

from math import pi

t = array([0, 10])
x = array([0, 0, 0, 0, 0, 0])

k = array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])

test = model_nocontrol(t,x,k)
