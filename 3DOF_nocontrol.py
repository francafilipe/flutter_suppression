# Importing dependencies
from scipy import *
from numpy import *

def model_3DOF(x,t):
    # Importing the model's parameters
    import parameters

    aero_Matrices = aerodynamics()

    return None