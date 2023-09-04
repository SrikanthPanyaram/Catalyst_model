import numpy as np

import parameter


def bc(ya,yb):
    """
    BC1: Inlet concenteration is equal to Ci
    BC2: At outlet dC/dr = 0 as no further change in concentration C
    """
    return np.array([ya[0] - parameter.Co ,yb[1]])

def bc_nonisothermal(ya,yb):
    """
    BC1: Inlet concenteration is equal to Ci
    BC2: At outlet dC/dr = 0 as no further change in concentration C
    BC3: Inlet temperature is equal to Ti
    BC4: At outlet dT/dr = 0 as no further change in Temperature T
    """
    return np.array([ya[0] - parameter.Co ,yb[1],ya[2] - parameter.To,yb[3]])


def bc_nonisothermal_autocatalytic(ya,yb):
    """
    BC1: Inlet concenteration is equal to Ci
    BC2: At outlet dC/dr = 0 as no further change in concentration C
    BC3: Inlet temperature is equal to Ti
    BC4: At outlet dT/dr = 0 as no further change in Temperature T
    """
    return np.array([ya[0] - 0.1,yb[1],ya[2] - parameter.To,yb[3]])