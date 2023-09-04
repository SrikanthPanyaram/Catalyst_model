import math

import numpy as np

import parameter

"""
Isothermal models
"""

"""
First order isothermal reaction 1D slab
"""
def nth_order_isothermal(x,y):
    """
    D(d2C/dr2) = k*C
    """
    return np.vstack((y[1], parameter.k/parameter.D*y[0]**parameter.n ))


def nth_order_non_isothermal(x,y):
    """
    D(d2C/dr2) = k*exp*(-E/RT)*C
    K(d2T/dr2) = k*exp(-E/RT)*C*delH
    """
    return np.vstack((y[1],
                      parameter.A/parameter.D*np.exp(-parameter.Ea/parameter.R/y[2])*y[0]**parameter.n,
                      y[3],
                      parameter.A/parameter.keff*parameter.delH*np.exp(-parameter.Ea/parameter.R/y[2])*y[0]))