"""
Catalyst reaction-diffusion profiles
Author: Srikanth Panyaram
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
from scipy.integrate import solve_bvp

"""
Import models from the model file
"""
import model

"""
Import parameters from parameters file
""" 
import parameter

"""
Import boundary conditions
"""
import boundary

"""
Solver
"""


#Intial guess
x = np.linspace(0.001,1,100)
y = np.ones((4,x.size))

phi2 = parameter.phi2

#Solver
solution = solve_bvp(model.nth_order_non_isothermal,boundary.bc_nonisothermal,x,y,verbose=2)



#Storing data in a dataframe
df = pd.DataFrame()
df['r'] = x
df['C'] = solution.sol(x)[0]
df['dC'] = solution.sol(x)[1]
df['T'] = solution.sol(x)[2]
df['dT'] = solution.sol(x)[3]

fig,axes = plt.subplots(2,2)

#plotting subplots
sns.lineplot(data=df,x = 'r', y= 'C', ax = axes[0,0])
sns.lineplot(data=df,x = 'r', y= 'dC', ax = axes[0,1])
sns.lineplot(data=df,x = 'r', y= 'T', ax = axes[1,0])
sns.lineplot(data=df,x = 'r', y= 'dT', ax = axes[1,1])


# sns.lineplot(data=solution.sol(x))
plt.show()