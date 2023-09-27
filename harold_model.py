import math
import model
import boundary

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.integrate import solve_bvp
from scipy.optimize import minimize

import parameter

#Intial guess
G = []
G.append(1)
G.append(5)

wg = 0.98
u = 1.0


#This module describes all the transport equations involved
def pelleteqns(x,y):
        #global G
        #dy=zeros(4,1);#arrays to store the diff eqns
        #Mass balance
        #dy(1)=y(2);
        #dy(2)=parameter.P[17]*parameter.P[2]*G(1)*y(2)+G(2)*y(1)*np.exp(parameter.P[14]*(1-1/y(3)))*(1-parameter.P[18]*y(1))-y(2)*parameter.P[18]*(y(2)-parameter.P[17]*G(1)*parameter.P[2]*y(1))/(1-parameter.P[18]*y(1));
        #Energy Balance
        #dy(3)=y(4);
        #dy(4)=-parameter.P[15]*G(2)*y(1)*np.exp(parameter.P[14]*(1-1/y(3)));
       
        # return np.vstack((y[1],
        #               parameter.A/parameter.D*np.exp(-parameter.Ea/parameter.R/y[2])*y[0]**parameter.n,
        #               y[3],
        #               parameter.A/parameter.keff*(-parameter.delH)*y[0]**parameter.n))#*np.exp(-parameter.Ea/parameter.R/y[2])*y[0]))
       
        
        return np.vstack((y[1],
                         parameter.P[17]*parameter.P[2]*G[0]*y[1] + G[1]*y[0]*np.exp(parameter.P[14]*(1-1/y[2]))*(1-parameter.P[18]*y[0])-y[1]*parameter.P[18]*(y[1]-parameter.P[17]*G[0]*parameter.P[2]*y[0])/(1-parameter.P[18]*y[0]),
                         y[3],
                         -parameter.P[15]*G[1]*y[0]*np.exp(parameter.P[14]*(1-1/y[2]))))
       
def bc_nonisothermal(ya,yb):
    """
    BC1: Inlet concenteration is equal to Ci
    BC2: At outlet dC/dr = 0 as no further change in concentration C
    BC3: Inlet temperature is equal to Ti
    BC4: At outlet dT/dr = 0 as no further change in Temperature T
    """
    return np.array([ya[0] - parameter.Co ,yb[1],ya[2] - 1,yb[3]]) 


def boundaryconditions1(ya,yb):
    return np.array([
                      ya[0]-(parameter.P[8]*np.exp((-parameter.P[11]*parameter.P[19]*parameter.P[3])/(parameter.P[6]*parameter.P[19]*wg)+parameter.P[5]*(1-parameter.P[9]/(parameter.P[19]*ya[2])))), #mole frac at interface
                      ya[1] + parameter.P[2]*(1-parameter.P[6]*u)*G[0], #mole frac gradient  at interface
                      yb[2] - 1, #temperatuer at non wetted end
                      yb[0] - 1, #mole fraction at non wetted end  
    ]
        
    )



def boundaryconditions(ya,yb):
#Boundary conditions
#ya[0]: Mole fraction at wetted end
#ya[1]: Mole fraction gradient
#ya[2]: Temperature at wetted end
#ya[3]: Temparture gradient at wetted end
    return np.array([ya[0] - (parameter.P[8]*np.exp((-parameter.P[11]*parameter.P[19]*parameter.P[3])/(parameter.P[6]*parameter.P[19]*ya[2])+parameter.P[5]*(1-parameter.P[9]/(parameter.P[19]*ya[1])))),#mole frac at interface
                     ya[1] + 0.1,#parameter.P[2]*(1-parameter.P[6]*ya[0])*G[0],#mole frac gradient  at interface
                     ya[2] - 1.1,# gas temperatuer at non wetted end 
                     #ya[3]-0.5,# gas tempaerature gradient #NEEDS TO be replaced with actual equation
                     ya[3] + 0.2#parameter.P[1]*G[0]*parameter.P[3]/(parameter.P[7]*parameter.P[19])*(parameter.P[9] + parameter.P[4]/(parameter.P[1]*parameter.P[3]) - (parameter.P[19]*ya[2] - np.exp(parameter.P[12]*parameter.P[1]))/(1-np.exp(parameter.P[12]*parameter.P[1])))
                    ])



#def function harold_model(G):
    
#P=[Ph,Pm,C,V,H,yba,k,uo,woL,pf,pc,p,z,gamma,beta,v,Y1,Y2,T]
#---1, 2, 3,4,5,6,  7,8, 9,  10, 1,2,3,4,    5,   6,7, 8, 9

#Initializing guess vlaues
#--------------------------------------------------------------------------
#options=optimset('MaxFunEvals',2000];

#UB=(pc*T*C/yba+H*woL]/(H+log(uo*yba));
#[xval,fval,exitflag]=fminsearch(@profile,0.98,options];
#wg=1;``


# def harold_model(G):
x = np.linspace(0.001,1,100)
y = np.ones((4,x.size))

wg = 1
u=parameter.P[8]*np.exp((-parameter.P[11]*parameter.P[19]*parameter.P[3])/(parameter.P[6]*parameter.P[19]*wg)+parameter.P[5]*(1-parameter.P[9]/(parameter.P[19]*wg)));
du=-parameter.P[2]*(1-parameter.P[6]*u)*G[0];
dwg=(parameter.P[1]*G[0]*parameter.P[3])/(parameter.P[7]*parameter.P[19])*(parameter.P[9]+parameter.P[4]/(parameter.P[1]*parameter.P[3])-(parameter.P[19]*wg-np.exp(parameter.P[12]*parameter.P[1]))/(1-np.exp(parameter.P[12]*parameter.P[1])));

#Initializing the guesses
y[0][:] = np.ones((1,x.size))#*u
y[1][:] = np.ones((1,x.size))#*du 
y[2][:] = np.ones((1,x.size))#*wg
y[3][:] = np.ones((1,x.size))#*dwg
#print(y[0])


#solution = solve_bvp(model.nth_order_non_isothermal ,boundary.bc_nonisothermal   boundaryconditions, x, y, verbose=2)
solution = solve_bvp(pelleteqns ,boundaryconditions1,x, y, verbose=2)

#Storing data in a dataframe
df = pd.DataFrame()
df['r'] = x
df['U'] = solution.sol(x)[0]
df['dU'] = solution.sol(x)[1]
df['W'] = solution.sol(x)[2]
df['dW'] = solution.sol(x)[3]

print(df.U) 
print(df.W)
yt = abs((df['U'].iloc[-1]**2 - 1) + (df['W'].iloc[-1]**2 -1)) 
print("The value of residual is", yt) 
    
    # return 100*abs((df['U'].iloc[-1]**2 - 1) + (df['W'].iloc[-1]**2 -1))
fig,axes = plt.subplots(2,2)

# plotting subplots
sns.lineplot(data=df,x = 'r', y= 'U', ax = axes[0,0])
sns.lineplot(data=df,x = 'r', y= 'dU', ax = axes[0,1])
sns.lineplot(data=df,x = 'r', y= 'W', ax = axes[1,0])
sns.lineplot(data=df,x = 'r', y= 'dW', ax = axes[1,1])
plt.show()



# result = minimize(harold_model,G)
# print(result)

