import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.integrate import solve_bvp

import parameter

#This module describes the assumed guesses and generates the boundary
#values
# def profile(g1):
#     #Solving model equations for partially wetted catalyst
    
#     wg=g1(1];
#     u=parameter.P[8]*np.exp((-parameter.P[11]*parameter.P[19]*parameter.P[3])/(parameter.P[6]*parameter.P[19]*wg)+parameter.P[5]*(1-parameter.P[9]/(parameter.P[19]*wg)));
#     du=-parameter.P[2]*(1-parameter.P[6]*u)*G(1);
#     dwg=(parameter.P[1]*G(1)*parameter.P[3])/(parameter.P[7]*parameter.P[19])*(parameter.P[9]+parameter.P[4]/(parameter.P[1]*parameter.P[3])-(parameter.P[19]*wg-exp(parameter.P[12]*parameter.P[1]))/(1-exp(parameter.P[12]*parameter.P[1])));
    
#     sol = ode15s(@pelleteqns,[0 1],[u du wg dwg]);
#     s=linspace(0,max(sol.x));
    
#     u=deval(sol,s);
    
#     error=sqrt((u(1,100)-1)^2+(u(3,100)-1)^2);
    
#     return error


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
        
        return np.vstack((y[1],
                         parameter.P[17]*y[0]*parameter.P[2]*G[0]*y[1]+G[1]*y[1]*np.exp(parameter.P[14]*(1-1/y[2]))*(1-parameter.P[18]*y[0])-y[1]*parameter.P[18]*(y[1]-parameter.P[17]*G[0]*parameter.P[2]*y[0])/(1-parameter.P[18]*y[0]),
                         y[3],
                         -parameter.P[15]*G[1]*y[1]*np.exp(parameter.P[14]*(1-1/y[3]))))
        


##Boundary conditions for vapor phase temperature and normalized mole
##fraction
def boundaryconditions(ya,yb):

    return np.array([ya[1]-(parameter.P[8]*np.exp((-parameter.P[11]*parameter.P[19]*parameter.P[3])/(parameter.P[6]*parameter.P[19]*wg)+parameter.P[5]*(1-parameter.P[9]/(parameter.P[19]*ya[3])))),#mole frac at interface
                     ya[2]+parameter.P[2]*(1-parameter.P[6]*u)*G[0],#mole frac gradient  at interface
                     yb[3]-1,# temperatuer at non wetted end
                     yb[1]-1,# mole fraction at non wetted end
                    ])



#def function harold_model(G):
    
#P=[Ph,Pm,C,V,H,yba,k,uo,woL,pf,pc,p,z,gamma,beta,v,Y1,Y2,T]
#---1, 2, 3,4,5,6,  7,8, 9,  10, 1,2,3,4,    5,   6,7, 8, 9

#Initializing guess vlaues
#--------------------------------------------------------------------------
#options=optimset('MaxFunEvals',2000];

#UB=(pc*T*C/yba+H*woL]/(H+log(uo*yba));
#[xval,fval,exitflag]=fminsearch(@profile,0.98,options];
#wg=1;
#Intial guess
G = []
G.append(100000)
G.append(0.1)


x = np.linspace(0.001,1,100)
y = np.ones((4,x.size))

wg = 1
u=parameter.P[8]*np.exp((-parameter.P[11]*parameter.P[19]*parameter.P[3])/(parameter.P[6]*parameter.P[19]*wg)+parameter.P[5]*(1-parameter.P[9]/(parameter.P[19]*wg)));
du=-parameter.P[2]*(1-parameter.P[6]*u)*G[0];
dwg=(parameter.P[1]*G[0]*parameter.P[3])/(parameter.P[7]*parameter.P[19])*(parameter.P[9]+parameter.P[4]/(parameter.P[1]*parameter.P[3])-(parameter.P[19]*wg-np.exp(parameter.P[12]*parameter.P[1]))/(1-np.exp(parameter.P[12]*parameter.P[1])));


solution = solve_bvp(pelleteqns,boundaryconditions, x, y, verbose=2)

#Storing data in a dataframe
df = pd.DataFrame()
df['r'] = x
df['C'] = parameter.Co*(1 - solution.sol(x)[0])
df['dC'] = solution.sol(x)[1]
df['T'] = solution.sol(x)[2]
df['dT'] = solution.sol(x)[3]

fig,axes = plt.subplots(2,2)

#plotting subplots
sns.lineplot(data=df,x = 'r', y= 'C', ax = axes[0,0])
sns.lineplot(data=df,x = 'r', y= 'dC', ax = axes[0,1])
sns.lineplot(data=df,x = 'r', y= 'T', ax = axes[1,0])
sns.lineplot(data=df,x = 'r', y= 'dT', ax = axes[1,1])
plt.show()
#sol = ode45(@pelleteqns,[0 1],[u du wg dwg]);

#s=linspace(0,max(sol.x));

# u=deval(sol,s);

# error=sqrt((u(1,100)-1)^2+(u(3,100)-1)^2);

# dwg=u(4,1);

# dwl=(parameter.P[1]*G(1)*parameter.P[3]*(parameter.P[19]*wg-1))/(1-exp(-parameter.P[12]*parameter.P[1]));

# Gcheck=(dwl-parameter.P[7]*parameter.P[19]*dwg)/(-parameter.P[4]+parameter.P[3]*parameter.P[1]*(parameter.P[9]+parameter.P[19]*wg));

#error=abs(Gcheck-G(1))/G(1);

#--------------------------------------------------------------------------

