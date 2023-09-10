
#Paramters

phi2    = 0.00005          #Thiele modulus
n       = 1                #order of the reaction
Co      = 5                #Initial concentration
To      = 300              #Initial Temperature
k       = 0.005            #Rate constant 1/s first orde reaction

D       = 0.001           #Diffusion coefficient
keff    = 1000             #Thermal diffusion coefficient
A       = 10          #Arrhenius coefficient
Ea      = 10000             #Activation energy
R       = 8.314            #Universal gas constant
delH    = -100            #heat of reaction

#Parameters for Harolds catalyst model
P = [];
P.append(0)#P[0] = 0
P.append(0.01)# P[1]=0.01#Ph-non dimensional heat transfer coeff in liquid
P.append(0.00001)# P[2]=0.00001#Pm-non dimensional mass tranfer coeff
P.append(10**(-5))# P[3]=10^-5#C-Coeff relating imbibation velocity and mass transfer flux at interface
P.append(10**(-7))# P[4]=1*10^-7#V-heat generation to removal ratio
P.append(5)# P[5]=5#H-non dimensional heat of vaporisation
P.append(0.1)# P[6]=0.1#yba-mole fraction of gas in bulk
P.append(1)# P[7]=1#k-ratio of gas to liquid filled thermal conductivities
P.append(8)# P[8]=1#uo-non dimensional vapor mole fraction
P.append(0.84)# P[9]=0.84#woL-non dimensional temperature
P.append(1)# P[10]=1#pf-dimensionless pressure of external liquid film
P.append(10)# P[11]=10#pc-dimenisonless capillary pressure
P.append(P[11] + P[10] -1)# P[12]=P[11]+P[10]-1#p=pc+pf-1%pressure difference
P.append(P[6]*P[2])# P[13]=P[6]*P[2]#z
P.append(5)# P[14]=5#gamma
P.append(0.05)# P[15]=0.05#beta
P.append(1)# P[16]=1#v
P.append(P[16]*P[6])# P[17]=P[16]*P[6]#Y1
P.append((1-P[16])*P[6])# P[18]=[1-P[16]]*P[6]#Y2
P.append(1.05)# P[19]=1.05#T