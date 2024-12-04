#Assignment 1: Flow through a converging nozzle with a finite reservoir

#Importing packages to use required functions.
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '16'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)

def plot(x, y, title, xlabel, ylabel, savefig):
    
    plt.plot(x, y, "-")
    #Giving a title to the graph.
    plt.title(title)
    #Labeling the x-axis.
    plt.ylabel(ylabel)
    #Labeling the y-axis.
    plt.xlabel(xlabel)
    #Setting grid properties.
    plt.grid(linestyle=':')
    plt.legend()
    plt.tight_layout()
    #Figure saved.
    plt.savefig(savefig, dpi=300, bbox_inches="tight")
    plt.show()
    
#Adiabatic reservoir drainage.

p0choked = np.linspace(1000000, 191000, 100000)
p0un = np.linspace(191000, 101000, 100000)
p0 = np.concatenate((p0choked, p0un))

rho0choked = p0choked**(1/1.4)
rho0un = p0un**(1/1.4)
rho0 = np.concatenate((rho0choked, rho0un))

mchoked = 0.040418*(p0choked**(6/7))
mun = 586.30465*np.sqrt(p0un**(2/7) - 26.9)
m = np.concatenate((mchoked, mun))

p0pb = p0/101000

t1 = []


for i in p0choked:
    t1.append((i**(-1/7) - 0.13895)*173.19)

for i in p0un:
    term1 = i**(1/7)/math.sqrt(26.9)
    term2 = math.sqrt((i**(2/7)/26.9) - 1)
    val = (term1**3)*term2 + 1.5*math.log(term1 + term2) + 1.5*term1*term2
    value = (4.94431 - val)/0.462996
    t1.append(value)

plot(t1, p0, "Reservoir Pressure VS Time (Adiabatic)", " Time -->", "$P_0$ -->", "PT_Adia")
plot(t1, rho0, "Reservoir Density VS Time (Adiabatic)", " Time -->", "Reservoir density -->", "rhoT_Adia")
plot(t1, m, "Mass Flow Rate VS Time (Adiabatic)", " Time -->", "Mass flow rate -->", "mT_Adia")
plot(p0pb, m, "Mass Flow Rate VS Pressure Ratio (Adiabatic)", " $P_0/P_b$ -->", "Mass flow rate -->", "mp0pb_AdiaT")


#Isothermal reservoir drainage.

p0choked = np.linspace(1000000, 191000, 100000)
p0un = np.linspace(191000, 101000, 100000)
p0 = np.concatenate((p0choked, p0un))

rho0choked = p0choked/114800
rho0un = p0un/114800
rho0 = np.concatenate((rho0choked, rho0un))

mchoked = 0.040418*(p0choked)/20
mun = 29.315*(p0un**(1/7))*np.sqrt(p0un**(2/7) - 26.9)
m = np.concatenate((mchoked, mun))

p0pb = p0/101000

t2 = []

for i in p0choked:
    t2.append(-math.log(i/1000000)*12.3707)
        
for i in p0un:
    term1 = (i**(2/7)-26.9033)
    e1 = 1 + 2*term1/(3*26.9033) + term1*term1/(5*26.9033*26.9033)
    e2 = math.sqrt(term1)*e1
    value = (e2 - 2.64469358895)/3.92089429
    val = 20.48 - value
    t2.append(val)

plot(t2, p0, "Reservoir Pressure VS Time (Isothermal)", " Time -->", "$P_0$ -->", "PT_IsoT")
plot(t2, rho0, "Reservoir Density VS Time (Isothermal)", " Time -->", "Reservoir density -->", "rhoT_IsoT")
plot(t2, m, "Mass Flow Rate VS Time (Isothermal)", " Time -->", "Mass flow rate -->", "mT_IsoT")
plot(p0pb, m, "Mass Flow Rate VS Pressure Ratio (Isothermal)", " $P_0/P_b$ -->", "Mass flow rate -->", "mp0pb_IsoT")


