#Assignment 2: Locate the shock in a nozzle

import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '16'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)

#Function for Newton Method
def Newton_method(f, df, x0, tol):
    
    #Set maximum number of iterations.
    n_max = 100000000
    
    while(abs(f(x0))>tol):
        #Compute the next better approximation using the formula.
        x1 = x0 - (f(x0)/df(x0))

        if abs(f(x1))<tol:
            return x1
        else:
            x0 = x1
            n_max -= 1
            
    #If no root is found, return None.
    return None


def plot_shock(AeAt, Pc1P0, Pc2P0):

    #Nozzle pressure ratio array.
    PbP0 = np.linspace(Pc2P0, Pc1P0)
    
    Me_array = np.empty(50)
    A1At = np.empty(50)
    
    x = sym.Symbol('x')
    
    for i in range(50):

        #To find exit Mach number.
        fMe = PbP0[i]*AeAt - 0.5787/(x*sym.sqrt(1 + 0.2*x*x))
        dfMe = fMe.diff(x)
        fMe = sym.lambdify(x, fMe)
        dfMe = sym.lambdify(x, dfMe)
        
        Me_array[i] = Newton_method(fMe, dfMe, 0.19, 0.001)
        
        #To find stagnation pressure ratios.
        p02p01 = PbP0[i]*((1+0.2*Me_array[i]*Me_array[i])**(3.5))

        #To find Mach number just before the shock.
        f = p02p01 - (((1.2*x*x)/(1 + 0.2*x*x))**3.5)*((6/(7*x*x - 1))**2.5)
        df = f.diff(x)
        f = sym.lambdify(x, f)
        df = sym.lambdify(x, df)

        M1 = Newton_method(f, df, 2.3, 0.001)

        #To find the Area ratio at the shock location.
        A1At[i] = ((1+0.2*M1*M1)**3)/(M1*1.728)

    plt.plot(PbP0, A1At)
    #Giving a title to the graph.
    plt.title("Shock Location: Area Ratio =" + str(AeAt))
    #Labeling the x-axis.
    plt.ylabel("$A_1$/$A_T$ -->")
    #Labeling the y-axis.
    plt.xlabel(" $P_b$/$P_{01}$ -->")
    #Setting grid properties.
    plt.grid(linestyle=':')
    plt.tight_layout()
    #Figure saved.
    plt.savefig("AreaRatio"+str(AeAt), dpi=300, bbox_inches="tight")
    plt.show()

    return PbP0, A1At
    
    
#Ae/At = 2:
AeAt = 2
Pc1P0 = 0.9372
Pc2P0 = 0.5132
PbP0_2, A1At_2 = plot_shock(AeAt, Pc1P0, Pc2P0)

#Ae/At = 3:
AeAt = 3
Pc1P0 = 0.9732
Pc2P0 = 0.3784
PbP0_2, A1At_2 = plot_shock(AeAt, Pc1P0, Pc2P0)

#Ae/At = 4:
AeAt = 4
Pc1P0 = 0.9851
Pc2P0 = 0.2956
PbP0_2, A1At_2 = plot_shock(AeAt, Pc1P0, Pc2P0)

#Ae/At = 5:
AeAt = 5
Pc1P0 = 0.9905
Pc2P0 = 0.2433
PbP0_2, A1At_2 = plot_shock(AeAt, Pc1P0, Pc2P0)

