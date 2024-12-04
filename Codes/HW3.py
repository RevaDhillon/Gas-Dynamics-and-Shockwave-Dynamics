#Homework 3: Shock Polars:

import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from Shock_Functions import *

#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '14'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)


def Transonic():
    
    theta1, y1 = ShockPolar(1.1)
    theta2, y2 = ShockPolar(1.2)
    theta3, y3 = ShockPolar(1.5)
    theta4, y4 = ShockPolar(2)

    plt.plot(theta1, y1, 'b', label = "M = " + str(1.1))
    plt.plot(-theta1, y1, 'b')

    plt.plot(theta2, y2, 'm', label = "M = " + str(1.2))
    plt.plot(-theta2, y2, 'm')

    plt.plot(theta3, y3, 'r', label = "M = " + str(1.5))
    plt.plot(-theta3, y3, 'r')

    plt.plot(theta4, y4, 'y', label = "M = " + str(2.0))
    plt.plot(-theta4, y4, 'y')

    plt.title("Shock Polar: Transonic")
    plt.xlabel("θ -->")
    plt.ylabel("$P_2/P_1$ -->")
    plt.legend(loc = (1.05, 0.25))
    plt.grid(linestyle = '--')
    plt.savefig("Transonic", dpi = 300, bbox_inches = 'tight')
    plt.show()

    return None


def HW2_RR():

    M0 = 3
    P0 = 10

    M1=1.9941316655645582
    P1=37.71257463082662
    t1 = 20
    
    M2 = 1.2037175733083607
    P2 = 107.1612355236284
    t2 = -20

    theta1, y1 = ShockPolar(M0)
    theta2, y2 = ShockPolar(M1)
    y2 = y2*(P1/P0)
    theta21 = t1 + theta2
    theta22 = t1 - theta2

    plt.plot(theta1, y1, 'b', label = "$P_1/P_0$")
    plt.plot(-theta1, y1, 'b')

    plt.plot(theta21, y2, 'm', label = "$P_2/P_0$")
    plt.plot(theta22, y2, 'm')

    plt.plot(t2 + t1, P2/P0, 'o')

    plt.title("Shock Polar: RR")
    plt.xlabel("θ -->")
    plt.ylabel("$P/P_0$ -->")
    plt.legend(loc = (1.05, 0.25))
    plt.grid(linestyle = '--')
    plt.savefig("RR", dpi = 300, bbox_inches = 'tight')
    plt.show()

    return None


def HW2_MR():

    M0 = 3
    P0 = 10

    M1 = 1.405933969110253
    P1 = 63.55884169515927
    t1 = 30
    
    M2 = 1.0471481037654962
    P2 = 100.88032736064744
    t2 = - 8.594

    M3 = 0.550915765313006
    P3 = 100.9993694037003
    t3 = 21.406

    theta1, y1 = ShockPolar(M0)
    
    theta2, y2 = ShockPolar(M1)
    y2 = y2*(P1/P0)
    theta21 = t1 + theta2
    theta22 = t1 - theta2

    plt.plot(theta1, y1, 'b', label = "$P_1/P_0$")
    plt.plot(-theta1, y1, 'b')

    plt.plot(theta21, y2, 'm', label = "$P_2/P_0$")
    plt.plot(theta22, y2, 'm')

    plt.plot(t2 + t1, P2/P0, '*', label = "$P_2/P_0$ Soln.")
    plt.plot(t3, P3/P0, 'x', label = "$P_3/P_0$ Soln.")

    plt.title("Shock Polar: MR")
    plt.xlabel("θ -->")
    plt.ylabel("$P/P_0$ -->")
    plt.legend(loc = (1.05, 0.25))
    plt.grid(linestyle = '--')
    plt.savefig("MR", dpi = 300, bbox_inches = 'tight')
    plt.show()

    return None


Transonic()
HW2_RR()
HW2_MR()


    
    

    
    









