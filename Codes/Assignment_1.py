#Assignment 1:
#Solve the RR and MR using two-shock and three-shock theory using Mach number relations.

import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from Shock_Functions import *

#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '16'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)

    
def RR(M1, P1, T1, p1, t):
    
    M2, P2P1, T2T1, p2p1 = M2_P2P1_T2T1_p2p1(M1, t, 'w')
    P2 = P2P1*P1
    T2 = T2T1*T1
    p2 = p2p1*p1

    print("M2:" + str(M2))
    print("P2:" + str(P2))
    print("T2:" + str(T2))
    print("p2:" + str(p2))

    return None


def MR(M0, P0, T0, p0, M1, P1, T1, p1, t):

    def funct(variables):
        (phi2, phi3, t2, t3, P2P1, P3P0) = variables
    
        eq1 = t3 + t2 - t
        eq2 = np.tan(t2)*np.tan(phi2) - ((M1*M1*np.sin(phi2)*np.sin(phi2) - 1)/(1.2*M1*M1 - M1*M1*np.sin(phi2)*np.sin(phi2) + 1))
        eq3 = np.tan(t3)*np.tan(phi3) - ((M0*M0*np.sin(phi3)*np.sin(phi3) - 1)/(1.2*M0*M0 - M0*M0*np.sin(phi3)*np.sin(phi3) + 1))
        eq4 = P2P1 - (7*M1*M1*np.sin(phi2)*np.sin(phi2) - 1)/6
        eq5 = P3P0 - (7*M0*M0*np.sin(phi3)*np.sin(phi3) - 1)/6
        eq6 = P2P1*P1 - P3P0*P0
    
        return [eq1, eq2, eq3, eq4, eq5, eq6]
        
    phi2, phi3, t2, t3, P2P1, P3P0 = fsolve(funct, (0.5, 1.57, 0.5, t-0.5, 10, 10))

    M2, P2P1, T2T1, p2p1 = M2_P2P1_T2T1_p2p1(M1, t2, 'w')
    P2 = P2P1*P1
    T2 = T2T1*T1
    p2 = p2p1*p1
    
    M3, P3P0, T3T0, p3p0 = M2_P2P1_T2T1_p2p1(M0, t3, 's')
    P3 = P3P0*P0
    T3 = T3T0*T0
    p3 = p3p0*p0
    
    print("M2:" + str(M2))
    print("P2:" + str(P2))
    print("T2:" + str(T2))
    print("p2:" + str(p2))
    print("θ_2:" + str(t2*180/np.pi) + "\n")
    
    print("M3:" + str(M3))
    print("P3:" + str(P3))
    print("T3:" + str(T3))
    print("p3:" + str(p3))
    print("θ_3:" + str(t3*180/np.pi))
    
    return t2, t3

    
M0 = float(input("Mach number:"))
print("Turn angle should be lesser than or equal to: " + str(theta_max(M0)*180/np.pi) + " degrees.")
t = float(input("Turn angle:"))*np.pi/180
P0 = float(input("Pressure:"))
T0 = float(input("Temperature:"))
p0 = float(input("Density:"))
print("\n")

if(t>theta_max(M0)):
    print("Detached shock")

else:
        
    M1, P1P0, T1T0, p1p0 = M2_P2P1_T2T1_p2p1(M0, t, 'w')
    P1 = P1P0*P0
    T1 = T1T0*T0
    p1 = p1p0*p0

    print("M1:" + str(M1))
    print("P1:" + str(P1))
    print("T1:" + str(T1))
    print("p1:" + str(p1) + "\n")

    check = theta_max(M1)

    if(check>=t):
        RR(M1, P1, T1, p1, t)

    else:
        t2, t3 = MR(M0, P0, T0, p0, M1, P1, T1, p1, t)
        














