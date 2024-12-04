
import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from Shock_Functions import *


#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '12'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)


def ShockPolar(M1):
    
    phi_list = np.linspace(np.arcsin(1/M1), np.pi/2)
    s = np.sin(phi_list)**2
    
    M2 = np.sqrt(((1 + 0.2*M1*M1)/(1.4*M1*M1*s - 0.2)) + (M1*M1*(1 - s)/(1 + 0.2*M1*M1*s)))
    p2p1 = 6*M1*M1*s/(M1*M1*s + 5)
    T2T1 = (5 + M1*M1*s)*(7*M1*M1*s - 1)/(36*M1*M1*s)
    P2P1 = (7*M1*M1*s - 1)/6

    y = P2P1
    a = 1+1.4*M1*M1
    b = (7*M1*M1 - 1)/6
    c = 1/6

    theta = np.arctan(((y-1)/(a - y))*np.sqrt((b - y)/(c + y)))
    theta2 = -theta
    theta = np.concatenate([np.flip(theta2), theta])
    y = np.concatenate([np.flip(y), y])
    T2T1 = np.concatenate([np.flip(T2T1), T2T1])
    p2p1 = np.concatenate([np.flip(p2p1), p2p1])
    M2 = np.concatenate([np.flip(M2), M2])
    theta = theta*180/np.pi

    return theta, y, T2T1, p2p1, M2


M = [1.2, 2.202, 3, 4, 5]

for i in range(5):
    theta, P2P1, T2T1, p2p1, M2 = ShockPolar(M[i])

    ts = theta_sonic(M[i])
    P2P1s = M2_P2P1(M[i], ts, 'w')[1]
    tm = theta_max(M[i])
    P2P1m = M2_P2P1(M[i], tm, 'w')[1]

    plot1 = plt.figure(1)
    plt.plot(theta, P2P1, label = "$M_1$ = "+str(M[i]))
    plt.plot(180*np.array([ts, tm])/np.pi, [P2P1s, P2P1m], '.')

    plot2 = plt.figure(2)
    plt.plot(theta, T2T1, label = "$M_1$ = "+str(M[i]))

    plot3 = plt.figure(3)
    plt.plot(theta, p2p1, label = "$M_1$ = "+str(M[i]))

    plot4 = plt.figure(4)
    plt.plot(theta, M2, label = "$M_1$ = "+str(M[i]))


plt.figure(1)
plt.title("Shock Polars: Pressure")
plt.ylabel("$ \\frac{P_2}{P_1} $ -->")
plt.xlabel("θ (Degrees) -->")
plt.grid(linestyle = "--")
plt.legend(bbox_to_anchor=(1, 0.5))
plt.savefig("Shock_Polars_P2P1", dpi = 300, bbox_inches = 'tight')

plt.figure(2)
plt.title("Shock Polars: Temperature")
plt.ylabel("$ \\frac{T_2}{T_1} $ -->")
plt.xlabel("θ (Degrees) -->")
plt.grid(linestyle = "--")
plt.legend(bbox_to_anchor=(1, 0.5))
plt.savefig("Shock_Polars_T2T1", dpi = 300, bbox_inches = 'tight')

plt.figure(3)
plt.title("Shock Polars: Density")
plt.ylabel("$ \\frac{p_2}{p_1} $ -->")
plt.xlabel("θ (Degrees) -->")
plt.grid(linestyle = "--")
plt.legend(bbox_to_anchor=(1, 0.5))
plt.savefig("Shock_Polars_density", dpi = 300, bbox_inches = 'tight')

plt.figure(4)
plt.title("Shock Polars: Downstream Mach number")
plt.ylabel("$ M_2 $ -->")
plt.xlabel("θ (Degrees) -->")
plt.grid(linestyle = "--")
plt.legend(bbox_to_anchor=(1, 0.5))
plt.savefig("Shock_Polars_M2", dpi = 300, bbox_inches = 'tight')
plt.show()













