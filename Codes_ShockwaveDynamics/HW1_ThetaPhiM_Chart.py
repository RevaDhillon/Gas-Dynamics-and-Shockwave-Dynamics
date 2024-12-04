#Homework 1: theta-phi-M chart:

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

#Input Quantities:
M = float(input("Mach number:"))

theta_max = theta_max(M)
theta_sonic = theta_sonic(M)

theta_list = np.linspace(0, theta_max, 500)
x1_list = np.empty(500)
x2_list = np.empty(500)

print("θ    ||    ɸ_weak    ||      ɸ_strong")

for i in range(500):
    theta = theta_list[i]
    x1_list[i], x2_list[i] = Shock_angles(theta, M)
    
    print("θ: " + str(theta*180/np.pi) + "\t" + "ɸ_weak: " + str(x1_list[i]*180/np.pi) + "\t" +  "ɸ_strong: " + str(x2_list[i]*180/np.pi)) 

plt.plot(theta_list*180/np.pi, x1_list*180/np.pi)
plt.plot(theta_list*180/np.pi, x2_list*180/np.pi)
plt.plot(theta_max*180/np.pi, (Shock_angles(theta_max, M)[0])*180/np.pi, 'o', label = "$θ_{max}$")
plt.plot([theta_sonic*180/np.pi, theta_sonic*180/np.pi], [Shock_angles(theta_sonic, M)[0]*180/np.pi, Shock_angles(theta_sonic, M)[1]*180/np.pi] , 'o', label = "$θ_{sonic}$")
plt.title("θ - ɸ - M Graph")
plt.xlabel("θ -->")
plt.ylabel("ɸ -->")
plt.legend()
plt.grid(linestyle = '--')
plt.savefig("thetaphiM", dpi = 300, bbox_inches = 'tight')
plt.show()


