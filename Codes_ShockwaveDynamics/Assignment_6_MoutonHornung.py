import sympy as sym
import numpy as np
from scipy.optimize import fsolve
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


theta1 = 20.726#float(input("Enter wedge angle in degrees:"))
theta1 = theta1*np.pi/180

M0 = 2.84#float(input("Enter inlet Mach number:"))
P0 = 287
p0 = 1
T0 = 1 

w =1.42#float(input("Enter length of the hypotenuse of the wedge w/H:"))
g = 1 - w*np.sin(theta1)#float(input("Enter height from line of symmetry to bottom corner of the wedge:"))
g_nd = g/w

M1, P1P0, T1T0, p1p0 = M2_P2P1_T2T1_p2p1(M0, theta1, 'w')
P1 = P1P0*P0
T1 = T1T0*T0
p1 = p1p0*p0


M2, theta2, P2, T2, p2, M3, delta, P3, T3, p3 = MR(M0, P0, T0, p0, M1, P1, T1, p1, theta1)
phi = (Shock_angles(theta2, M1))[0] - theta1
alpha = (Shock_angles(theta1, M0))[0]


Ar = (1/M3)*(((5+M3*M3)/6)**3)

vM2 = PM_angle(M2)


v_M2 = vM2 + delta
M_2 = M_expansion(v_M2)
mu_2 = np.arcsin(1/M_2)

vM1 = PM_angle(M1)

print(theta1 - delta)

def funct(variables):
    (M, delta_e, vM) = variables
    
    eq1 = -vM+ (np.sqrt(6)*np.arctan(np.sqrt((M*M -1)/6)) - np.arctan(np.sqrt(M*M -1)))
    eq2 = -vM + vM1 + delta+ delta_e
    eq3 = -1+np.tan(theta1 - delta - delta_e)*np.tan(theta1 - delta - delta_e + phi)*(2.4*M*M/(2*(M*M*np.sin(theta1 - delta - delta_e + phi)*np.sin(theta1 - delta - delta_e + phi)) - 1) - 1)
    return [eq1, eq2, eq3]

M_1, delta_e, vM_1 = fsolve(funct, (M1, theta1 - delta, vM1))
mu_1 = np.arcsin(1/M_1)


A = np.array([[np.sin(alpha), Ar, 0, 0, 0],
              [0, 1, np.sin(delta + mu_2), np.sin(theta1 + mu_1), 0],
              [-np.cos(alpha), (1-Ar)/np.tan(delta), np.cos(delta + mu_2), np.cos(theta1 + mu_1), 0],
              [-np.cos(alpha), 0, 0, np.cos(theta1 + mu_1), -np.cos(phi)],
              [-np.sin(alpha), 0, 0, np.sin(theta1 + mu_1), np.sin(phi)]])

b = np.array([g_nd + np.sin(theta1), g_nd, -np.cos(theta1), -np.cos(theta1), -np.sin(theta1)])
x = np.linalg.inv(A)
x = np.matmul(x, b)

s = Ar*x[1]*w


print("The Mach stem height Hm/H = %.2f" % s)
print(M3)

##plt.plot(theta1_list, P2P0, label = "$M_1$ = "+str(M1))
##plt.title("Shock Polars: Downstream Mach number")
##plt.ylabel("$ M_2 $ -->")
##plt.xlabel("θ (Degrees) -->")
##plt.grid(linestyle = "--")
##plt.legend(bbox_to_anchor=(1, 0.5))
##plt.savefig("Shock_Polars_M2", dpi = 300, bbox_inches = 'tight')
##plt.show()











