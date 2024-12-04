#Assignment 3: Using shock polar solutions, find the flow deflection angle beta at the trailing edge of a symmetric diamond airfoil for angles of attack (0, 5, 10, -5, -10). theta=10 degrees.

import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
from Shock_Functions import *


#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '12'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)


M0 = 4.96
theta0, PP0_M0 = ShockPolar(M0)


def theta_u2(theta_l):
        
        Ml, P1P0_l = M2_P2P1(M0, theta_l, 'w')
        theta1_l, PP1_M1_l = ShockPolar(Ml)
        theta1_l = theta1_l + theta_l*180/np.pi
        P1P0_l_list = PP1_M1_l*P1P0_l

        a0 = 1+1.4*M0*M0
        b0 = (7*M0*M0 - 1)/6
        c0 = 1/6

        al = 1+1.4*Ml*Ml
        bl = (7*Ml*Ml - 1)/6
        cl = 1/6

        k = P1P0_l


        #For point A:
        def funct(variables):
                (y, theta) = variables
            
                eq1 = np.tan(theta) + (((y - 1)/(a0 - y))*np.sqrt((b0 - y)/(c0 + y))) 
                eq2 = np.tan(theta - theta_l) - (((y/k - 1)/(al - y/k))*np.sqrt((bl - y/k)/(cl + y/k))) 
                
                return [eq1, eq2]

               
        PaP0, theta_A = fsolve(funct, (28.535, 0.0001))
        ##

        #theta_u = 25*np.pi/180

        def fun(variables):
                (y, theta_u2, Mu, phi) = variables

                au = 1+1.4*Mu*Mu
                bu = (7*Mu*Mu - 1)/6
                cu = 1/6

                eq1 = Mu - np.sqrt(((1 + 0.2*M0*M0)/(1.4*M0*M0*(np.sin(phi))**2 - 0.2)) + (M0*M0*(1 - (np.sin(phi))**2)/(1 + 0.2*M0*M0*(np.sin(phi))**2)))
                eq2 = np.tan(theta_A - theta_u2) + (((PaP0/y - 1)/(au - PaP0/y))*np.sqrt((bu - PaP0/y)/(cu + PaP0/y)))
                eq3 = y - (7*(M0**2)*((np.sin(phi))**2) - 1)/6
                eq4 = np.tan(theta_u2)*np.tan(phi) - ((M0*M0*np.sin(phi)*np.sin(phi) - 1)/(1.2*M0*M0 - M0*M0*np.sin(phi)*np.sin(phi) + 1))
                
                return [eq1, eq2, eq3, eq4]

               
        PuP0, theta_u, M1_u, phi = fsolve(fun, (20, 0.5, 1.2, 1.2))


        return(theta_u)


def theta_u1(theta_l):

        Ml, PlP0 = M2_P2P1(M0, theta_l, 'w')
        theta1_l, PP1_M1_l = ShockPolar(Ml)
        theta1_l = theta1_l + theta_l*180/np.pi
        P1P0_l_list = PP1_M1_l*PlP0

        a0 = 1+1.4*M0*M0
        b0 = (7*M0*M0 - 1)/6
        c0 = 1/6

        al = 1+1.4*Ml*Ml
        bl = (7*Ml*Ml - 1)/6
        cl = 1/6
        

        def fun1(variables):
                (y, theta_u1, Mu, phi, theta_B) = variables

                au = 1+1.4*Mu*Mu
                bu = (7*Mu*Mu - 1)/6
                cu = 1/6
                PuP0 = (7*(M0**2)*((np.sin(phi))**2) - 1)/6

                eq1 = Mu - np.sqrt(((1 + 0.2*M0*M0)/(1.4*M0*M0*(np.sin(phi))**2 - 0.2)) + (M0*M0*(1 - (np.sin(phi))**2)/(1 + 0.2*M0*M0*(np.sin(phi))**2)))
                
                eq2 = np.tan(theta_B - theta_u1) + (((y/PuP0 - 1)/(au - y/PuP0))*np.sqrt((bu - y/PuP0)/(cu + y/PuP0)))
                
                eq3 = np.tan(theta_B - theta_l) - (((y/PlP0 - 1)/(al - y/PlP0))*np.sqrt((bl - y/PlP0)/(cl + y/PlP0)))
                
                eq4 = np.tan(theta_u1)*np.tan(phi) - ((M0*M0*np.sin(phi)*np.sin(phi) - 1)/(1.2*M0*M0 - M0*M0*np.sin(phi)*np.sin(phi) + 1))

                eq5 =(theta_max(Ml)+ theta_max(Mu)) - (theta_u1 - theta_l) - 0.05
                
                return [eq1, eq2, eq3, eq4, eq5]

       
        PbP0, theta_u1, Mu, phi, theta_B = fsolve(fun1, (50, 0.53, 2.2, 0.73, 0.09))


##        M1_u, P1P0_u = M2_P2P1(M0, theta_u1, 'w')
##        theta1_u, PP1_M1_u = ShockPolar(M1_u)
##        theta1_u= theta1_u + theta_u1*180/np.pi
##        P1P0_u_list = PP1_M1_u*P1P0_u
##
##        plt.plot(theta1_u, P1P0_u_list)
##        plt.plot(theta1_l, P1P0_l_list)
##        plt.plot(theta0, PP0_M0)
##        plt.show()
        
        return theta_u1


theta_l_list = -np.linspace(20, 33)*np.pi/180
theta_u2_list = np.empty(len(theta_l_list))
theta_u1_list = np.empty(len(theta_l_list))

for i in range(len(theta_l_list)):

        theta_u2_list[i] = theta_u2(theta_l_list[i])
        theta_u1_list[i] = theta_u1(theta_l_list[i])


plt.plot(-(180/np.pi)*theta_l_list, (180/np.pi)*theta_u2_list, 'o-', label = "von-Neumann transition line")
plt.plot(-(180/np.pi)*theta_l_list, (180/np.pi)*theta_u1_list, 'o-', label = "Detachment transition line")

plt.title("Dual domain solution: M = 4.96")
plt.ylabel("θ_upper (Degrees)-->")
plt.xlabel("- θ_lower (Degrees) -->")
plt.grid(linestyle = "--")
plt.legend()#bbox_to_anchor=(1, 0.5))
#plt.savefig("Shock_Polars_M2", dpi = 300, bbox_inches = 'tight')
plt.show()













