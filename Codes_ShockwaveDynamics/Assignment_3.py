#Assignment 3: Using shock polar solutions, find the flow deflection angle beta at the trailing edge of a symmetric diamond airfoil for angles of attack (0, 5, 10, -5, -10). theta=10 degrees.

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

M0 = 3
theta = 10*np.pi/180
alpha = -10*np.pi/180

#Shock polar for M0:
theta_polar1, P1P0_list = ShockPolar(M0)
ts = theta_sonic(M0)
P2P1s = M2_P2P1(M0, ts, 'w')[1]
tm = theta_max(M0)
P2P1m = M2_P2P1(M0, tm, 'w')[1]
plt.plot(180*np.array([ts, tm])/np.pi, [P2P1s, P2P1m], '.')

#Solutions in 1 and 2:
theta1 = -theta + alpha
theta2 = theta + alpha

M1, P1P0, T1T0, p1p0= M2_P2P1_T2T1_p2p1(M0, theta1, 'w')
M2, P2P0, T2T0, p2p0= M2_P2P1_T2T1_p2p1(M0, theta2, 'w')

#Solving for expansion fans:
v_M1 = PM_angle(M1)
v_M2 = PM_angle(M2)

t_list = np.linspace(0, 21*np.pi/180, 500)

P3P0_list = np.zeros(500)
theta3_list = theta1 + t_list
theta3 = 2*theta + theta1
v3 = 2*theta + v_M1
M3 = M_expansion(v3)
P3P1 = P2P1_expansion(M1, M3)
P3P0 = P3P1*P1P0

P4P0_list = np.zeros(500)
theta4_list = theta2 - t_list
theta4 = theta2 - 2*theta
v4 = 2*theta + v_M2
M4 = M_expansion(v4)
P4P2 = P2P1_expansion(M2, M4)
P4P0 = P4P2*P2P0

for i in range(len(t_list)):
    v = t_list[i] + v_M1
    M3 = M_expansion(v)
    P3P1 = P2P1_expansion(M1, M3)
    P3P0_list[i] = P3P1*P1P0

for i in range(len(t_list)):
    v = t_list[i] + v_M2
    M4 = M_expansion(v)
    P4P2 = P2P1_expansion(M2, M4)
    P4P0_list[i] = P4P2*P2P0


if(P3P0<= P4P0):
    #Initial guess.
    
    t = 2*theta
    t6 = 0

    #Values corresponding to initial guess.
    M5, P5P3 = M2_P2P1(M3, t, 'w')
    P5P0 = P5P3*P3P0

    M6 = M4
    P6P0 = P4P0

    if(P6P0<P5P0):
        
        theta5_list, P5P3_list = ShockPolar(M3)
        P5P0_list = P5P3_list*P3P0
        theta5_list = -theta5_list + (alpha + theta)*180/np.pi
        
        ts = theta_sonic(M3)
        P2P1s = M2_P2P1(M3, ts, 'w')[1]*P3P0
        tm = theta_max(M3)
        P2P1m = M2_P2P1(M3, tm, 'w')[1]*P3P0
        plt.plot(180*np.array([ts, tm])/np.pi  + (alpha + theta)*180/np.pi, [P2P1s, P2P1m], '.')

        theta6_list, P6P4_list = ShockPolar(M4)
        P6P0_list = P6P4_list*P4P0
        theta6_list = theta6_list + (alpha - theta)*180/np.pi

        ts = theta_sonic(M4)
        P2P1s = M2_P2P1(M4, ts, 'w')[1]*P4P0
        tm = theta_max(M4)
        P2P1m = M2_P2P1(M4, tm, 'w')[1]*P4P0
        plt.plot(180*np.array([ts, tm])/np.pi + (alpha - theta)*180/np.pi, [P2P1s, P2P1m], '.')
        

    else:
        
        theta5_list, P5P3_list = ShockPolar(M3)
        P5P0_list = P5P3_list*P3P0
        theta5_list = theta5_list + (alpha + theta)*180/np.pi

        ts = theta_sonic(M3)
        P2P1s = M2_P2P1(M3, ts, 'w')[1]*P3P0
        tm = theta_max(M3)
        P2P1m = M2_P2P1(M3, tm, 'w')[1]*P3P0
        plt.plot(180*np.array([ts, tm])/np.pi + (alpha + theta)*180/np.pi, [P2P1s, P2P1m], '.')
        
        P6P0_list = np.zeros(500)
        theta6_list = ((alpha - theta) - t_list)*180/np.pi
        v_M4 = PM_angle(M4)

        for i in range(len(t_list)):
            v = t_list[i] + v_M4
            M6 = M_expansion(v)
            P6P4 = P2P1_expansion(M4, M6)
            P6P0_list[i] = P6P4*P4P0
            
if(P3P0> P4P0):

    t = 0
    t6 = 2*theta

    M5 = M3
    P5P0 = P3P0

    M6, P6P4 = M2_P2P1(M4, t6, 'w')
    P6P0 = P6P4*P4P0

    if(P6P0>P5P0):
        
        theta5_list, P5P3_list = ShockPolar(M3)
        P5P0_list = P5P3_list*P3P0
        theta5_list = theta5_list + (alpha + theta)*180/np.pi

        ts = theta_sonic(M3)
        P2P1s = M2_P2P1(M3, ts, 'w')[1]*P3P0
        tm = theta_max(M3)
        P2P1m = M2_P2P1(M3, tm, 'w')[1]*P3P0
        plt.plot(180*np.array([ts, tm])/np.pi + (alpha + theta)*180/np.pi, [P2P1s, P2P1m], '.')

        theta6_list, P6P4_list = ShockPolar(M4)
        P6P0_list = P6P4_list*P4P0
        theta6_list = theta6_list + (alpha - theta)*180/np.pi

        ts = theta_sonic(M4)
        P2P1s = M2_P2P1(M4, ts, 'w')[1]*P4P0
        tm = theta_max(M4)
        P2P1m = M2_P2P1(M4, tm, 'w')[1]*P4P0
        plt.plot(180*np.array([ts, tm])/np.pi + (alpha - theta)*180/np.pi, [P2P1s, P2P1m], '.')

    else:
        
        theta6_list, P6P4_list = ShockPolar(M4)
        P6P0_list = P6P4_list*P4P0
        theta6_list = theta6_list + (alpha - theta)*180/np.pi ###

        ts = theta_sonic(M4)
        P2P1s = M2_P2P1(M4, ts, 'w')[1]*P4P0
        tm = theta_max(M4)
        P2P1m = M2_P2P1(M4, tm, 'w')[1]*P4P0
        plt.plot(180*np.array([ts, tm])/np.pi + (alpha - theta)*180/np.pi, [P2P1s, P2P1m], '.')
        
        P5P0_list = np.zeros(500)
        theta5_list =((alpha + theta) + t_list)*180/np.pi ###
        v_M3 = PM_angle(M3)

        for i in range(len(t_list)):
            v = t_list[i] + v_M3
            M5 = M_expansion(v)
            P5P3 = P2P1_expansion(M3, M5)
            P5P0_list[i] = P5P3*P3P0


    
plt.plot(theta_polar1, P1P0_list, label = "$M_0$")
plt.plot(theta3_list*180/np.pi, P3P0_list, label = "$M_1$")
plt.plot(theta4_list*180/np.pi, P4P0_list, label = "$M_2$")
plt.plot(theta5_list, P5P0_list, label = "$M_3$")
plt.plot(theta6_list, P6P0_list, label = "$M_4$")

plt.plot([0.6965],[0.98755], 'ob', label = "Solution β")
#plt.plot(M0_list, tN*180/np.pi, 'm', label = "$θ_N$")
plt.title("Diamond Airfoil: α = -10°")
plt.ylabel("P/P0 -->")
plt.xlabel("θ (Degrees) -->")
plt.legend(bbox_to_anchor=(1, 0.5))
plt.grid(linestyle = '--')
plt.savefig("Assign3", dpi = 300, bbox_inches = 'tight')
plt.show()













