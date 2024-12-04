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

 #float(input("Enter design Mach number:"))


def Hm(PbPe, NPR, M0):

    M0 = 3
    P0 = 287
    p0 = 1
    T0 = 1

    Pb = PbPe*P0

    print("Pe:" + str(P0))
    print("M0:" + str(M0))
    
    print("Pb:" + str(Pb))
    
    term3 = NPR**(2/7) - 1
    M1 = math.sqrt(5*term3)
    P1 = Pb
    T1 = ((1 + 0.2*M0*M0)/(1 + 0.2*M1*M1))*T0


    mu_1 = np.arcsin(1/M1)
    phi1 = np.arcsin(1/M0)
    theta = PM_angle(M1) - PM_angle(M0)

    H = np.tan(theta + mu_1)*(np.tan(theta) + np.tan(phi1))/(np.tan(phi1)*(np.tan(theta + mu_1) - np.tan(theta)))

    v0 = PM_angle(M1) + theta
    M0 = M_expansion(v0)
    P0 = P2P1_expansion(M1, M0)*Pb
    T0 = ((1 + 0.2*M1*M1)/(1 + 0.2*M0*M0))*T1
    p0 = P0/(287*T0)

    print("P_0:" + str(P0))
    print("M_0:" + str(M0))
    

    PbPe = Pb/P0

    #1
    phi1 = np.arcsin(np.sqrt(((6*PbPe + 1)/(7*M0*M0))))
    s =  (np.sin(phi1))**2
    M1 = np.sqrt(((1 + 0.2*M0*M0)/(1.4*M0*M0*s - 0.2)) + (M0*M0*(1 - s)/(1 + 0.2*M0*M0*s)))
    P1P0 = PbPe
    print(M0, M1)
    print(P1P0)
    p1p0 = 6*M0*M0*s/(M0*M0*s + 5)
    T1T0 = (5 + M0*M0*s)*(7*M0*M0*s - 1)/(36*M0*M0*s)
    P1 = P1P0*P0
    T1 = T1T0*T0
    p1 = p1p0*p0
    theta1 = theta_phi_M(M0, phi1)


    #2, 3
    M2, theta2, P2, T2, p2, M3, theta3, P3, T3, p3 = MR(M0, P0, T0, p0, M1, P1, T1, p1, theta1)
    phi2, a = Shock_angles(theta2, M1)
    a, phi3 = Shock_angles(theta3, M0)

    #G
    Mg, PgP0, TgT0, pgp0 = M2_P2P1_T2T1_p2p1(M0, 0, 's')
    Pg = PgP0*P0
    Tg = TgT0*T0
    pg = pgp0*p0
    phi_g = np.pi*0.5

    #Avg.
    p_av = 0.5*(p3 + pg)
    a3 = math.sqrt(1.4*287*T3)
    ag = math.sqrt(1.4*287*Tg)
    a_av = 0.5*(a3+ag)
    u_av = (p3*M3*a3*math.cos(theta3) + pg*Mg*ag)/(2*p_av)
    M_av = u_av/a_av


    #HmHs:
    HmHs = (((5+M_av*M_av)/6)**3)/M_av

    #C, C', D:
    vM1 = PM_angle(M1)
    vM2 = PM_angle(M2)

    def functs(variables):
        (vMd, Md) = variables
        
        eq1 = -vMd+ (np.sqrt(6)*np.arctan(np.sqrt((Md*Md-1)/6)) - np.arctan(np.sqrt(Md*Md -1)))
        eq2 = -vMd + vM2 + theta3
        
        return [eq1, eq2]

    vMd, Md = fsolve(functs, (vM1, M1))

    #Geometry:
    mu_D = math.asin(1/Md)
    mu_2 = math.asin(1/M2)

    def funct(variables):
        (xT, Hm, xD, yD, xF, yF, xE, Hs) = variables
        
        eq1 = xT*math.tan(phi1) - (H-Hm)
        eq2 = xD*math.tan(theta1) - (H-yD)
        eq3 = (xD - xT)*math.tan(phi2 - theta1) - (yD - Hm)
        eq4 = (xF - xT)*math.tan(theta3) - (Hm - yF)
        eq5 = (xE - xD)*math.tan(mu_D) - (yD - Hs)
        eq6 = (xF - xD)*math.tan(mu_2 + theta3) - (yF - yD)
        eq7 = (xE - xF)*math.tan(theta3) - (2 + math.tan(theta3)*math.tan(theta3))*(yF - Hs)
        eq8 = Hm/Hs - HmHs
        
        return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8]


    xT, Hm, xD, yD, xF, yF, xE, Hs = fsolve(funct, (8.465, 1.5653, 11.3213, 1.70466, 10.3975, 1.2666, 14.2496, 0.9725))
    print(P2, P3)

    return Hm, theta3


####

M0 = 3
NPR = 195
P0Pe = (1 + 0.2*M0*M0)**3.5
PbPe = P0Pe/NPR
Hm(PbPe, NPR, M0)
##Ar = np.linspace(3.5, 20)
##
##    
##P0Pe = (1 + 0.2*M0*M0)**3.5
##NPR = P0Pe*np.linspace(5.2, 9, 10)
####print(Hm(P0Pe/6.20))
##
##PbPe = P0Pe/NPR
##
##Hm_list = np.empty(len(NPR))
##theta3list = np.empty(len(NPR))
##
##for i in range(len(NPR)):
##    Hm_list[i], theta3list[i] = Hm(PbPe[i], NPR[i], 10)
##
##plt.plot(theta3list, Hm_list, 'o-', label = "A/A* = 4.23")
##
######
##
##M0 = 2.75
##NPR = P0Pe*np.linspace(5, 9.2, 10)
##Ar = np.linspace(3.5, 20)
##
##    
##P0Pe = (1 + 0.2*M0*M0)**3.5
####
####print(Hm(P0Pe/6.20))
##
##PbPe = P0Pe/NPR
##print(PbPe)
##
##Hm_list = np.empty(len(NPR))
##theta3list = np.empty(len(NPR))
##
##for i in range(len(NPR)):
##    Hm_list[i], theta3list[i] = Hm(PbPe[i], NPR[i], M0)
##
##plt.plot(theta3list, Hm_list, 's-', label = "A/A* = 3.34")
##
#######
##
##M0 = 4
##NPR = P0Pe*np.linspace(6, 14, 10)
##Ar = np.linspace(3.5, 20)
##
##    
##P0Pe = (1 + 0.2*M0*M0)**3.5
####
####print(Hm(P0Pe/6.20))
##
##PbPe = P0Pe/NPR
##print(PbPe)
##
##Hm_list = np.empty(len(NPR))
##theta3list = np.empty(len(NPR))
##
##for i in range(len(NPR)):
##    Hm_list[i], theta3list[i] = Hm(PbPe[i], NPR[i], M0)
##
##plt.plot(theta3list, Hm_list, '*-', label = "A/A* = 10.719")
##
##
##plt.title("Underexpanded regime - Mach stem height")
##plt.ylabel("$ \\frac{H_m}{H} $ -->")
##plt.xlabel("θ3 -->")
##plt.grid(linestyle = "--")
##plt.legend()
###plt.savefig("Shock_Polars_M2", dpi = 300, bbox_inches = 'tight')
##plt.show()
##
##
##
##
##


