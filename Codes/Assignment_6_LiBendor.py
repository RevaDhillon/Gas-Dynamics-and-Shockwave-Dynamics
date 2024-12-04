
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


#thetaw =float(input("Enter wedge angle in degrees:"))#20.726

def Hmreturn(thetaw):
    thetaw = thetaw*np.pi/180

    H = 1

    M0 =2.84#float(input("Enter inlet Mach number:")) #2.84
    P0 = 287
    p0 = 1
    T0 = 1 

    w = 1.42 #float(input("Enter length of the hypotenuse of the wedge w/H:")) # 1.42
    g = 1 - w*np.sin(thetaw)#float(input("Enter height from line of symmetry to bottom corner of the wedge:"))
    g_nd = g/w

    #1
    theta1 = thetaw
    M1, P1P0, T1T0, p1p0 = M2_P2P1_T2T1_p2p1(M0, theta1, 'w')
    P1 = P1P0*P0
    T1 = T1T0*T0
    p1 = p1p0*p0
    phi1, a = Shock_angles(theta1, M0)


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
        (Mc, alpha, vMc, Pc, theta_c, P_c, M_c, Pd, Md, vMd) = variables
        
        eq1 = -vMc+ (np.sqrt(6)*np.arctan(np.sqrt((Mc*Mc-1)/6)) - np.arctan(np.sqrt(Mc*Mc -1)))
        eq2 = -vMc + vM1 + thetaw - alpha
        eq3 = (Pc/P1) - ((5 + M1*M1)/(5 + Mc*Mc))**3.5
        eq4 = alpha - theta_c
        eq5 = (P_c/Pc) - M2_P2P1(Mc, theta_c, 'w')[1]
        eq6 = M_c - M2_P2P1(Mc, theta_c, 'w')[0]
        eq7 = P_c - Pd
        eq8 = (Pd/P2) - ((5 + M2*M2)/(5 + Md*Md))**3.5
        eq9 = -vMd + (np.sqrt(6)*np.arctan(np.sqrt((Md*Md-1)/6)) - np.arctan(np.sqrt(Md*Md -1)))
        eq10 = vMd - vM2 -theta3
        
        return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10]

    Mc, alpha, vMc, Pc, theta_c, P_c, M_c, Pd, Md, vMd = fsolve(functs, (M1, 0.01, vM1, P1, theta2, P2, M2, P2, M2, vM2))

    phiC, a = Shock_angles(theta_c, Mc)

    #Geometry:

    xR = w*math.cos(thetaw)
    yR = H - w*math.sin(thetaw)
    mu_B = math.asin(1/M1)
    mu_C = math.asin(1/Mc)
    mu_F = math.asin(1/M2)
    mu_D = math.asin(1/Md)
    del_bbc = phi2 - thetaw 
    del_cbc = phiC - alpha
    del_ccd = -math.asin(1/M_c)
    del_dcd = - mu_D
    del_bbd = - theta3
    del_dbd = 0
    del_ffe = - theta3
    del_efe = 0


    def tan_funct(delta1, delta2):
        
        term1 = 2*math.tan(delta1) + math.tan(delta2 - delta1)
        term2 = 2 - math.tan(delta1)*math.tan(delta2 - delta1)
        term3 = math.atan(term1/term2)

        return term3
        

    def funct(variables):
        (yB, xB, yC, xC, yF, xF, yE, Hs, yD, xD, xE, yT, xT, Hm) = variables
        
        eq1 = yB - yR+(xB - xR)*math.tan(mu_B + thetaw)
        eq2 = yC - yR+(xC - xR)*math.tan(mu_C + alpha)
        eq3 = yF - yB +(xF - xB)*math.tan(mu_F + theta3)
        eq4 = yE - Hs  
        eq5 = yE - yD +(xE - xD)*math.tan(mu_D)
        eq6 = yB - yT - (xB - xT)*math.tan(phi2 - thetaw)
        eq7 = Hm - yT  
        eq8 = xT - (H - Hm)/math.tan(phi1)
        eq9 = yF - yT +(xF - xT)*math.tan(theta3)

        eq10 = yB - yC -(xB-xC)*tan_funct(del_bbc, del_cbc)
        eq11 = yC - yD - (xC-xD)*tan_funct(del_ccd, del_dcd)
        eq12 = yB-yD -(xB-xD)*tan_funct(del_bbd, del_dbd)
        eq13 = yF-yE-(xF-xE)*tan_funct(del_ffe, del_efe)
        eq14 = Hm/Hs - HmHs
        
        return [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14]


    yB, xB, yC, xC, yF, xF, yE, Hs, yD, xD, xE, yT, xT, Hm = fsolve(funct, (0.4535,1.3855,0.45416,1.3867, 0.1315,1.5448, 0.127, 0.127,0.4534,1.3876, 1.6692, 0.1737, 0.9607, 0.1737))

    return Hm

##theta_list = np.linspace(19.15, 20.25, 51)
##Hm_l = np.empty(51)
##for i in range(51):
theta = input("Enter:")
Hm = Hmreturn(theta)
print("The Mach stem height Hm/H = {0:.2f}".format(Hm))



##plt.plot(theta_list, Hm_l, 'o-')
##plt.title("Mach stem height variation with wedge angle:")
##plt.ylabel("$ Hm/H $ -->")
##plt.xlabel("θ (Degrees) -->")
##plt.grid(linestyle = "--")
##plt.legend(bbox_to_anchor=(1, 0.5))
####plt.savefig("Shock_Polars_M2", dpi = 300, bbox_inches = 'tight')
##plt.show()











