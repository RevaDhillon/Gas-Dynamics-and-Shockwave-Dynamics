#AS6060: Functions for shocks:

import sympy as sym
from scipy.optimize import fsolve
import numpy as np
import math


def Shock_angles(theta, M):
    
    A = np.tan(theta)*(1 + 1.2*M*M)
    B = M*M - 1
    C = A - np.tan(theta)*M*M

    leave, x1, x2 = np.roots([1, A, -B, C])

    if(x1 == 0):
        x1 = np.pi/2
    else:
        x1 = np.arctan(1/x1)
        
    if(x2 == 0):
        x2 = np.pi/2
    else:
        x2 = np.arctan(1/x2)

    return x1, x2


#To return M2, P2/P1, etc. after an oblique shock.
def M2_P2P1_T2T1_p2p1(M1, t, p):
    
    phi_w, phi_s = Shock_angles(t, M1)

    if(p == 'w'):
        s = (np.sin(phi_w))**2
    else:
        s = (np.sin(phi_s))**2
    
    M2 = np.sqrt(((1 + 0.2*M1*M1)/(1.4*M1*M1*s - 0.2)) + (M1*M1*(1 - s)/(1 + 0.2*M1*M1*s)))
    P2P1 = (7*(M1**2)*s - 1)/6
    p2p1 = 6*M1*M1*s/(M1*M1*s + 5)
    T2T1 = (5 + M1*M1*s)*(7*M1*M1*s - 1)/(36*M1*M1*s)
    
    return M2, P2P1, T2T1, p2p1


#To return M2, P2/P1 after an oblique shock.
def M2_P2P1(M1, t, p):
    
    phi_w, phi_s = Shock_angles(t, M1)

    if(p == 'w'):
        s = (np.sin(phi_w))**2
    else:
        s = (np.sin(phi_s))**2
    
    M2 = np.sqrt(((1 + 0.2*M1*M1)/(1.4*M1*M1*s - 0.2)) + (M1*M1*(1 - s)/(1 + 0.2*M1*M1*s)))

    P2P1 = (7*(M1**2)*s - 1)/6
    
    return M2, P2P1


#Find theta_max:
def theta_max(M):
    
    sin2_phi_max = (1/(1.4*M*M))*(0.6*M*M -1 + np.sqrt((2.4)*((2.4)*M*M*M*M/16 + 0.2*M*M + 1)))
    phi_max = np.arcsin(np.sqrt(sin2_phi_max))
    N = M*M*sin2_phi_max - 1
    tan_theta_max = N/((np.tan(phi_max))*(1.2*M*M - N))
    theta_max = np.arctan(tan_theta_max)

    return theta_max


#Find theta_sonic:
def theta_sonic(M):
    
    sin2_phi_sonic = (1/(1.4*M*M))*(0.6*M*M -0.4 + np.sqrt((2.4)*((2.4)*M*M*M*M/16 - 0.2*M*M + 0.65)))
    phi_sonic = np.arcsin(np.sqrt(sin2_phi_sonic))
    N = M*M*sin2_phi_sonic - 1
    tan_theta_sonic = N/((np.tan(phi_sonic))*(1.2*M*M - N))
    theta_sonic = np.arctan(tan_theta_sonic)

    return theta_sonic


def theta_phi_M(M, phi):
    N = M*M*np.sin(phi)*np.sin(phi) - 1
    tan_theta = N/(np.tan(phi)*(1.2*M*M - N))
    theta = np.arctan(tan_theta)
    return theta


def ShockPolar(M):
    
    phi_list = np.linspace(np.arcsin(1/M), np.pi/2)
    P2P1 = (7*M*M*(np.sin(phi_list)**2) - 1)/6

    y = P2P1
    a = 1+1.4*M*M
    b = (7*M*M - 1)/6
    c = 1/6

    theta = np.arctan(((y-1)/(a - y))*np.sqrt((b - y)/(c + y)))
    theta2 = -theta
    theta = np.concatenate([np.flip(theta2), theta])
    y = np.concatenate([np.flip(y), y])
    theta = theta*180/np.pi

    return theta, y


def Intersections(theta_u, theta_l, PlP0, PuP0, Ml, Mu):
    
    au = 1+1.4*Mu*Mu
    bu = (7*Mu*Mu - 1)/6
    cu = 1/6

    al = 1+1.4*Ml*Ml
    bl = (7*Ml*Ml - 1)/6
    cl = 1/6

    lhs = theta_u - theta_l    

    y = sym.symbols('y')
    f = sym.Eq(sym.atan(((y/(PuP0)-1)/(au - y/(PuP0)))*sym.sqrt((bu - y/(PuP0))/(cu + y/(PuP0)))) + sym.atan(((y/(PlP0)-1)/(al - y/(PlP0)))*sym.sqrt((bl - y/(PlP0))/(cl + y/(PlP0)))), lhs)
    res = sym.solve(f, y)

    return res


def MR(M0, P0, T0, p0, M1, P1, T1, p1, t):

    def funct(variables):
        (phi2, phi3, t2, t3, P2, P3) = variables
    
        eq1 = np.tan(t2)*np.tan(phi2) - ((M1*M1*np.sin(phi2)*np.sin(phi2) - 1)/(1.2*M1*M1 - M1*M1*np.sin(phi2)*np.sin(phi2) + 1))
        eq2 = np.tan(t3)*np.tan(phi3) - ((M0*M0*np.sin(phi3)*np.sin(phi3) - 1)/(1.2*M0*M0 - M0*M0*np.sin(phi3)*np.sin(phi3) + 1))
        eq3 = abs(t3) + abs(t2) - t
        eq4 = P2/P1 - (7*M1*M1*np.sin(phi2)*np.sin(phi2) - 1)/6
        eq5 = P3/P0 - (7*M0*M0*np.sin(phi3)*np.sin(phi3) - 1)/6
        eq6 = P2 - P3
    
        return [eq1, eq2, eq3, eq4, eq5, eq6]
        
    phi2, phi3, t2, t3, P2, P3 = fsolve(funct, (np.arcsin(1/M1), 1.57, 0.0015, t-0.0015, P1*0.3, P1*0.3))

    
    M2, P2P1, T2T1, p2p1 = M2_P2P1_T2T1_p2p1(M1, t2, 'w')
    P2 = P2P1*P1
    T2 = T2T1*T1
    p2 = p2p1*p1

    
    M3, P3P0, T3T0, p3p0 = M2_P2P1_T2T1_p2p1(M0, t3, 's')
    P3 = P3P0*P0
    T3 = T3T0*T0
    p3 = p3p0*p0
    
    return M2, t2, P2, T2, p2, M3, t3, P3, T3, p3
    

#To return PM angle given M:
def PM_angle(M):
    v = np.sqrt(6)*np.arctan(math.sqrt((M*M -1)/6)) - math.atan(math.sqrt(M*M -1))
    return v


#To return M given the PM angle.
def M_expansion(v):

    def funct(x):
        return sym.sqrt(6)*sym.atan(sym.sqrt((x[0]*x[0] -1)/6)) - sym.atan(sym.sqrt(x[0]*x[0] -1)) - v

    M = float(fsolve(funct, [2]))

    return float(M)


#To return P2P1 across an expansion fan.
def P2P1_expansion(M1, M2):
    P2P1 = ((1 + 0.2*M1*M1)/(1 + 0.2*M2*M2))**3.5
    return P2P1







    
    

