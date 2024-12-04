
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



def Shock_angles(theta, M):
    
    A = np.tan(theta)*(1 + 4*M*M/3)
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


#To return M2, P2/P1 after an oblique shock.
def M2_P2P1(M1, t, p):
    
    phi_w, phi_s = Shock_angles(t, M1)

    if(p == 'w'):
        s = (np.sin(phi_w))**2
    else:
        s = (np.sin(phi_s))**2
    
    M2 = np.sqrt(((3 + M1*M1)/(5*M1*M1*s - 1)) + (3*M1*M1*(1 - s)/(3 + M1*M1*s)))

    P2P1 = (5*(M1**2)*s - 1)/4
    
    return M2, P2P1


#Find theta_max:
def theta_max(M):
    
    sin2_phi_max = (3/(5*M*M))*(2*M*M/3 -1 + np.sqrt((8/3)*((8/3)*M*M*M*M/16 + (1/3)*M*M + 1)))
    phi_max = np.arcsin(np.sqrt(sin2_phi_max))
    N = M*M*sin2_phi_max - 1
    tan_theta_max = N/((np.tan(phi_max))*((4/3)*M*M - N))
    theta_max = np.arctan(tan_theta_max)

    return theta_max


#Find theta_sonic:
def theta_sonic(M):
    
    sin2_phi_sonic = (1/(5*M*M))*(2*M*M - 1 + 2*np.sqrt(M*M*M*M - M*M + 4))
    phi_sonic = np.arcsin(np.sqrt(sin2_phi_sonic))
    
    N = M*M*sin2_phi_sonic - 1
    tan_theta_sonic = N/((np.tan(phi_sonic))*((4/3)*M*M - N))
    theta_sonic = np.arctan(tan_theta_sonic)

    return theta_sonic

def ShockPolar(M):
    
    phi_list = np.linspace(np.arcsin(1/M), np.pi/2)
    P2P1 = (5*M*M*(np.sin(phi_list)**2) - 1)/4

    y = P2P1
    a = 1+(5/3)*M*M
    b = (5*M*M - 1)/4
    c = 0.25

    theta = np.arctan(((y-1)/(a - y))*np.sqrt((b - y)/(c + y)))
    theta2 = -theta
    theta = np.concatenate([np.flip(theta2), theta])
    y = np.concatenate([np.flip(y), y])
    theta = theta*180/np.pi

    return theta, y

#To return PM angle given M:
def PM_angle(M):
    v = 2*np.arctan(math.sqrt((M*M -1)/4)) - math.atan(math.sqrt(M*M -1))
    return v


#To return M given the PM angle.
def M_expansion(v):

    def funct(x):
        return 2*sym.atan(sym.sqrt((x[0]*x[0] -1)/4)) - sym.atan(sym.sqrt(x[0]*x[0] -1)) - v

    M = float(fsolve(funct, [1.2]))

    return float(M)


#To return P2P1 across an expansion fan.
def P2P1_expansion(M1, M2):
    P2P1 = ((1 + M1*M1/3)/(1 + M2*M2/3))**2.5
    return P2P1


#M0:

M0 = 1.8

theta0, P1P0_list = ShockPolar(M0)
plt.plot(theta0, P1P0_list, label = "$M_0$ = "+str(M0))

ts = abs(theta_sonic(M0))
P2P1s = M2_P2P1(M0, ts, 'w')[1]
tm = abs(theta_max(M0))
P2P1m = M2_P2P1(M0, tm, 'w')[1]


t = 15.9*np.pi/180

M1, P1P0 = M2_P2P1(M0, t, 'w')

theta1, P2P1_list = ShockPolar(M1)
theta1_list = t*180/np.pi + theta1
P2P0 = P1P0*P2P1_list 

plt.plot(theta1_list, P2P0, label = "$M_{1_{R}}$ = "+str(int(M1*100)/100))
plt.plot([0,0], [1,4])
plt.plot(180*np.array([ts, tm])/np.pi, [P2P1s, P2P1m], '.')

#M1:

ts = abs(theta_sonic(M1))
P2P1s = M2_P2P1(M1, ts, 'w')[1]
tm = abs(theta_max(M1))
P2P1m = M2_P2P1(M1, tm, 'w')[1]

plt.plot((180*np.array([ts, tm])/np.pi + t*180/np.pi), P1P0*np.array([P2P1s, P2P1m]), '.')


P2P0_list = np.empty(50)

t_list= np.linspace(0, 0.001)

for i in range(len(t_list)):
    v = t_list[i] 
    M2 = M_expansion(v)
    P2P1 = P2P1_expansion(1, M2)
    P2P0_list[i] = P2P1*P2P1s*P1P0

#plt.plot((t_list + ts + t)*180/np.pi, P2P0_list, '-')


t = -15.8*np.pi/180

M1, P1P0 = M2_P2P1(M0, t, 'w')

theta1, P2P1_list = ShockPolar(M1)
theta1_list = t*180/np.pi + theta1
P2P0 = P1P0*P2P1_list 

plt.plot(theta1_list, P2P0, label = "$M_{1_{L}}$ = "+str(int(M1*100)/100))
plt.plot([0,0], [1,4])

#M1:

ts = abs(theta_sonic(M1))
P2P1s = M2_P2P1(M1, ts, 'w')[1]
tm = abs(theta_max(M1))
P2P1m = M2_P2P1(M1, tm, 'w')[1]

plt.plot(-(180*np.array([ts, tm])/np.pi) + t*180/np.pi, P1P0*np.array([P2P1s, P2P1m]), '.')


##P2P0_list = np.empty(50)
##
##t_list= np.linspace(0, 0.001)
##
##for i in range(len(t_list)):
##    v = t_list[i] 
##    M2 = M_expansion(v)
##    P2P1 = P2P1_expansion(1, M2)
##    P2P0_list[i] = P2P1*P2P1s*P1P0
##
##plt.plot((t_list + ts + t)*180/np.pi, P2P0_list, '-')
##


##plt.figure(1)
plt.title("Shock Polars: M = " + str(M0))
plt.ylabel("$ \\frac{P}{P_0} $ -->")
plt.xlabel("θ (Degrees) -->")
plt.grid(linestyle = "--")
plt.legend()#bbox_to_anchor=(1, 0.5))
#plt.savefig("Shock_Polars_P2P1", dpi = 300, bbox_inches = 'tight')
plt.show()











