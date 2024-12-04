#Assignment 4: Flow over a compression ramp

import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math


#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '16'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)


#Function for Newton Method
def Newton_method(f, df, x0, tol):
    
    #Set maximum number of iterations.
    n_max = 100000000
    
    while(abs(f(x0))>tol):
        #Compute the next better approximation using the formula.
        x1 = x0 - (f(x0)/df(x0))

        if abs(f(x1))<tol:
            return x1
        else:
            x0 = x1
            n_max -= 1
            
    #If no root is found, return None.
    return None


#To return sin^2(phi):
def sin_2_phi(a,b,c):

    #Creating the function to find the root for.
    x = sym.Symbol('x')
    f = (a*a + 1)*b*b*(x**3) - (b*b +2*b + 2*a*b*c)*x*x + (c*c + 2*b + 1)*x -1
    df = f.diff(x)
    f = sym.lambdify(x, f)
    df = sym.lambdify(x, df)

    #To obtain sin^2(phi):
    sin2phi = Newton_method(f, df, 0.5, 0.001)

    return sin2phi


#To return M2, P2/P1 and sin^2(ɸ) after an oblique shock.
def M2_P2P1(M1, t):
    a = math.tan(t)
    b = M1**2
    c = 1.2*a*b + a
    s = sin_2_phi(a,b,c)
    M2 = math.sqrt(((1 + 0.2*M1*M1)/(1.4*M1*M1*s - 0.2)) + (M1*M1*(1 - s)/(1 + 0.2*M1*M1*s)))
    P2P1 = (7*(M1**2)*s - 1)/6
    return M2,P2P1, s


#To return PM angle given M:
def PM_angle(M):
    v = math.sqrt(6)*math.atan(math.sqrt((M*M -1)/6)) - math.atan(math.sqrt(M*M -1))
    return v


#To return M given the PM angle.
def M_expansion(v):
    
    #Creating the function to find the root for.
    x = sym.Symbol('x')
    f = sym.sqrt(6)*sym.atan(sym.sqrt((x*x -1)/6)) - sym.atan(sym.sqrt(x*x -1)) - v
    df = f.diff(x)
    f = sym.lambdify(x, f)
    df = sym.lambdify(x, df)

    M = Newton_method(f, df, 1.2, 0.001)

    return M


#To return P2P1 across an expansion fan.
def P2P1_expansion(M1, M2):
    P2P1 = ((1 + 0.2*M1*M1)/(1 + 0.2*M2*M2))**3.5
    return P2P1


#Known quantities.
M0 = 3
v0 = PM_angle(M0)


#Generating the compression ramp:
x = np.linspace(0,10,10)
l = len(x)
theta = np.pi/9
dtheta = (np.pi/9)/l
dx = x[1] - x[0]
y = np.empty(l)

y[0] = 0
for i in range(1,l):
    y[i] = np.tan(i*dtheta)*dx + y[i-1]

plt.plot(x, y, '-b')
#Extending ramp to abscissa = 25.
x_0 = 25
plt.plot([x[-1], x_0], [y[-1], y[-1] + math.tan(theta)*(x_0 - x[-1])], '-b')


#Characteristic waves.
x_1 = 21.7 #Approximate point of intersection abscissa.

v = np.empty(l+1)
M = np.empty(l+1)
u = np.empty(l+1)

v[0] = v0
M[0] = M0
u[0] = math.asin(1/M0)

for i in range(1,l+1):
    v[i] = v[0] - i*dtheta
    M[i] = M_expansion(v[i])
    u[i] = math.asin(1/M[i])

P2P1 = P2P1_expansion(M[0], M[-1])
    
#Plotting the characteristic waves.
for i in range(l):
    y_char = y[i] + math.tan(u[i])*(x_1 - x[i])
    plt.plot([x[i], x_1], [y[i], y_char])
    

#Plotting the oblique shock:
turn_angle = theta
M3, P3P1, sin_2_phi_val = M2_P2P1(M0, turn_angle)
phi = math.asin(math.sqrt(sin_2_phi_val))

x1 = x_1
x2 = x_0
y1 = y[0] + math.tan(u[0])*(x_1 - x[0])
y2 = y1 + math.tan(phi)*(x2 - x1)

plt.plot([x1, x2], [y1, y2], '-r',linewidth = 3, label = 'Oblique shock')

#In this case, the flow prefers to maintain flow tangency over pressure balance across the slip plane.
#Hence, the slip plane is parallel to the ramp.

#Plotting the slip plane:
plt.plot([x1, x_0], [y1, y1 + math.tan(theta)*(x_0 - x1)], '--b', label = 'Slip plane')    

plt.title("Flow over a Compression Ramp")
plt.xlabel("x -->")
plt.ylabel("y -->")
plt.grid(':')
plt.legend()
plt.tight_layout()
plt.savefig("Compression_Ramp", dpi = 300, bbox_inches = 'tight')
plt.show()

print("Turn angle made by the compression ramp θ:" + str(theta*180/np.pi) + "°")
print("Angle made by the slip plane 𝞪:" + str(theta*180/np.pi) + "°")
print("Pressure Ratio P₃/P₁:" + str(P3P1))
print("Pressure Ratio P₂/P₁:" + str(P2P1))
print("M₃:" + str(M3))
print("M₂:" + str(M[-1]))





