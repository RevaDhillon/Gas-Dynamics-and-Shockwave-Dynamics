#Assignment 3: Double wedge problem

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


#To return M2 and P2/P1 after an oblique shock.
def M2_P2P1(M1, t):
    a = math.tan(t)
    b = M1**2
    c = 1.2*a*b + a
    s = sin_2_phi(a,b,c)
    M2 = math.sqrt(((1 + 0.2*M1*M1)/(1.4*M1*M1*s - 0.2)) + (M1*M1*(1 - s)/(1 + 0.2*M1*M1*s)))
    P2P1 = (7*(M1**2)*s - 1)/6
    return M2,P2P1


#Function to return alpha:
def alpha(tol, alpha0, tu, tl, M2, P2P1, M3, P3P1):
    alpha = alpha0
    n_max = 10000 #Maximum number of iterations.
    err = 1
    while(abs(err)>0.001 and n_max>0):

        #Setting the turn angles as per the guess.
        t4 = tu + alpha
        t5 =  tl - alpha

        #1: For region 4, corresponding to this guess:
        M4, P4P2 = M2_P2P1(M2, t4)
        P4P1 = P4P2 * P2P1

        #2: For region 5, corresponding to this guess:
        M5, P5P3 = M2_P2P1(M3, t5)
        P5P1 = P5P3 * P3P1

        #Checking the difference between the pressure ratios. Should be 0 ideally.
        err = P4P1 - P5P1

        #Making a better guess for alpha.
        if(err>0):
            alpha -= 0.001
            n_max -= 1

        else:
            alpha += 0.001
            n_max -= 1
        
    return M4, P4P1, M5, P5P1, alpha


#Initializing the values for the specific case.

#To radians.
theta_u = np.pi/18
theta_l = 15*np.pi/180
M1 = 3

#Values computed: Non - iterative.
M2, P2P1 = M2_P2P1(M1, theta_u)

M3, P3P1 = M2_P2P1(M1, theta_l)

M4, P4P1, M5, P5P1, alpha = alpha(0.01, 0,theta_u, theta_l, M2, P2P1, M3, P3P1)

print("Angle made by the slip plane 𝞪:" + str(alpha*180/math.pi) + "°")
print("Pressure Ratio P₅/P₁:" + str(P5P1))
print("Pressure Ratio P₄/P₁:" + str(P4P1))
print("M₄:" + str(M4))
print("M₅:" + str(M5))
