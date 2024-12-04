#Assignment 5: Flow over a diamond airfoil

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


#Function to get the beta values.
def Beta_array(alpha_array, theta, M_inf, tol):

    l = len(alpha_array)
    beta = np.empty(l)
        
    #Obtaining beta for each alpha.
    for i in range(l):
        alpha = alpha_array[i]

        #Obtaining quantities in regions 1, 2, 3, 4:
        if(alpha >= theta):

            #Region 0 - 1:
            v_inf = PM_angle(M_inf)
            v_1 = alpha - theta + v_inf
            M1 = M_expansion(v_1)
            P1P_inf = P2P1_expansion(M_inf, M1)

            #Region 1 - 3:
            v_3 = 2*theta + v_1
            M3 = M_expansion(v_3)
            P3P1 = P2P1_expansion(M1, M3)
            P3P_inf = P3P1 * P1P_inf

            #Region 0 - 2:
            t = alpha + theta
            M2, P2P_inf = M2_P2P1(M_inf, t)

            #Region 2 - 4:
            v_2 = PM_angle(M2)
            v_4 = 2*theta + v_2
            M4 = M_expansion(v_4)
            P4P2 = P2P1_expansion(M2, M4)
            P4P_inf = P4P2 * P2P_inf
            

        else:

            #Region 0 - 1:
            t = theta - alpha 
            M1, P1P_inf = M2_P2P1(M_inf, t)

            #Region 1 - 3:
            v_1 = PM_angle(M1)
            v_3 = 2*theta + v_1
            M3 = M_expansion(v_3)
            P3P1 = P2P1_expansion(M1, M3)
            P3P_inf = P3P1 * P1P_inf

            #Region 0 - 2:
            if(alpha > -theta):
                t = alpha + theta
                M2, P2P_inf = M2_P2P1(M_inf, t)
            else:
                v_inf = PM_angle(M_inf)
                v_2 = v_inf - (alpha + theta)
                M2 = M_expansion(v_2)
                P2P_inf = P2P1_expansion(M_inf, M2)

            #Region 2 - 4:
            v_2 = PM_angle(M2)
            v_4 = 2*theta + v_2
            M4 = M_expansion(v_4)
            P4P2 = P2P1_expansion(M2, M4)
            P4P_inf = P4P2 * P2P_inf


        #Use the quantities to get pressures in regions 5 and 6.
        #Compare these pressures to obtain beta.
            
        if(P3P_inf < P4P_inf):
            #Initial guess.
            beta_guess = theta - alpha
            t = 2*theta
            t6 = 0

            #Values corresponding to initial guess.
            M5, P5P3 = M2_P2P1(M3, t)
            P5P_inf = P5P3*P3P_inf

            M6 = M4
            P6P_inf = P4P_inf

            #Loops to make a better guess and obtain beta.
            
            if(P5P_inf - P6P_inf <0):
                
                while(abs(P5P_inf - P6P_inf)>tol):
                    
                    beta_guess = beta_guess + 0.001
                    t = t + 0.001
                    v_6 = v_4 + 0.001

                    M5, P5P3 = M2_P2P1(M3, t)
                    P5P_inf = P5P3*P3P_inf

                    M6 = M_expansion(v_6)
                    P6P4 = P2P1_expansion(M4, M6)
                    P6P_inf = P6P4*P4P_inf

                beta_final = beta_guess
                
            if(P5P_inf - P6P_inf >0):
                
                while(abs(P5P_inf - P6P_inf)>tol):
                    
                    beta_guess = beta_guess - 0.001
                    t = t - 0.001
                    t6 = t6 + 0.001

                    M5, P5P3 = M2_P2P1(M3, t)
                    P5P_inf = P5P3*P3P_inf

                    M6, P6P4 = M2_P2P1(M4, t6)
                    P6P_inf = P6P4*P4P_inf

                beta_final = beta_guess

            
        if(P3P_inf - P4P_inf  == 0):
            
                beta_final = alpha
                    
               
        if(P3P_inf > P4P_inf):
            
            beta_guess = -(theta + alpha)
            t = 0
            t6 = 2*theta

            M5 = M3
            P5P_inf = P3P_inf

            M6, P6P4 = M2_P2P1(M4, t6)
            P6P_inf = P6P4*P4P_inf

            if(P5P_inf - P6P_inf <0):
                while(abs(P5P_inf - P6P_inf)>tol):
                    beta_guess = beta_guess + 0.001
                    t = t + 0.001
                    t6 = t6 - 0.001

                    M5, P5P3 = M2_P2P1(M3, t)
                    P5P_inf = P5P3*P3P_inf

                    M6, P6P4 = M2_P2P1(M4, t6)
                    P6P_inf = P6P4*P4P_inf

                beta_final = beta_guess
                
            if(P5P_inf - P6P_inf >0):
                while(abs(P5P_inf - P6P_inf)>tol):
                    beta_guess = beta_guess - 0.001
                    t6 = t6 + 0.001
                    v_5 = v_3 + 0.001

                    M5 = M_expansion(v_5)
                    P5P3 = P2P1_expansion(M3, M5)
                    P5P_inf = P5P3*P3P_inf
                    
                    M6, P6P4 = M2_P2P1(M4, t6)
                    P6P_inf = P6P4*P4P_inf

                beta_final = beta_guess

        #Store the values of beta obtained.    
        beta[i] = beta_final 

    beta = beta*180/np.pi
    return beta

#Cases:

#Setting quantities.
alpha_array = np.linspace(-10*np.pi/180, 10*np.pi/180, 5)
theta = np.pi/18
tol = 0.01

beta3 = Beta_array(alpha_array, theta, 3, tol) #For M = 3
#beta4 = Beta_array(alpha_array, theta, 4, tol) #For M = 4

#Converting to degrees.
alpha_array = alpha_array*180/np.pi

print(alpha_array)
print(beta3)

#Plotting and saving the graph obtained.
plt.plot(alpha_array, beta3, label = "M = 3")
#plt.plot(alpha_array, beta4, label = "M = 4")
plt.title("Diamond airfoil: α vs. β curve")
plt.xlabel("α (Degrees) -->")
plt.ylabel("β (Degrees) -->")
plt.grid('--')
plt.legend()
plt.tight_layout()
plt.savefig("DiamondAirfoil_M34", dpi = 300, bbox_inches = 'tight')
plt.show()







