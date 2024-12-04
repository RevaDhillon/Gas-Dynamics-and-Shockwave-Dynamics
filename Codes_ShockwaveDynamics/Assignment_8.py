import sympy as sym
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from Shock_Functions import *

#Selecting the font size and type to be used in the figure.
font={ 'family' : 'Times New Roman',
       'weight' : 'normal',
       'size' : '14'}

#Setting the selected font properties.          
matplotlib.rc('font', **font)

    
M0 = np.array([1.015, 1.124, 1.234, 1.343, 1.453, 1.562, 1.672, 1.781, 1.89, 2, 2, 2.02631579, 2.05263158, 2.07894737, 2.10526316, 2.13157895, 2.15789474, 2.18421053, 2.21052632, 2.23684211, 2.26315789, 2.28947368, 2.31578947, 2.34210526, 2.36842105, 2.39473684, 2.42105263, 2.44736842, 2.47368421, 2.5, 2.5, 2.75862069,  3.01724138,  3.27586207,  3.53448276,  3.79310345, 4.05172414,  4.31034483,  4.56896552,  4.82758621,  5.0862069,   5.34482759,   5.60344828,  5.86206897,  6.12068966,  6.37931034,  6.63793103,  6.89655172, 7.15517241,  7.4137931,   7.67241379,  7.93103448,  8.18965517,  8.44827586,  8.70689655,  8.96551724,  9.22413793,  9.48275862,  9.74137931, 10])
tN = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.27057217, 0.27498609, 0.27923453, 0.28326741, 0.28713017, 0.29082132, 0.29433952, 0.29773404, 0.30090332, 0.30399933, 0.30692013, 0.30692013, 0.32960996, 0.34410763, 0.35333738, 0.35910567, 0.36259535,  0.36441125, 0.36520985, 0.36519739, 0.36470276, 0.36392643, 0.36280598, 0.3616048,  0.3602711,  0.35887881, 0.35748763, 0.35607111, 0.35466815, 0.3532352,  0.35187413, 0.35060738, 0.34930104, 0.34812306, 0.34693358, 0.34574329, 0.34463844, 0.34362738, 0.34256267, 0.34160477, 0.34068164])
tD = np.array([0.00076641, 0.01700157, 0.04074034, 0.06751891, 0.09540399, 0.12333589, 0.15059823, 0.17679197, 0.20167243, 0.22512357, 0.22512357, 0.23052778, 0.23586004, 0.24107905, 0.24622688, 0.25126099, 0.25622511, 0.26112064, 0.26590315, 0.27061863, 0.27526856, 0.27980674, 0.28423288, 0.28859562, 0.2928965,  0.29713705, 0.30126828, 0.30529011, 0.30925402, 0.31316152, 0.31316152, 0.34773401, 0.37621479, 0.39976674, 0.41933358, 0.43568586, 0.44943148, 0.46109908, 0.47108635, 0.47959838, 0.48701178, 0.49338705, 0.49897503, 0.50393658, 0.50826926, 0.51217983, 0.51556756, 0.51869746, 0.52145325, 0.52393861, 0.52617766, 0.52819084, 0.53007211, 0.53176013, 0.53334459, 0.53475987, 0.53601548, 0.53727437, 0.53831212, 0.53936736] )

#Properly expanded:
P0Pe_PerExp = (1 + 0.2*M0*M0)**3.5

#Overexpanded:
PbPe_overMR = np.empty(35)
PbPe_overRR = np.empty(35)

#Underexpanded:
P0Pb_underMR = np.empty(35)
P0Pb_underRR = np.empty(35)

for i in range(35):
    
    phi_w, phi_s = Shock_angles(tN[i], M0[i])
    s = (np.sin(phi_w))**2
    PbPe_overMR[i] = (7*(M0[i]**2)*s - 1)/6

    phi_w, phi_s = Shock_angles(tD[i], M0[i])
    s = (np.sin(phi_w))**2
    PbPe_overRR[i] = (7*(M0[i]**2)*s - 1)/6

P0Pb_overMR = P0Pe_PerExp[:35]/PbPe_overMR
P0Pb_overRR = P0Pe_PerExp[:35]/PbPe_overRR

for i in range(35):

    P0Pe = P0Pe_PerExp[i]
    PbP_0_MR = PbPe_overMR[i]
    PbP_0_RR = PbPe_overRR[i]
    Me = M0[i]
    
    def functMR(variables):
        (P0, Pe, M_0, M1, Pb, P_0) = variables
            
        eq1 = P0/Pe - P0Pe
        eq2 = P0Pe - (1 + 0.2*Me*Me)**3.5
        eq3 = PbP_0_MR - ((1 + 0.2*M_0*M_0)/(1 + 0.2*M1*M1))**3.5
        eq4 = Pb/P_0 - PbP_0_MR
        eq5 = P_0/Pe - ((1 + 0.2*Me*Me)/(1 + 0.2*M_0*M_0))**3.5 
        eq6 = Pb/Pe - ((1 + 0.2*Me*Me)/(1 + 0.2*M1*M1))**3.5

        return [eq1, eq2, eq3, eq4, eq5, eq6]

    P0, Pe, M_0, M1, Pb, P_0 = fsolve(functMR, (100, 100/P0Pe, Me+1.5, Me+1.5, 0.25*100/P0Pe, 0.1*100/P0Pe))

    P0Pb_underMR[i] = P0/Pb

    def functRR(variables):
        (P0, Pe, M_0, M1, Pb, P_0) = variables
            
        eq1 = P0/Pe - P0Pe
        eq2 = P0Pe - (1 + 0.2*Me*Me)**3.5
        eq3 = PbP_0_RR - ((1 + 0.2*M_0*M_0)/(1 + 0.2*M1*M1))**3.5
        eq4 = Pb/P_0 - PbP_0_RR
        eq5 = P_0/Pe - ((1 + 0.2*Me*Me)/(1 + 0.2*M_0*M_0))**3.5 
        eq6 = Pb/Pe - ((1 + 0.2*Me*Me)/(1 + 0.2*M1*M1))**3.5

        return [eq1, eq2, eq3, eq4, eq5, eq6]

    P0, Pe, M_0, M1, Pb, P_0 = fsolve(functRR, (100, 100/P0Pe, Me+1.5, Me+1.5, 0.25*100/P0Pe, 0.1*100/P0Pe))

    P0Pb_underRR[i] = P0/Pb
    

plt.plot(M0[:35], P0Pe_PerExp[:35],  '.-', label = " Fully expanded jet.")
plt.plot(M0[20:35], P0Pb_overMR[20:],  '.-', label = "Overexpanded jet: von-Neumann criterion.")
plt.plot(M0[:35], P0Pb_overRR,  '.-', label = "Overexpanded jet: Detachment criterion.")
plt.plot(M0[20:35], P0Pb_underMR[20:],  '.-', label = "Underexpanded jet: von-Neumann criterion.")
plt.plot(M0[:35], P0Pb_underRR,  '.-', label = "Underexpanded jet: Detachment criterion.")

plt.title("CD Nozzle Flow")
plt.xlabel("M -->")
plt.ylabel("P0/Pb -->")
plt.legend()
plt.grid(linestyle = '--')
plt.savefig("Transition", dpi = 300, bbox_inches = 'tight')
plt.show()













