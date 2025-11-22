import numpy as np
import matplotlib.pyplot as plt
from Formulas import f

def main():
    #Propeller Properties
    S = 511
    mass = 700 *1000#klb
    mass = f.lb__kg(mass,1)
    CL_M = 1.16
    dCL_alpha = 0.1
    V = f.Mach_To_V_Calc(0.8,228.71)

    #Calculations
    U_de = 18 #m/s from graph in course Note
    W_S = mass*9.81/S
    print(W_S)
    K = 0.8 - (5.25/(W_S)**(3/4))
    print(K)
    n = 1 + dCL_alpha *(rhp)
if __name__ == "__main__":
    main() 