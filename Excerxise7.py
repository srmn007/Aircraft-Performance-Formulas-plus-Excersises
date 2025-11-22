import numpy as np
import matplotlib.pyplot as plt
from Formulas import f

def main():
    #Propeller Properties
    S = 511
    mass = 700 *1000#klb
    mass = f.lb__kg(mass,1)
    Thrust_A = 4*200 *1000
    CL_Max_Take_off = 1.584
    mu = 0.025
    T_FL0 = 288.16
    P_FL0 = 101325
    rho = f.rho_Calc(P_FL0,T_FL0)
    #Calculation
    V_LOF = 1.2 *np.sqrt((2*mass*9.81)/(rho*S*CL_Max_Take_off))
    print(f"Lift of speed: {V_LOF:.2f}")
    time_during_ground_roll = 45 #s from Appendix H
    acceleration = ((mass * (V_LOF- f.knots_To_m_s(10))**2)/2)*(1/(Thrust_A - mu*mass*9.81 -()))
    print(f"Acceleration : {acceleration} m/s^2")

if __name__ == "__main__":
    main() 