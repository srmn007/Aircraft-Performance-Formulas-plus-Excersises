import numpy as np
import matplotlib.pyplot as plt
from Formulas import f

def main():
    #Propeller Properties
    span_prop = 11.4
    S_prop = 19.2
    mass_prop = 2359
    e = 0.9
    eff_prop = 0.83
    Thrust_A = 2*186*1000 #W
    Pa_prop = Thrust_A* eff_prop
    CD_0 = 0.025
    CL_M = 1.16
    
    #Calculations
    AR = f.AR_Calc(span_prop,S_prop)
    k = f.Induced_Drag_Factor_Calc(AR,e)
    
    #Sea Level
    T_FL0 = 288.16
    P_FL0 = 101325
    Rho_FL0 = f.rho_Calc(P_FL0,T_FL0)

    V_Range = np.arange(5,110,5)
    Cl_prop_range = f.V_To_CL_Calc(mass_prop,V_Range,Rho_FL0,S_prop)
    Cd_prop_range = f.Cd_Calc(CD_0,k,Cl_prop_range)
    Thrust_R = f.Thrust_Calc(CD_0,Rho_FL0,V_Range,S_prop)
    gam = 0
    #First Iteration gam =0
    R_C = ((Pa_prop)/(mass_prop*9.81)) - V_Range*((Thrust_R/(mass_prop*9.81))+(2*k*mass_prop*9.81)/(Rho_FL0*S_prop*(V_Range**2)))
    V_Hor = np.sqrt(V_Range**2 -(R_C)**2)
   
    gam1 = (np.arcsin(R_C/V_Range))
    
    R_C = ((Pa_prop)/(mass_prop*9.81)) - V_Range*((Thrust_R/(mass_prop*9.81))+(2*k*mass_prop*9.81*(np.cos(gam1)**2))/(Rho_FL0*S_prop*(V_Range**2)))
    V_Hor2 = np.sqrt(V_Range**2 -(R_C)**2)
    gam2 = np.rad2deg(np.arcsin(R_C/V_Hor2))
    idx =np.nanargmax(gam2)

    plt.figure(figsize=(7,6))

    plt.plot(V_Hor2,R_C, linewidth=2)
    plt.plot(V_Hor2[idx],R_C[idx],"o",color="red",linewidth=5)

    plt.xlabel("Speed Horizontal (m/s)")
    plt.ylabel("Climb rate (m/s)")
    plt.title("Vertical vs Horizontal Speed Profile")
    plt.grid(True, linestyle=":")
    plt.tight_layout()
    plt.show()
if __name__ == "__main__":
    main() 