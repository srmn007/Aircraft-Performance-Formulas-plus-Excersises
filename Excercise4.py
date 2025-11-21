import numpy as np
import matplotlib.pyplot as plt
from Formulas import f

def main():
    #Propeller Properties
    span = 20
    S = 65
    mass_t = 25883 #Mass aircraft +Fuel
    mass_f = 6754
    mass_aircraft = mass_t - mass_f 
    e = 0.9
    CD_0 = 0.029
    CL_M = 1.16
    Masses = [mass_t,22506,19129]
    #Calculations
    AR = f.AR_Calc(span,S)
    k = f.Induced_Drag_Factor_Calc(AR,e)
    TSFC = 0.70 #lb/hr.lb
    TSFC_N = 0.70* 2.835e-5
    #Sea Level
    T_FL0 = 288.16
    P_FL0 = 101325
    Rho_FL0 = f.rho_Calc(P_FL0,T_FL0)

    print("-----------At Sea Level ------------\n Max Range")
    for i in Masses:
        V_opt_FL0 , CL_CD_OPT , CL_OPT = f.Optimal_Val_Max_J_range(i*9.81,S,Rho_FL0,CD_0,k)
        V_stall = f.CL_To_V_Calc(i,Rho_FL0,S,CL_M)
        print(f"Optimal speed for {i} kg:{V_opt_FL0:.2f} m/s")
        print(f"Mach for {i} kg:{f.V_To_Mach_Calc(V_opt_FL0,T_FL0):.2f}")
        print(f"Thrust for {i} kg:{f.lbf__N(f.Thrust_Calc(f.Cd_Calc(CD_0,k,CL_OPT),Rho_FL0,V_opt_FL0,S)/1000,-1):.2f} klbf")
        Range = (V_opt_FL0/(9.81*TSFC_N))*(CL_OPT/f.Cd_Calc(CD_0,k,CL_OPT)) *np.log(i/(i-mass_f))
        print(f"Range for {i} kg:{Range/1000:.2f} km \n")

    print(f"CL {CL_OPT:.3f}")
    print(f"CD {f.Cd_Calc(CD_0,k,CL_OPT):.3f}")

    print("-----------At Sea Level ------------\nMax Endurance")

    for i in Masses:
        V_opt_FL0 , CL_CD_OPT , CL_OPT = f.Optimal_Val_Min_T(i*9.81,S,Rho_FL0,CD_0,k)
        V_stall = f.CL_To_V_Calc(i,Rho_FL0,S,CL_M)
        print(f"Optimal speed for {i} kg:{V_opt_FL0:.2f} m/s")
        print(f"Mach for {i} kg:{f.V_To_Mach_Calc(V_opt_FL0,T_FL0):.2f}")
        print(f"Thrust for {i} kg:{f.lbf__N(f.Thrust_Calc(f.Cd_Calc(CD_0,k,CL_OPT),Rho_FL0,V_opt_FL0,S)/1000,-1):.2f} klbf")
        Endurance = (1/(9.81*TSFC_N))*(CL_OPT/f.Cd_Calc(CD_0,k,CL_OPT)) *np.log(i/(i-mass_f))
        print(f"Endurance for {i} kg:{Endurance/3600:.2f} hr \n")

    print(f"CL {CL_OPT:.3f}")
    print(f"CD {f.Cd_Calc(CD_0,k,CL_OPT):.3f}")
    #20 000 ft
    T_FL20 = 248.53
    P_FL20 = f.mbar_To_Pa(465.6)
    Rho_FL20 = f.rho_Calc(P_FL20,T_FL20)
    
    print("-----------At 20 000 ft ------------\nMax Range")
   
    TSFC_N = 0.6* 2.835e-5
    for i in Masses:
        V_opt_FL20 , CL_CD_OPT , CL_OPT = f.Optimal_Val_Max_J_range(i*9.81,S,Rho_FL20,CD_0,k)
        V_stall = f.CL_To_V_Calc(i,Rho_FL20,S,CL_M)
        print(f"Optimal speed for {i} kg:{V_opt_FL20:.2f} m/s")
        print(f"Mach for {i} kg:{f.V_To_Mach_Calc(V_opt_FL20,T_FL20 ):.2f}")
        print(f"Thrust for {i} kg:{f.lbf__N(f.Thrust_Calc(f.Cd_Calc(CD_0,k,CL_OPT),Rho_FL20,V_opt_FL20,S)/1000,-1):.2f} klbf")
        Range = (V_opt_FL20/(9.81*TSFC_N))*(CL_OPT/f.Cd_Calc(CD_0,k,CL_OPT)) *np.log(i/(i-mass_f))
        print(f"Range for {i} kg:{Range/1000:.2f} km \n")

    print(f"CL {CL_OPT:.3f}")
    print(f"CD {f.Cd_Calc(CD_0,k,CL_OPT):.3f}")

    print("-----------At 20 000 ft ------------\nMax Endurance")

    for i in Masses:
        V_opt_FL20 , CL_CD_OPT , CL_OPT = f.Optimal_Val_Min_T(i*9.81,S,Rho_FL20,CD_0,k)
        V_stall = f.CL_To_V_Calc(i,Rho_FL20,S,CL_M)
        print(f"Optimal speed for {i} kg:{V_opt_FL20:.2f} m/s")
        print(f"Mach for {i} kg:{f.V_To_Mach_Calc(V_opt_FL20,T_FL20 ):.2f}")
        print(f"Thrust for {i} kg:{f.lbf__N(f.Thrust_Calc(f.Cd_Calc(CD_0,k,CL_OPT),Rho_FL20,V_opt_FL20,S)/1000,-1):.2f} klbf")
        Endurance = (1/(9.81*TSFC_N))*(CL_OPT/f.Cd_Calc(CD_0,k,CL_OPT)) *np.log(i/(i-mass_f))
        print(f"Endurance for {i} kg:{Endurance/3600:.2f} hr \n")

    print(f"CL {CL_OPT:.3f}")
    print(f"CD {f.Cd_Calc(CD_0,k,CL_OPT):.3f}")
if __name__ == "__main__":
    main() 