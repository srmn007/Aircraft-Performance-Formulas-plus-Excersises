import numpy as np
import matplotlib.pyplot as plt
from Formulas import f 

def main():
    #Given 
    span = 20 #m
    S = 65 #m^2
    Mass = 25833
    e  = 0.9
    Cl_Max =1.16
    AR = f.AR_Calc(span,S)
    k = f.Induced_Drag_Factor_Calc(AR,e)
    cd = f.Cd_Calc(0.029,k,Cl_Max)

    #@20 000 ft
    P_FL20 = f.mbar_To_Pa(465.6) 
    T_FL20 = 248.53 
    rho_FL20 = f.rho_Calc(P_FL20,T_FL20)
    V_Stall = f.CL_To_V_Calc(Mass,rho_FL20,S,Cl_Max)

    print(f"V Stall = {V_Stall:.2f}")
    
    Speed_range =np.arange(20,600,1)
        
    Cl_FL20 = f.V_To_CL_Calc(Mass,Speed_range,rho_FL20,S)

    Cd_FL20 = f.Cd_Calc(0.029,k,Cl_FL20)

    T_FL20 = f.Thrust_Calc(Cd_FL20,rho_FL20,Speed_range,S)

    T_A_FL20 = f.Thrust_0_Thrust_FL(22,rho_FL20)

    #To find Max en Minimum speed, look for the intersections of FA, with FR, 
    
    V_Min_intersect, V_max = f.V_Int_FA_FR(f.lbf__N(T_A_FL20*1000,1),Mass*9.81,S,0.029,k,rho_FL20)

    print(f"The Plane must no go below the stall speed: V stall = {V_Stall} m/s")
    
    plt.figure()
    plt.plot(Speed_range,f.lbf__N(T_FL20/1000,-1),label='thrust required')
    plt.vlines(V_Stall, 0, np.max(f.lbf__N(T_FL20/1000,-1)),color='red',label='V_stall')
    plt.vlines(V_Min_intersect, 0, np.max(f.lbf__N(T_FL20/1000,-1)),color='black',label='V_min_intersect')
    plt.vlines(np.round(V_max,2), 0, np.max(f.lbf__N(T_FL20/1000,-1)),color='green',label='V_Max')
    plt.hlines(T_A_FL20,np.min(Speed_range),np.max(Speed_range),color='yellow',label='Thrust available')
    plt.title("At 20 000 ft")
    plt.xlabel("TAS [m/s]")
    plt.ylabel("Thrust klbf")
    plt.legend()
    plt.grid()
    plt.show()

    #@ 30 000 ft
    P_FL30 = f.mbar_To_Pa(300.9) 
    T_FL30 = 228.71 
    rho_FL30 = f.rho_Calc(P_FL30,T_FL30)
    
    V_Stall = f.CL_To_V_Calc(Mass,rho_FL30,S,Cl_Max)

    print(f"V Stall = {V_Stall:.2f}")
            
    Cl_FL30 = f.V_To_CL_Calc(Mass,Speed_range,rho_FL30,S)

    Cd_FL30 = f.Cd_Calc(0.029,k,Cl_FL30)

    T_FL30 = f.Thrust_Calc(Cd_FL30,rho_FL30,Speed_range,S)

    T_A_FL30 = f.Thrust_0_Thrust_FL(22,rho_FL30)

    #To find Max en Minimum speed, look for the intersections of FA, with FR, 
    V_Min_intersect, V_max = f.V_Int_FA_FR(f.lbf__N(T_A_FL30*1000,1),Mass*9.81,S,0.029,k,rho_FL30)

    print(f"The Plane must no go below the stall speed: V stall = {V_Stall} m/s")
    
    plt.figure()
    plt.plot(Speed_range,f.lbf__N(T_FL30/1000,-1),label='thrust required')
    plt.vlines(V_Stall, 0, np.max(f.lbf__N(T_FL30/1000,-1)),color='red',label='V_stall')
    plt.vlines(V_Min_intersect, 0, np.max(f.lbf__N(T_FL30/1000,-1)),color='black',label='V_min_intersect')
    plt.vlines(V_max, 0, np.max(f.lbf__N(T_FL30/1000,-1)),color='green',label='V_Max')
    plt.hlines(T_A_FL30,np.min(Speed_range),np.max(Speed_range),color='yellow',label='Thrust available')
    plt.title("At 30 000 ft")
    plt.xlabel("TAS [m/s]")
    plt.ylabel("Thrust klbf")
    plt.legend()
    plt.grid()
    plt.show()

    # 2)
    delta_T = 0.5
    v1 = 126.4
    FA = (rho_FL20/1.2250)*26*4448.221615
    t1 =0
    while v1 < 158:
        CL_I = f.V_To_CL_Calc(Mass,v1,rho_FL20,S)
        CD_I = f.Cd_Calc(0.029,k,CL_I)
        Fr = f.Thrust_Calc(CD_I,rho_FL20,v1,S)
        Delta_V = delta_T *((FA-Fr)/(Mass))
        v1 = v1 + Delta_V
        t1 = t1 + delta_T

    print(f"Time to Accelerate: {t1} s")
    
if __name__ == "__main__":
    main()