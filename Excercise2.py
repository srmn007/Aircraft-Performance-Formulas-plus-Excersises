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
    FA = f.lbf__N(22,1) #From Graph interpolated, 22 lbf for N1 = 100%

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

      
    plt.figure()

    plt.plot(Speed_range,T_FL20/4448.221615,label='thrust required')
    plt.vlines(V_Stall, 0, np.max(T_FL20/4448.221615),color='red',label='V_stall')
    plt.hlines(T_A_FL20,np.min(Speed_range),np.max(Speed_range),color='yellow',label='Thrust available')
    plt.title("At 20 000 ft")
    plt.xlabel("TAS [m/s]")
    plt.ylabel("Thrust klbf")
    plt.legend()
    plt.grid()
    plt.show()
    T_A_FL20 = f.lbf__N(T_A_FL20*1000,1) #In the graph it is given in klbf
    V_Max_FL20 = f.Thrust_To_V(T_A_FL20,Cd_FL20,rho_FL20,S)
    
    #@ 30 000 ft
    P_FL30 = 300.9 *100
    T_FL30 = 228.71 
    rho_FL30 = rho(P_FL30,T_FL30)
    V_Stall = CL_To_V(Mass,rho_FL30,S,Cl_Max)
    
    Cl_FL30 = V_To_CL(Mass,Speed_range,rho_FL30,S)

    Cd_FL30 = DragCoeff(0.029,k,Cl_FL30)

    T_FL30 = Thrust(Cd_FL30,rho_FL30,Speed_range,S)
    T_A_FL30 = 22 * (rho_FL30/1.225)
    plt.figure()

    plt.plot(Speed_range,T_FL30/4448.221615)
    plt.vlines(V_Stall, 0, np.max(T_FL30/4448.221615),color='red',label='V_stall')
    plt.hlines(22,np.min(Speed_range),np.max(Speed_range),color='yellow',label='Thrust available')
    plt.title("At 30 000 ft")
    plt.xlabel("TAS [m/s]")
    plt.grid()
    plt.legend()
    plt.ylabel("Thrust klbf")
    plt.show()

    # 2)
    delta_T = 0.5
    v1 = 126.4
    FA = (rho_FL20/1.2250)*26*4448.221615
    t1 =0
    while v1 < 158:
        CL_I = V_To_CL(Mass,v1,rho_FL20,S)
        CD_I = DragCoeff(0.029,k,CL_I)
        Fr = Thrust(CD_I,rho_FL20,v1,S)
        Delta_V = delta_T *((FA-Fr)/(Mass))
        v1 = v1 + Delta_V
        t1 = t1 + delta_T
        print(t1)
    
if __name__ == "__main__":
    main()