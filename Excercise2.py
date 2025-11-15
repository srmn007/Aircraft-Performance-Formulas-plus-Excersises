import numpy as np
import matplotlib.pyplot as plt

g = 9.81          # Gravity
R = 287.0         # Gas constant for air
gamma = 1.4       # Specific heat ratio

#----------------------------------------------------------
# TAS → Lift Coefficient
def V_To_CL_Calc(Mass, V, rho, S):
    return (Mass * g) / (0.5 * rho * V**2 * S)
#----------------------------------------------------------
# Mach → TAS
def Mach_To_V_Calc(M, rho, T):
    return M * np.sqrt(R * 1.4 * T)
#----------------------------------------------------------
# Thrust = Drag
def Thrust_Calc(Cd, rho, V, S):
    return Cd * (0.5 * rho * V**2 * S)
#----------------------------------------------------------
# Lift + Drag → Aerodynamic resultant
def Aero_Force_Resultant(L, D):
    Res = np.sqrt(L**2 + D**2)
    Theta = np.degrees(np.atan(D / L))
    return Res, Theta
#----------------------------------------------------------
# Drag coefficient from CD0 + k CL²
def Cd_Calc(CD_0, k, Cl):
    return CD_0 + k * (Cl**2)
#----------------------------------------------------------
# ----------------------------------------
# Optimal Values; Jet-Max Range
# ----------------------------------------
def Optimal_Val_Max_J_range(L, S, rho, CD_0, k):
    v_min = np.sqrt((L * 2 * np.sqrt(3 * k)) / (S * rho * np.sqrt(CD_0)))
    CL_CD_OPT = 3 / (4 * (3 * k * CD_0**3)**0.25)
    CL_Opt = np.sqrt(CD_0 / (3 * k))
    return v_min, CL_CD_OPT, CL_Opt

#----------------------------------------------------------
# Optimal Values Minimal Thrust; Prop-Max Range, Jet-Max Endurance
# ----------------------------------------
def Optimal_Val_Min_T(L, S, rho, CD_0, k):
    v_min = np.sqrt((L * 2 * np.sqrt(k)) / (S * rho * np.sqrt(CD_0)))
    CL_CD_OPT = 1 / np.sqrt(4 * k * CD_0)
    CL_Opt = np.sqrt(CD_0 / k)
    return v_min, CL_CD_OPT, CL_Opt
#----------------------------------------------------------
# Optimal Values for Minimal Power; Prop – Max Endurance
def Optimal_Val_Min_P(L, S, rho, CD_0, k):
    v_min = np.sqrt((L * 2 * np.sqrt(k)) / (S * rho * np.sqrt(3 * CD_0)))
    CL_CD_OPT = 0.25 * (27 / (k**3 * CD_0))**0.25
    CL_OPT = np.sqrt((3 * CD_0) / k)
    return v_min, CL_CD_OPT, CL_OPT
#----------------------------------------------------------
# Density from ISA
def rho_Calc(P, T):
    return P / (R * T)
#----------------------------------------------------------
# Aspect ratio
def AR_Calc(span, S):
    return (span**2) / S   
#----------------------------------------------------------
# Induced drag factor
def Induced_Drag_Factor_Calc(AR, e):
    return 1 / (np.pi * e * AR)
#----------------------------------------------------------
# CL → TAS
def CL_To_V_Calc(Mass, rho, S, Cl):
    return np.sqrt(Mass * g / (0.5 * rho * Cl * S))

#----------------------------------------------------------
# TAS → EAS
def EAS(TAS, rho):
    return TAS * np.sqrt(rho / 1.225)

def main():
        span = 20 #m
        S = 65 #m^2
        Mass = 25833
        e =0.9
        Cl_Max =1.16
        k =1/(np.pi*e*AR(span,S))

        cd = DragCoeff(0.029,k,Cl_Max)
    #@ Sea Level

    #@20 000 ft
        P_FL20 = 465.6 *100
        T_FL20 = 248.53 
        rho_FL20 = rho(P_FL20,T_FL20)
        V_Stall = CL_To_V(Mass,rho_FL20,S,Cl_Max)

        Speed_range =np.arange(20,600,1)
        
        Cl_FL20 = V_To_CL(Mass,Speed_range,rho_FL20,S)

        Cd_FL20 = DragCoeff(0.029,k,Cl_FL20)

        T_FL20 = Thrust(Cd_FL20,rho_FL20,Speed_range,S)
        T_A_FL20 = 22 *rho_FL20/1.2250

      
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