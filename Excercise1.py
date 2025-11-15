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
    
    #BOEING 777_200    
        S_B = 427.8 #m^2
        Mass_B = 304000 #kg
        P = 196.8 #mbar
        P_pa = P*100 #Pa
        T = 216.83
        L = Mass_B * 9.81
        #Calculations
        rho_E = rho_Calc(P_pa,T)
        V_B = Mach_To_V_Calc(0.82,rho_E,216.83)
        CL_B = V_To_CL_Calc(Mass_B,V_B,rho_E,S_B)
        CD_B = Cd_Calc(0.018,0.045,CL_B)

        FR_B = Thrust_Calc(CD_B,rho_E,V_B,S_B)
        F_B,theta = Aero_Force_Resultant(L,FR_B)
        ratio_CL_CD = Mass_B*9.81/FR_B
        v_min_B, CL_CD_OPT_B, CL_Opt_B = Optimal_Val_Min_T(L,S_B,rho_E,0.045,0.018)

    #F_16
        S_F16 = 27.87
        Mass_F16 = 11000
        L_F16 =Mass_F16*9.81
        CL_F16 = V_To_CL_Calc(Mass_F16,V_B,rho_E,S_F16)
        CD_F16 = Cd_Calc(0.018,0.045,CL_F16)

        FR_F16 = Thrust_Calc(CD_F16,rho_E,V_B,S_F16)
      
        F_F16,theta_F16 = Aero_Force_Resultant(L_F16,FR_F16)
        ratio_CL_CD = Mass_F16*9.81/FR_F16
    #2)
        CL_CD_MAX_B = 1/np.sqrt(4*0.045*0.018)    
        CL_CD_MAX_F16 = 1/np.sqrt(4*0.0117*0.0169)

    #4)
        
        v_min_B, CL_CD_OPT_B, CL_Opt_B = Optimal_Val_Min_T(L,S_B,rho_E,0.045,0.018)
        Cl_Min_Drag = V_To_CL_Calc(Mass_B,v_min_B,rho_E,)
    
    
if __name__ == "__main__":
    main()