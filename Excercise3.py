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
        #Given
        span = 20
        S = 65
        Mass = 25883 
        e = 0.9
        AspR =AR(span,S)
        k = inducedDragCoeff(AspR,e)
        Cd_0 = 0.029
        Speed_Range = np.arange(30,400,1)
        rho_FL0 = 1.2250

        #Calculations
        Cl_3 = V_To_CL(Mass,Speed_Range,rho_FL0,S)

        Cd_3 = DragCoeff(Cd_0,k,Cl_3)

        Fr = Thrust(Cd_3,rho_FL0,Speed_Range,S)

        Pr = Fr*Speed_Range

        V_Min_T = np.sqrt((Mass*9.81*2*np.sqrt(k/Cd_0))/(S*(rho_FL0)))
        Fr_min = Mass*9.81/(np.max(Cl_3/Cd_3))

        V_Min_P = np.sqrt((Mass*9.81*2*np.sqrt(k/(3*Cd_0)))/(S*(rho_FL0)))
        Pr_min = Mass*9.81/(np.max((Cl_3**3/2)/Cd_3))

        print("the minimum thrust required",Fr_min)
        print("the minimum power required",Pr_min)
        print("the minimum thrust required velocity",V_Min_T)
        print("the minimum power required velocity",V_Min_P)
        ratio_V_minFR_MinPR = V_Min_P/V_Min_T
        print("the ratio of speed at minimal power on the speed at minimal thrust;",ratio_V_minFR_MinPR)
        print("the ratio of the thrust required at minimal power on the thrust required at minimal thrust;",)

        fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # adjust figsize as needed

        # Left subplot: Required Thrust
        axs[0].plot(Speed_Range * 1.94384, Fr / 1000, color='blue')
        axs[0].set_xlabel('TAS [knots]')
        axs[0].set_ylabel('Required Thrust [kN]')
        axs[0].grid(True)
        axs[0].set_title('Required Thrust vs TAS')

        # Right subplot: Required Power
        axs[1].plot(Speed_Range * 1.94384, Pr / 1e6, color='red')
        axs[1].set_xlabel('TAS [knots]')
        axs[1].set_ylabel('Required Power [MW]')
        axs[1].grid(True)
        axs[1].set_title('Required Power vs TAS')

        # Adjust layout
        plt.tight_layout()
        plt.show()    

if __name__ == "__main__":
    main()