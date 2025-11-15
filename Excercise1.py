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
    
    # Boeing 777-200LR
    S_B = 427.8
    m_B = 304000
    CD0_B = 0.018
    k_B = 0.045

    # F16
    S_F = 27.87
    m_F = 11000
    CD0_F = 0.0169
    k_F = 0.117

    # Atmosphere FL390
    T = 218.8       # K
    P = 18700       # Pa
    rho = rho_Calc(P, T)

    M = 0.82
    V = Mach_To_V_Calc(M, rho, T)

    # ----------------------------------------------------------
    # 1. Aerodynamic forces at FL390
    # ----------------------------------------------------------
    CL_B = V_To_CL_Calc(m_B, V, rho, S_B)
    CD_B = Cd_Calc(CD0_B, k_B, CL_B)
    D_B = Thrust_Calc(CD_B, rho, V, S_B)
    L_B = m_B * g
    Res_B, Theta_B = Aero_Force_Resultant(L_B, D_B)

    CL_F = V_To_CL_Calc(m_F, V, rho, S_F)
    CD_F = Cd_Calc(CD0_F, k_F, CL_F)
    D_F = Thrust_Calc(CD_F, rho, V, S_F)
    L_F = m_F * g
    Res_F, Theta_F = Aero_Force_Resultant(L_F, D_F)

    print("=== 1. Aerodynamic forces at FL390 ===")
    print(f"B777 L/D = {L_B/D_B:.2f}, Resultant={Res_B:.0f} N, Angle={Theta_B:.2f}°")
    print(f"F16  L/D = {L_F/D_F:.2f}, Resultant={Res_F:.0f} N, Angle={Theta_F:.2f}°")

    # ----------------------------------------------------------
    # 2. Max L/D
    # ----------------------------------------------------------
    _, LDmax_B, _ = Optimal_Val_Max_J_range(m_B*g, S_B, rho, CD0_B, k_B)
    _, LDmax_F, _ = Optimal_Val_Max_J_range(m_F*g, S_F, rho, CD0_F, k_F)

    print("\n=== 2. Maximum L/D ===")
    print(f"B777 max L/D ≈ {LDmax_B:.2f}")
    print(f"F16  max L/D ≈ {LDmax_F:.2f}")

    # ----------------------------------------------------------
    # 3. Drag vs TAS & EAS
    # ----------------------------------------------------------
    V_range = np.linspace(80, 350, 300)

    def drag_curve(V, S, m, CD0, k):
        CL = V_To_CL_Calc(m, V, rho, S)
        CD = Cd_Calc(CD0, k, CL)
        D = Thrust_Calc(CD, rho, V, S)
        return D

    D_B_curve = drag_curve(V_range, S_B, m_B, CD0_B, k_B)
    D_F_curve = drag_curve(V_range, S_F, m_F, CD0_F, k_F)

    EAS_range = EAS(V_range, rho)

    plt.figure(figsize=(10,6))
    plt.plot(V_range, D_B_curve, label="B777 Drag (TAS)")
    plt.plot(V_range, D_F_curve, label="F-16 Drag (TAS)")
    plt.xlabel("TAS (m/s)")
    plt.ylabel("Drag (N)")
    plt.legend()
    plt.grid()
    plt.title("Drag vs TAS at FL390")
    plt.show()

    plt.figure(figsize=(10,6))
    plt.plot(EAS_range, D_B_curve, label="B777 Drag (EAS)")
    plt.plot(EAS_range, D_F_curve, label="F-16 Drag (EAS)")
    plt.xlabel("EAS (m/s)")
    plt.ylabel("Drag (N)")
    plt.legend()
    plt.grid()
    plt.title("Drag vs EAS at FL390")
    plt.show()

    # ----------------------------------------------------------
    # 4. Airspeed that minimizes thrust --> CL_opt
    # ----------------------------------------------------------
    _, _, CLopt_B = Optimal_Val_Min_T(m_B*g, S_B, rho, CD0_B, k_B)
    _, _, CLopt_F = Optimal_Val_Min_T(m_F*g, S_F, rho, CD0_F, k_F)

    V_opt_B = CL_To_V_Calc(m_B, rho, S_B, CLopt_B)
    V_opt_F = CL_To_V_Calc(m_F, rho, S_F, CLopt_F)

    print("\n=== 4. Airspeed that minimizes Thrust ===")
    print(f"B777: V_opt = {V_opt_B:.1f} m/s, EAS = {EAS(V_opt_B, rho):.1f} m/s")
    print(f"F16 : V_opt = {V_opt_F:.1f} m/s, EAS = {EAS(V_opt_F, rho):.1f} m/s")
    
    
if __name__ == "__main__":
    main()