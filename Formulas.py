import numpy as np

class f:
    g = 9.81
    R = 287.0
    gamma = 1.4

    @staticmethod
    def V_To_CL_Calc(Mass, V, rho, S):
        return (Mass * f.g) / (0.5 * rho * V**2 * S)

    @staticmethod
    def Mach_To_V_Calc(M, rho, T):
        return M * np.sqrt(f.R * f.gamma * T)

    @staticmethod
    def Thrust_Calc(Cd, rho, V, S):
        return Cd * (0.5 * rho * V**2 * S)

    @staticmethod
    def Thrust_To_V(Thrust, Cd, rho, S):
        return np.sqrt((2 * Thrust) / (S * Cd))

    @staticmethod
    def Aero_Force_Resultant(L, D):
        Res = np.sqrt(L**2 + D**2)
        Theta = np.degrees(np.atan(D / L))
        return Res, Theta

    @staticmethod
    def Cd_Calc(CD_0, k, Cl):
        return CD_0 + k * (Cl**2)

    @staticmethod
    def Optimal_Val_Max_J_range(L, S, rho, CD_0, k):
        v_min = np.sqrt((L * 2 * np.sqrt(3 * k)) / (S * rho * np.sqrt(CD_0)))
        CL_CD_OPT = 3 / (4 * (3 * k * CD_0**3)**0.25)
        CL_Opt = np.sqrt(CD_0 / (3 * k))
        return v_min, CL_CD_OPT, CL_Opt

    @staticmethod
    def Optimal_Val_Min_T(L, S, rho, CD_0, k):
        v_min = np.sqrt((L * 2 * np.sqrt(k)) / (S * rho * np.sqrt(CD_0)))
        CL_CD_OPT = 1 / np.sqrt(4 * k * CD_0)
        CL_Opt = np.sqrt(CD_0 / k)
        return v_min, CL_CD_OPT, CL_Opt

    @staticmethod
    def Optimal_Val_Min_P(L, S, rho, CD_0, k):
        v_min = np.sqrt((L * 2 * np.sqrt(k)) / (S * rho * np.sqrt(3 * CD_0)))
        CL_CD_OPT = 0.25 * (27 / (k**3 * CD_0))**0.25
        CL_OPT = np.sqrt((3 * CD_0) / k)
        return v_min, CL_CD_OPT, CL_OPT

    @staticmethod
    def rho_Calc(P, T):
        return P / (f.R * T)

    @staticmethod
    def AR_Calc(span, S):
        return (span**2) / S

    @staticmethod
    def Induced_Drag_Factor_Calc(AR, e):
        return 1 / (np.pi * e * AR)

    @staticmethod
    def CL_To_V_Calc(Mass, rho, S, Cl):
        return np.sqrt(Mass * f.g / (0.5 * rho * Cl * S))

    @staticmethod
    def EAS(TAS, rho):
        return TAS * np.sqrt(rho / 1.225)

    @staticmethod
    def CL_To_Lift_calc(Cl, V, rho, S):
        return 0.5 * Cl * rho * S * (V**2)

    @staticmethod
    def Weight_To_Lift_CalC(Mass, gam):
        return Mass * f.g * np.cos(np.deg2rad(gam))

    @staticmethod
    def FA_Calc(Thrust, Mass, gam):
        return Thrust + Mass * f.g * np.sin(np.deg2rad(gam))

    @staticmethod
    def lbf__N(Force, x):
        return Force * (4.4482216153**x)

    @staticmethod
    def mbar_To_Pa(mbar):
        return mbar * 100

    @staticmethod
    def Thrust_0_Thrust_FL(Thrust, rho):
        return Thrust * np.sqrt(rho / 1.225)
