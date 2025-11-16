import numpy as np

class f:
    g = 9.81          # Gravity
    R = 287.0         # Gas constant for air
    gamma = 1.4       # Specific heat ratio

    #----------------------------------------------------------
    # TAS → Lift Coefficient
    @staticmethod
    def V_To_CL_Calc(Mass, V, rho, S,gam=0):
        return (Mass * f.g * np.cos(np.deg2rad(gam))) / (0.5 * rho * V**2 * S)

    #----------------------------------------------------------
    # Mach → TAS
    @staticmethod
    def Mach_To_V_Calc(M, T):
        return M * np.sqrt(f.R * f.gamma * T)

    #----------------------------------------------------------
    # Thrust = Drag
    @staticmethod
    def Thrust_Calc(Cd, rho, V, S):
        return Cd * (0.5 * rho * V**2 * S)

    #----------------------------------------------------------
    # Thrust --> V
    @staticmethod
    def Thrust_To_V(Thrust, Cd, rho, S):
        return np.sqrt((2 * Thrust) / (S * Cd))

    #----------------------------------------------------------
    # Lift + Drag → Aerodynamic resultant
    @staticmethod
    def Aero_Force_Resultant(L, D):
        Res = np.sqrt(L**2 + D**2)
        Theta = np.degrees(np.atan(D / L))
        return Res, Theta

    #----------------------------------------------------------
    # Drag coefficient from CD0 + k CL²
    @staticmethod
    def Cd_Calc(CD_0, k, Cl):
        return CD_0 + k * (Cl**2)
    #----------------------------------------------------------

    #----------------------------------------------------------
    # Optimal Values; Jet-Max Range
    @staticmethod
    def Optimal_Val_Max_J_range(L, S, rho, CD_0, k):
        v_min = np.sqrt((L * 2 * np.sqrt(3 * k)) / (S * rho * np.sqrt(CD_0)))
        CL_CD_OPT = 3 / (4 * (3 * k * CD_0**3)**0.25)
        CL_Opt = np.sqrt(CD_0 / (3 * k))
        return v_min, CL_CD_OPT, CL_Opt

    #----------------------------------------------------------
    # Optimal Values Minimal Thrust; Prop-Max Range, Jet-Max Endurance
    @staticmethod
    def Optimal_Val_Min_T(L, S, rho, CD_0, k):
        v_min = np.sqrt((L * 2 * np.sqrt(k)) / (S * rho * np.sqrt(CD_0)))
        CL_CD_OPT = 1 / np.sqrt(4 * k * CD_0)
        CL_Opt = np.sqrt(CD_0 / k)
        return v_min, CL_CD_OPT, CL_Opt

    #----------------------------------------------------------
    # Optimal Values for Minimal Power; Prop – Max Endurance
    @staticmethod
    def Optimal_Val_Min_P(L, S, rho, CD_0, k):
        v_min = np.sqrt((L * 2 * np.sqrt(k)) / (S * rho * np.sqrt(3 * CD_0)))
        CL_CD_OPT = 0.25 * (27 / ((k**3) * CD_0))**0.25
        CL_OPT = np.sqrt((3 * CD_0) / k)
        return v_min, CL_CD_OPT, CL_OPT

    #----------------------------------------------------------
    # Density from ISA
    @staticmethod
    def rho_Calc(P, T):
        return P / (f.R * T)

    #----------------------------------------------------------
    # Aspect ratio
    @staticmethod
    def AR_Calc(span, S):
        return (span**2) / S

    #----------------------------------------------------------
    # Induced drag factor
    @staticmethod
    def Induced_Drag_Factor_Calc(AR, e):
        return 1 / (np.pi * e * AR)

    #----------------------------------------------------------
    # CL → TAS, gam == Flight Path Angle in degress
    @staticmethod
    def CL_To_V_Calc(Mass, rho, S, Cl,gam=0):
        
        return np.sqrt((Mass * f.g * np.cos(np.deg2rad(gam))) / (0.5 * rho * Cl * S))
        
    #----------------------------------------------------------
    # TAS → EAS
    @staticmethod
    def EAS(TAS, rho):
        return TAS * np.sqrt(rho / 1.225)

    #----------------------------------------------------------
    # CL -> Lift
    @staticmethod
    def CL_To_Lift_calc(Cl, V, rho, S):
        return 0.5 * Cl * rho * S * (V**2)

    #----------------------------------------------------------
    # Gam is to take into account the flight path angle gamma
    # W --> L
    @staticmethod
    def Weight_To_Lift_CalC(Mass, gam=0):
        return Mass * f.g * np.cos(np.deg2rad(gam))

    #----------------------------------------------------------
    # FR --> FA (Available Thrust)
    @staticmethod
    def FA_Calc(Thrust, Mass, gam):
        return Thrust + Mass * f.g * np.sin(np.deg2rad(gam))
     #----------------------------------------------------------
    # Returns V, intersections between FR and FA SLF
    def V_Int_FA_FR(FA,W,S,CD_0,k,rho):
        V2 = np.sqrt(((FA * W)/(W *S) + (W/S) * np.sqrt((FA/W)**2 -4 *CD_0*k) ) / (rho*CD_0))
        V1 = np.sqrt(((FA * W)/(W *S) - (W/S) * np.sqrt((FA/W)**2 -4 *CD_0*k) ) / (rho*CD_0))
        return V1,V2    
    #----------------------------------------------------------
    # Unit Conversions
    # lbf <--> N , x = 1 : lbf -> N
    @staticmethod
    def lbf__N(Force, x):
        return Force * (4.4482216153**x)

    #----------------------------------------------------------
    # mbar --> Pa
    @staticmethod
    def mbar_To_Pa(mbar):
        return mbar * 100

    #----------------------------------------------------------
    # Thrust @FL00  --> Thrust @FLXX
    @staticmethod
    def Thrust_0_Thrust_FL(Thrust, rho):
        return Thrust * np.sqrt(rho / 1.225)
    #----------------------------------------------------------
    # For intersection Between FA and FR, adjust tolerance, to find two values 
    # far from each other 
    @staticmethod
    def Find_Match_FA_FR(Fr, FA, tol):
        return np.where(np.abs(Fr - FA) <= tol)[0]