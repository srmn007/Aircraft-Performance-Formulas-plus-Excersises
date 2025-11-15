#--------------------------
#TAS To Lift Coefficient
def V_To_CL_Calc(Mass,V,rho,S):
    CL = (Mass*9.81) /(0.5*rho*(V**2)*S)
    return CL
#--------------------------
#Mach To TAS Coefficient
def Mach_To_V_Calc(M,rho,T):
    #Assuming Gamma = 1.4
    V = M* np.sqrt(287*1.4*T)
    return V
#--------------------------
#Thrust Calculation, Note: Thrust = Drag
def Thrust_Calc(Cd,rho,V,S):
    Thrust =Cd*(0.5*rho*(V**2)*S)
    return Thrust
#--------------------------
#Resultant Force Calculations, from lift and drag
def Aero_Force_Resultant(L,D):
    #Drag Is the Thrust
    R = np.sqrt(L**2+D**2)
    Theta = np.rad2deg(np.atan(D/L))
    return R, Theta
#--------------------------
#Drag Coefficient ns from Parasite drag coefficient ( CD_0) and induced drag factor and CL
def Cd_Calc(CD_0,k,Cl):
    Drag = CD_0 +k*(Cl**2)
    return Drag
#----------------------------
#Optimal Values; Jet-Max Range
def Optimal_Val_Min_T(L,S,rho,CD_0,k):

    v_min = np.sqrt( (L * 2 * np.sqrt(3 * k)) /(S * rho * np.sqrt(CD_0)) )
    CL_CD_OPT = 3 / (4 * np.power(3 *k * CD_0**3), 0.25)
    CL_Opt = np.sqrt(CD_0/ (3 * k))

    return v_min, CL_CD_OPT, CL_Opt
#----------------------------
#Optimal Values Minimal; Thrust, Prop-Max Range, Jet-Max Endurance
def Optimal_Val_Min_T(L,S,rho,CD_0,k):

    v_min = np.sqrt( (L * 2 * np.sqrt(k)) /(S * rho * np.sqrt(CD_0)) )
    CL_CD_OPT = 1/np.sqrt(4 * k * CD_0) 
    CL_Opt = np.sqrt(CD_0/k)
    return v_min, CL_CD_OPT, CL_Opt
#------------------------------------------
#Optimal Values for Minimal Power, Prop -Max Endurance
def Optimal_Val_Min_P(L,S,rho,CD_0,k):
    v_min = np.sqrt( (L * 2 * np.sqrt(k)) /(S * rho * np.sqrt(3 * CD_0)) )
    CL_CD_OPT = 0.25 * np.power( 27/ (k**3 * CD_0 ),0.25)
    CL_OPT = np.sqrt((3 *CD_0)/( k ))
    return v_min, CL_CD_OPT , CL_OPT
#------------------------------------------
#Density Calculation from ISA
def rho_Calc(P,T):
    density = P/(287*T)
    return density
#------------------------------------------
#Aspect Ratio Calculations
def AR_Calc(span,S):
    AR =(2*(span**2))/S
    return AR
#------------------------------------------
#Induced Drag factor
def Induced_Drag_Factor_Calc(AR,e):
    k =1/(np.pi*e*AR)
    return k
#------------------------------------------
#Lift coefficient to TAS
def CL_To_V_Calc(Mass,rho,S,Cl):
    V = np.sqrt(Mass*9.81/(0.5*rho*Cl*S))
    return V
#------------------------------------------
#TAS --> EAS
def EAS(TAS,rho):
    return TAS*np.sqrt(rho/1.2250)