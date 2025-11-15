import numpy as np
import matplotlib.pyplot as plt
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
def main():
    
    def Exercise1():
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
        print(FR_B)
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
    def Exercise2():
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
    def Exercise3():
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
    Exercise1()
if __name__ == "__main__":
    main()