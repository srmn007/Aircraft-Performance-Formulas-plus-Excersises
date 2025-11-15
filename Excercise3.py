import numpy as np
import matplotlib.pyplot as plt
from Formulas import f


def main():
        #Given
        span = 20
        S = 65
        Mass = 25883 
        e = 0.9
        AspR =f.AR_Calc(span,S)
        k = f.Induced_Drag_Factor_Calc(AspR,e)
        Cd_0 = 0.029
        Speed_Range = np.arange(30,400,1)
        rho_FL0 = 1.2250

        #Calculations
        Cl_3 = f.V_To_CL_Calc(Mass,Speed_Range,rho_FL0,S)
        Cd_3 = f.Cd_Calc(Cd_0,k,Cl_3)
        Fr = f.Thrust_Calc(Cd_3,rho_FL0,Speed_Range,S)
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