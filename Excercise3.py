import numpy as np
import matplotlib.pyplot as plt
from Formulas import f


def main():
        #Given
        span = 20
        S = 65
        
        e = 0.9
        AspR =f.AR_Calc(span,S)
        k = f.Induced_Drag_Factor_Calc(AspR,e)
        Cd_0 = 0.029
        Speed_Range = np.arange(30,400,1)
        rho_FL0 = 1.2250
        def EX1():
            Mass = 25883 
            #Calculations
            Cl_3 = f.V_To_CL_Calc(Mass,Speed_Range,rho_FL0,S)
        
            Cd_3 = f.Cd_Calc(Cd_0,k,Cl_3)
            Fr = f.Thrust_Calc(Cd_3,rho_FL0,Speed_Range,S)
            Pr = Fr*Speed_Range

            Fr_min = min(Fr)
            Pr_min = min(Pr)
    
            #Calculate Minimal Speed, Table in course Note
            V_Min_T,CL_CD_Min_T,CL_Min_T = f.Optimal_Val_Min_T(Mass*9.81,S,rho_FL0,Cd_0,k) 
            V_Min_P,CL_CD_Min_P,CL_Min_P = f.Optimal_Val_Min_P(Mass*9.81,S,rho_FL0,Cd_0,k)

            print(f"the minimum thrust required {Fr_min:.3f} N")
            print(f"the minimum power required {Pr_min:.3f} W")
            print(f"the minimum thrust required velocity {V_Min_T:.3f} m/s")
            print(f"the minimum power required velocity {V_Min_P:.3f} m/s")
            ratio_V_minFR_MinPR = V_Min_P/V_Min_T
            print(f"the ratio of speed at minimal power on the speed at minimal thrust {ratio_V_minFR_MinPR:.3f}")
            Fr_At_Min_Pr = Pr_min/V_Min_P
            print(f"the ratio of the thrust required at minimal power on the thrust required at minimal thrust: {Fr_At_Min_Pr/Fr_min:.3f}")
            Pr_At_Min_Fr = Fr_min * V_Min_T
            print(f"the ratio of the power required at minimal power on the power required at minimal thrust: {Pr_min / Pr_At_Min_Fr:.3f}")
            print(f"the ratio of the lift coefficient at minimal power on the lift coefficient at minimal thrust: {CL_Min_P / CL_Min_T} is equal to square root of 3")
            print(f"the ratio of the drag coefficient at minimal power on the drag coefficient at minimal thrust: {(CL_Min_P/CL_CD_Min_P )/ (CL_Min_T/CL_CD_Min_T):.0f}")
            print(f"the ratio of the lift-to-drag ratio at minimal power on the lift-to-drag ratio at minimal thrust: {CL_CD_Min_P / CL_CD_Min_T} is equal to square root of 3/2: {np.sqrt(3)/2}")
            print(f"the ratio of the drag coefficient on the parasite drag coefficient at minimal thrust: {(CL_Min_T/CL_CD_Min_T) / Cd_0:.3f}")
            print(f"the ratio of the drag coefficient on the parasite drag coefficient at minimal power.: {(CL_Min_P/CL_CD_Min_P ) / Cd_0:.3f}")


            fig, axs = plt.subplots(1, 2, figsize=(12, 5))  
            # Left subplot: Required Thrust
            axs[0].plot(Speed_Range , Fr / 1000, color='blue')
            axs[0].set_xlabel('TAS [knots]')
            axs[0].set_ylabel('Required Thrust [kN]')
            axs[0].grid(True)
            axs[0].set_title('Required Thrust vs TAS')

            # Right subplot: Required Power
            axs[1].plot(Speed_Range , Pr / 1e6, color='red')
            axs[1].set_xlabel('TAS [knots]')
            axs[1].set_ylabel('Required Power [MW]')
            axs[1].grid(True)
            axs[1].set_title('Required Power vs TAS')

            # Adjust layout
            plt.tight_layout()
            plt.show()    
    #2)
        def Weights():
            Weights = np.arange(15000,30000,50)
            for Mass in Weights:
                
                Speed_Range = np.arange(30,400,1)
                rho_FL0 = 1.2250

                #Calculations
                Cl_3 = f.V_To_CL_Calc(Mass,Speed_Range,rho_FL0,S)
            
                Cd_3 = f.Cd_Calc(Cd_0,k,Cl_3)
                Fr = f.Thrust_Calc(Cd_3,rho_FL0,Speed_Range,S)
                Pr = Fr*Speed_Range

                Fr_min = min(Fr)
                Pr_min = min(Pr)
        
                #Calculate Minimal Speed, Table in course Note
                V_Min_T,CL_CD_Min_T,CL_Min_T = f.Optimal_Val_Min_T(Mass*9.81,S,rho_FL0,Cd_0,k) 
                V_Min_P,CL_CD_Min_P,CL_Min_P = f.Optimal_Val_Min_P(Mass*9.81,S,rho_FL0,Cd_0,k)

                print(f"the minimum thrust required {Fr_min:.3f} N")
                print(f"the minimum power required {Pr_min:.3f} W")
                print(f"the minimum thrust required velocity {V_Min_T:.3f} m/s")
                print(f"the minimum power required velocity {V_Min_P:.3f} m/s")
                ratio_V_minFR_MinPR = V_Min_P/V_Min_T
                print(f"the ratio of speed at minimal power on the speed at minimal thrust {ratio_V_minFR_MinPR:.3f}")
                Fr_At_Min_Pr = Pr_min/V_Min_P
                print(f"the ratio of the thrust required at minimal power on the thrust required at minimal thrust: {Fr_At_Min_Pr/Fr_min:.3f}")
                Pr_At_Min_Fr = Fr_min * V_Min_T
                print(f"the ratio of the power required at minimal power on the power required at minimal thrust: {Pr_min / Pr_At_Min_Fr:.3f}")
                print(f"the ratio of the lift coefficient at minimal power on the lift coefficient at minimal thrust: {CL_Min_P / CL_Min_T} is equal to square root of 3")
                print(f"the ratio of the drag coefficient at minimal power on the drag coefficient at minimal thrust: {(CL_Min_P/CL_CD_Min_P )/ (CL_Min_T/CL_CD_Min_T):.0f}")
                print(f"the ratio of the lift-to-drag ratio at minimal power on the lift-to-drag ratio at minimal thrust: {CL_CD_Min_P / CL_CD_Min_T} is equal to square root of 3/2: {np.sqrt(3)/2}")
                print(f"the ratio of the drag coefficient on the parasite drag coefficient at minimal thrust: {(CL_Min_T/CL_CD_Min_T) / Cd_0:.3f}")
                print(f"the ratio of the drag coefficient on the parasite drag coefficient at minimal power.: {(CL_Min_P/CL_CD_Min_P ) / Cd_0:.3f}")
        Weights()

if __name__ == "__main__":
    main()