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
    # TAS → Mach
    @staticmethod
    def V_To_Mach_Calc(V, T):
        return V / np.sqrt(f.R * f.gamma * T)
    #----------------------------------------------------------
    # Thrust_required = Drag
    @staticmethod
    def Thrust_Calc(Cd, rho, V, S):
        return Cd * (0.5 * rho * V**2 * S)

    #----------------------------------------------------------
    # Thrust --> V
    @staticmethod
    def Thrust_To_V(Thrust, Cd, rho, S):

        return np.sqrt((2 * Thrust) / (rho* S * Cd))

    #----------------------------------------------------------
    # Lift + Drag → Aerodynamic resultant
    @staticmethod
    def Aero_Force_Resultant(L, D):
        """
        Lift + Drag → Aerodynamic resultant
        """
        Res = np.sqrt(L**2 + D**2)
        Theta = np.degrees(np.atan(D / L))
        return Res, Theta
    #----------------------------------------------------------
    @staticmethod
    def V_Stall_In_Climb(W,gam,rho,S,CL_max):
        """
            Parameters:
            W : Weight (N)
            gam ; flight path Angle (deg)
            rho: density
            CL_max : Max lift coefficient
            
            Returns 
            V stall (m/s)
        """
        return np.sqrt((2*W*np.cos(np.deg2rad(gam)))/(rho*S*CL_max))
    #----------------------------------------------------------
    # Drag coefficient from CD0 + k CL²
    @staticmethod
    def Cd_Calc(CD_0, k, Cl):
        return CD_0 + k * (Cl**2)
     #----------------------------------------------------------
    @staticmethod
    def Rate_of_Climb(Fa,mass,v,CD_0,k,gam,rho,S):
        """
            Computes the aircraft rate of climb (ROC).

            Parameters
            ----------
            Fa : float
                Available thrust [N].
            mass : float
                Aircraft mass [kg].
            v : float
                True airspeed [m/s].
            FR : float
                Required drag force (total drag) for CD_0 [N].
            k : float
                Induced drag factor (1/(pi*e*AR)).
            gam : float
                Flight path angle gamma [deg].
            rho : float
                Air density [kg/m^3].
            S : float
                Wing reference area [m^2].

            Returns
            -------
            float
                Rate of climb [m/s].
            float
                Horizontal velocity [m/s].
            float
                Flight path angle [deg].           
            """
        R_C= ((Fa*v)/(mass*9.81)) - v*(((CD_0*0.5*rho*S*v**2)/(mass*9.81)) + (2*k*mass*9.81*(np.cos(np.deg2rad(gam))**2))/(rho*S*(v**2)))
        V_hor = np.sqrt(v**2 - R_C**2)
        gamma = np.rad2deg(np.asin(R_C/v))
        return R_C , V_hor, gamma
    #----------------------------------------------------------
    @staticmethod
    def calculate_jet_range(W,S,TSFC_N,CD_0, k, rho, mass_f):
        """
        Computes jet maximum range using the given parameters.

        Parameters:
            W : Weight of total mass
            TSFC_N (float)    : Thrust-specific fuel consumption (per Newton)
            CD_0 (float)      : Zero-lift drag coefficient
            k (float)         : Induced drag factor
             
            mass_f (float)    : Fuel mass burned

        Returns:
            float: Aircraft range in m
        """
        Mass_t = W/f.g
        V_opt_FL0,_,CL_OPT = f.Optimal_Val_Max_J_range(W, S, rho, CD_0, k)
        return (V_opt_FL0/(9.81*TSFC_N))*(CL_OPT/f.Cd_Calc(CD_0,k,CL_OPT)) *np.log((Mass_t)/(Mass_t-mass_f))
    #----------------------------------------------------------
    @staticmethod
    def calculate_jet_Endurance(W,S,TSFC_N,CD_0, k, rho, mass_f):
        """
        Computes jet maximum Endurance using the given parameters.

        Parameters:
            W : Weight of total mass
            TSFC_N (float)    : Thrust-specific fuel consumption (per Newton)
            CD_0 (float)      : Zero-lift drag coefficient
            k (float)         : Induced drag factor
             
            mass_f (float)    : Fuel mass burned

        Returns:
            float: Aircraft endurance in s
        """
        Mass_t = W/f.g
        _,_,CL_OPT = f.Optimal_Val_Max_J_range(W, S, rho, CD_0, k)
        return (1/(9.81*TSFC_N))*(CL_OPT/f.Cd_Calc(CD_0,k,CL_OPT)) *np.log((Mass_t)/(Mass_t-mass_f))
   
    #----------------------------------------------------------
    # Optimal Values; Jet-Max Range
    @staticmethod
    def Optimal_Val_Max_J_range(W, S, rho, CD_0, k):
        """Optimal Values; Jet-Max Range"""
        v_min = np.sqrt(((W * 2) / (S * rho)) * np.sqrt(3*k/CD_0))
        CL_CD_OPT = 3 / (4 * (3 * k * CD_0**3)**0.25)
        CL_Opt = np.sqrt(CD_0 / (3 * k))
        return v_min, CL_CD_OPT, CL_Opt

    #----------------------------------------------------------
    # Optimal Values Minimal Thrust; Prop-Max Range, Jet-Max Endurance
    @staticmethod
    def Optimal_Val_Min_T(W, S, rho, CD_0, k):
        """
        Optimal Values Minimal Thrust; Prop-Max Range, Jet-Max Endurance
        """
        v_min = np.sqrt((W * 2 * np.sqrt(k)) / (S * rho * np.sqrt(CD_0)))
        CL_CD_OPT = 1 / np.sqrt(4 * k * CD_0)
        CL_Opt = np.sqrt(CD_0 / k)
        return v_min, CL_CD_OPT, CL_Opt

    #----------------------------------------------------------
    # Optimal Values for Minimal Power; Prop – Max Endurance
    @staticmethod
    def Optimal_Val_Min_P(L, S, rho, CD_0, k):
        """Optimal Values for Minimal Power; Prop - Max Endurance"""
        v_min = np.sqrt((L * 2 * np.sqrt(k)) / (S * rho * np.sqrt(3 * CD_0)))
        CL_CD_OPT = 0.25 * (27 / ((k**3) * CD_0))**0.25
        CL_OPT = np.sqrt((3 * CD_0) / k)
        return v_min, CL_CD_OPT, CL_OPT

     #----------------------------------------------------------
    @staticmethod
    def liftoff_Speed(W,rho,S,CL_max_TO):
        """
        Lift of speed:
        Parameters: 
        W :Aircraft weight [Newtons or pounds-force]. 
        rho : float Air density [kg/m³]. 
        S : float Wing reference area [m² ].
        CL_max_TO : float Maximum lift coefficient during takeoff (dimensionless).
        """

        return 1.2 *np.sqrt((2*W)/(rho*S*CL_max_TO))
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
    # CL → TAS, gam == Flight Path Angle in degrees
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
        """Gam is to take into account the flight path angle gamma
             W --> L"""
        return Mass * f.g * np.cos(np.deg2rad(gam))

    #----------------------------------------------------------
    # FR --> FA (Available Thrust)
    @staticmethod
    def FA_Calc(Thrust, Mass, gam):
        return Thrust + Mass * f.g * np.sin(np.deg2rad(gam))
     #----------------------------------------------------------
    # Returns V, intersections between FR and FA SLF
    def V_Int_FA_FR(FA,W,S,CD_0,k,rho):
        """Returns V, intersections between FR and FA SLF"""
        V2 = np.sqrt(((FA * W)/(W *S) + (W/S) * np.sqrt((FA/W)**2 -4 *CD_0*k) ) / (rho*CD_0))
        V1 = np.sqrt(((FA * W)/(W *S) - (W/S) * np.sqrt((FA/W)**2 -4 *CD_0*k) ) / (rho*CD_0))
        return V1,V2    
    #----------------------------------------------------------
    # Unit Conversions
    # lbf <--> N , x = 1 : lbf -> N
    @staticmethod
    def lbf__N(Force, x):
        return Force * (4.4482216153**x)
    # lb <--> kg , x = 1 : lb -> kg
    @staticmethod
    def lb__kg(mass, x):
        return mass * (0.45359237**x)
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
    # knots --> m/s
    @staticmethod
    def knots_To_m_s(v):
        return v * 0.514444
    #----------------------------------------------------------
    # For intersection Between FA and FR, adjust tolerance, to find two values 
    # far from each other 
    @staticmethod
    def Find_Match_FA_FR(Fr, FA, tol):
        return np.where(np.abs(Fr - FA) <= tol)[0]
    @staticmethod
    def get_isa_density(p, T):
    
        return p / (f.R * T)

    @staticmethod
    def get_atmosphere_properties(altitude, units='m'):
        """
        Returns (pressure_Pa, temperature_K, density_kgm3)
        for a given altitude using ISA troposphere + lower stratosphere.

        altitude: input altitude
        units: 'm' (default) or 'ft'
        """
        T_SL = 288.15      # Sea-level temperature (K)
        P_SL = 101325      # Sea-level pressure (Pa)
        L   = 0.0065       # Lapse rate (K/m) (troposphere)
        R   = 287.053      # Gas constant for dry air (J/kg·K)
        g   = 9.80665 
        # Convert to meters
        if units.lower() == 'ft':
            h = altitude * 0.3048
        else:
            h = altitude

        # --- Troposphere (0–11 km) ---
        if h <= 11000:
            T = T_SL - L * h
            p = P_SL * (T / T_SL) ** (g / (R * L))
            rho = f.get_isa_density(p, T)

        # --- Lower Stratosphere (11–20 km) ---
        else:
            # Values at 11 km
            T_11 = T_SL - L * 11000
            p_11 = P_SL * (T_11 / T_SL) ** (g / (R * L))

            # Isothermal layer
            T = T_11
            p = p_11 * np.exp(-g * (h - 11000) / (R * T))
            rho = f.get_isa_density(p, T)

        return p, T, rho
    
    
    @staticmethod
    def ground_roll(W, V_lo, vw, FA, mu_r, phi, Cd, Cl, rho, S):
        """
        Compute the ground-roll distance, average acceleration, and time.

        Parameters
        ----------
        W : float
            Aircraft weight (N).
        V_lo : float
            Liftoff speed (m/s).
        vw : float
            Headwind (+) or tailwind (–) component (m/s).
        FA : float
            Available thrust (N).
        mu_r : float
            Rolling friction coefficient.
        phi : float
            Runway upslope angle (rad). (positive = uphill)
        Cd : float
            Drag coefficient in takeoff configuration.
        Cl : float
            Lift coefficient used for ground roll.
        rho : float
            Air density (kg/m^3).
        S : float
            Wing reference area (m^2).

        Returns
        -------
        sg : float
            Ground roll distance (m).
        ag : float
            Average acceleration during ground roll (m/s²).
        tg : float
            Time to liftoff speed (s).
        """
        q_lo = 0.5 * rho * S * (V_lo / np.sqrt(2))**2

        sg = (W * (V_lo - vw)**2) / (
            2 * 9.81 * ((FA - mu_r * W) - W * phi - (Cd - mu_r * Cl) * q_lo)
        )

        ag = (V_lo - vw)**2 / (2 * sg)
        tg = (V_lo - vw) / ag

        return sg, ag, tg
    @staticmethod
    def rotation_distance(V_lo, t_rot):
        """
        Distance covered during rotation phase.

        Parameters
        ----------
        V_lo : float
            Liftoff speed (m/s).
        t_rot : float
            Rotation duration (s).

        Returns
        -------
        float
            Rotation distance (m).
        """
        return V_lo * t_rot
    @staticmethod
    def delta_CL_transition(V_lo, V_stall, CL_max_TO):
        """
        Compute the increase in lift coefficient during transition from rotation
        to climb.

        Parameters
        ----------
        V_lo : float
            Liftoff speed (m/s).
        V_stall : float
            Stall speed at climb configuration (m/s).
        CL_max_TO : float
            Maximum lift coefficient in takeoff configuration.

        Returns
        -------
        float
            ΔCL during the transition.
        """
        v_ratio = V_lo**2 / V_stall**2
        return 0.5 * (v_ratio - 1) * (CL_max_TO * (1/v_ratio - 0.53) + 0.38)
    @staticmethod
    def transition_arc(W, S, rho, dCL, FA, Fr, V_lo):
        """
        Compute transition arc radius, climb angle, arc distance, height gained,
        and time to complete the transition.

        Parameters
        ----------
        W : float
            Aircraft weight (N).
        S : float
            Wing area (m^2).
        rho : float
            Air density (kg/m^3).
        dCL : float
            Increment in lift coefficient during transition.
        FA : float
            Available thrust (N).
        Fr : float
            Drag force at liftoff speed (N).
        V_lo : float
            Liftoff speed (m/s).

        Returns
        -------
        Radius : float
            Radius of transition arc (m).
        gamma : float
            Climb angle (rad).
        s_trans : float
            Horizontal distance covered in transition (m).
        h_trans : float
            Height gained in transition (m).
        t_trans : float
            Time to complete transition (s).
        """
        Radius = (2 * W) / (S * rho * 9.81 * dCL)
        gamma = np.asin((FA - Fr) / W)

        s_trans = Radius * np.sin(gamma)
        h_trans = Radius - Radius * np.cos(gamma)
        t_trans = (gamma * Radius) / V_lo

        return Radius, gamma, s_trans, h_trans, t_trans
    @staticmethod
    def climb_to_obstacle(h_obstacle, h_transition, gamma, V_lo, tg):
        """
        Compute the climb distance & time needed after transition arc
        to clear a given obstacle height.

        Parameters
        ----------
        h_obstacle : float
            Obstacle height to clear (m).
        h_transition : float
            Height after end of transition (m).
        gamma : float
            Climb angle (rad).
        V_lo : float
            Airspeed during climb (m/s).
        tg : float
            Ground roll time (s), used in FAA/AIAA climb time models.

        Returns
        -------
        s_climb : float
            Additional ground distance to clear obstacle (m).
        t_climb : float
            Time to clear obstacle (s).
        """
        if h_obstacle <= h_transition:
            return 0.0, 0.0

        s_climb = (h_obstacle - h_transition) / (tg * gamma)
        t_climb = s_climb / (np.cos(gamma) * V_lo)

        return s_climb, t_climb


    # ---------------------------------------------------------------------
    # 7. Total Takeoff Distance and Time
    # ---------------------------------------------------------------------
    @staticmethod
    def total_takeoff(sg, sr, s_trans, s_climb, tg, t_rot, t_trans, t_climb):
        """
        Aggregate all distances and times to compute total takeoff performance.

        Parameters
        ----------
        sg : float
            Ground roll distance (m).
        sr : float
            Rotation distance (m).
        s_trans : float
            Transition arc horizontal distance (m).
        s_climb : float
            Additional climb distance to clear obstacle (m).
        tg : float
            Ground roll time (s).
        t_rot : float
            Rotation time (s).
        t_trans : float
            Transition time (s).
        t_climb : float
            Additional climb time (s).

        Returns
        -------
        S_total : float
            Total takeoff distance (m).
        T_total : float
            Total takeoff time (s).
        """
        S_total = sg + sr + s_trans + s_climb
        T_total = tg + t_rot + t_trans + t_climb
        return S_total, T_total