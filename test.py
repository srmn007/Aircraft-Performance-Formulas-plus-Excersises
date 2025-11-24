import numpy as np
import matplotlib.pyplot as plt

# Constants and given parameters
g = 9.81  # m/s²
b_prop = 11.4  # m
S_prop = 19.2  # m²
mass_prop = 2359  # kg
W_prop = mass_prop * g  # N

engine_power_each = 186000  # W (186 kW)
num_engines = 2
P_avail_sea_level = engine_power_each * num_engines  # total engine power W
eta_p = 0.83  # propeller efficiency

CD0_prop = 0.025
e_prop = 0.9
CL_max_prop = 1.16
AR_prop = b_prop**2 / S_prop
k_prop = 1 / (np.pi * AR_prop * e_prop)

# ISA Standard atmosphere densities (approximate)
rho_SL = 1.225  # sea level, kg/m^3
rho_10k = 1.112  # approx 10,000 ft
rho_20k = 0.653  # approx 20,000 ft

altitudes = [rho_SL, rho_10k, rho_20k]
altitude_names = ['Sea Level (0 ft)', '10,000 ft', '20,000 ft']

# Velocity range (horizontal velocity), m/s
V = np.arange(5, 110, 5)

def thrust_available(P_avail, V):
    # Thrust = Power / Velocity, avoid division by zero
    return np.where(V > 0, P_avail / V, 0)

def lift_coefficient(W, rho, V):
    return 2 * W / (rho * V**2 * S_prop)

def drag_coefficient(CD0, k, CL):
    return CD0 + k * CL**2

def drag_force(rho, V, S, CD):
    return 0.5 * rho * V**2 * S * CD

def rate_of_climb(Thrust, Drag, V, W):
    # Excess power available = (Thrust - Drag) * Velocity
    excess_power = (Thrust - Drag) * V
    # Rate of climb = Excess power / Weight
    ROC = excess_power / W
    return ROC

# Store results to plot
results = {}

for i, rho in enumerate(altitudes):
    P_avail = P_avail_sea_level * eta_p  # max power available adjusted by propeller efficiency
    thrust = thrust_available(P_avail, V)
    
    CL = lift_coefficient(W_prop, rho, V)
    # To avoid unrealistic lift values near stall, clip CL to max
    CL = np.clip(CL, 0, CL_max_prop)
    
    CD = drag_coefficient(CD0_prop, k_prop, CL)
    D = drag_force(rho, V, S_prop, CD)
    
    ROC = rate_of_climb(thrust, D, V, W_prop)
    
    # clip ROC to zero for negative excess power (no climb possible)
    ROC = np.clip(ROC, 0, None)
    
    results[altitude_names[i]] = {'V': V, 'ROC': ROC}

# Plot Horizontal Speed vs Rate of Climb
plt.figure(figsize=(10, 6))

for alt_name, data in results.items():
    plt.plot(data['V'], data['ROC'], label=alt_name)

plt.xlabel('Horizontal Speed (m/s)')
plt.ylabel('Rate of Climb (m/s)')
plt.title('Rate of Climb vs Horizontal Speed for Propeller Aircraft at Different Altitudes')
plt.legend()
plt.grid(True)
plt.show()
