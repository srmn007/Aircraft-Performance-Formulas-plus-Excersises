import numpy as np

# -----------------------------
# Constants & atmosphere
# -----------------------------
g = 9.80665
rho0 = 1.225  # sea level density

# FL390 (ISA)
rho_FL390 = 0.31641
a_FL390 = 295.07     # speed of sound (m/s)
M = 0.82
V = M * a_FL390      # TAS at FL390

# Dynamic pressure
q_FL390 = 0.5 * rho_FL390 * V**2

# -----------------------------
# Aircraft data
# -----------------------------
# Boeing 777-200LR
S_777 = 427.8
m_777 = 304000
W_777 = m_777 * g

CD0_777 = 0.018
k_777   = 0.045

# F-16
S_f16 = 27.87
m_f16 = 11000
W_f16 = m_f16 * g

CD0_f16 = 0.0169
k_f16   = 0.117

# -----------------------------
# Helper functions
# -----------------------------
def compute_CL(W, rho, V, S):
    return 2 * W / (rho * V**2 * S)

def compute_CD(CD0, k, CL):
    return CD0 + k * CL**2

def aerodynamic_forces(W, S, CD0, k, rho, V):
    q = 0.5 * rho * V**2
    CL = compute_CL(W, rho, V, S)
    CD = compute_CD(CD0, k, CL)
    L = W
    D = q * S * CD
    R = np.sqrt(L**2 + D**2)
    theta = np.degrees(np.arctan(D/L))
    LD = L / D
    return CL, CD, L, D, R, theta, LD

def max_LD(CD0, k):
    CL_opt = np.sqrt(CD0 / k)
    LD_max = 1 / (2 * np.sqrt(CD0 * k))
    return CL_opt, LD_max

def speed_min_drag(W, S, CD0, k, rho, rho0):
    # Minimum-drag TAS
    V_min = ((4 * k * W**2) / (rho**2 * S**2 * CD0))**0.25
    # Convert to EAS
    EAS = V_min * np.sqrt(rho / rho0)
    return V_min, EAS

# -----------------------------
# 1) Aerodynamic forces at FL390, M=0.82
# -----------------------------
print("=== 1) Aerodynamic Forces at FL390 ===")

aircraft = [
    ("Boeing 777-200LR", W_777, S_777, CD0_777, k_777),
    ("F-16",             W_f16, S_f16, CD0_f16, k_f16)
]

for name, W, S, CD0, k in aircraft:
    CL, CD, L, D, R, theta, LD = aerodynamic_forces(W, S, CD0, k, rho_FL390, V)
    print(f"\n{name}:")
    print(f" CL = {CL:.5f}")
    print(f" CD = {CD:.5f}")
    print(f" L  = {L:.2f} N")
    print(f" D  = {D:.2f} N")
    print(f" L/D = {LD:.3f}")
    print(f" Resultant magnitude = {R:.2f} N")
    print(f" Resultant angle θ = {theta:.3f} deg")

# -----------------------------
# 2) Maximum L/D
# -----------------------------
print("\n=== 2) Maximum L/D ===")
for name, W, S, CD0, k in aircraft:
    CL_opt, LD_max = max_LD(CD0, k)
    print(f"{name}: CL_opt = {CL_opt:.5f},  L/D_max = {LD_max:.3f}")

# -----------------------------
# 4) Speeds that minimize drag
# -----------------------------
print("\n=== 4) Minimum-Drag Speeds (FL390 & SL) ===")

for name, W, S, CD0, k in aircraft:
    Vmin_FL390, EAS_FL390 = speed_min_drag(W, S, CD0, k, rho_FL390, rho0)
    Vmin_SL, EAS_SL = speed_min_drag(W, S, CD0, k, rho0, rho0)

    print(f"\n{name}:")
    print(f" V_min TAS at FL390 = {Vmin_FL390:.3f} m/s")
    print(f" EAS at that speed  = {EAS_FL390:.3f} m/s")
    print(f" V_min TAS at SL    = {Vmin_SL:.3f} m/s")