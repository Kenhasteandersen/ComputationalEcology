import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#
# Define the default set of parameters as a class:
#
class Parameters:
    def __init__(self, D=3, u=1, depth=100, dz=1):
        self.D = D           # Diffusivity (m2/day)
        self.u = u           # Settling velocity (m/day)
        self.depth = depth   # Depth (meter)
        self.dz = dz         # Grid spacing (m)
#
# Run the model with a given set of parameters:
#
def advection_diffusion(param):
    z = np.arange(param.dz / 2, param.depth, param.dz)
    param.nGrid = len(z)  # No. of grid cells

    # Initialization:
    P0 = np.exp(-((z - param.depth / 2)**2) / 5)

    # Run model
    sol = solve_ivp(pmodel_deriv, [0, 200], P0, args=(param,), method='RK45')
    t = sol.t
    P = sol.y

    return t, z, P
#
# The right-hand-side of the equations in the model:
#
def pmodel_deriv(t, P, param):
    # Advective fluxes
    Jadv = np.zeros(param.nGrid + 1)
    ix = np.arange(1, param.nGrid)
    Jadv[ix] = param.u * P[ix - 1]
    Jadv[0] = 0  # No input from the surface
    Jadv[param.nGrid] = 0  # Closed bottom

    # Diffusive fluxes:
    Jdiff = np.zeros(param.nGrid + 1)
    Jdiff[ix] = -param.D * (P[ix] - P[ix - 1]) / param.dz
    Jdiff[0] = 0  # No flux at the surface...
    Jdiff[param.nGrid] = 0  # ...or the bottom

    # Rate-of-change due to advection and diffusion:
    J = Jadv + Jdiff
    dPdt = -(J[1:param.nGrid + 1] - J[:param.nGrid]) / param.dz

    return dPdt
    
#
# Plotting a solution as a surface plot
#
def plot_solution_matrix(t, z, P):
    plt.figure(0)
    plt.imshow(P, extent=[t[0], t[-1], z[-1], z[0]], aspect='auto', cmap='viridis')
    plt.colorbar(label='Plankton concentration')
    plt.xlabel('Time (days)')
    plt.ylabel('Depth (m)')
    plt.title('Advection-Diffusion Solution Matrix')
#
# Do a grid sensitivity
#
def gridsensitivity(param, dz_values):
    plt.figure(1)
    for dz in dz_values:
        param.dz = dz
        t, z, P = advection_diffusion(param)
        plt.plot(P[:, -1], -z, 'o-', label=f'dz={dz} m')

    plt.legend()
    plt.xlabel('Concentration (X/m^3)')
    plt.ylabel('Depth (m)')
    plt.title('Plankton profiles')

# Parameters
param = Parameters()

# Run an example
t, z, P = advection_diffusion(param)
plot_solution_matrix(t, z, P)

# Doing grid sensitivity
dz_values = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
gridsensitivity(param, dz_values)

plt.show()
