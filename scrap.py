"""
Scratch test: inclined-pipe Poiseuille flow.

Run from the repository root:
    python scrap.py

Assumptions:
- Steady, fully developed, laminar flow of a Newtonian fluid.
- Circular pipe with constant radius and inclination.
- theta is measured from horizontal; positive theta means the pipe rises
  in the positive flow direction.
- P_in and P_out are static pressures at the inlet and outlet.

The gravity-adjusted driving pressure gradient is:
    G = (P_in - P_out) / L - rho * g * sin(theta)

For G > 0, flow is in the positive axial direction. For G < 0, flow
is in the negative axial direction.
"""

import math

import numpy as np
import matplotlib.pyplot as plt


# Pipe and fluid properties (SI units)
P_in = 20_000.0       # Pa
P_out = 10_000.0      # Pa
L = 1.0               # m
R = 0.01              # m
mu = 0.001            # Pa.s, water near room temperature
rho = 1000.0          # kg/m^3
g = 9.81              # m/s^2
theta_deg = 30.0      # degrees above horizontal

# Radial sampling; include the wall at r = R.
n_r = 101
r = np.linspace(0.0, R, n_r)
theta = math.radians(theta_deg)

# Pressure-driven gradient minus the hydrostatic gradient along the pipe.
pressure_gradient = (P_in - P_out) / L
gravity_gradient = rho * g * math.sin(theta)
effective_gradient = pressure_gradient - gravity_gradient

# Signed axial velocity profile, positive along the pipe's chosen axis.
u = effective_gradient * (R**2 - r**2) / (4.0 * mu)

# Analytical volumetric flow rate for a circular pipe.
Q_analytical = math.pi * R**4 * effective_gradient / (8.0 * mu)

# Integrate the sampled profile as a numerical cross-check:
# Q = integral_0^R u(r) * 2*pi*r dr.
Q_numerical = float(np.trapezoid(2.0 * math.pi * r * u, r))
relative_error = abs(Q_numerical - Q_analytical) / max(
    abs(Q_analytical), np.finfo(float).eps
)

print("Inclined-pipe Poiseuille test")
print(f"  Inclination:             {theta_deg:.2f} deg")
print(f"  Pressure gradient:       {pressure_gradient:.3f} Pa/m")
print(f"  Gravity contribution:    {gravity_gradient:.3f} Pa/m")
print(f"  Effective gradient:      {effective_gradient:.3f} Pa/m")
print(f"  Centreline velocity:     {u[0]:.6g} m/s")
print(f"  Wall velocity:           {u[-1]:.6g} m/s")
print(f"  Analytical flow rate:    {Q_analytical:.6g} m^3/s")
print(f"  Integrated flow rate:    {Q_numerical:.6g} m^3/s")
print(f"  Relative integration error: {relative_error:.3e}")

# Basic checks: no-slip at the wall and agreement with analytical flow.
assert np.isclose(u[-1], 0.0, atol=1e-12), "No-slip wall condition failed."
assert np.isclose(Q_numerical, Q_analytical, rtol=1e-3), (
    "Numerical integration does not agree with analytical flow rate."
)

# Plot velocity profile and the direction of flow.
fig, ax = plt.subplots()
ax.plot(u, r, label="Positive radius")
ax.plot(u, -r, label="Negative radius")
ax.axvline(0.0, color="black", linewidth=0.8)
ax.set_xlabel("Axial velocity u (m/s)")
ax.set_ylabel("Radius coordinate (m)")
ax.set_title(f"Inclined pipe velocity profile (theta = {theta_deg:g} deg)")
ax.grid(True)
ax.legend()
fig.tight_layout()
plt.show()
