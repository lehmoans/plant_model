import matplotlib.pyplot as plt
import numpy as np

"""
velocity field for a steady, fully developed, laminar Newtonian fluid in a pipe. 
 - pressure driven flow with gravity effects given by theta (angle of pipe with respect to horizontal)
"""

r_num_nodes = 10
z_num_nodes= 10

#qtys
P_in = 10 #Pa
P_out = 1 #Pa
L = 1 #m
R = 0.1 #m
dz = L / z_num_nodes
dr = R / r_num_nodes

mu = 0.001 #Pa·s

density = 2 #kg/m^3
g = 9.81 #m/s^2
theta = 0.5 #degrees
theta_rad = np.radians(theta)

P_drop_length = (P_in - P_out)/(L)
z = L
u = lambda r: (1/(4*mu))*(P_drop_length - density*g*np.sin(theta_rad) )*(R**2 - r**2)
u_p = lambda r: (1/(4*mu))*(P_drop_length )*(R**2 - r**2)

r_nodes = np.linspace(0, R, r_num_nodes-1)
z_nodes = np.linspace(0, L, z_num_nodes-1)
u_values = [u(r) for r in r_nodes]
u_p_values = [u_p(r) for r in r_nodes]

""" 
abs_diff = [abs(u_values[i]-u_p_values[i]) for i in range(len(u_values))]
plt.plot(abs_diff,r_nodes, linewidth=2)
plt.show() 
"""

plt.plot(u_values,r_nodes, linewidth=2, label = "Pressure + Gravity")
plt.plot(u_values,-r_nodes, markersize=2)

plt.plot(u_p_values,r_nodes, linewidth=2, label = "Pressure Only")
plt.plot(u_p_values,-r_nodes, markersize=2)

plt.ylabel("Radius r (m)")
plt.xlabel("Axial velocity u (m/s)")
plt.title("Analytical Soln: Pipe Poiseuille Velocity Profile with gravity effects")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()



