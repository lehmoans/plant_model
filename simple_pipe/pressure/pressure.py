import matplotlib.pyplot as plt
import numpy as np

"""
velocity field for a steady, fully developed, laminar Newtonian fluid in a pipe. 
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

P_drop_length = (P_in - P_out)/(L)

u = lambda r: (1/(4*mu))*(P_drop_length)*(R**2 - r**2)

r_nodes = np.linspace(0, R, r_num_nodes-1)
z_nodes = np.linspace(0, L, z_num_nodes-1)
u_values = [u(r) for r in r_nodes]

plt.plot(u_values,r_nodes, linewidth=2)
plt.plot(u_values,-r_nodes, markersize=2)

plt.ylabel("Radius r (m)")
plt.xlabel("Axial velocity u (m/s)")
plt.title("Analytical Soln: Pipe Poiseuille Velocity Profile")
plt.grid(True)
plt.tight_layout()
plt.show()



