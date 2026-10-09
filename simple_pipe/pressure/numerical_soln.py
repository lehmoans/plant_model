import matplotlib.pyplot as plt
import numpy as np

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

dp_dz = (P_in - P_out)/(L)




r_nodes = np.linspace(0, R, r_num_nodes)
z_nodes = np.linspace(0, L, z_num_nodes)
u_values = [u(r) for r in r_nodes]

plt.plot(u_values,r_nodes, linewidth=2)
plt.plot(u_values,-r_nodes, markersize=2)

plt.ylabel("Radius r (m)")
plt.xlabel("Axial velocity u (m/s)")
plt.title("Analytical Soln: Pipe Poiseuille Velocity Profile")
plt.grid(True)
plt.tight_layout()
plt.show()



