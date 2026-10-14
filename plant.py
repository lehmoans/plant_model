from plant_model.root import Root
mu = 0.001  # Viscosity (Pa·s)
density = 1000  # Density (kg/m³)


rootA = Root(R=0.1, L=2, mu=mu, density=density, theta=60)
rootB = Root(R=0.1, L=2, mu=mu, density=density, theta=80)
rootC = Root(R=0.1, L=2, mu=mu, density=density, theta=60)

#junction
rootC_Q = rootA.get_flow_rate(P_in=10, P_out=1) + rootB.get_flow_rate(P_in=10, P_out=1)
