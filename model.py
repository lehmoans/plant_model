import numpy as np
import matplotlib.pyplot as plt
import math

mu = 1e-3
rho = 1000
g = 9.81

P_root_tip = 0
P_leaf_tip = 0

# Nodes: pressure is known at boundaries; otherwise None.
nodes = {
    #above ground
    1:{
        1: {
            "m": {"p": None, "z": 0.0},
        },
        2:{
            "n1": {"p": None, "z": 0.0},
            "n2": {"p": None, "z": 0.0},
        },
        3: {
            "n11": {"p": P_leaf_tip, "z": 0.0},
                        "n12": {"p": P_leaf_tip, "z": 0.0},
                    
                        "n21": {"p": P_leaf_tip, "z": 0.0},
                        "n22": {"p": P_leaf_tip, "z": 0.0},
        },
    },
    #ground
    0: {
         "g": {"p": None, "z": 0.0},
    },
    #below ground
    -1:{
        1:{
             "k":   {"p": None,  "z": None},
        },
        2:{
            "j1":   {"p": None,  "z": None},
            "j2": {"p": None,   "z": None},
        },
        3:{
            "j11": {"p": P_root_tip,   "z": None},
            "j12": {"p": P_root_tip,   "z": None},
            "j21": {"p": P_root_tip,   "z": None},
            "j22": {"p": P_root_tip,   "z": None},
        },
        },
    }

        # Branch: (node_a, node_b, length_m, radius_m, theta_deg)
        # theta is measured from horizontal, positive upward from a to b.
branches = [
            ("n1", "n11",   0.10, 0.0020, 45),
            ("n1", "n12",   0.10, 0.0020, 90),
        
            ("n2", "n21",   0.10, 0.0020, 90),
            ("n2", "n22",   0.10, 0.0020, 135),
        
            ("n2", "m",   0.10, 0.0020, -25),
            ("n1", "m",   0.10, 0.0020, -130),
        
            ("m", "g",   0.10, 0.0020, -90),
        
            ("g", "k",   0.10, 0.0020, 270),
            
            ("k",   "j1",   0.15, 0.0015, 225),
            ("k",   "j2", 0.20, 0.0010, -45),
        
            ("j1",   "j11", 0.1, 0.0010, -135),
            ("j1",   "j12", 0.1, 0.0010, -45),
        
            ("j2",   "j21", 0.1, 0.0010, -135),
            ("j2",   "j22", 0.1, 0.0010, -45),
   
]

########### DRAWING PLANT ################

# 1. Calculate coordinates independently of branch order.
# Seed one node; all other positions follow from branch geometry.
positions = {"g": (0.0, 0.0)}
pending = branches.copy()

while pending:
    progressed = False
    remaining = []

    for a, b, length, radius, angle_deg in pending:
        theta = math.radians(angle_deg)
        dx = length * math.cos(theta)
        dz = length * math.sin(theta)

        if a in positions and b not in positions:
            xa, za = positions[a]
            positions[b] = (xa + dx, za + dz)
            progressed = True

        elif b in positions and a not in positions:
            xb, zb = positions[b]
            positions[a] = (xb - dx, zb - dz)
            progressed = True

        elif a not in positions and b not in positions:
            remaining.append(
                (a, b, length, radius, angle_deg)
            )

    pending = remaining

    if not progressed:
        break

if pending:
    raise ValueError(
        f"Could not position all nodes: {pending}"
    )

# 2. Plot branches, scaling line width by radius.
max_radius = max(branch[3] for branch in branches)

fig, ax = plt.subplots(figsize=(10, 8))

for a, b, length, radius, angle_deg in branches:
    xa, za = positions[a]
    xb, zb = positions[b]

    width = 1.5 + (radius / max_radius) * 5.5

    ax.plot(
        [xa, xb],
        [za, zb],
        linewidth=width,
        solid_capstyle="round",
    )

# 3. Plot and label nodes.
for node, (x, z) in positions.items():
    ax.scatter(x, z, s=28, zorder=3)
    ax.annotate(
        node,
        (x, z),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=9,
    )

ax.set_title("Tree Hydraulic Network")
ax.set_xlabel("Horizontal position (m)")
ax.set_ylabel("Vertical position (m)")
ax.set_aspect("equal", adjustable="datalim")
ax.grid(True, alpha=0.25)

plt.tight_layout()
#plt.show()

# Set elevations from the root along the tree.
# This assumes each node has one parent branch.
nodes["g"]["z"] = 0.0

for a, b, L, r, theta_deg in branches:
    if nodes[a]["z"] is not None:
        nodes[b]["z"] = (
            nodes[a]["z"] + L * np.sin(np.deg2rad(theta_deg))
        )


# Hydraulic resistance of each branch
R = [
    8 * mu * L / (np.pi * r**4)
    for a, b, L, r, theta in branches
]

# Known hydraulic potentials (BCs)
phi = {
    name: node["p"] + rho * g * node["z"]
    for name, node in nodes.items()
    if node["p"] is not None
}

unknown = [name for name in nodes if nodes[name]["p"] is None]
idx = {name: i for i, name in enumerate(unknown)}
n = len(unknown)

A = np.zeros((n, n))
b_vec = np.zeros(n)

# Build conservation equations at each internal junction
for (a, b, L, r, theta), resistance in zip(branches, R):
    G = 1 / resistance

    for node, other in [(a, b), (b, a)]:
        #checks for forward pair (downstream) and backward pair (upstream)
        if node not in idx:
            continue

        i = idx[node]
        A[i, i] += G #diag: G values on nodes that are linked to current unknown node

        if other in idx:
            A[i, idx[other]] -= G #off-diag: G values on nodes that are linked to other unknown nodes
        else:
            b_vec[i] += G * phi[other] #constants

solution = np.linalg.solve(A, b_vec)

for name, value in zip(unknown, solution):
    phi[name] = value

# Recover pressures and calculate flows
print("Node results:")
for name, node in nodes.items():
    p = phi[name] - rho * g * node["z"]
    print(f"{name}: z={node['z']:.4f} m, p={p:.2f} Pa")

Q_values = []
print("\nBranch results:")
for (a, b, L, r, theta), resistance in zip(branches, R):
    Q = (phi[a] - phi[b]) / resistance
    print(f"{a} -> {b}: theta={theta}°, Q={Q:.6e} m^3/s")
    Q_values.append(Q)

linked_nodes
for node in nodes:
    if nodes[node]["p"] is None:
        print(node)
        for a, b, _ , _ , _ in branches:
            if node == a:


        (node_a, node_b, length_m, radius_m, theta_deg)



#plotting results

# Q_values[i] corresponds to branches[i]
Q_values = np.asarray(Q_values, dtype=float)

fig, ax = plt.subplots(figsize=(10, 8))

max_Q = np.max(np.abs(Q_values))
if max_Q == 0:
    max_Q = 1.0  # avoid division by zero

for i, (a, b, length, radius, theta_deg) in enumerate(branches):
    xa, za = positions[a]
    xb, zb = positions[b]
    Q = Q_values[i]

    # Thicker lines indicate larger flow magnitude
    width = 1.0 + 5.0 * abs(Q) / max_Q

    ax.plot(
        [xa, xb], [za, zb],
        linewidth=width,
        color="tab:blue" if Q >= 0 else "tab:orange",
        solid_capstyle="round",
    )

    # Arrow points in the actual flow direction
    dx = xb - xa
    dz = zb - za

    if Q < 0:
        xa, za = xb, zb
        dx, dz = -dx, -dz

    # Place arrow near the middle of the branch
    fraction = 0.55
    ax.annotate(
        "",
        xy=(xa + fraction * dx, za + fraction * dz),
        xytext=(xa + (fraction - 0.12) * dx,
                za + (fraction - 0.12) * dz),
        arrowprops=dict(
            arrowstyle="->",
            color="black",
            lw=1.2,
        ),
    )

    # Label flow rate near branch midpoint
    xm = (positions[a][0] + positions[b][0]) / 2
    zm = (positions[a][1] + positions[b][1]) / 2

    ax.text(
        xm, zm,
        f"{Q:.2e}",
        fontsize=8,
        ha="center",
        va="center",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
    )

# Plot and label nodes
for node, (x, z) in positions.items():
    ax.scatter(x, z, s=20, color="black", zorder=3)
    ax.annotate(
        node, (x, z),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=8,
    )

ax.set_title("Flow Rates Through the Plant Network")
ax.set_xlabel("Horizontal position (m)")
ax.set_ylabel("Vertical position (m)")
ax.set_aspect("equal", adjustable="datalim")
ax.grid(True, alpha=0.25)

plt.tight_layout()
plt.show()

"""
- adjust the code for the new node dictionary

- mass conservtion at each node
- deriving min energy needed to drive water upward flow
    - energy dissipation due to resistence
        - find total resistence
    - energy lost to overcoming gravitational effects
        - find total gravitation effects-
"""
