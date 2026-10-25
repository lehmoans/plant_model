"""Steady laminar flow through a branched plant hydraulic network.

Each branch is modeled as a Poiseuille resistance. The solved variable is
hydraulic potential phi = p + rho*g*z, and positive Q follows the branch's
declared direction (node_a -> node_b).
"""
from __future__ import annotations

import math
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

MU = 1e-3       # Dynamic viscosity [Pa s]
RHO = 1000.0    # Fluid density [kg/m^3]
G = 9.81        # Gravitational acceleration [m/s^2]

# BCs
P_ROOT_TIP = 0.0
P_LEAF_TIP = 0.0

#above ground, where L = 1
CHANNEL_R1 = 2
CHANNEL_R2 = 1
CHANNEL_R3 = 0.5

#below ground where R = 0.5
CHANNEL_LEN1 = 2
CHANNEL_LEN2 = 1
CHANNEL_LEN3  = 0.5


# Nodes are grouped by region/order. A known pressure marks a boundary node;
# p=None marks an internal junction whose potential is solved for.
NODE_TREE = {
    # Above ground
    1: {
        1: {"m": {"p": None, "z": None}},
        2: {
            "n1": {"p": None, "z": None},
            "n2": {"p": None, "z": None},
        },
        3: {
            "n11": {"p": P_LEAF_TIP, "z": None},
            "n12": {"p": P_LEAF_TIP, "z": None},
            "n21": {"p": P_LEAF_TIP, "z": None},
            "n22": {"p": P_LEAF_TIP, "z": None},
        },
    },
    # Ground
    0: {"g": {"p": None, "z": None}},
    # Below ground
    -1: {
        1: {"k": {"p": None, "z": None}},
        2: {
            "j1": {"p": None, "z": None},
            "j2": {"p": None, "z": None},
        },
        3: {
            "j11": {"p": P_ROOT_TIP, "z": None},
            "j12": {"p": P_ROOT_TIP, "z": None},
            "j21": {"p": P_ROOT_TIP, "z": None},
            "j22": {"p": P_ROOT_TIP, "z": None},
        },
    },
}

# Branch: (node_a, node_b, length_m, radius_m, theta_deg)
# theta is measured from horizontal, positive upward from a to b.
BRANCHES = [
    ("n1", "n11", CHANNEL_LEN1, CHANNEL_R3, 45),
    ("n1", "n12", CHANNEL_LEN1, CHANNEL_R3, 90),
    ("n2", "n21", CHANNEL_LEN1, CHANNEL_R3, 90),
    ("n2", "n22", CHANNEL_LEN1, CHANNEL_R3, 135),
    ("n2", "m", CHANNEL_LEN1, CHANNEL_R2, -25),
    ("n1", "m", CHANNEL_LEN1, CHANNEL_R2, -155),
    ("m", "g", CHANNEL_LEN1, CHANNEL_R1, -90),
    ("g", "k", CHANNEL_LEN1, CHANNEL_R3, 270),
    ("k", "j1", CHANNEL_LEN2, CHANNEL_R3, -145),
    ("k", "j2", CHANNEL_LEN2, CHANNEL_R3, -35),
    ("j1", "j11", CHANNEL_LEN3, CHANNEL_R3, -135),
    ("j1", "j12", CHANNEL_LEN3, CHANNEL_R3, -45),
    ("j2", "j21", CHANNEL_LEN3, CHANNEL_R3, -135),
    ("j2", "j22", CHANNEL_LEN3, CHANNEL_R3, -45),
]




def flatten_nodes(node_tree: dict[Any, Any]) -> dict[str, dict[str, float | None]]:
    """Flatten the nested node tree to {node_name: {p, z}}."""
    flattened: dict[str, dict[str, float | None]] = {}

    def visit(branch: dict[Any, Any]) -> None:
        for key, value in branch.items():
            if not isinstance(value, dict):
                continue
            if "p" in value and "z" in value:
                name = str(key)
                if name in flattened:
                    raise ValueError(f"Duplicate node name: {name}")
                flattened[name] = {"p": value["p"], "z": value["z"]}
            else:
                visit(value)

    visit(node_tree)
    return flattened


def calculate_positions(
    branches: list[tuple[str, str, float, float, float]],
    seed_node: str = "g",
    seed_position: tuple[float, float] = (0.0, 0.0),
) -> dict[str, tuple[float, float]]:
    """Calculate (x, z) coordinates from branch lengths and angles."""
    positions = {seed_node: seed_position}
    pending = list(branches)

    while pending:
        progressed = False
        remaining = []

        for a, b, length, _radius, angle_deg in pending:
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
            elif a in positions and b in positions:
                # This edge has already been geometrically determined.
                continue
            else:
                remaining.append((a, b, length, _radius, angle_deg))

        pending = remaining
        if not progressed and pending:
            raise ValueError(
                "Could not position all nodes. Check that every branch is "
                f"connected to seed node {seed_node!r}: {pending}"
            )

    return positions


def solve_network(
    node_tree: dict[Any, Any] = NODE_TREE,
    branches: list[tuple[str, str, float, float, float]] = BRANCHES,
    mu: float = MU,
    rho: float = RHO,
    gravity: float = G,
) -> dict[str, Any]:
    """Solve node hydraulic potentials and branch flow rates."""
    node_data = flatten_nodes(node_tree)
    positions = calculate_positions(branches)

    branch_nodes = {name for a, b, *_ in branches for name in (a, b)}
    missing = branch_nodes - set(node_data)
    if missing:
        raise ValueError(f"Branches reference undefined nodes: {sorted(missing)}")

    unconnected = set(node_data) - branch_nodes
    if unconnected:
        raise ValueError(f"Nodes are not connected to any branch: {sorted(unconnected)}")

    # Use geometry-derived elevations consistently in the gravity term.
    for name, (_x, z) in positions.items():
        node_data[name]["z"] = z

    resistance = np.asarray(
        [8.0 * mu * length / (math.pi * radius**4)
         for _a, _b, length, radius, _theta in branches],
        dtype=float,
    )
    if np.any(resistance <= 0.0) or not np.all(np.isfinite(resistance)):
        raise ValueError("Every branch must have positive, finite resistance.")

    phi: dict[str, float] = {
        name: float(node["p"]) + rho * gravity * float(node["z"])
        for name, node in node_data.items()
        if node["p"] is not None
    }
    unknown = [name for name, node in node_data.items() if node["p"] is None]
    index = {name: i for i, name in enumerate(unknown)}

    matrix = np.zeros((len(unknown), len(unknown)), dtype=float)
    rhs = np.zeros(len(unknown), dtype=float)

    # Kirchhoff-style flow balance: sum of signed branch flows = 0
    # at each unknown/internal junction.
    for (a, b, _length, _radius, _theta), r_branch in zip(branches, resistance):
        conductance = 1.0 / r_branch
        for node, other in ((a, b), (b, a)):
            if node not in index:
                continue
            i = index[node]
            matrix[i, i] += conductance
            if other in index:
                matrix[i, index[other]] -= conductance
            else:
                rhs[i] += conductance * phi[other]

    if unknown:
        try:
            solution = np.linalg.solve(matrix, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "Could not solve the network. Check boundary conditions and connectivity."
            ) from exc
        phi.update({name: float(value) for name, value in zip(unknown, solution)})

    pressures = {
        name: phi[name] - rho * gravity * float(node["z"])
        for name, node in node_data.items()
    }
    flows = np.asarray(
        [(phi[a] - phi[b]) / r_branch
         for (a, b, _length, _radius, _theta), r_branch in zip(branches, resistance)],
        dtype=float,
    )

    return {
        "nodes": node_data,
        "positions": positions,
        "phi": phi,
        "pressures": pressures,
        "resistances": resistance,
        "matrix": matrix.copy(),
        "rhs": rhs.copy(),
        "unknowns": list(unknown),
        "flows": flows,
        "branches": branches,
    }


def check_mass_conservation(
    node_data: dict[str, dict[str, float | None]],
    branches: list[tuple[str, str, float, float, float]],
    flows: np.ndarray,
    rtol: float = 1e-8,
    atol: float = 1e-12,
) -> list[dict[str, Any]]:
    """Return signed flow-balance residuals; positive net means net outflow."""
    net_outflow = {name: 0.0 for name in node_data}
    sum_abs_flow = {name: 0.0 for name in node_data}

    for (a, b, *_), q in zip(branches, flows):
        # q > 0 means flow follows a -> b.
        net_outflow[a] += float(q)
        net_outflow[b] -= float(q)
        sum_abs_flow[a] += abs(float(q))
        sum_abs_flow[b] += abs(float(q))

    report = []
    for name, node in node_data.items():
        residual = net_outflow[name]
        scale = sum_abs_flow[name]
        is_internal = node["p"] is None
        tolerance = atol + rtol * scale
        report.append({
            "node": name,
            "net_outflow": residual,
            "relative_residual": abs(residual) / scale if scale > 0.0 else 0.0,
            "is_internal": is_internal,
            "passed": abs(residual) <= tolerance if is_internal else None,
            "tolerance": tolerance,
        })
    return report


def flow_direction_color(
    a: str,
    b: str,
    q: float,
    positions: dict[str, tuple[float, float]],
    flow_tol: float = 1e-15,
    vertical_tol: float = 1e-12,
) -> str:
    """Green for upward flow, red for downward flow, grey for zero/horizontal."""
    if abs(q) <= flow_tol:
        return "0.55"

    source, target = (a, b) if q > 0.0 else (b, a)
    dz = positions[target][1] - positions[source][1]
    if dz > vertical_tol:
        return "tab:green"
    if dz < -vertical_tol:
        return "tab:red"
    return "0.55"


def plot_flow_results(results: dict[str, Any]) -> None:
    """Plot the network with branch colour showing actual vertical flow direction."""
    branches = results["branches"]
    positions = results["positions"]
    flows = results["flows"]

    fig, ax = plt.subplots(figsize=(10, 8))
    max_flow = float(np.max(np.abs(flows))) if len(flows) else 0.0
    if max_flow == 0.0:
        max_flow = 1.0

    for (a, b, _length, _radius, _theta), q in zip(branches, flows):
        xa, za = positions[a]
        xb, zb = positions[b]
        colour = flow_direction_color(a, b, float(q), positions)
        width = 1.0 + 5.0 * abs(float(q)) / max_flow
        ax.plot([xa, xb], [za, zb], linewidth=width, color=colour,
                solid_capstyle="round", zorder=1)

        # Arrow follows the actual flow, not necessarily the listed branch direction.
        if q >= 0.0:
            source, target = (xa, za), (xb, zb)
        else:
            source, target = (xb, zb), (xa, za)
        fraction = 0.55
        back_fraction = 0.43
        ax.annotate(
            "",
            xy=(source[0] + fraction * (target[0] - source[0]),
                source[1] + fraction * (target[1] - source[1])),
            xytext=(source[0] + back_fraction * (target[0] - source[0]),
                    source[1] + back_fraction * (target[1] - source[1])),
            arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.1},
            zorder=2,
        )

        ax.text(
            (xa + xb) / 2.0, (za + zb) / 2.0, f"{q:.2e}",
            fontsize=8, ha="center", va="center",
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
            zorder=3,
        )

    for name, (x, z) in positions.items():
        ax.scatter(x, z, s=20, color="black", zorder=4)
        ax.annotate(name, (x, z), xytext=(4, 4),
                    textcoords="offset points", fontsize=8)

    ax.legend(handles=[
        Patch(color="tab:green", label="Upward flow"),
        Patch(color="tab:red", label="Downward flow"),
        Patch(color="0.55", label="Horizontal / near-zero flow"),
    ], loc="best")
    ax.set_title("Flow Rates Through the Plant Network")
    ax.set_xlabel("Horizontal position (m)")
    ax.set_ylabel("Vertical position (m)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.show()



def plot_hydraulic_resistance_network(results: dict[str, Any]) -> None:
    """Plot the network with branch colour and width representing resistance."""
    branches = results["branches"]
    positions = results["positions"]
    values = np.asarray(results["resistances"], dtype=float)
    if values.size == 0 or np.any(values <= 0.0):
        raise ValueError("All branches must have positive hydraulic resistance.")

    logs = np.log10(values)
    low, high = float(logs.min()), float(logs.max())
    norm = plt.Normalize(low, high if high > low else low + 1.0)
    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(10, 8))

    for (a, b, *_), value, log_value in zip(branches, values, logs):
        xa, za = positions[a]
        xb, zb = positions[b]
        fraction = 0.0 if high == low else (log_value - low) / (high - low)
        ax.plot([xa, xb], [za, zb], color=cmap(norm(log_value)),
                linewidth=1.2 + 3.0 * fraction, solid_capstyle="round", zorder=1)

    top_indices = np.argsort(values)[-min(5, len(values)):]
    for idx in top_indices:
        a, b, *_ = branches[int(idx)]
        xa, za = positions[a]
        xb, zb = positions[b]
        ax.annotate(f"{a}–{b}", ((xa + xb) / 2, (za + zb) / 2),
                    xytext=(3, 5), textcoords="offset points", fontsize=8,
                    bbox={"boxstyle": "round,pad=0.2", "fc": "white",
                          "alpha": 0.85, "ec": "none"})

    _draw_nodes(ax, positions)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    fig.colorbar(scalar, ax=ax, pad=0.02).set_label(
        r"Hydraulic resistance $R_i$ (Pa s/m$^3$), log scale"
    )
    ax.set_title("Plant Network — Hydraulic Resistance by Branch")
    ax.set_xlabel("Horizontal position (m)")
    ax.set_ylabel("Vertical position (m)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()


def plot_hydraulic_resistance_ranking(results: dict[str, Any]) -> None:
    """Rank branches from lowest to highest hydraulic resistance."""
    branches = results["branches"]
    values = np.asarray(results["resistances"], dtype=float)
    labels = [f"{a} → {b}" for a, b, *_ in branches]
    order = np.argsort(values)
    values = values[order]
    labels = [labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, max(5.0, 0.34 * len(labels))))
    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=plt.get_cmap("viridis")(
        np.linspace(0.15, 0.9, len(labels))))
    ax.set_yticks(y, labels)
    ax.set_xlabel(r"Hydraulic resistance $R_i$ (Pa s/m$^3$)")
    ax.set_title("Branch Hydraulic Resistance — Lowest to Highest")
    ax.set_xscale("log")
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    for bar, value in zip(bars, values):
        ax.annotate(f"{value:.2e}", (value, bar.get_y() + bar.get_height()/2),
                    xytext=(5, 0), textcoords="offset points",
                    va="center", fontsize=8)
    fig.tight_layout()



def main() -> None:
    results = solve_network()
    nodes = results["nodes"]
    branches = results["branches"]
    flows = results["flows"]

    print("Node results:")
    for name in nodes:
        print(
            f"{name}: z={nodes[name]['z']:.4f} m, "
            f"p={results['pressures'][name]:.6g} Pa, "
            f"phi={results['phi'][name]:.6g} Pa"
        )

    print("\nBranch flow results (positive follows listed a -> b):")
    for (a, b, _length, _radius, theta), q in zip(branches, flows):
        print(f"{a} -> {b}: theta={theta:g} deg, Q={q:.6e} m^3/s")

    print("\nMass conservation check (internal junctions only):")
    report = check_mass_conservation(nodes, branches, flows)
    failed = []
    for item in report:
        if item["is_internal"]:
            status = "PASS" if item["passed"] else "FAIL"
            print(
                f"{item['node']}: {status}, net outflow={item['net_outflow']:.3e} "
                f"m^3/s, relative residual={item['relative_residual']:.3e}"
            )
            if not item["passed"]:
                failed.append(item["node"])

    boundary_net = sum(
        item["net_outflow"] for item in report if not item["is_internal"]
    )
    print(f"Boundary flow imbalance: {boundary_net:.3e} m^3/s")
    if failed:
        print(f"WARNING: mass conservation failed at: {', '.join(failed)}")
    else:
        print("All internal junctions satisfy mass conservation within tolerance.")

    plot_flow_results(results)
    plot_hydraulic_resistance_network(results)
    plot_hydraulic_resistance_ranking(results)
    plt.show()


if __name__ == "__main__":
    main()
