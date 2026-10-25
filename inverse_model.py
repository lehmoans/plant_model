"""Inverse flow / energy model using the network defined in model.py.

Root-tip pressure is fixed at P_ROOT_TIP. A minimum outflow Q_min is prescribed
at each leaf tip, and the leaf-tip pressures are solved (not fixed). Internal
junctions satisfy mass conservation. The required steady mechanical power is
computed as viscous dissipation plus the net rate of gravitational potential
energy increase.

Numerical stability: solve for hydraulic potential relative to the common root
potential, rather than subtracting large absolute potentials to obtain tiny
flows. This changes only the reference level, not the physical solution.

Important boundary-condition note:
A passive resistance network cannot independently prescribe both pressure and
flow at every leaf. This model prescribes leaf flows and solves leaf pressures.
If root and leaf pressures must both be fixed to P0, an explicit pump/active
source model is needed.

Run from the repository root:
    python inverse_model.py

Q_min is in m^3/s per leaf tip. Geometry and branch properties come from model.py.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from model import (
    BRANCHES,
    G,
    MU,
    NODE_TREE,
    P_ROOT_TIP,
    RHO,
    calculate_positions,
    flatten_nodes,
)

#deriving root/leaf tips as per NODE_TREE structure
for region in NODE_TREE:
    if region == 1:
        max_branch_iter = max(k for k,v in NODE_TREE[region].items())
        LEAF_TIPS = [k for k,v in NODE_TREE[region][max_branch_iter].items()]
    elif region == -1:
        max_branch_iter = max(k for k,v in NODE_TREE[region].items())
        ROOT_TIPS = [k for k,v in NODE_TREE[region][max_branch_iter].items()]


# Minimum required flow at EACH leaf tip [m^3/s].
# Change this value to explore different uptake requirements.
Q_MIN = 1.0e-9

# Plot styling. Line widths are in Matplotlib points.
MIN_LINE_WIDTH = 1.2
FLOW_WIDTH_PER_QMIN = MIN_LINE_WIDTH
FRICTION_MIN_WIDTH = MIN_LINE_WIDTH
FLOW_TOL = 1.0e-18
VERTICAL_TOL = 1.0e-12


def solve_for_minimum_leaf_flows(
    q_min: float | Mapping[str, float] = Q_MIN,
) -> dict[str, Any]:
    """Solve hydraulic potentials and branch flows for prescribed leaf demands.

    Root tips remain pressure boundaries at P_ROOT_TIP, as specified in
    model.py. Each leaf tip is a flow sink: its net inflow must equal Q_min.
    If q_min is a mapping, provide an individual demand for each leaf tip.
    """
    node_data = flatten_nodes(NODE_TREE)
    positions = calculate_positions(BRANCHES, seed_node="g", seed_position=(0.0, 0.0))

    missing = (set(ROOT_TIPS) | set(LEAF_TIPS)) - set(node_data)
    if missing:
        raise ValueError(f"Configured root/leaf tips missing from NODE_TREE: {sorted(missing)}")

    demands = (
        {tip: float(q_min) for tip in LEAF_TIPS}
        if isinstance(q_min, (int, float))
        else {tip: float(q_min[tip]) for tip in LEAF_TIPS}
    )
    if any(value < 0.0 for value in demands.values()):
        raise ValueError("Q_min demands must be non-negative.")

    for name, (_x, z) in positions.items():
        node_data[name]["z"] = z

    # Fix root potentials to zero in a shifted reference frame. All root tips
    # share the same elevation in this geometry, so their physical P_ROOT_TIP
    # gives one common reference potential. Solving relative potentials avoids
    # cancellation when tiny flows are calculated from large absolute values.
    root_elevations = [float(node_data[name]["z"]) for name in ROOT_TIPS]
    if not np.allclose(root_elevations, root_elevations[0], rtol=0.0, atol=VERTICAL_TOL):
        raise ValueError(
            "Root tips do not share one elevation; use individual fixed root "
            "potentials and a numerically stable reference shift."
        )
    reference_phi = P_ROOT_TIP + RHO * G * root_elevations[0]
    fixed_relative_phi = {name: 0.0 for name in ROOT_TIPS}
    unknowns = [name for name in node_data if name not in fixed_relative_phi]
    index = {name: i for i, name in enumerate(unknowns)}

    branch_resistances = []
    for a, b, length, radius, _theta in BRANCHES:
        if length <= 0.0 or radius <= 0.0:
            raise ValueError(f"Invalid dimensions for branch {a!r}->{b!r}.")
        branch_resistances.append(8.0 * MU * length / (math.pi * radius**4))

    matrix = np.zeros((len(unknowns), len(unknowns)), dtype=float)
    rhs = np.zeros(len(unknowns), dtype=float)

    # Nodal balance in hydraulic-potential form:
    # sum_j (phi_i - phi_j)/R_ij = -Q_demand at a leaf sink, 0 otherwise.
    # Root relative potentials are zero, so their RHS contribution is zero.
    for (a, b, *_), resistance in zip(BRANCHES, branch_resistances):
        conductance = 1.0 / resistance
        for node, other in ((a, b), (b, a)):
            if node not in index:
                continue
            i = index[node]
            matrix[i, i] += conductance
            if other in index:
                matrix[i, index[other]] -= conductance
            else:
                rhs[i] += conductance * fixed_relative_phi[other]

    for tip, demand in demands.items():
        rhs[index[tip]] -= demand

    try:
        solution = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Could not solve the flow-demand problem. Check network connectivity "
            "and ensure each leaf is hydraulically connected to a fixed-pressure root."
        ) from exc

    relative_phi = dict(fixed_relative_phi)
    relative_phi.update({name: float(value) for name, value in zip(unknowns, solution)})
    # Restore the physical hydraulic-potential reference only after solving.
    phi = {name: value + reference_phi for name, value in relative_phi.items()}
    pressures = {
        name: phi[name] - RHO * G * float(node_data[name]["z"])
        for name in node_data
    }

    branch_results = []
    total_viscous_power = 0.0
    total_gravity_power = 0.0
    for (a, b, length, radius, theta), resistance in zip(BRANCHES, branch_resistances):
        # Positive Q follows the declared BRANCHES direction a -> b.
        q = (relative_phi[a] - relative_phi[b]) / resistance
        delta_z = positions[b][1] - positions[a][1]
        viscous_power = resistance * q**2
        gravity_power = RHO * G * q * delta_z
        branch_results.append({
            "start": a,
            "end": b,
            "length": length,
            "radius": radius,
            "angle_deg": theta,
            "flow_rate": float(q),
            "resistance": resistance,
            "pressure_drop_viscous": resistance * q,
            "delta_z": delta_z,
            "viscous_power": viscous_power,
            "gravity_power": gravity_power,
        })
        total_viscous_power += viscous_power
        total_gravity_power += gravity_power

    total_root_inflow = 0.0
    for item in branch_results:
        a, b, q = item["start"], item["end"], item["flow_rate"]
        if a in ROOT_TIPS:
            total_root_inflow += q
        if b in ROOT_TIPS:
            total_root_inflow -= q

    total_leaf_outflow = sum(demands.values())
    total_power = total_viscous_power + total_gravity_power
    mass_balance_error = total_root_inflow - total_leaf_outflow
    mass_balance_scale = max(abs(total_root_inflow), abs(total_leaf_outflow), 1e-30)
    mass_balance_relative_error = abs(mass_balance_error) / mass_balance_scale

    return {
        "q_min_by_leaf": demands,
        "root_tips": ROOT_TIPS,
        "leaf_tips": LEAF_TIPS,
        "root_inflow": total_root_inflow,
        "total_leaf_demand": total_leaf_outflow,
        "mass_balance_error": mass_balance_error,
        "mass_balance_relative_error": mass_balance_relative_error,
        "nodes": node_data,
        "positions": positions,
        "hydraulic_potential": phi,
        "relative_hydraulic_potential": relative_phi,
        "reference_hydraulic_potential": reference_phi,
        "pressures": pressures,
        "matrix": matrix,
        "rhs": rhs,
        "unknowns": unknowns,
        "branches": branch_results,
        "total_viscous_power": total_viscous_power,
        "total_gravity_power": total_gravity_power,
        "net_mechanical_power": total_power,
    }


def _flow_colour(branch: dict[str, Any], positions: dict[str, tuple[float, float]]) -> str:
    """Return green for upward flow, red for downward flow, grey for horizontal/zero."""
    q = branch["flow_rate"]
    if abs(q) <= FLOW_TOL:
        return "0.55"

    start, end = branch["start"], branch["end"]
    # Positive flow follows start -> end; negative flow follows end -> start.
    source, target = (start, end) if q > 0.0 else (end, start)
    dz = positions[target][1] - positions[source][1]
    if dz > VERTICAL_TOL:
        return "tab:green"
    if dz < -VERTICAL_TOL:
        return "tab:red"
    return "0.55"


def _draw_nodes(ax, positions: dict[str, tuple[float, float]]) -> None:
    """Draw and label the shared network geometry."""
    for name, (x, z) in positions.items():
        ax.scatter(x, z, s=20, color="black", zorder=3)
        ax.annotate(name, (x, z), xytext=(4, 4), textcoords="offset points", fontsize=8)


def _finish_network_axes(ax, title: str) -> None:
    ax.set_title(title)
    ax.set_xlabel("Horizontal position (m)")
    ax.set_ylabel("Vertical position (m)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25)


def plot_flow_structure(results: dict[str, Any]) -> None:
    """Plot branch flow; width is scaled by |Q|/Q_min and colour by vertical direction.

    The reference flow is the smallest positive prescribed leaf demand. Thus a
    branch carrying one Q_min has the minimum configured width, and branches
    carrying more flow are drawn proportionally thicker.
    """
    branches = results["branches"]
    positions = results["positions"]
    positive_demands = [q for q in results["q_min_by_leaf"].values() if q > 0.0]
    reference_flow = min(positive_demands) if positive_demands else 0.0

    fig, ax = plt.subplots(figsize=(10, 8))
    for branch in branches:
        a, b = branch["start"], branch["end"]
        xa, za = positions[a]
        xb, zb = positions[b]
        q_abs = abs(branch["flow_rate"])
        if reference_flow > 0.0:
            width = FLOW_WIDTH_PER_QMIN * q_abs / reference_flow
        else:
            width = MIN_LINE_WIDTH
        width = max(width, 0.25)
        ax.plot(
            [xa, xb], [za, zb],
            linewidth=width,
            color=_flow_colour(branch, positions),
            solid_capstyle="round",
            zorder=1,
        )

        # Arrow direction follows the solved flow.
        source, target = ((xa, za), (xb, zb)) if branch["flow_rate"] >= 0.0 else ((xb, zb), (xa, za))
        ax.annotate(
            "",
            xy=(source[0] + 0.58 * (target[0] - source[0]),
                source[1] + 0.58 * (target[1] - source[1])),
            xytext=(source[0] + 0.42 * (target[0] - source[0]),
                    source[1] + 0.42 * (target[1] - source[1])),
            arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.0},
            zorder=2,
        )

    _draw_nodes(ax, positions)
    ax.legend(handles=[
        Line2D([0], [0], color="tab:green", lw=2, label="Upward flow"),
        Line2D([0], [0], color="tab:red", lw=2, label="Downward flow"),
        Line2D([0], [0], color="0.55", lw=2, label="Horizontal / near-zero flow"),
    ], loc="best")
    _finish_network_axes(ax, "Plant Network — Flow Rate")
    fig.tight_layout()


def plot_hydraulic_resistance_network(results: dict[str, Any]) -> None:
    """Plot branch resistance with a log-scaled colour map and line widths."""
    branches = results["branches"]
    positions = results["positions"]
    values = np.asarray([b["resistance"] for b in branches], dtype=float)
    if len(values) == 0 or np.any(values <= 0.0):
        raise ValueError("All branches must have positive hydraulic resistance.")
    lo, hi = float(values.min()), float(values.max())
    log_values = np.log10(values)
    norm = plt.Normalize(float(log_values.min()), float(log_values.max()) if hi > lo else float(log_values.min()) + 1.0)
    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(10, 8))
    for branch, lv in zip(branches, log_values):
        a, b = branch["start"], branch["end"]
        xa, za = positions[a]; xb, zb = positions[b]
        fraction = 0.0 if hi == lo else (lv - log_values.min()) / (log_values.max() - log_values.min())
        ax.plot([xa, xb], [za, zb], color=cmap(norm(lv)), linewidth=1.2 + 3.0 * fraction,
                solid_capstyle="round", zorder=1)
    for branch in sorted(branches, key=lambda b: b["resistance"], reverse=True)[:5]:
        a, b = branch["start"], branch["end"]
        xa, za = positions[a]; xb, zb = positions[b]
        ax.annotate(f"{a}–{b}", ((xa + xb) / 2, (za + zb) / 2), xytext=(3, 5),
                    textcoords="offset points", fontsize=8,
                    bbox={"boxstyle": "round,pad=0.2", "fc": "white", "alpha": 0.85, "ec": "none"})
    _draw_nodes(ax, positions)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    fig.colorbar(scalar, ax=ax, pad=0.02).set_label(r"Hydraulic resistance $R_i$ (Pa s/m$^3$), log scale")
    _finish_network_axes(ax, "Plant Network — Hydraulic Resistance by Branch")
    fig.tight_layout()


def plot_hydraulic_resistance_ranking(results: dict[str, Any]) -> None:
    """Rank all branches by hydraulic resistance."""
    branches = sorted(results["branches"], key=lambda b: b["resistance"])
    labels = [f'{b["start"]} → {b["end"]}' for b in branches]
    values = np.asarray([b["resistance"] for b in branches], dtype=float)
    fig, ax = plt.subplots(figsize=(10, max(5.0, 0.34 * len(branches))))
    y = np.arange(len(branches))
    bars = ax.barh(y, values, color=plt.get_cmap("viridis")(np.linspace(0.15, 0.9, len(branches))))
    ax.set_yticks(y, labels)
    ax.set_xlabel(r"Hydraulic resistance $R_i$ (Pa s/m$^3$)")
    ax.set_title("Branch Hydraulic Resistance — Lowest to Highest")
    ax.set_xscale("log")
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    for bar, value in zip(bars, values):
        ax.annotate(f"{value:.2e}", (value, bar.get_y() + bar.get_height() / 2),
                    xytext=(5, 0), textcoords="offset points", va="center", fontsize=8)
    fig.tight_layout()


def plot_results(results: dict[str, Any]) -> None:
    """Show the flow plot and both hydraulic-resistance plots."""
    plot_flow_structure(results)
    plot_hydraulic_resistance_network(results)
    plot_hydraulic_resistance_ranking(results)
    plt.show()


def print_report(results: dict[str, Any], duration: float = 1.0) -> None:
    """Print required leaf flows, solved pressures, and energy over duration."""
    if duration < 0.0:
        raise ValueError("duration must be non-negative.")

    print("Inverse plant water-uptake model")
    print(f"Required total leaf uptake : {results['total_leaf_demand']:.6e} m^3/s")
    print(f"Calculated root inflow     : {results['root_inflow']:.6e} m^3/s")
    print(
        f"Mass-balance relative error: "
        f"{results['mass_balance_relative_error']:.3e}"
    )
    print("\nLeaf-tip demand and solved pressure (leaf pressure is solved, not fixed):")
    for tip in results["leaf_tips"]:
        print(
            f"  {tip}: Q_min={results['q_min_by_leaf'][tip]:.6e} m^3/s, "
            f"p={results['pressures'][tip]:.6e} Pa"
        )

    print("\nEnergy rate:")
    print(f"Viscous dissipation power  : {results['total_viscous_power']:.6e} W")
    print(f"Signed gravitational power : {results['total_gravity_power']:.6e} W")
    print(f"Net mechanical power       : {results['net_mechanical_power']:.6e} W")
    print("  (viscous dissipation + net gravitational potential-energy rate)")
    print(f"Duration                   : {duration:.6g} s")
    print(f"Viscous energy             : {results['total_viscous_power'] * duration:.6e} J")
    print(f"Signed gravity energy      : {results['total_gravity_power'] * duration:.6e} J")
    print(f"Net mechanical energy      : {results['net_mechanical_power'] * duration:.6e} J")

    print("\nPer-branch flow and power:")
    for item in results["branches"]:
        direction = "a -> b" if item["flow_rate"] >= 0.0 else "b -> a"
        print(
            f"  {item['start']} -> {item['end']}: "
            f"Q={item['flow_rate']:.6e} m^3/s ({direction}), "
            f"P_visc={item['viscous_power']:.6e} W, "
            f"P_gravity={item['gravity_power']:.6e} W"
        )


def main() -> None:
    results = solve_for_minimum_leaf_flows(Q_MIN)
    print_report(results, duration=1.0)
    plot_results(results)


if __name__ == "__main__":
    main()
