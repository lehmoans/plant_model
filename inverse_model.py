"""Energy accounting using the channel network defined in model.py.

This module imports MU, RHO, G, BRANCHES and solve_network from model.py,
so the geometry and boundary conditions have one source of truth.

Important: model.solve_network() solves flow from the existing pressure
boundary conditions (P_ROOT_TIP and P_LEAF_TIP). It does not impose a
user-specified minimum flow at each tip. This module calculates the energy
associated with those solved flows.
"""

from __future__ import annotations

from math import pi
from typing import Any

from model import BRANCHES, G, MU, RHO, solve_network


def calculate_energy(results: dict[str, Any] | None = None,
                     duration: float = 1.0) -> dict[str, Any]:
    """Calculate viscous dissipation and gravitational power for model.py's network.

    Args:
        results: Optional output of model.solve_network(). If omitted, the
            network is solved using model.py's current boundary conditions.
        duration: Time interval in seconds over which to integrate power.

    Returns:
        A dictionary with per-branch details and total power/energy in SI units.

    For each branch:
        R = 8 μ L / (π r⁴)
        P_viscous = R Q²
        P_gravity = ρ g Q_signed Δz

    Q_signed is positive in the declared BRANCHES direction and negative
    when the actual flow is opposite. Positive gravitational power means
    the actual flow moves net uphill along that branch; negative means
    gravity assists the flow.
    """
    if duration < 0:
        raise ValueError("duration must be non-negative.")

    if results is None:
        results = solve_network()

    branches = results["branches"]
    flows = results["flows"]
    positions = results["positions"]

    branch_results = []
    total_viscous_power = 0.0
    total_gravity_power = 0.0

    for (a, b, length, radius, theta), q in zip(branches, flows):
        if length <= 0 or radius <= 0:
            raise ValueError(
                f"Branch {a!r}->{b!r} must have positive length and radius."
            )

        resistance = 8.0 * MU * length / (pi * radius**4)
        z_a = positions[a][1]
        z_b = positions[b][1]
        delta_z = z_b - z_a

        q_signed = float(q)
        viscous_power = resistance * q_signed**2
        gravity_power = RHO * G * q_signed * delta_z

        branch_results.append({
            "start": a,
            "end": b,
            "length": length,
            "radius": radius,
            "angle_deg": theta,
            "flow_rate": q_signed,
            "resistance": resistance,
            "pressure_drop_viscous": resistance * q_signed,
            "delta_z": delta_z,
            "viscous_power": viscous_power,
            "gravity_power": gravity_power,
        })
        total_viscous_power += viscous_power
        total_gravity_power += gravity_power

    # Viscous energy is always dissipated. Gravitational energy is signed:
    # it can assist flow or represent energy that must be supplied uphill.
    # Do not clip negative gravitational contributions; retain the physical
    # signed energy balance.
    total_power = total_viscous_power + total_gravity_power

    return {
        "branches": branch_results,
        "root_tip_pressures": {
            name: results["pressures"][name]
            for name in results["nodes"]
            if name.startswith("j") and name[1:].isdigit()
        },
        "total_viscous_power": total_viscous_power,
        "total_gravity_power": total_gravity_power,
        "net_mechanical_power": total_power,
        "duration": duration,
        "viscous_energy": total_viscous_power * duration,
        "gravity_energy": total_gravity_power * duration,
        "net_mechanical_energy": total_power * duration,
    }


def print_report(energy: dict[str, Any]) -> None:
    """Print the total energy and each channel's contribution."""
    print("Plant network energy report")
    print(f"Viscous dissipation power : {energy['total_viscous_power']:.6e} W")
    print(f"Signed gravity power      : {energy['total_gravity_power']:.6e} W")
    print(f"Net mechanical power      : {energy['net_mechanical_power']:.6e} W")
    print(f"Integration duration      : {energy['duration']:.6g} s")
    print(f"Viscous energy             : {energy['viscous_energy']:.6e} J")
    print(f"Signed gravity energy      : {energy['gravity_energy']:.6e} J")
    print(f"Net mechanical energy     : {energy['net_mechanical_energy']:.6e} J")

    print("\nPer-branch results (Q sign follows BRANCHES direction):")
    for item in energy["branches"]:
        direction = "listed direction" if item["flow_rate"] >= 0 else "reverse direction"
        print(
            f"  {item['start']} -> {item['end']}: "
            f"Q={item['flow_rate']:.6e} m^3/s ({direction}), "
            f"R={item['resistance']:.6e} Pa s/m^3, "
            f"P_visc={item['viscous_power']:.6e} W, "
            f"P_gravity={item['gravity_power']:.6e} W"
        )


def main() -> None:
    # Uses the actual NODE_TREE, BRANCHES, channel dimensions, angles and
    # pressure boundary conditions already defined in model.py.
    network_results = solve_network()
    energy = calculate_energy(network_results, duration=1.0)
    print_report(energy)


if __name__ == "__main__":
    main()
