"""Inverse hydraulic model for water uptake through a branched channel network.

The model takes a channel network and prescribed leaf-tip flow demands, works
backwards to calculate flow in every upstream channel, then estimates the
mechanical power/energy needed to overcome viscous resistance and elevation.

SI units throughout:
    radius [m], length [m], elevation [m], flow [m^3/s],
    viscosity [Pa s], density [kg/m^3], time [s], energy [J].

Network format
--------------
channels is an iterable of mappings, one per channel:
    {"start": "root", "end": "junction_1", "length": 0.1,
     "radius": 1e-3, "z_start": 0.0, "z_end": 0.1}

Each channel is directed from upstream/root toward downstream/tip. The
network is assumed to be a tree (no loops), with one or more root nodes and
one or more terminal tip nodes. Tip demands map tip node names to positive
outward flow rates, e.g. {"leaf_a": 1e-9, "leaf_b": 2e-9}.

This file intentionally does not import the project's other modules, so it
can be tested independently. Convert the existing channel objects to the
small mapping format above at the call site.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from math import pi
from typing import Hashable, Iterable, Mapping, Any


Node = Hashable


@dataclass(frozen=True)
class ChannelResult:
    """Flow and energy contribution for one channel."""

    start: Node
    end: Node
    flow_rate: float
    resistance: float
    pressure_loss: float
    viscous_power: float
    gravity_power: float


def hydraulic_resistance(length: float, radius: float,
                         viscosity: float = 1.0e-3) -> float:
    """Hagen–Poiseuille resistance R = 8 μ L / (π r⁴), in Pa s/m³."""
    if length <= 0:
        raise ValueError(f"Channel length must be > 0; got {length!r}")
    if radius <= 0:
        raise ValueError(f"Channel radius must be > 0; got {radius!r}")
    if viscosity <= 0:
        raise ValueError(f"Dynamic viscosity must be > 0; got {viscosity!r}")
    return 8.0 * viscosity * length / (pi * radius**4)


def _channel_value(channel: Mapping[str, Any], key: str,
                   default: Any = None) -> Any:
    if key in channel:
        return channel[key]
    return default


def calculate_inverse_model(
    channels: Iterable[Mapping[str, Any]],
    tip_demands: Mapping[Node, float],
    *,
    viscosity: float = 1.0e-3,
    density: float = 1000.0,
    gravity: float = 9.81,
    duration: float = 1.0,
) -> dict[str, Any]:
    """Calculate channel flows and the energy required for the prescribed uptake.

    Flow is calculated by conservation of mass: each channel carries the sum
    of all tip demands downstream of its end node. This is valid for a
    directed, acyclic tree with no internal sources/sinks.

    Viscous hydraulic power is sum(R_i Q_i²). Gravity power is calculated
    channel-by-channel as sum(ρ g Q_i Δz_i), where Δz is end elevation minus
    start elevation. A positive value means net uphill transport. Negative
    contributions represent downhill segments and may offset gravitational
    potential energy in an ideal energy balance.

    Total energy = total power * duration. This is the required ideal
    mechanical energy for the specified steady flow, not the plant's
    biochemical/metabolic energy cost or an efficiency-adjusted estimate.
    """
    channels = [dict(channel) for channel in channels]
    if not channels:
        raise ValueError("At least one channel is required.")
    if not tip_demands:
        raise ValueError("Provide at least one prescribed tip flow demand.")
    if viscosity <= 0 or density <= 0 or gravity < 0 or duration < 0:
        raise ValueError("Viscosity/density must be positive; gravity/duration non-negative.")
    for tip, flow in tip_demands.items():
        if flow < 0:
            raise ValueError(f"Tip demand for {tip!r} cannot be negative.")

    # Validate and build directed adjacency / incoming-edge count.
    outgoing: dict[Node, list[int]] = defaultdict(list)
    incoming_count: dict[Node, int] = defaultdict(int)
    nodes: set[Node] = set()
    for i, ch in enumerate(channels):
        if "start" not in ch or "end" not in ch:
            raise KeyError("Each channel must define 'start' and 'end'.")
        start, end = ch["start"], ch["end"]
        if start == end:
            raise ValueError(f"Self-loop channel at node {start!r}.")
        if "length" not in ch or "radius" not in ch:
            raise KeyError("Each channel must define 'length' and 'radius'.")
        outgoing[start].append(i)
        incoming_count[end] += 1
        nodes.update((start, end))

    unknown_tips = set(tip_demands) - nodes
    if unknown_tips:
        raise ValueError(f"Tip nodes not found in channels: {sorted(map(str, unknown_tips))}")

    # Topological sort detects loops and ensures children are processed first
    # when calculating downstream demand.
    indegree = {node: incoming_count[node] for node in nodes}
    queue = deque(node for node in nodes if indegree[node] == 0)
    topo: list[Node] = []
    while queue:
        node = queue.popleft()
        topo.append(node)
        for edge_idx in outgoing.get(node, []):
            child = channels[edge_idx]["end"]
            indegree[child] -= 1
            if indegree[child] == 0:
                queue.append(child)
    if len(topo) != len(nodes):
        raise ValueError("Network contains a directed cycle; this implementation expects a tree/DAG.")

    # Every node's demand is its own prescribed terminal demand plus the
    # demands of all children. Reverse topological order propagates flow back
    # toward roots, automatically summing demands across branches.
    downstream_demand = {node: float(tip_demands.get(node, 0.0)) for node in nodes}
    for node in reversed(topo):
        for edge_idx in outgoing.get(node, []):
            child = channels[edge_idx]["end"]
            downstream_demand[node] += downstream_demand[child]

    # A channel flow equals the total demand downstream of its end node.
    results: list[ChannelResult] = []
    total_viscous_power = 0.0
    total_gravity_power = 0.0
    total_root_flow = sum(downstream_demand[node] for node in nodes
                          if incoming_count[node] == 0)
    for ch in channels:
        start, end = ch["start"], ch["end"]
        q = downstream_demand[end]
        resistance = hydraulic_resistance(
            float(ch["length"]), float(ch["radius"]), viscosity
        )
        dz = float(_channel_value(ch, "z_end", 0.0)) - float(
            _channel_value(ch, "z_start", 0.0)
        )
        viscous_power = resistance * q * q
        gravity_power = density * gravity * q * dz
        results.append(ChannelResult(
            start=start,
            end=end,
            flow_rate=q,
            resistance=resistance,
            pressure_loss=resistance * q,
            viscous_power=viscous_power,
            gravity_power=gravity_power,
        ))
        total_viscous_power += viscous_power
        total_gravity_power += gravity_power

    total_power = total_viscous_power + total_gravity_power
    # Negative total mechanical power would mean net gravitational assistance
    # exceeds viscous dissipation. Report it rather than silently clipping it.
    return {
        "channels": results,
        "tip_demands": dict(tip_demands),
        "root_flow_rate": total_root_flow,
        "total_viscous_power": total_viscous_power,
        "total_gravity_power": total_gravity_power,
        "total_power": total_power,
        "duration": duration,
        "viscous_energy": total_viscous_power * duration,
        "gravity_energy": total_gravity_power * duration,
        "total_energy": total_power * duration,
        "units": {
            "flow_rate": "m^3/s",
            "resistance": "Pa s/m^3",
            "pressure_loss": "Pa",
            "power": "W",
            "energy": "J",
        },
    }


def print_report(result: Mapping[str, Any]) -> None:
    """Print a compact human-readable report."""
    print("Inverse water-uptake model")
    print(f"Total root uptake: {result['root_flow_rate']:.6e} m^3/s")
    print(f"Viscous power:     {result['total_viscous_power']:.6e} W")
    print(f"Gravity power:     {result['total_gravity_power']:.6e} W")
    print(f"Total power:       {result['total_power']:.6e} W")
    print(f"Duration:          {result['duration']:.6g} s")
    print(f"Viscous energy:    {result['viscous_energy']:.6e} J")
    print(f"Gravity energy:    {result['gravity_energy']:.6e} J")
    print(f"Total energy:      {result['total_energy']:.6e} J")
    print("\nPer-channel flows:")
    for item in result["channels"]:
        print(
            f"  {item.start!r} -> {item.end!r}: "
            f"Q={item.flow_rate:.6e} m^3/s, "
            f"R={item.resistance:.6e} Pa s/m^3, "
            f"P_visc={item.viscous_power:.6e} W, "
            f"P_gravity={item.gravity_power:.6e} W"
        )


if __name__ == "__main__":
    # Small illustrative branched network. Replace with channels from your
    # model and the minimum flow demands you want at each leaf tip.
    example_channels = [
        {"start": "root", "end": "j1", "length": 0.10, "radius": 1.0e-3,
         "z_start": 0.0, "z_end": 0.10},
        {"start": "j1", "end": "leaf_a", "length": 0.20, "radius": 5.0e-4,
         "z_start": 0.10, "z_end": 0.30},
        {"start": "j1", "end": "leaf_b", "length": 0.25, "radius": 4.0e-4,
         "z_start": 0.10, "z_end": 0.35},
    ]
    example_demands = {"leaf_a": 1.0e-9, "leaf_b": 1.0e-9}
    report = calculate_inverse_model(
        example_channels, example_demands, duration=3600.0
    )
    print_report(report)
