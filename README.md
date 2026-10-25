# Pipe Flow and Plant Water-Transport Models

This repository explores fluid mechanics through progressively more complex, code-based models. The broad aim is to understand the governing equations, implement them transparently, and validate results against analytical relationships before extending the models toward more realistic systems.

The two main directories address different scales of the same problem: **`simple_pipe/` isolates the physics of flow through an individual conduit; `plant_model/` connects many conduits into a branched plant-inspired transport network.**

## Project areas

### `simple_pipe/` — Fundamental pipe-flow physics

**Idea:** Build a reliable baseline for understanding how pressure, viscous resistance, pipe geometry, and elevation affect steady water flow through a single cylindrical pipe.

**Goals**
- Implement and understand analytical solutions for steady, incompressible, laminar flow of a Newtonian fluid.
- Calculate the pressure-driven velocity profile and volumetric flow rate using the Hagen–Poiseuille relationship.
- Account for pipe dimensions and, in the relevant cases, elevation and gravity.
- Check units, signs, boundary conditions, and limiting cases so each result has a clear physical interpretation.
- Use these simple cases as reference tests for later numerical or computational fluid-dynamics (CFD) work.

This directory is intended to establish trustworthy, interpretable building blocks. It is not intended to be a full general-purpose pipe-network or CFD solver.

### `plant_model/` — Branched plant water-transport network

**Idea:** Represent a simplified plant water-transport system as a network of connected cylindrical conduits, from roots through the stem to branches and leaf tips. Treat the network first as a hydraulic problem, then investigate what its predicted flow implies for energy expenditure.

**Goals**
- Represent conduits as links between boundary nodes and junctions, with each segment having properties such as length, radius, and orientation.
- Calculate each segment's hydraulic resistance, initially using the laminar Hagen–Poiseuille model.
- Solve for node pressures or hydraulic potentials and the flow rate in each individual segment.
- Enforce mass conservation at junctions: flow into a junction must balance flow out.
- Account for elevation and gravity consistently with the chosen pressure/potential convention.
- Estimate hydraulic power associated with viscous pressure losses and gravitational lifting, then develop a careful model of the energy required to sustain transport.
- Keep hydraulic energy requirements distinct from the biological mechanisms that supply that energy.
- Build toward plant-inspired questions such as how network structure and conduit dimensions influence flow distribution and energy use.

The initial scope is deliberately simplified: water is treated as an incompressible Newtonian fluid, and xylem-like water transport is the primary focus. More complex biological effects, transpiration-driven boundary conditions, and possible phloem comparisons can be introduced later as separate extensions.

## How the projects relate

1. **Understand one conduit** in `simple_pipe/`: derive and check the pressure–flow relationship, including elevation where appropriate.
2. **Connect conduits into a network** in `plant_model/`: apply segment relationships together with junction conservation to find the distribution of flow.
3. **Investigate system-level behaviour:** calculate the network's hydraulic losses and gravitational requirements, and explore how geometry and branching affect the total.
4. **Validate each step:** compare code against hand calculations and analytical cases, then add complexity only when the simpler model is understood and checked.

The aim is not merely to produce a working result. The equations, assumptions, numerical implementation, and validation should remain understandable enough that each model can be explained and independently checked.
