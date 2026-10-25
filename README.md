# Plant Hydraulic Network Models

A code-based, reduced-order model of water transport through a plant-like branched network. The project began with single-pipe Poiseuille flow and has progressed to a network solver and an inverse flow/energy calculation.

The two main scripts answer complementary questions:

- **Forward model — `model.py`:** Given prescribed pressures at the root and leaf tips, what junction pressures and branch flow rates result?
- **Inverse model — `inverse_model.py`:** Given fixed root-tip pressure and a target flow at each leaf tip, what junction/leaf pressures and branch flows are required, and what steady mechanical power is associated with viscous losses and gravity?

> **Boundary-condition distinction:** The current inverse model prescribes leaf flow rates and solves the leaf pressures. It does not simultaneously enforce the leaf-tip pressure values in `model.py`. A passive resistance network cannot independently prescribe both pressure and flow at the same terminal. To enforce both root and leaf pressures while also imposing a flow target, an explicit active source/pump model is needed.

## Contents

- [Model scope](#model-scope)
- [Repository files](#repository-files)
- [Physical assumptions](#physical-assumptions)
- [Network geometry and data](#network-geometry-and-data)
- [Governing equations](#governing-equations)
- [Forward model](#forward-model)
- [Inverse model](#inverse-model)
- [Numerical stability](#numerical-stability)
- [Visualizations](#visualizations)
- [Running the models](#running-the-models)
- [Validation and limitations](#validation-and-limitations)
- [Next steps](#next-steps)

## Model scope

The goal is to understand the mechanical requirements for moving water through a branched, plant-inspired network. The implementation focuses on:

1. Single-pipe laminar flow as the starting point.
2. A branched network represented by circular channels with specified lengths and radii.
3. Gravity represented through elevation and hydraulic potential.
4. Junction mass conservation.
5. Forward prediction of flow from pressure boundary conditions.
6. Inverse calculation for prescribed leaf-tip uptake.
7. Branch resistance, viscous power, gravitational potential-energy rate, and total mechanical energy over a chosen duration.
8. Visual inspection of the network, solved flows, and branch hydraulic resistance.

This is a **reduced-order hydraulic model**, not a full 3D CFD solver and not yet a complete biological model of xylem transport.

## Repository files

| File | Purpose |
|---|---|
| `model.py` | Network geometry, physical parameters, forward pressure/flow solver, and associated checks/plots. |
| `inverse_model.py` | Inverse flow-demand solve, branch-flow and power calculations, report, and plots. |
| `plantmodel.md` | Longer-form background notes on hydraulic transport, gravity, energy, and the development plan. |
| `simple_pipe/pressure/analytical_soln.py` | Analytical single-pipe Poiseuille reference solution. |
| `simple_pipe/pressure/numerical_soln.py` | Earlier numerical-solution work; check its current implementation before treating it as a validated numerical solver. |

## Physical assumptions

The current branched-network models use these assumptions:

- Steady-state flow.
- Water is incompressible and Newtonian.
- Each branch is approximated as a rigid, circular pipe.
- Flow in each branch follows the Hagen–Poiseuille laminar-flow resistance relation.
- Branch resistance is constant for a given geometry and viscosity.
- Internal junctions conserve mass; no water is created or lost at a junction.
- Elevations are calculated from the branch lengths and angles in the geometry.
- Gravity is included through hydraulic potential.
- The network does not yet model changing channel radius, elastic walls, transient storage, cavitation, solute transport, or detailed biological processes.

These assumptions define the model's current range of applicability. In particular, Poiseuille resistance is appropriate for fully developed laminar flow in a circular tube; it is not a universal pressure-loss model for every possible flow regime or geometry.

## Network geometry and data

The network is defined in `model.py`.

### Physical parameters

Current default values:

| Parameter | Default | Units |
|---|---:|---|
| Dynamic viscosity, `MU` | `1e-3` | Pa·s |
| Density, `RHO` | `1000` | kg/m³ |
| Gravitational acceleration, `G` | `9.81` | m/s² |
| Root-tip gauge pressure, `P_ROOT_TIP` | `0.0` | Pa |
| Leaf-tip gauge pressure in forward model, `P_LEAF_TIP` | `0.0` | Pa |

Lengths and radii are supplied in metres. The constants for channel lengths and radii are used by the branch definitions in `BRANCHES`.

### Nodes and branches

- `NODE_TREE` stores node names and their pressure/elevation metadata.
- A node with a known pressure is a pressure-boundary node in the forward model.
- A node with `p=None` is an internal junction whose hydraulic potential is solved.
- `BRANCHES` stores each connection as `(node_a, node_b, length_m, radius_m, theta_deg)`.
- Positive branch flow is defined in the declared direction `node_a -> node_b`; a negative value means the actual flow goes in the reverse direction.
- `calculate_positions()` reconstructs node positions from the branch lengths and angles, using `g` as the reference point.

The geometry is a 2D schematic of the plant-inspired network. The calculated vertical coordinate is used as elevation in the gravity term. Any changes to branch lengths or angles can therefore change both resistance and elevation head.

## Governing equations

### 1. Hydraulic resistance

For steady, fully developed laminar flow in a rigid circular pipe, the Hagen–Poiseuille relation gives

$
R_i = \frac{8\mu L_i}{\pi r_i^4}
$

where:

- $R_i$ is branch hydraulic resistance in Pa·s/m³;
- $\mu$ is dynamic viscosity in Pa·s;
- $L_i$ is branch length in metres;
- $r_i$ is internal radius in metres.

Resistance grows linearly with length and with the inverse fourth power of radius. Consequently, small changes in radius can have a large effect on the network solution.

### 2. Hydraulic potential and gravity

The network uses hydraulic potential

$
\phi_i = p_i + \rho g z_i
$

where $p_i$ is pressure, $z_i$ is elevation, $\rho$ is fluid density, and $g$ is gravitational acceleration.

For a branch directed from node $a$ to node $b$, the signed flow is

$
Q_{a\to b} = \frac{\phi_a-\phi_b}{R_{ab}}
$

Positive flow follows the branch's declared direction. Negative flow means flow is in the opposite direction.

Equivalently, the potential difference accounts for both pressure and elevation:

$
\phi_a-\phi_b = (p_a-p_b) - \rho g(z_b-z_a)
$

This is why a pressure difference alone is not sufficient to describe the driving potential when the endpoints have different elevations.

### 3. Mass conservation at junctions

For each internal junction $i$, steady flow requires

$
\sum_j Q_{j\to i} = 0
$

or equivalently, the signed sum of flows leaving the junction is zero. The solvers assemble a linear conductance matrix from each branch resistance and enforce this condition at unknown-pressure junctions.

### 4. Viscous power

The rate at which mechanical energy is dissipated by viscous resistance in branch $i$ is

$
P_{\mathrm{viscous},i}=R_i Q_i^2
$

and for the whole network,

$
P_{\mathrm{viscous}}=\sum_i R_i Q_i^2
$

This contribution is non-negative regardless of the direction of flow.

### 5. Gravitational potential-energy rate

For signed flow $Q_i$ in the declared direction from $a$ to $b$, with $\Delta z_i=z_b-z_a$,

$
P_{\mathrm{gravity},i}=\rho g Q_i\Delta z_i
$

The signed network contribution is

$
P_{\mathrm{gravity}}=\sum_i \rho g Q_i\Delta z_i
$

Upward transport contributes positively; downward transport can contribute negatively. Using signed branch flows makes the expression consistent when the actual flow opposes the declared branch direction.

### 6. Total mechanical power and energy

The inverse model reports the steady mechanical power requirement as

$
P_{\mathrm{net}}=P_{\mathrm{viscous}}+P_{\mathrm{gravity}}
$

For a constant steady power over duration $T$,

$
E=P T
$

The report prints viscous, gravitational, and net mechanical power and the corresponding energy over the selected duration.

**Interpretation:** this is a mechanical hydraulic-energy accounting for the current idealized network. It is not a prediction of biochemical energy expenditure or metabolic efficiency. Those would require a biological source/transport model and additional assumptions.

## Forward model

Run `model.py` to solve the network with pressure boundary conditions.

### Inputs

- Network geometry and branch dimensions.
- Fluid properties (`MU`, `RHO`, `G`).
- Root-tip and leaf-tip pressures, as specified in `NODE_TREE` through `P_ROOT_TIP` and `P_LEAF_TIP`.

### Method

1. Flatten the nested node tree and validate that nodes and branches agree.
2. Calculate node positions and elevations from branch geometry.
3. Calculate the Poiseuille resistance of each branch.
4. Convert prescribed boundary pressures to hydraulic potentials.
5. Assemble the linear system from junction mass conservation.
6. Solve the unknown junction potentials.
7. Recover pressure and calculate signed branch flows from potential differences.

### Outputs

The forward solver returns node data, positions, hydraulic potentials, pressures, branch resistances, the assembled matrix and right-hand side, unknown-node ordering, branch flows, and branch definitions. The script also contains mass-conservation checking and visualization/reporting utilities.

The forward model is the natural reference case for checking network equations because the boundary pressures are prescribed and the resulting flows are determined by the network.

## Inverse model

Run `inverse_model.py` to prescribe leaf uptake and solve the corresponding network state.

### Current boundary conditions

- Root-tip pressure is fixed at `P_ROOT_TIP` from `model.py`.
- Each leaf tip is assigned a prescribed demand using `Q_MIN`.
- The default is `Q_MIN = 1.0e-9 m^3/s` **per leaf tip**.
- The four leaf tips are `n11`, `n12`, `n21`, and `n22`.
- The four root tips are `j11`, `j12`, `j21`, and `j22`.
- Internal junction potentials and leaf-tip pressures are solved rather than fixed.

The total prescribed leaf uptake at the default setting is $4\times10^{-9}\,\mathrm{m^3/s}$.

### Method

1. Reuse the geometry, fluid constants, node tree, and branch definitions from `model.py`.
2. Calculate branch Poiseuille resistances.
3. Shift the hydraulic-potential reference so the common root potential is zero for the linear solve.
4. Apply fixed root potentials and impose the leaf-tip flow demands in the nodal balance equations.
5. Solve for unknown relative hydraulic potentials.
6. Calculate each signed branch flow using the relative potential difference divided by branch resistance.
7. Recover physical hydraulic potentials and pressures for reporting.
8. Calculate branch viscous power and signed gravitational power, then sum them over the network.
9. Calculate root inflow and compare it with the total leaf demand as a mass-balance check.

If root tips do not share the same elevation, the current inverse implementation raises an error rather than silently applying a single shared reference potential.

### Outputs

The inverse solver returns prescribed leaf demands, root/leaf tip names, calculated root inflow, total leaf demand, mass-balance error and relative error, node positions, absolute and relative hydraulic potentials, solved pressures, matrix and right-hand side, unknown ordering, per-branch results, total viscous power, total gravity power, and net mechanical power.

The console report includes solved leaf pressures, flow and power per branch, total uptake, mass-balance error, and energy over a chosen duration.

### Boundary-condition limitation

The forward model currently marks both root tips and leaf tips as pressure boundaries. The inverse model instead uses root-tip pressure boundaries and prescribed leaf flows. This difference is intentional: in a passive network, specifying both pressure and flow independently at the same terminal generally overconstrains the problem.

If the intended biological scenario requires both root and leaf pressures to remain at a common reference pressure while a nonzero target flow is maintained, the model will need explicit active pressure sources/pumps or another specified active mechanism. A zero gauge pressure does not mean zero hydraulic potential at every elevation, because $\phi=p+\rho gz$.

## Numerical stability

The inverse model solves in hydraulic potential relative to the common root potential rather than solving with large absolute potentials and subtracting them to calculate tiny flow rates.

This matters because

$
Q=\frac{\phi_a-\phi_b}{R}
$

can be sensitive to floating-point cancellation when $\phi_a$ and $\phi_b$ are large and nearly equal. Defining

$
\phi_i'=\phi_i-\phi_{\mathrm{ref}}
$

preserves potential differences mathematically:

$
\phi_a'-\phi_b'=\phi_a-\phi_b
$

and can improve numerical behavior when the equations are solved consistently in the shifted reference. This is a change of reference, not a change in physical pressure or gravity.

The inverse model also reports mass-balance error so numerical results can be checked rather than judged by plots alone. Increasing `Q_MIN` is not a substitute for a numerically stable formulation; it changes the physical demand being simulated.

## Visualizations

The inverse model produces three figures:

1. **Flow structure:** branch width represents flow magnitude; branch color indicates upward, downward, or horizontal/near-zero flow. Arrows show the solved direction.
2. **Hydraulic resistance network:** each branch is colored by resistance using a logarithmic color mapping; line thickness also increases with resistance. Only the five highest-resistance branches receive text labels.
3. **Ranked hydraulic resistance:** a horizontal bar chart lists all branches from lowest to highest resistance on a logarithmic axis, with numeric values annotated.

The resistance visualizations use the computed branch resistance $R_i=8\mu L_i/(\pi r_i^4)$, not a separate dimensionless friction factor. High resistance identifies a branch that strongly opposes flow for a given potential difference; it does not by itself establish that the branch dissipates the most power. Viscous power also depends on the flow through the branch, $R_iQ_i^2$.

## Running the models

Use Python 3 with NumPy and Matplotlib installed. From the repository root:

```bash
python model.py
python inverse_model.py
```

If using a virtual environment, activate it first and install the required dependencies, for example:

```bash
python -m pip install numpy matplotlib
```

Change `Q_MIN` in `inverse_model.py` to explore other prescribed leaf-flow requirements. The current solver is steady-state; the reported energy for a duration assumes the calculated power remains constant over that duration.

## Validation and limitations

### Checks implemented or used during development

- The forward network solver validates that branch nodes exist, nodes are connected to branches, and branch resistances are positive and finite.
- Junction equations enforce steady mass conservation.
- The inverse model reports root inflow versus total prescribed leaf demand and a relative mass-balance error.
- The forward model has been hand-checked against calculations during development.
- The inverse implementation has been developed with a shifted-potential solve to reduce cancellation when computing small flows.

### Checks still recommended

- Run both scripts locally and inspect all figures and console output.
- Independently test junction residuals, not only the total root-versus-leaf balance.
- Compare simple one-pipe and small-network cases with hand calculations.
- Test sensitivity to radius, length, viscosity, elevation, and `Q_MIN`.
- Verify expected scaling for this linear fixed-resistance model: scaling all prescribed flows by a factor $c$ scales flows and viscous power by $c$ and $c^2$, respectively; gravitational power scales linearly with flow for fixed geometry.
- Check energy accounting against the boundary mechanical-energy balance, especially when changing the topology, elevations, or boundary conditions.
- Add automated regression tests for matrix assembly, flow signs, conservation, and the resistance calculation.

No claim is made here that a full automated test suite has passed. The model should be treated as a developing engineering tool whose equations, numerical behavior, and assumptions need continued verification.

### Current limitations

- Steady-state rather than transient flow.
- Rigid channels with constant geometry.
- Poiseuille resistance in every branch.
- No pump/source terms in the current inverse model.
- No explicit biochemical energy or active transport mechanism.
- Geometry and elevations are derived from the branch list; inconsistent cycles or geometry should be checked when the network is edited.
- The inverse reference shift currently assumes all root tips share an elevation.

## Next steps

1. Run the forward and inverse scripts and save representative outputs.
2. Add automated tests for nodal conservation and branch-flow calculations.
3. Validate energy accounting independently against a mechanical-energy balance.
4. Decide how active sources should be represented if both root and leaf pressures must be prescribed alongside a target uptake.
5. Extend the model only after the existing hydraulic calculations and boundary conditions are independently validated.
6. Consider more complex flow physics after the reduced-order network is reliable.

