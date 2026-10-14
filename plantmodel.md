# Plant Model — Hydraulic Transport and Energy

## Purpose

Develop a reduced-order model of water transport through a plant-like network, starting from Poiseuille flow in a single pipe and extending to branches, hydraulic energy, and plant water-potential boundaries.

The model is built in stages so each part can be understood and checked before adding biological complexity.

## Model scope and assumptions

The initial hydraulic model assumes:

- steady, incompressible flow;
- water treated as a Newtonian fluid with constant density $\rho$ and dynamic viscosity $\mu$;
- rigid, circular pipe segments;
- fully developed laminar flow in each segment;
- negligible kinetic-energy changes at junctions;
- no leaks or fluid sources inside a pipe segment.

These assumptions are useful for a first model, but real xylem has complex geometry, compliant tissues, changing conductance, and possible embolism. Do not interpret the initial model as a complete physiological model.

## 1. Symbols and units

| Symbol | Meaning | SI unit |
|---|---|---|
| $Q_i$ | Signed volumetric flow in segment $i$ | $\mathrm{m^3/s}$ |
| $p_j$ | Pressure at node $j$ | Pa |
| $z_j$ | Elevation of node $j$, positive upward | m |
| $L_i$ | Segment length | m |
| $r_i$ | Segment internal radius | m |
| $\mu$ | Dynamic viscosity | $\mathrm{Pa\,s}$ |
| $\rho$ | Fluid density | $\mathrm{kg/m^3}$ |
| $g$ | Gravitational acceleration | $\mathrm{m/s^2}$ |
| $\mathcal R_i$ | Hydraulic resistance of segment | $\mathrm{Pa\,s/m^3}$ |
| $P$ | Power | W |
| $E$ | Energy | J |

Use SI units internally. Convert millilitres per second to cubic metres per second using $1\,\mathrm{mL/s}=10^{-6}\,\mathrm{m^3/s}$.

## 2. Single-pipe Poiseuille model

For steady, fully developed laminar flow in a circular pipe of radius $r$, the axial velocity profile is

$
u(r')=-\frac{1}{4\mu}\frac{dp}{ds}(R^2-r'^2),
$

where $s$ points from inlet to outlet and $r'$ is the radial coordinate. Equivalently, this matches the sign convention in the existing derivation:

$
u(r')=\frac{1}{4\mu}\frac{dp}{ds}(r'^2-R^2).
$

Integrating the velocity profile over the cross-section gives the volumetric flow rate:

$
Q=\int_A u\,dA
=-\frac{\pi R^4}{8\mu}\frac{dp}{ds}
$

for a horizontal pipe. For a pipe whose elevation varies along its length, the axial momentum balance gives

$
Q=-\frac{\pi R^4}{8\mu}
\left(\frac{dp}{ds}+\rho g\frac{dz}{ds}\right).
$

Thus,

$
\frac{dp}{ds}=-\frac{8\mu Q}{\pi R^4}-\rho g\frac{dz}{ds}.
$

Define the Poiseuille hydraulic resistance for a constant-radius segment:

$
\boxed{\mathcal R_i=\frac{8\mu L_i}{\pi r_i^4}}
$

The frictional pressure drop in the positive flow direction is $Q_i\mathcal R_i$.

## 3. Pressure relation for one network segment

Orient each segment from endpoint node $a$ to endpoint node $b$. Define $Q_i>0$ as flow from $a$ to $b$. Let $\Delta z_i=z_b-z_a$.

Integrating the axial equation over the segment gives

$
\boxed{p_a-p_b=Q_i\mathcal R_i+\rho g(z_b-z_a)}
$

or

$
\boxed{Q_i=\frac{p_a-p_b-\rho g(z_b-z_a)}{\mathcal R_i}}.
$

A negative $Q_i$ means the actual flow is opposite to the chosen orientation. The sign of the elevation term follows from defining elevation as positive upward.

For a series path with a single flow rate $Q$, summing segment pressure drops gives

$
p_{\mathrm{in}}-p_{\mathrm{out}}
=Q\sum_i\mathcal R_i+\rho g(z_{\mathrm{out}}-z_{\mathrm{in}}).
$

Therefore,

$
\boxed{
Q=\frac{p_{\mathrm{in}}-p_{\mathrm{out}}-\rho g(z_{\mathrm{out}}-z_{\mathrm{in}})}
{\sum_i\mathcal R_i}
}
$

This series formula is for one unbranched path. Do not apply one global $Q$ to all segments in a branched network.

## 4. Branched network: mass conservation

Represent the network as nodes connected by segments. Each segment has a resistance and a chosen orientation.

At a steady junction with no storage or fluid source,

$
\boxed{\sum Q_{\mathrm{in}}=\sum Q_{\mathrm{out}}}
$

For a node-based signed formulation, define each edge flow as positive from its chosen start node $a$ to end node $b$. At an internal node $j$, the sum of signed flows leaving the node is zero:

$
\boxed{\sum_{i\in\mathrm{incident}(j)} s_{ji}Q_i=0}
$

where $s_{ji}=+1$ if edge $i$ is oriented outward from node $j$, and $s_{ji}=-1$ if it is oriented into node $j$.

For one inlet and two outlets:

$
\boxed{Q_{\mathrm{in}}=Q_1+Q_2}.
$

Each branch has its own flow rate. The shared junction has one pressure value.

## 5. A convenient network variable: hydraulic potential

Define pressure plus gravitational potential per unit volume at node $j$:

$
\boxed{\phi_j=p_j+\rho gz_j}
$

Then the segment equation becomes

$
\boxed{Q_i=\frac{\phi_a-\phi_b}{\mathcal R_i}}.
$

This is algebraically equivalent to the pressure/elevation equation above. It can make the network equations simpler because gravity is incorporated into the node potential.

If endpoint pressures are specified, calculate their potentials using their elevations. Solve for the unknown internal-node potentials, then recover pressure from

$
\boxed{p_j=\phi_j-\rho gz_j}.
$

### Linear system for unknown internal nodes

For each internal node $j$, substitute the edge flow expressions into mass conservation:

$
\sum_{k\text{ connected to }j}\frac{\phi_j-\phi_k}{\mathcal R_{jk}}=0.
$

Separate known boundary-node potentials from unknown internal potentials. This produces a linear system

$
\boxed{A\boldsymbol{\phi}=\boldsymbol b}.
$

Here, the diagonal coefficient for node $j$ is the sum of conductances $1/\mathcal R_{jk}$ connected to that node; an off-diagonal coefficient for connected nodes $j,k$ is $-1/\mathcal R_{jk}$. Contributions from known boundary potentials go into $\boldsymbol b$.

After solving for node potentials:

1. recover node pressures;
2. calculate each edge flow;
3. check mass balance at every internal node.

The system needs enough boundary conditions to be uniquely solvable. For the first solver, use specified pressures at the inlet and outlets. A prescribed-flow boundary can be added later.

## 6. Hydraulic power and energy

### 6.1 Viscous dissipation

For segment $i$, the signed pressure loss due to friction in the direction of its flow is $Q_i\mathcal R_i$. The power dissipated is nonnegative:

$
\boxed{P_{\mathrm{friction},i}=Q_i^2\mathcal R_i}
$

Total viscous dissipation across the network is

$
\boxed{P_{\mathrm{friction,total}}=\sum_i Q_i^2\mathcal R_i}.
$

Units: $\mathrm{m^3/s}\times\mathrm{Pa}=\mathrm W$. If flow and conditions are constant over time $t$,

$
\boxed{E_{\mathrm{friction}}=P_{\mathrm{friction,total}}t}.
$

If flow changes with time,

$
E_{\mathrm{friction}}=\int_0^t P_{\mathrm{friction,total}}(\tau)\,d\tau.
$

### 6.2 Gravitational energy

For a volume $V$ of water moved through elevation difference $\Delta z$,

$
\boxed{\Delta E_g=\rho Vg\Delta z}
$

For a steady flow rate $Q$,

$
\boxed{P_g=\rho gQ\Delta z}
$

for the chosen direction along that segment. This term is positive for an upward elevation change and negative for a downward one. In a network, use signed segment flows and elevations consistently.

Gravity may assist or oppose flow; gravitational energy change is not automatically energy actively expended by the plant.

### 6.3 Pressure work and interpretation

Pressure power crossing a boundary is related to pressure and volumetric flow by

$
P_p=pQ.
$

For a complete control-volume energy balance, account consistently for pressure work at inlets/outlets, gravitational potential-energy change, kinetic-energy change if significant, and viscous dissipation. Do not simply add every pressure drop and every gravitational term as positive quantities: the terms have signs and may represent energy supplied, stored, or dissipated.

The initial project's main hydraulic outputs should be:
- flow rate in every segment;
- pressure at every node;
- viscous power dissipated in every segment and in total;
- gravitational potential-energy change;
- the assumptions and boundary conditions used.

## 7. Connect the hydraulic network to plant water transport

Pressure alone does not capture all drivers of water movement in plants. Water potential combines several contributions. A useful form is

$
\boxed{\Psi_w=\Psi_p+\Psi_s+\Psi_m+\Psi_g}
$

where:
- $\Psi_p$: pressure potential;
- $\Psi_s$: solute/osmotic potential;
- $\Psi_m$: matric potential, relevant to binding of water in porous materials and tissues;
- $\Psi_g$: gravitational potential, commonly $\rho gz$ when expressed in pressure units.

The relevant terms depend on the tissue and model scale. Avoid counting gravity twice: if using $\Psi_g$ in the total water potential, do not add a separate gravity term again in the same driving-potential difference.

For a simple hydraulic pathway with conductance $K$,

$
\boxed{Q=K(\Psi_{\mathrm{upstream}}-\Psi_{\mathrm{downstream}})}
$

or, with hydraulic resistance $R_{\mathrm{hyd}}=1/K$,

$
\boxed{Q=\frac{\Delta\Psi_w}{R_{\mathrm{hyd}}}}.
$

For a pipe-only model, this is equivalent to the pressure-plus-gravity formulation. For living tissues, effective conductance and additional water-potential terms may be needed.

### Transpiration boundary condition

At steady state with negligible change in internal water storage,

$
\boxed{Q_{\mathrm{root uptake}}\approx E_{\mathrm{transpiration}}}
$

where transpiration is expressed as a volumetric flow rate. During transient conditions, water storage changes matter:

$
\boxed{Q_{\mathrm{root uptake}}-E_{\mathrm{transpiration}}=\frac{dV_{\mathrm{stored}}}{dt}}.
$

A later model can prescribe transpiration as a demand at leaf nodes or model it from environmental conditions. Begin with a prescribed rate rather than adding stomatal physics immediately.

## 8. Biological energy: what can and cannot be inferred

The hydraulic model can estimate:
- the viscous energy dissipated while water moves;
- the gravitational potential-energy change of transported water;
- the hydraulic power associated with specified flows and water-potential differences.

It cannot, by itself, determine total metabolic energy expenditure or ATP consumption. Plant xylem transport is primarily driven by water-potential gradients associated with transpiration; root ion transport and other cellular processes can consume metabolic energy, but their costs require additional physiological data or a separately justified model.

For a specified biological process, metabolic energy over time could be estimated from a measured metabolic power:

$
E_{\mathrm{metabolic}}=\int_0^t P_{\mathrm{metabolic}}(\tau)\,d\tau.
$

Do not equate $P_{\mathrm{friction,total}}$ with $P_{\mathrm{metabolic}}$, and do not apply a universal conversion factor between them.

## 9. Verification plan

Complete these checks in order:

1. **Horizontal single pipe:** compare with analytical Poiseuille flow:
   $Q=\pi R^4(p_{\mathrm{in}}-p_{\mathrm{out}})/(8\mu L)$.
2. **Inclined single pipe:** verify the gravity sign and hydrostatic contribution.
3. **Series segments:** verify that summed resistance reproduces the equivalent single-pipe result for constant radius.
4. **One inlet, two outlets:** verify $Q_{\mathrm{in}}=Q_1+Q_2$ and equal junction pressure for both branches.
5. **General network:** verify mass balance at every internal node and check that every edge flow agrees with its endpoint potential difference.
6. **Energy:** confirm each $Q_i^2\mathcal R_i$ is nonnegative and units are watts.
7. **Plant comparison:** when possible, compare model outputs with published plant hydraulic measurements under matching conditions.

Analytical checks verify the implementation for idealized cases. Published measurements are needed to validate predictions of real plant behaviour.

## 10. Recommended implementation order

Keep the implementation focused. Do not build every phase at once.

- [ ] Preserve the existing single-pipe analytical solution.
- [ ] Verify the horizontal and inclined single-pipe cases.
- [ ] Implement a small network with one inlet and two outlets.
- [ ] Solve unknown internal node potentials using mass conservation.
- [ ] Calculate pressure and flow on every segment.
- [ ] Add frictional power and gravitational energy outputs.
- [ ] Add root/leaf water-potential or transpiration boundary conditions.
- [ ] Only then investigate metabolic energy using separate physiological evidence.

**First milestone:** a branched hydraulic network that calculates flow and pressure correctly. Biological energy estimation is a later layer, not a prerequisite for completing the hydraulic solver.
