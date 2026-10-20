# Current Projects

## Purpose

Build a code-first engineering portfolio for graduate roles or student teams in **motorsport and rocketry**, with emphasis on fluid mechanics, simulation, numerical methods, verification, and engineering decisions.

## Constraints and working principles

- **Code-based only:** do not depend on ANSYS Fluent, OpenFOAM, or other existing CFD packages.
- **Limited hardware:** design models to run on a modest laptop. Start with small grids and reduced-order models; profile before optimising.
- **Evidence over complexity:** each project should have a clear physical question, documented assumptions, automated checks, and comparison against analytical or published reference results.
- **Finish before expanding:** complete a defensible minimum version before adding more physics.
- **Keep projects independently testable:** shared utilities are fine, but avoid turning the portfolio into one oversized codebase.

## Portfolio at a glance

| Project | Main capability | Priority |
|---|---|---|
| Plant flow network model | Steady hydraulic networks, gravity, conservation | Finish and validate |
| Pulsatile membrane pump | Transient pressure/flow and compliant boundaries | Next reduced-order model |
| Centrifugal heart pump | Rotating machinery and pump performance | Follow-on pump model |
| Code-based 2D CFD solver | Numerical solution of incompressible flow | Flagship simulation project |
| Compressible nozzle solver | Compressible flow, shocks, propulsion fundamentals | Later / rocketry specialisation |
| AutoFluent workflow project | Reproducible simulation orchestration and automation | Maintain as a supporting software project |

---

## 1. Plant flow network model

**Status:** Current project; approaching completion.

**Purpose:** Model steady water flow through a branched network, accounting for hydraulic resistance, pressure, and elevation.

### It should do

- Represent a branched network using a clear node/edge structure.
- Calculate node elevations from branch geometry, with geometry and topology kept explicit.
- Calculate Poiseuille resistance for each idealised cylindrical channel.
- Solve unknown node pressures or hydraulic potentials for specified boundary conditions.
- Calculate flow rate and direction through each branch.
- Check mass conservation at every internal junction.
- Plot the network with flow directions and distinguish upward/downward flow consistently.
- Report the assumptions and physical limits of the model.

### Verification and tests

- Compare a single horizontal pipe against the analytical Poiseuille solution.
- Check hydrostatic/elevation effects with zero-flow equilibrium cases.
- Check that setting gravity to zero removes elevation-dependent effects.
- Check symmetric branches for symmetric flow under symmetric conditions.
- Verify that the signed flow sum at each internal node is approximately zero.
- Test invalid radii, lengths, disconnected nodes, and inconsistent inputs.
- Where geometry is defined by branch angles, detect and report inconsistent geometric closure rather than silently accepting it.

### Deliverables

- Readable source code and automated tests.
- A short derivation of the governing equations.
- A plot of the network and branch flows.
- A validation note explaining error tolerances, assumptions, and limitations.

### Scope boundary

Do not add plant physiology, transpiration, cavitation, or a large dynamic model until the steady network model is verified. Equal pressure boundary conditions alone do not represent the full biological mechanism that pulls water through a plant.

**Definition of done:** benchmark cases pass, internal-node mass balance is demonstrated, and the README explains what the model can and cannot predict.

---

## 2. Pulsatile membrane pump

**Status:** Planned next reduced-order model.

**Purpose:** Simulate a pump whose flexible membrane changes chamber volume over time, producing pulsatile pressure and flow. Use it to extend the plant model's steady network ideas into transient fluid-system dynamics.

### It should do

- Prescribe membrane displacement or chamber volume as a function of time.
- Model chamber pressure and flow through inlet/outlet restrictions.
- Include a clear pressure–flow relation for each restriction, starting with a simple resistance model.
- Enforce volume/mass conservation as the chamber fills and empties.
- Represent inlet/outlet valves if needed, including their one-way behaviour.
- Produce time histories of membrane displacement, chamber pressure, inlet flow, and outlet flow.
- Explore the effect of actuation frequency, stroke volume, resistance, compliance, and valve timing.
- Compare the simulated mean flow and pulsation with simple analytical or limiting cases.

### Verification and tests

- Confirm conservation of volume over a full cycle, accounting for any net accumulation.
- Test zero actuation and zero-resistance/large-resistance limits where mathematically appropriate.
- Check periodic steady state by comparing successive cycles.
- Reduce to a simpler analytical model and compare results.
- Test time-step sensitivity and demonstrate that conclusions do not depend strongly on the chosen step.

### Deliverables

- A small, runnable model using Python/NumPy or another chosen language.
- Plots of pressure, flow, and membrane motion over several cycles.
- A parameter study showing how frequency and system resistance affect delivered flow and pressure ripple.
- A short explanation of the model's assumptions and limitations.

### Scope boundary

Start with a **lumped-parameter model** rather than a moving-mesh CFD or full fluid–structure interaction solver. A reduced-order model is computationally cheap and makes the coupling and conservation laws easier to verify.

**Engineering relevance:** transient flow systems, pressure ripple, compliant components, and fluid-system dynamics transfer to motorsport cooling/fuel systems and to rocket feed-system analysis. This is a fundamentals connection, not a claim that a membrane pump directly represents rocket turbopumps.

**Definition of done:** the model conserves volume, reaches repeatable periodic behaviour under defined conditions, and produces interpretable parameter-study results.

---

## 3. Centrifugal heart pump

**Status:** Planned / develop after the membrane pump.

**Purpose:** Model a centrifugal blood pump as a rotating-fluid-machinery problem. Emphasise the transferable engineering principles of head/pressure rise, flow rate, rotational speed, and losses.

### It should do

- Define a simplified impeller geometry and operating assumptions.
- Begin with a reduced-order model using angular momentum / Euler turbomachinery principles and velocity triangles where applicable.
- Estimate pressure rise or head as a function of rotational speed and flow rate.
- Include clearly identified efficiency or loss assumptions rather than presenting an ideal model as a real pump.
- Generate performance curves: head or pressure rise versus flow, and estimated efficiency versus operating point if the model supports it.
- Investigate how speed and flow rate shift the operating point.
- Compare predictions with published pump data or established pump-theory relationships when suitable data is available.

### Verification and tests

- Check units and dimensional consistency.
- Verify expected scaling trends with rotational speed under stated assumptions.
- Compare ideal limits with the corresponding analytical equations.
- Perform sensitivity analysis for slip, losses, geometry, and efficiency assumptions.
- Clearly separate calculated values from empirical inputs.

### Deliverables

- Documented equations and assumptions.
- A runnable model and performance plots.
- A comparison against a suitable reference dataset or published performance curve.
- A discussion of where a simplified model becomes unreliable.

### Scope boundary

Do not begin with detailed 3D impeller CFD. Build the analytical/reduced-order performance model first. Blood compatibility, shear exposure, haemolysis, and medical-device safety are specialised concerns and should not be claimed as validated by a simple hydraulic model.

**Engineering relevance:** rotating machinery, energy transfer, pressure generation, pump maps, and loss modelling connect to motorsport oil/coolant pumps and rocket turbomachinery fundamentals.

**Definition of done:** the model produces defensible performance trends, clearly states its loss assumptions, and is checked against theory or reference data.

---

## 4. Code-based 2D incompressible CFD solver

**Status:** Main future numerical-methods project.

**Purpose:** Implement a small flow solver from first principles. Demonstrate numerical understanding and verification without depending on commercial or third-party CFD software.

### Recommended progression

1. Implement and test the discrete operators and boundary conditions.
2. Solve a simple benchmark such as 2D lid-driven cavity flow.
3. Demonstrate convergence and compare velocity profiles with a published benchmark.
4. Extend to 2D flow past a circular cylinder at a suitable laminar Reynolds number.
5. Measure the vortex-shedding frequency and Strouhal number, then compare with reference data at a matching Reynolds number.

### It should do

- Solve the incompressible Navier–Stokes equations in two dimensions.
- State the discretisation, pressure–velocity coupling method, and boundary conditions explicitly.
- Begin with a small structured grid suitable for a modest laptop.
- Track residuals or other meaningful convergence indicators.
- Monitor mass conservation and relevant physical quantities.
- Support configurable grid size, time step, Reynolds number, and output interval.
- Save reproducible simulation settings and results.
- Visualise velocity, pressure, streamlines, and/or vorticity using ordinary plotting libraries.
- Separate the numerical solver from plotting and case configuration.

### Verification and tests

- Unit-test discrete derivative/divergence operators and boundary conditions.
- Check simple manufactured or analytical solutions where feasible.
- Compare a standard benchmark with published reference results.
- Conduct grid-refinement and time-step sensitivity studies.
- Check global and local mass conservation.
- Demonstrate that the vortex-shedding frequency is stable enough under refinement to support the stated conclusion.
- Profile runtime and memory use on the actual available machine.

### Deliverables

- Source code with automated tests.
- A derivation or numerical-methods note.
- Benchmark comparisons and convergence plots.
- A reproducible case configuration.
- A short engineering report explaining what the result means and where the model is limited.

### Scope boundary

Do not start with 3D turbulence, complex vehicle geometries, moving meshes, or full fluid–structure interaction. A small, verified 2D solver is the first milestone. A plausible-looking contour plot is not proof of correctness.

**Engineering relevance:** broadly applicable to aerodynamic wakes, separation, numerical modelling, and flow simulation. It is the portfolio's flagship demonstration of numerical implementation.

**Definition of done:** at least one benchmark is verified and documented, with a clear convergence study, before the solver is extended to more complex cases.

---

## 5. Compressible nozzle solver

**Status:** Later project; particularly useful for rocketry.

**Purpose:** Build a separate, lightweight solver for compressible flow, progressing from one-dimensional theory to numerical conservation laws.

### Recommended progression

1. Implement isentropic relations and choking calculations.
2. Solve quasi-one-dimensional converging–diverging nozzle flow.
3. Implement a 1D compressible Euler solver and test it on a shock tube.
4. Compare numerical results with the analytical shock-tube solution.
5. Explore nozzle operating regimes and back-pressure effects within the limits of the model.

### It should do

- Calculate Mach number, pressure, density, and temperature under clearly stated assumptions.
- Predict the choking condition and mass-flow behaviour.
- Solve the quasi-1D area–Mach relation.
- If implementing the Euler solver, use a conservative finite-volume method with a suitable numerical flux.
- Handle shocks without relying on a non-conservative discretisation.
- Include stability controls such as a CFL-based time-step restriction.
- Compare numerical predictions against analytical or published benchmarks.

### Verification and tests

- Recover known isentropic relations.
- Verify choking and nozzle trends against theory.
- Compare shock-tube density, pressure, and velocity profiles against the analytical solution.
- Test conservation and grid sensitivity.
- State assumptions such as ideal gas, inviscid flow, and whether heat transfer or chemical reactions are omitted.

### Deliverables

- A small, runnable codebase.
- Analytical-versus-numerical comparison plots.
- A validation report for the shock tube and nozzle cases.
- A clear explanation of the model's applicability to propulsion problems.

### Scope boundary

Do not claim to predict a real rocket engine's complete performance. A 1D ideal-gas model does not resolve viscous losses, combustion chemistry, injector behaviour, cooling, or three-dimensional nozzle flow.

**Engineering relevance:** compressible flow, shocks, conservation laws, and nozzle fundamentals directly support a rocketry-focused portfolio.

**Definition of done:** the analytical nozzle model and numerical shock-tube case are checked against established solutions before adding further complexity.

---

## 6. AutoFluent simulation workflow project

**Status:** Supporting software/automation project; separate from the code-based CFD solver.

**Purpose:** Demonstrate that simulation workflows can be configured, executed, checked, and reported reproducibly. The project may automate an existing solver when available, but the portfolio must not depend on access to Fluent.

### It should do

- Keep environment configuration separate from case/model configuration.
- Support a mock execution path that works without Fluent or HPC access.
- Record input configuration, run status, errors, and output locations.
- Include tests for successful execution and common failure paths.
- Provide a clear interface for adding solver-specific adapters later.
- Where a solver exposes relevant data, support post-processing and convergence checks without claiming a check was performed when the required data are absent.

### Verification and deliverables

- A clean setup guide.
- Tests runnable on a machine without Fluent.
- A reproducible mock pipeline example.
- Example run metadata and a concise result report.
- Clear distinction between implemented functionality, mock behaviour, and planned integrations.

### Scope boundary

Do not spend time adding abstractions or solver integrations before the mock pipeline and its tests are dependable. This project demonstrates software engineering and workflow automation; it does not replace verification of the physics solver itself.

**Definition of done:** a new user can run the mock pipeline, understand its outputs, and run the tests without proprietary software.

---

## Suggested order of work

1. **Close the plant model:** validate equations, geometry, mass balance, and limiting cases.
2. **Build the pulsatile membrane pump:** keep it lumped-parameter and transient.
3. **Develop the centrifugal pump model:** focus on theory, performance curves, and reference comparisons.
4. **Make the 2D incompressible solver the main technical project:** start with a benchmark and do not progress until it is verified.
5. **Add the compressible nozzle solver** if rocketry remains a priority.
6. **Maintain AutoFluent as a supporting tool**, keeping its scope constrained.

This order is a guide, not a requirement to complete every project before applying for roles. Prefer one rigorous, finished project over several ambitious unfinished ones.

## Definition of portfolio quality

For every project, aim to include:

- [ ] A concise physical problem statement.
- [ ] Governing equations and assumptions.
- [ ] Reproducible code and configuration.
- [ ] Automated tests for important invariants.
- [ ] Analytical or published benchmark comparisons.
- [ ] Sensitivity or convergence analysis where appropriate.
- [ ] Plots with units, labels, and interpretable scales.
- [ ] A limitations section that does not overstate accuracy.
- [ ] A short engineering conclusion explaining what the model enables you to decide.

## Career positioning

Present the portfolio as **computational modelling and fluid-systems engineering**. For motorsport applications, foreground flow modelling, wakes, cooling/oil-system applications, verification, and design trade-offs. For rocketry applications, foreground compressible flow, nozzle theory, transient feed-system fundamentals, and rotating machinery.

The goal is not to claim professional-grade CFD from a laptop. It is to demonstrate that you can derive a model, implement it, test it against evidence, understand its limitations, and use its results to make a defensible engineering argument.
