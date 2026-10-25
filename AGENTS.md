# AGENTS.md — Team Architecture & Engineering Guidelines

> **PROJECT**: `urban-flood-simulator`  
> **LEVEL**: First-Year Engineering Mathematics / Numerical Methods  
> **PURPOSE**: Small, explainable 2D urban waterlogging simulator.  
> **DISCLAIMER**: **Educational only, not a certified flood tool.**

---

## 1. Core Mathematical Model & Physics Rules

All code in this repository models 2D shallow ponding and overland cell-to-cell water balance using discrete finite-volume principles.

### Fundamental Equation
For every grid cell \((i, j)\):
$$\frac{\partial h}{\partial t} = R - I - D - \nabla \cdot \mathbf{q}$$

Discretized in time and space with explicit forward Euler:
$$h^{\text{new}}_{i,j} = h_{i,j} + \Delta t \left( R_{i,j} - I_{i,j} - D_{i,j} - \text{div}(\mathbf{q})_{i,j} \right)$$

Where:
- $h$ is ponded water depth (m).
- $z$ is bed elevation / ground topography (m).
- $\eta = h + z$ is total water surface elevation (hydraulic head, m).
- $R$ is rainfall intensity rate (m/s).
- $I$ is soil infiltration rate (m/s), bounded by soil infiltration capacity $I_{\max}$ and available water.
- $D$ is storm drain extraction rate (m/s), bounded by local drain capacity $D_{\max}$ and available water.
- $\mathbf{q} = (q_x, q_y)$ is the 2D discharge flux per unit width ($\text{m}^2/\text{s}$).

### Conservation and Flux Rules
1. **Shared Cell-Edge Fluxes**: Fluxes between neighboring cells $(i, j)$ and $(i, j+1)$ (or $(i+1, j)$) must be computed at cell faces and shared with equal and opposite signs ($q_{\text{face}}$ leaving one cell is identically $q_{\text{face}}$ entering the adjacent cell). This guarantees strict numerical mass conservation: $\sum \text{div}(\mathbf{q}) = 0$ on closed boundaries.
2. **Water Surface Gradient**: Flux across face $e$ between cell 1 and cell 2 is driven by water surface slope:
   $$S_e = -\frac{\eta_2 - \eta_1}{\Delta x} = -\frac{(h_2 + z_2) - (h_1 + z_1)}{\Delta x}$$
   Flow moves from higher water surface elevation to lower water surface elevation.
3. **No Negative Depths**: Numerical stability limits must enforce $h \ge 0$ at all times. If total sink terms $(I + D + \text{outflow}) \cdot \Delta t > h + R \Delta t$, fluxes must be limited proportionally.

---

## 2. Universal Code Constraints

1. **Strict SI Units Only**:
   - Depth $h$: meters ($\text{m}$)
   - Bed elevation $z$: meters ($\text{m}$)
   - Grid cell spacing $dx$: meters ($\text{m}$)
   - Time step $dt$: seconds ($\text{s}$)
   - Simulation time $t, t_{\text{end}}$: seconds ($\text{s}$)
   - Rates $R, I, D$: meters per second ($\text{m/s}$)
   - Cumulative volumes: cubic meters ($\text{m}^3$)
   - Grid shape: 2D arrays of shape `(rows, cols)`
   - Floating point precision: `np.float64` (`numpy.float64`) throughout.
2. **Type Hints Everywhere**:
   - Every function and method must declare comprehensive parameter types and return type hints (e.g., `np.ndarray`, `float`, `int`, `Tuple[...]`, `Optional[...]`).
3. **Maths Docstrings on Every Function**:
   - Every single function must have a docstring documenting:
     - The mathematical formula / physical definition implemented.
     - Units of all inputs and return values.
     - Discrete numerical approximation method used.
4. **Contract Discipline & Ownership**:
   - **Never change the contract without explicit team consensus.**
   - Each team member edits **ONLY** the files they own.
   - When introducing new helper methods, keep them private (`_helper`) within your module.

---

## 3. Team File Ownership

| Module / File | Primary Owner | Scope & Responsibility |
| :--- | :--- | :--- |
| `floodsim/config.py` | Lead Architect | `SimConfig`, `SimResult` dataclasses and validation |
| `floodsim/terrain.py` | Geomatics Specialist | Topography generation (`bowl`, `slope`, `random`), land cover, drain layouts |
| `floodsim/rainfall.py` | Hydrology Specialist | Synthetic storm hyetographs (constant, triangular, SI conversion) |
| `floodsim/operators.py`| Numerical Analyst | Discrete differential operators: $\nabla z$, $\text{div}(\mathbf{q})$, cell face fluxes |
| `floodsim/solver.py` | Numerical Analyst | 1D ODE reference solver, 2D time stepper, simulation loop |
| `floodsim/scenarios.py`| Urban Planner | Standard scenario definitions (baseline, improved drains, green infrastructure) |
| `floodsim/metrics.py` | Data Analyst | Water volume conservation, peak depth maps, flood duration hotspots |
| `floodsim/viz.py` | Visualization Dev | Matplotlib figures: depth colormaps, hydrographs, scenario comparisons |
| `floodsim/validation.py`| QA / Validation Lead | Analytical verification, mass balance tolerance checks |
| `floodsim/sensitivity.py`| QA / Validation Lead | Sensitivity runs on roughness, infiltration capacity, drain sizes |
| `app/dashboard.py` | Full-Stack UI Dev | Streamlit interactive web dashboard for real-time demonstration |
| `tests/` | QA Engineer | Pytest test suite, contract validation, regression tests |
| `docs/` | Documentation Lead| Educational derivations, mathematical documentation, user tutorials |

---

## 4. Coding Conventions

- PEP 8 compliant.
- No side effects on module import.
- Pure functions preferred where possible.
- NumPy vectorized operations preferred over raw Python loops for 2D cell updates.
