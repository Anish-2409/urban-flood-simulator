# Urban Flood Simulator — Engineering Mathematics Technical Report

> **DISCLAIMER: EDUCATIONAL USE ONLY**  
> This simulator was developed exclusively for first-year engineering mathematics and numerical methods instruction.  
> **The outputs are scenario estimates and comparative policy illustrations, not reliable street-level forecasts or certified flood risk assessments.** This software must not be used for emergency planning, infrastructure certification, or regulatory flood hazard zoning.

---

## 1. Problem Statement

Modern urban catchments face escalating pluvial flood and waterlogging risks due to two converging factors:
1. **Urban Imperviousness**: Replacement of natural permeable soils and vegetated cover with asphalt, concrete, and rooftops prevents natural infiltration, converting up to 80–90% of rainfall directly into overland runoff.
2. **Extreme Convective Storms**: Short-duration, high-intensity precipitation events routinely exceed conventional municipal storm sewer design capacities (typically engineered for 5- to 10-year return intervals).

When rainfall intensity surpasses both local soil infiltration capacity and sewer inlet intake rates, excess water ponds in topographic depressions and street corridors, causing pedestrian hazards, vehicular stalling, basement flooding, and economic disruption.

The objective of this project is to develop an explainable, computationally transparent 2D discrete finite-volume numerical model that enables engineering students to evaluate:
- The fundamental continuity equation governing shallow surface ponding.
- The hydrological effectiveness of conventional grey infrastructure (expanded storm drains) versus nature-based green infrastructure (permeable pavements and Sponge City retrofits).

---

## 2. Conceptual Physical Model

The simulation domain is discretized into a uniform Cartesian 2D grid of square cells of size $\Delta x \times \Delta x$ ($\text{m}^2$). The physical state at each cell $(i, j)$ is characterized by:
- **Bed Topography $z(x, y)$** ($\text{m}$): Fixed ground elevation relative to a vertical datum.
- **Water Depth $h(x, y, t)$** ($\text{m}$): Thickness of the standing surface water layer.
- **Hydraulic Head $\eta(x, y, t) = h + z$** ($\text{m}$): Total water surface elevation driving overland gravity flow.

```
          Water Surface Elevation eta = h + z
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  v
      |                                         | | Water Depth h
      |               Ponded Water              | |
======+=========================================+== v =================== Ground Surface z
      |                                         |
      |             Subsurface / Soil           |
      +-----------------------------------------+
```

Water mass enters the surface via precipitation $R(t)$ and exits through two vertical sinks:
- Soil infiltration $I(x, y, t)$ bounded by maximum infiltration capacity $I_{\max}$.
- Storm sewer drainage $D(x, y, t)$ bounded by inlet capture capacity $D_{\max}$.

Excess water is routed horizontally between adjacent 4-connected cells along steepest water surface slopes.

---

## 3. Mathematical Equations and Discretization

### 3.1 Depth-Integrated Continuity Equation
The conservation of water mass over a 2D surface element is governed by the shallow water continuity equation:

$$\frac{\partial h}{\partial t} = R(t) - I(x, y, t) - D(x, y, t) - \nabla \cdot \mathbf{q}$$

where:
- $h$ is ponded water depth ($\text{m}$).
- $R$ is rainfall intensity ($\text{m/s}$).
- $I$ is soil infiltration rate ($\text{m/s}$).
- $D$ is storm drain intake rate ($\text{m/s}$).
- $\mathbf{q} = (q_x, q_y)$ is horizontal discharge flux per unit width ($\text{m}^2/\text{s}$).
- $\nabla \cdot \mathbf{q} = \frac{\partial q_x}{\partial x} + \frac{\partial q_y}{\partial y}$ is the flux divergence ($\text{m/s}$).

### 3.2 Time Integration via Explicit Forward Euler
Discretizing time into uniform increments $\Delta t$:

$$h^{n+1}_{i, j} = h^n_{i, j} + \Delta t \left( R^n - I^n_{i, j} - D^n_{i, j} - \text{div}(\mathbf{q})^n_{i, j} \right)$$

### 3.3 Inter-Cell Fluxes and Strict Mass Conservation
To ensure machine-precision conservation of water volume, fluxes are evaluated at cell faces and shared identically between neighboring cells:

$$\text{div}(\mathbf{q})_{i, j} = \frac{q_{x, i, j+1/2} - q_{x, i, j-1/2}}{\Delta x} + \frac{q_{y, i+1/2, j} - q_{y, i-1/2, j}}{\Delta x}$$

The discharge across vertical face $(i, j+1/2)$ between cell $(i, j)$ and $(i, j+1)$ is driven by hydraulic head slope:

$$q_{x, i, j+1/2} = k_{\text{flow}} \cdot h_{\text{edge}} \cdot \left( \frac{\eta_{i, j} - \eta_{i, j+1}}{\Delta x} \right)$$

where:
- $h_{\text{edge}} = h_{i, j}$ if $\eta_{i, j} > \eta_{i, j+1}$ else $h_{i, j+1}$ (upwind donor depth).
- $k_{\text{flow}}$ is the dimensionless conveyance scaling parameter ($k_{\text{flow}} \le 0.25$ for 2D forward Euler numerical stability).

Because $q_{x, i, j+1/2}$ leaving cell $(i, j)$ is identically the flux entering cell $(i, j+1)$, summing $\sum_{i,j} \text{div}(\mathbf{q}) \Delta x^2$ over the closed domain cancels all interior face fluxes in pairs:

$$\sum_{i, j} \text{div}(\mathbf{q})_{i, j} \Delta x^2 = 0 \quad \text{(Closed Boundary)}$$

### 3.4 Physical Sink Bounds and Depth Non-Negativity
To prevent unphysical negative depths ($h < 0$), sinks and outflows are sequentially limited by available water:

1. **Infiltration**: $I_{\text{actual}} = \min(I_{\max} \Delta t, \, h^n + R \Delta t)$
2. **Drainage**: $D_{\text{actual}} = \min(D_{\max} \Delta t, \, h^n + R \Delta t - I_{\text{actual}})$
3. **Adaptive Outflow Limiter**: For cell $i$, total proposed outgoing flux $\sum_k F_{\text{out}, i \to k}$ is scaled by:
   $$\alpha_i = \min\left(1.0, \, \frac{h_{\text{avail}, i}}{\sum_k F_{\text{out}, i \to k}}\right)$$
   guaranteeing $h^{n+1} \ge 0$ unconditionally.

---

## 4. Scenario Results and Comparative Analysis

Three standard policy interventions were simulated on a $30 \times 30$ grid ($\Delta x = 2.0\text{ m}$, area $3600\text{ m}^2$) with a central bowl depression topography under a design storm of $50\text{ mm/hr}$ for 30 minutes ($90.07\text{ m}^3$ total rainfall):

| Metric | Baseline (80% Paved) | Improved Drains (2x Capacity) | Green Infrastructure (Sponge City) |
| :--- | :---: | :---: | :---: |
| **Total Infiltration Volume ($m^3$)** | 6.51 | 6.51 | **38.42** |
| **Total Drained Volume ($m^3$)** | 0.40 | **0.80** | 0.35 |
| **Final Surface Ponding ($m^3$)** | 83.16 | 82.76 | **51.30** |
| **Peak Inundation Depth ($cm$)** | 20.44 | 20.12 | **14.28** |
| **Peak Surface Volume ($m^3$)** | 83.16 | 82.76 | **51.30** |
| **Mass Balance Error ($m^3$)** | $< 10^{-12}$ | $< 10^{-12}$ | $< 10^{-12}$ |

### Key Findings
1. **Grey Infrastructure (Improved Drains)**:
   - Doubling drain capacity and adding inlets directly at low-lying cells effectively increases local drainage rate and lowers peak surface volume.
   - However, point drains are localized sinks; overland runoff must first travel across the surface to reach the inlets, meaning peripheral cells still experience surface runoff.
2. **Green Infrastructure (Sponge City Retrofits)**:
   - Increasing permeable surface area to 60% of paved cells increases soil absorption by nearly 6x (from $6.5\text{ m}^3$ to $38.4\text{ m}^3$).
   - Peak inundation depth is reduced by **30%** (from $20.4\text{ cm}$ to $14.3\text{ cm}$), mitigating flood risk catch-wide at the source before water concentrates in depressions.

---

## 5. Sensitivity Analysis

A parametric one-at-a-time (OAT) sweep varied physical inputs by $\pm 50\%$:
- **Rainfall Intensity**: Demonstrates nearly linear scaling with peak water depth and total ponding. A $50\%$ increase in storm intensity increases peak inundation depth proportionally.
- **Infiltration Capacity**: High sensitivity in permeable areas; retrofitting high-traffic impervious surfaces provides substantial attenuation.
- **Drainage Capacity**: Exhibits diminishing returns if inlet density is low or if inlets become water-starved during short-duration peak surges.

---

## 6. Engineering Limitations & Educational Disclaimers

> [!WARNING]
> This model is an educational demonstration tool and must not be used for life-critical engineering design or regulatory flood studies.

The model incorporates deliberate simplifying approximations:
1. **Simplified Flow Physics**: Overland routing uses a diffusive head-slope approximation rather than full 2D Saint-Venant (shallow water) momentum equations. Advective acceleration ($\mathbf{u} \cdot \nabla \mathbf{u}$), bed shear stress tensors (Manning's $n$), and turbulent dissipation are not modeled.
2. **No Pipe Network Hydraulics**: Storm sewers are represented as sink boundary conditions with constant intake rates. Pressurized pipe flow, surcharging, backwater effects, and manhole overflow are omitted.
3. **Unsaturated Zone Simplification**: Soil absorption uses an effective bucket capacity rate rather than solving non-linear 1D Richards' equation or Green-Ampt infiltration with soil suction head and variable moisture content.
4. **Coarse Resolution**: Sub-meter curb heights, building walls, roadside gutters, and micro-topographic features are averaged across the cell area $\Delta x^2$.

---

## 7. References

1. Chow, V. T., Maidment, D. R., & Mays, L. W. (1988). *Applied Hydrology*. McGraw-Hill.
2. LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*. Cambridge University Press.
3. Butler, D., Digman, C. J., Makropoulos, C., & Davies, J. W. (2018). *Urban Drainage* (4th ed.). CRC Press.
4. Bates, P. D., Horritt, M. S., & Fewtrell, T. J. (2010). A simple inertial formulation of the shallow water equations for efficient two-dimensional flood inundation modelling. *Journal of Hydrology*, 387(1-2), 33-45.
