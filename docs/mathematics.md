# Mathematical Formulation of the 2D Waterlogging Model

> **Educational Project**: First-Year Engineering Mathematics  
> **Notice**: Educational only; not a certified or regulatory flood risk modeling tool.

---

## 1. Governing Equation: The Cell Water Balance

At every spatial location $(x, y)$ on an urban ground surface, the depth of ponded water $h(x, y, t)$ (in meters) satisfies the local continuity equation:

$$\frac{\partial h}{\partial t} + \nabla \cdot \mathbf{q} = R(t) - I(x, y, t) - D(x, y, t)$$

Rearranging for depth rate of change:

$$\frac{\partial h}{\partial t} = R(t) - I(x, y, t) - D(x, y, t) - \nabla \cdot \mathbf{q}$$

### Physical Terms and Units (Strict SI)
| Symbol | Quantity | Units | Description |
| :--- | :--- | :--- | :--- |
| $h$ | Water depth | $\text{m}$ | Height of ponded water column above ground bed |
| $z$ | Ground elevation | $\text{m}$ | Topographic bed elevation above datum |
| $\eta = h + z$ | Water surface elevation | $\text{m}$ | Hydraulic head driving gravity flow |
| $R$ | Rainfall intensity | $\text{m/s}$ | Volumetric water input rate per unit ground area |
| $I$ | Soil infiltration rate | $\text{m/s}$ | Rate of water absorbed into permeable soil |
| $D$ | Drain extraction rate | $\text{m/s}$ | Rate of water collected into storm sewers |
| $\mathbf{q} = (q_x, q_y)$ | Unit discharge flux | $\text{m}^2/\text{s}$ | Volumetric flow rate per unit width of cell face |
| $\nabla \cdot \mathbf{q}$ | Flux divergence | $\text{m/s}$ | Net lateral outflow per unit ground area |

---

## 2. Discrete Finite-Volume Formulation

The urban landscape is discretized into a uniform 2D Cartesian grid of square cells, each of size $\Delta x \times \Delta x$ (area $\Delta A = \Delta x^2$).

Integrating the continuity equation over cell $(i, j)$ from time $t$ to $t + \Delta t$ via explicit Forward Euler:

$$h^{\text{new}}_{i, j} = h_{i, j} + \Delta t \left( R_{i, j} - I_{i, j} - D_{i, j} - \text{div}(\mathbf{q})_{i, j} \right)$$

### Shared Cell-Edge Fluxes (Mass Conservation)
To guarantee strict numerical mass conservation, inter-cell fluxes are evaluated at cell faces and shared identically between neighbors:

$$\text{div}(\mathbf{q})_{i, j} = \frac{q_{x, i, j+1/2} - q_{x, i, j-1/2}}{\Delta x} + \frac{q_{y, i+1/2, j} - q_{y, i-1/2, j}}{\Delta x}$$

The flux leaving cell $(i, j)$ into $(i, j+1)$ across face $(i, j+1/2)$ is:

$$q_{x, i, j+1/2} = k_{\text{flow}} \cdot h_{\text{edge}} \cdot S_{x, i, j+1/2}$$

where:
- Hydraulic slope: $S_{x, i, j+1/2} = -\frac{\eta_{i, j+1} - \eta_{i, j}}{\Delta x} = \frac{\eta_{i, j} - \eta_{i, j+1}}{\Delta x}$
- Upwind depth: $h_{\text{edge}} = h_{i, j}$ if $\eta_{i, j} > \eta_{i, j+1}$ else $h_{i, j+1}$
- If $h_{\text{edge}} = 0$, $q = 0$ (dry cells do not transfer water).

Summing $\sum_{i,j} \text{div}(\mathbf{q})_{i, j}$ over a domain with closed boundaries cancels every interior face flux identically in pairs, preserving total water mass to machine precision.

---

## 3. Physical Sink Constraints & Depth Non-Negativity

Water depth $h$ cannot be negative. The sink terms $(I, D)$ and outflows are physically constrained by available water:

1. **Infiltration Limit**:
   $$I_{i, j} = \min\left(I_{\max, i, j},\, \frac{h_{i, j}}{\Delta t} + R_{i, j}\right)$$
2. **Drainage Limit**:
   $$D_{i, j} = \min\left(D_{\max, i, j},\, \max\left(0,\, \frac{h_{i, j}}{\Delta t} + R_{i, j} - I_{i, j}\right)\right)$$
3. **Depth Non-Negativity**:
   $$h^{\text{new}}_{i, j} = \max\left(0,\, h_{i, j} + \Delta t(R - I - D - \text{div}(\mathbf{q}))\right)$$

---

## 4. Stability Condition (CFL)

Explicit time integration requires the Courant-Friedrichs-Lewy (CFL) condition:

$$\Delta t \le \frac{\Delta x}{2 \cdot \max |\mathbf{u}|}$$

For typical urban shallow ponding with $\Delta x = 2\,\text{m}$ and flow velocities $u \sim 0.1 - 0.5\,\text{m/s}$, $\Delta t \approx 1.0\,\text{s}$ provides stable, non-oscillatory time integration.
