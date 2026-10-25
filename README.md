# Urban Flood Simulator (`urban-flood-simulator`)

> **First-Year Engineering Mathematics & Numerical Methods**  
> **A Small, Explainable 2D Urban Waterlogging Simulator**  
> ⚠️ **DISCLAIMER: Educational only, not a certified flood-warning tool.**  
> *This software is strictly an academic demonstration model for teaching discrete finite-volume methods and water balance principles. It is NOT certified, calibrated, or suitable for real-world emergency planning, life-safety warnings, or municipal flood hazard mitigation.*

---

## 🌊 Project Summary

`urban-flood-simulator` is a lightweight, mathematically transparent 2D numerical simulator designed for first-year engineering students to explore urban hydrology, discrete differential equations, and computational water balance.

The simulator models surface water accumulation (ponding) across a 2D Cartesian grid under rainfall forcing, soil infiltration, storm drain extraction, and gravity-driven overland routing.

### Key Engineering Concepts Explored
- **2D Discrete Water Balance**: Explicit forward Euler time integration on finite volumes:
  $$\frac{\partial h}{\partial t} = R - I - D - \nabla \cdot \mathbf{q}$$
  $$h^{\text{new}}_{i,j} = h_{i,j} + \Delta t \left( R_{i,j} - I_{i,j} - D_{i,j} - \text{div}(\mathbf{q})_{i,j} \right)$$
- **Topographic Water Routing**: Overland flux $\mathbf{q} = (q_x, q_y)$ driven by water-surface elevation ($\eta = h + z$).
- **Strict SI Physical Units**:
  - Depths ($h$) and bed elevations ($z$): meters ($\text{m}$)
  - Grid cell spacing ($dx$): meters ($\text{m}$)
  - Time step ($dt$) and simulation duration ($t_{\text{end}}$): seconds ($\text{s}$)
  - Fluxes and capacities ($R, I, D$): meters per second ($\text{m/s}$)
  - Cumulative water volumes: cubic meters ($\text{m}^3$)
- **Mass Conservation**: Strict mass balance verification via shared cell-edge fluxes with opposing signs ($\sum \text{div}(\mathbf{q}) = 0$ on closed boundaries).
- **Policy & Civil Interventions**: Fair comparison of Grey (drain capacity upgrades) vs. Green (permeable pavements / Sponge City) infrastructure.

---

## 📁 Repository Structure

```text
urban-flood-simulator/
├── AGENTS.md                  # Architectural contract, team rules, and ownership boundaries
├── README.md                  # Project overview, installation, usage, and folder structure
├── requirements.txt           # Python dependencies (numpy, matplotlib, pandas, pytest, streamlit)
├── main.py                    # CLI runner executing baseline, grey, and green scenarios
├── .gitignore                 # Standard Python gitignore rules
├── floodsim/                  # Core simulation engine package
│   ├── __init__.py            # Package initialization (exports SimConfig, SimResult)
│   ├── config.py              # SimConfig and SimResult dataclasses with validation
│   ├── terrain.py             # Topography generation (bowl, slope, random), land cover, drains
│   ├── rainfall.py            # Rainfall hyetographs: constant & triangular design storms
│   ├── scenarios.py           # Pre-configured policy scenarios (baseline, improved drains, green)
│   ├── operators.py           # Discrete spatial differential operators: gradient, divergence
│   ├── solver.py              # 0D ODE reference bucket, 2D time stepper, simulation loop
│   ├── metrics.py             # Volume conservation, boundary losses, peak depth, hotspot duration
│   ├── viz.py                 # Matplotlib visualization routines: depth colormap, hydrographs
│   ├── validation.py          # Mass balance tolerance verification and analytical benchmarks
│   └── sensitivity.py         # Parametric sensitivity sweeps
├── app/
│   ├── __init__.py
│   └── dashboard.py           # Interactive Streamlit dashboard for real-time demonstration
├── tests/
│   ├── __init__.py
│   ├── test_contract.py       # Dataclass structure and contract tests
│   └── test_inputs.py         # Input generation tests: shapes, non-negative rates, mm/hr conversion
└── docs/
    └── mathematics.md         # Full mathematical derivations, numerical notes, and formulas
```

---

## 🛠️ Installation

### Prerequisites
- Python 3.10+ (tested up to Python 3.15)
- Standard virtual environment tools (`venv` or `conda`)

### Step-by-Step Setup

1. **Navigate to project directory**:
   ```bash
   cd c:\Projects\urban-flood-simulator
   ```

2. **Create and activate a virtual environment**:
   - On Windows (PowerShell):
     ```powershell
     python -m venv .venv
     .venv\Scripts\Activate.ps1
     ```
   - On Windows (Command Prompt):
     ```cmd
     python -m venv .venv
     .venv\Scripts\activate.bat
     ```
   - On Linux / macOS:
     ```bash
     python -m venv .venv
     source .venv/bin/activate
     ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

---

## 🚀 How to Run

### 1. Run the Main CLI Scenario Suite
Execute all three standard scenarios (`baseline`, `improved_drains`, `green_infrastructure`) and call visualization routines:
```bash
python main.py
```
*Note: Uses `try/except` guards to run cleanly whether backend modules are stubs or fully implemented.*

### 2. Run the Pytest Test Suite
Verify input shapes, rate constraints, physical unit conversions, and architectural contracts:
```bash
pytest tests/ -v
```

To run individual test modules:
```bash
pytest tests/test_inputs.py -v
pytest tests/test_contract.py -v
```

### 3. Launch the Streamlit Web Dashboard
Launch the interactive web-based dashboard to visualize topography, storm profiles, and scenario comparisons:
```bash
streamlit run app/dashboard.py
```

---

## 🧩 Core Input Components

### 1. Terrain Generation (`floodsim/terrain.py`)
- `make_terrain(rows=30, cols=30, kind="bowl", seed=0)`:
  - `"bowl"`: Parabolic depression with low centre; water ponds naturally in the interior.
  - `"slope"`: Tilted plane with uniform gradient from upstream (row 0) to downstream.
  - `"random"`: Continuous smoothed micro-topography using 2D binomial/Gaussian filtering with seed reproducibility.
  - Returns: 2D `np.float64` array of bed elevations in meters ($\text{m}$).
- `make_land_cover(shape, paved_fraction, seed=0)`:
  - Assigns urban land cover and returns maximum infiltration capacity $I_{\max}$ in $\text{m/s}$.
  - Impermeable paved cells: $\sim 0$ to $10^{-7}\text{ m/s}$.
  - Permeable green/soil cells: $\sim 10^{-6}$ to $5 \times 10^{-6}\text{ m/s}$.
- `place_drains(shape, positions, capacity)`:
  - Generates drain intake grid $D_{\max}$ in $\text{m/s}$ with specified capacity at target inlet coordinates.

### 2. Rainfall Hyetographs (`floodsim/rainfall.py`)
- **Unit Conversion**:
  $$1\text{ mm/hr} = \frac{10^{-3}\text{ m}}{3600\text{ s}} \approx 2.78 \times 10^{-7}\text{ m/s}$$
- `constant_storm(intensity_mm_hr, duration_min, dt, t_end)`:
  - Rectangular storm profile with constant intensity during the storm and zero afterwards.
  - Returns 1D array of length $N = \lceil t_{\text{end}} / \Delta t \rceil$ with rates in $\text{m/s}$.
- `triangular_storm(peak_mm_hr, duration_min, dt, t_end)`:
  - Symmetric triangular design storm peaking at $t_{\text{half}} = t_{\text{dur}} / 2$.

### 3. Standard Scenarios (`floodsim/scenarios.py`)
All scenarios share identical 30x30 bowl terrain and identical rainfall forcing for strictly fair comparison:
1. `baseline()`: 80% paved impervious urban domain with conventional central drainage inlets.
2. `improved_drains()`: Same terrain and land cover, but doubles drain intake capacity ($2\times$) and expands inlet count.
3. `green_infrastructure()`: Same terrain and drainage, but retrofits ~60% of paved cells with permeable green surfaces ($I \sim 10^{-6} \text{ to } 5 \times 10^{-6}\text{ m/s}$).

---

## ⚠️ Scope Disclaimer

> **ACADEMIC AND EDUCATIONAL USE ONLY**  
> This simulator was developed solely as an educational teaching aid for engineering students studying numerical differential equations, conservation laws, and computational fluid mechanics. It is **not** a certified hydraulic modeling package (such as US Army Corps HEC-RAS, EPA SWMM, or TUFLOW). Do not use this software for flood risk zoning, real-world emergency evacuation planning, structure design, or life-safety flood hazard assessments.
