# Urban Flood Simulator (`urban-flood-simulator`)

> **First-Year Engineering Mathematics Project**  
> **A Small, Explainable 2D Urban Waterlogging Simulator**  
> ⚠️ **Educational only, not a certified flood tool.**

---

## 🌊 Overview

`urban-flood-simulator` is an educational numerical simulation framework designed for engineering students to explore:
- 2D finite-volume water balance on uniform Cartesian grids.
- Explicit time-stepping schemes (Forward Euler) and numerical stability (CFL conditions).
- Topographic water routing driven by water-surface elevation ($\eta = h + z$).
- Urban drainage physics: rainfall deposition, soil infiltration, and storm drain capacity.
- Mass conservation verification and policy interventions (Grey vs. Green infrastructure).

---

## 📐 Mathematical Formulation

At every grid cell $(i, j)$ with dimensions $\Delta x \times \Delta x$ (area $\Delta x^2$):

$$h^{\text{new}}_{i,j} = h_{i,j} + \Delta t \left( R_{i,j} - I_{i,j} - D_{i,j} - \text{div}(\mathbf{q})_{i,j} \right)$$

### Strict Physical Units (SI Only)
- **Depth $h$**: meters ($\text{m}$)
- **Bed elevation $z$**: meters ($\text{m}$)
- **Grid cell spacing $dx$**: meters ($\text{m}$)
- **Time step $dt$**: seconds ($\text{s}$)
- **Rates ($R, I, D$)**: meters per second ($\text{m/s}$)
- **Volume $V$**: cubic meters ($\text{m}^3$)
- **Data type**: `np.float64` 2D arrays of shape `(rows, cols)`

---

## 📁 Repository Structure

```text
urban-flood-simulator/
├── AGENTS.md                  # Team code rules, ownership, and contract guidelines
├── README.md                  # Project overview and instructions
├── requirements.txt           # Dependencies: numpy, matplotlib, pandas, pytest, streamlit
├── main.py                    # CLI entrypoint demonstrating contract instantiation
├── .gitignore                 # Python git ignore rules
├── floodsim/                  # Core simulation package
│   ├── __init__.py            # Package init (exposes SimConfig, SimResult, __version__)
│   ├── config.py              # SimConfig and SimResult dataclasses
│   ├── terrain.py             # DEM generation (bowl, slope, random), land cover, drains
│   ├── rainfall.py            # Hyetographs: constant & triangular design storms
│   ├── scenarios.py           # Pre-configured baseline, grey, and green scenarios
│   ├── operators.py           # Discrete differential operators: gradient, divergence
│   ├── solver.py              # 0D ODE single cell, 2D time step, simulation loop
│   ├── metrics.py             # Volume conservation, peak depth, hotspot duration
│   ├── viz.py                 # Matplotlib figures: depth map, timeseries hydrograph
│   ├── validation.py          # Mass conservation tests and analytical benchmarks
│   └── sensitivity.py         # Parametric sweeps (k_flow, roughness, infiltration)
├── app/
│   ├── __init__.py
│   └── dashboard.py           # Streamlit interactive exploratory dashboard
├── tests/
│   ├── __init__.py
│   └── test_contract.py       # Pytest suite validating dataclasses and contracts
└── docs/
    └── mathematics.md         # Detailed mathematical derivations and numerical notes
```

---

## 📜 Development Guidelines & Contract

1. **SI Units Exclusively**: No conversion factors inside numerical routines; all inputs are pre-converted to meters, seconds, m/s, and m³.
2. **Mass Conservation**: Fluxes between neighboring cells use shared cell-edge fluxes with opposing signs so $\sum \text{div}(\mathbf{q}) = 0$ on closed domains.
3. **Non-Negativity**: Water depth $h$ is strictly non-negative ($h \ge 0$).
4. **Docstrings & Types**: Every function documents the mathematical formula, inputs/outputs with physical units, and includes full Python type hints.
5. **Team Ownership**: Each team member edits **only** the files they own as outlined in [AGENTS.md](AGENTS.md). Contracts are strictly immutable without team consensus.

---

## 🚀 Getting Started

### 1. Installation
```bash
python -m venv .venv
# On Windows:
.venv\Scripts\activate
# On Linux/macOS:
source .venv/bin/activate

pip install -r requirements.txt
```

### 2. Verify Architectural Contracts
```bash
pytest tests/ -v
```

### 3. Run Sample CLI Entrypoint
```bash
python main.py
```

### 4. Launch Streamlit Web Dashboard
```bash
streamlit run app/dashboard.py
```
