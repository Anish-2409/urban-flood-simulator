# 🌊 Urban Flood & Drainage Flow Simulator

**A small, explainable Python simulator that shows where rain water collects in a neighbourhood, and tests whether bigger drains or green ground reduce the flooding.**

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Dashboard-Streamlit-ff4b4b)](https://urban-flood-simulator.streamlit.app/)
[![Tests](https://img.shields.io/badge/Tests-pytest-green)](#-testing-and-validation)
[![Use](https://img.shields.io/badge/Use-Educational%20only-orange)](#-limitations-and-responsible-use)

> **First-year CSE Engineering Mathematics project**
> Differential equations, vector calculus, integrals and numerical methods, applied to one real problem: urban waterlogging.

### 👉 [Try the live dashboard](https://urban-flood-simulator.streamlit.app/)

> ⚠️ **Educational use only.** This is a learning model. It is **not** a flood-warning system and must not be used for emergency decisions or real infrastructure design. See [Limitations](#-limitations-and-responsible-use).

---

## 📖 Table of contents

1. [The problem](#-the-problem)
2. [How it works](#-how-it-works)
3. [The mathematics](#-the-mathematics)
4. [Scenarios compared](#-scenarios-compared)
5. [Example result](#-example-result)
6. [Getting started](#-getting-started)
7. [Project structure](#-project-structure)
8. [Testing and validation](#-testing-and-validation)
9. [Limitations and responsible use](#-limitations-and-responsible-use)
10. [Team](#-team)
11. [References](#-references)

---

## 🌧️ The problem

During heavy rain, water reaches the streets faster than drains can remove it. Low ground, paved surfaces that do not soak up water, and small or blocked drains make waterlogging worse.

**Our engineering question:** given rainfall, terrain and a simple drainage system, where does water collect, how does its depth change over time, and which fixes might help?

## ⚙️ How it works

The neighbourhood is split into a **30 × 30 grid of small squares**, like a chocolate bar. The computer tracks the water depth in every square and repeats four steps in small time steps:

1. **Rain** adds water to each square.
2. **Infiltration**: soil and grass soak some water up, while paved ground soaks up almost none.
3. **Drains** remove water, up to a fixed capacity.
4. **Flow**: water moves downhill to lower neighbouring squares. Each edge shares one flow value, so water leaving one square exactly enters the next.

After thousands of steps, the grid shows where the flood spots are.

## 🧮 The mathematics

The core rule is a water balance (the continuity equation):

```
∂h/∂t + ∇·q = R − I − D
```

where `h` is water depth, `q` is the water flow, `R` is rainfall, `I` is infiltration and `D` is drain removal. On the grid, one time step looks like this:

```
h_new[i,j] = h[i,j] + dt * ( R − I − D − div(q)[i,j] )
```

| Concept from the course | What it does in this project | Where |
|---|---|---|
| **ODE** | Water balance for one square: `dH/dt = R − I − D` | `floodsim/solver.py` |
| **PDE** | The same balance across the whole map and over time | `floodsim/solver.py` |
| **Partial derivatives and gradient** | Water-surface slope decides which way water flows (downhill) | `floodsim/operators.py` |
| **Divergence** | Net flow leaving a square, used in the water balance | `floodsim/operators.py` |
| **Double integral** | Total water volume: `V ≈ Σ h · dx²` | `floodsim/metrics.py` |
| **Line integral** | Flow rate crossing the map boundary | `floodsim/metrics.py` |
| **Numerical methods** | Forward Euler time stepping and finite differences on a grid | `floodsim/solver.py` |

All quantities use SI units: depth and elevation in metres, time in seconds, rates in m/s, volumes in m³.

The full derivations are in [`docs/`](docs/).

## 🏙️ Scenarios compared

All three scenarios use the **same terrain and the same rainfall**, so the comparison is fair.

| Scenario | Idea |
|---|---|
| **Baseline** | A mostly paved area with a few small drains |
| **Improved drains** | The same area with larger drain capacity and more inlets |
| **Green infrastructure** | The same area, but part of the paving is replaced by soil or grass that soaks up water |

## 📊 Example result

The dashboard prints a water balance for each run. One example run:

| Item | Volume | Share of rain |
|---|---|---|
| Total rain | 126.1 m³ | 100% |
| Soaked into the soil | 29.4 m³ | 23.3% |
| Taken by drains | 1.4 m³ | 1.1% |
| Left through the boundary | 0.0 m³ | 0% |
| **Still on the surface at the end** | **95.3 m³** | **75.5%** |
| Peak depth anywhere on the map | 19.95 cm | |
| **Water-balance error** | **about 4 × 10⁻¹² m³** | effectively zero |

The error is practically zero because the solver never creates or loses water. Switch scenarios in the dashboard to see how the "still on the surface" number changes.

<!--
Add screenshots after saving them in docs/images/, then remove this comment:

![Depth map](docs/images/depth_map.png)
![Scenario comparison](docs/images/scenario_comparison.png)
-->

## 🚀 Getting started

### Option 1: Use the live dashboard

Open **https://urban-flood-simulator.streamlit.app/**. If the app is asleep, click the button to wake it up and wait a few seconds.

### Option 2: Run it on your own computer

You need **Python 3.10 or newer** (3.12 is recommended) and Git.

```bash
# 1. Get the code
git clone https://github.com/Anish-2409/urban-flood-simulator.git
cd urban-flood-simulator

# 2. Create and activate a virtual environment
python -m venv .venv
.venv\Scripts\activate            # Windows (PowerShell)
# source .venv/bin/activate       # macOS / Linux

# 3. Install the packages
pip install -r requirements.txt
```

If PowerShell says that running scripts is disabled, run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once and try again.

### Run it

```bash
# Interactive dashboard (opens in your browser at http://localhost:8501)
streamlit run app/dashboard.py

# Command-line run of all three scenarios
python main.py

# Water-balance sanity check
python check_sanity.py

# All automated tests
pytest
```

In the dashboard you can change the rainfall intensity, storm duration, drain capacity and paved fraction, switch between the three scenarios, and view the depth map, the volume chart and the water balance table.

## 📁 Project structure

```
urban-flood-simulator/
├── app/
│   └── dashboard.py        # Streamlit web dashboard
├── floodsim/               # The simulation package
│   ├── config.py           #   SimConfig and SimResult data classes
│   ├── terrain.py          #   Terrain, land cover and drain layout
│   ├── rainfall.py         #   Constant and triangular design storms
│   ├── scenarios.py        #   Baseline, improved drains, green infrastructure
│   ├── operators.py        #   Gradient and divergence
│   ├── solver.py           #   Water balance and time stepping
│   ├── metrics.py          #   Volume, boundary flow, hotspot duration
│   ├── viz.py              #   Plots
│   ├── validation.py       #   Water-balance checks
│   └── sensitivity.py      #   Rainfall / infiltration / drain sensitivity
├── tests/                  # pytest test suite
├── docs/                   # Mathematics notes and report
├── AGENTS.md               # Shared rules and file ownership for the team
├── check_sanity.py         # Quick water-balance check
├── main.py                 # Command-line runner
└── requirements.txt        # Python packages
```

## ✅ Testing and validation

We do not just check that the program runs. We check that it is **physically correct**:

- **Water balance:** rain in − infiltration − drains − boundary outflow = change in stored water (tested for all three scenarios).
- **No rain, no water:** a dry map with no rainfall stays dry.
- **No negative depth:** water depth never drops below zero.
- **Closed boundary:** with no drains and no infiltration, the water stored equals the rain that fell.
- **Bowl terrain:** water collects in the low centre.
- **Better drains help:** the improved-drains scenario never has a higher peak volume than the baseline.
- **ODE check:** the single-square solver is compared with its exact analytic solution.

Run them with `pytest`. At the time of writing, all 35 tests pass.

## ⚠️ Limitations and responsible use

- The terrain is generated, not measured. Real curbs, road crowns, underpasses and small channels are not modelled.
- Flow uses a simplified rule, not the full shallow-water equations.
- Real drain performance depends on pipe layout, blockages and water levels downstream, which are simplified here.
- The model has **not** been calibrated or validated against a real flood event.
- Results are **scenario estimates, not street-level forecasts**.
- **Do not use this software** for emergency decisions, evacuation planning, or certified infrastructure design.

**Possible next steps:** use real elevation data (Copernicus DEM) and rainfall data (IMD, NASA GPM IMERG), compare with an established model such as EPA SWMM, and add a cost-versus-benefit ranking of the interventions.

## 👥 Team

First-year CSE team of four. Each member owned a separate part of the code.

| Name | GitHub | Role |
|---|---|---|
| Anish P | [@Anish-2409](https://github.com/Anish-2409) | Inputs and scenarios: terrain, rainfall, scenario setup, project lead |
| Vishwa  | [@Vishwa-0403](https://github.com/Vishwa-0403) | Core solver: operators and time stepping |
| Chaitanya DeshPande | [@chaitanyadeshpande2007-lang](https://github.com/chaitanyadeshpande2007-lang) | Metrics and visualisation |
| Rahul Rathi | [@rahulrathi14022007](https://github.com/rahulrathi14022007) | Testing, validation, sensitivity analysis and dashboard |

We worked with Git branches and pull requests, and used AI coding assistants to help write code, which we reviewed and tested ourselves.

## 📚 References

- EPA, [Storm Water Management Model (SWMM)](https://www.epa.gov/water-research/storm-water-management-model-swmm)
- US Army Corps of Engineers, [HEC-RAS Hydraulic Equations](https://www.hec.usace.army.mil/confluence/rasdocs/hecras/latest/technical-reference/hydraulic-equations)
- Guo et al. (2021), [Urban surface water flood modelling: a review](https://hess.copernicus.org/articles/25/2843-2021.html)
- [IMD rainfall, National Water Data Portal](https://nwdp.nwic.gov.in/dataset/rainfall-daily-imd)
- [NASA GPM IMERG](https://gpm.nasa.gov/data/imerg)
- [Copernicus DEM](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)
- [NRSC thematic datasets / Bhuvan](https://www.nrsc.gov.in/nrscnew/Dataproducts_Thematic_overview.php)
- [Copernicus Emergency Management Service](https://emergency.copernicus.eu/data/)

---

*Educational project. Not a certified flood-warning or engineering design tool.*
