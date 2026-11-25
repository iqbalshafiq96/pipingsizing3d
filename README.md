# Piping Hydraulics & Temperature Drop Builder (3D)

Lay out a piping network in 3D (East/West/North/South/Up/Down + sloped runs) in a single
`index.html`, then calculate hydraulics and temperature drop with Python, hosted on Streamlit.

```
repository/
└── pipingsizing/
    ├── app.py           # Streamlit host - bridges index.html <-> Python, report, 3D true-scale plot, downloads
    ├── index.html       # Interactive 3D builder: inputs, canvas, validation, results, profile chart
    ├── hydraulics.py    # Units, pipe data, fluids, friction factor, Darcy-Weisbach, Bernoulli, 3D coordinates
    ├── thermal.py       # Heat loss through insulation, Joule-Thomson / friction heating
    ├── requirements.txt
    └── README.md
```

## Run locally
```bash
cd pipingsizing
pip install -r requirements.txt
streamlit run app.py
```

## Deploy on Streamlit Community Cloud
1. Push the repo to GitHub.
2. New app -> pick the repo -> **Main file path:** `pipingsizing/app.py`.

## 3D builder
- Axes: **E = +X, N = +Y, Z = elevation (Up)**. The inlet sits at (0, 0, inlet elevation).
- Pipes follow the current flow direction. **Elbows** turn to any of the 4 perpendicular directions,
  **tees** branch to any perpendicular direction.
- Horizontal pipes take an optional **Δz** (rise + / fall −) for sloped lines; vertical pipes use Δz = ±L.
- Canvas: drag = rotate, Shift/right-drag = pan, wheel = zoom, arrow keys = rotate/tilt.
  Presets: Iso, Plan, Front (looking North), Side (looking West).
- **Schematic lengths** (log-scaled, readable) or **True scale** (metres; fittings drawn with a small gap).
- Streamlit adds a true-scale Plotly 3D route coloured by pressure; every segment row reports its end E/N/Z.

## How the pieces talk
`app.py` registers `index.html` as a Streamlit custom component (`declare_component(path=<this folder>)`),
no npm build. `index.html` sends the network JSON on **Calculate**, `hydraulics.solve()` (which calls
`thermal.py` at every pipe step) returns results, and the canvas re-colours by pressure.
Opening `index.html` directly still lets you build/export a network; calculation needs Python.

## Rules enforced
- Size mismatch (e.g. 6" -> 4" without a reducer) -> red alert, Add/Calculate blocked.
- Every branch end must be closed by an **Outlet**; outlet flows must close the mass balance (one may be blank).
- |Δz| must be smaller than the pipe length.
- **ΔP not feasible**: vacuum, below vapour pressure (liquids), choked gas (Mach ≥ 1), or below an outlet's minimum P.

## Method & assumptions
- Single-phase, steady flow; tree networks with specified outlet flows (no loops).
- Darcy f: 64/Re laminar, Colebrook-White turbulent, linear blend 2300-4000.
- Fittings: Crane TP-410 `K = n·f_T`; reducers by Crane contraction/expansion formulae; fittings are zero-length.
- Pipes integrated in up to 200 steps (400 for gas): density, velocity, temperature, elevation update along the line.
- Heat loss: Gnielinski inside film, wall, insulation, outside (Churchill-Chu / Churchill-Bernstein + radiation) or buried.
- Built-in fluid properties are engineering approximations - verify critical cases with CoolProp or a simulator.
- Old 2D network files (right/left) are converted automatically (right = East, left = West).
