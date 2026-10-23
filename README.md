# 1D Finite Element Solver for Stepped Axial Bars

A Python finite element solver for axial bar structures made of any number of elements in series.
It performs **static analysis** (displacements, stresses, reactions) and **free-vibration (modal) analysis**.
The universal solver asks for all input values first and prints results after.

## Files

| File | Purpose |
|---|---|
| `bar_solver_universal.py` | Interactive solver for any stepped bar. Asks for all inputs, then prints results |
| `bar_analysis_1d.py` | Worked example: three-element stepped steel bar with fixed values |
| `bar_analysis_1d.ipynb` | Same worked example as a notebook, one cell per step |
| `example_input.txt` | Ready-made answers for the worked example |
| `requirements.txt` | Python packages needed |

## Features

| Feature | Supported |
|---|---|
| Any number of bar elements in series | Yes |
| Different length, area, Young's modulus and density per element | Yes |
| Input units: m / cm / mm, m2 / cm2 / mm2, Pa / MPa / GPa, N / kN | Yes |
| Point loads at any node | Yes |
| Any set of fixed nodes | Yes |
| Displacements, reactions, equilibrium check | Yes |
| Element strain, stress (tension +, compression -) and axial force | Yes |
| Natural frequencies and mode shapes (consistent mass matrix) | Yes |
| Distributed loads, thermal loads, 2D trusses, beams | Not yet |

## Method

| Step | What the code does |
|---|---|
| 1 | Element stiffness: k = (EA/L) [[1, -1], [-1, 1]] |
| 2 | Consistent element mass: m = (rho A L / 6) [[2, 1], [1, 2]] |
| 3 | Assemble global K and M matrices |
| 4 | Apply boundary conditions by removing the fixed degrees of freedom |
| 5 | Solve K U = F for nodal displacements |
| 6 | Reactions R = K U - F and equilibrium check |
| 7 | Element strain = (u_j - u_i) / L, stress = E x strain, force = stress x A |
| 8 | Solve (K - w^2 M) phi = 0 as a generalized eigenvalue problem |

## How to run

```bash
pip install -r requirements.txt
python bar_solver_universal.py
```

To run the worked example without typing the inputs:

```bash
python bar_solver_universal.py < example_input.txt
```

## Input order of the universal solver

| Order | Question |
|---|---|
| 1 | Number of elements |
| 2 | Units for length, area, Young's modulus, force |
| 3 | Same Young's modulus for all elements? |
| 4 | Length and area of each element |
| 5 | Loads: node number and force (empty line to finish) |
| 6 | Fixed node numbers, for example `1,4` |
| 7 | Modal analysis? Density |

## Worked example

Three-element stepped steel bar, E = 207 GPa, density 7850 kg/m3, P = 2000 N at node 3, nodes 1 and 4 fixed.

| Element | Length | Area |
|---|---|---|
| 1 | 20 cm | 5 cm2 |
| 2 | 30 cm | 3 cm2 |
| 3 | 40 cm | 1 cm2 |

| Quantity | Result |
|---|---|
| u2 | 2.863 x 10^-6 m |
| u3 | 1.002 x 10^-5 m |
| R1 / R4 | -1481.48 N / -518.52 N |
| Equilibrium (sum of R + P) | ~0 N |

| Element | Stress (MPa) | Axial force (N) |
|---|---|---|
| 1 | +2.963 (tension) | 1481.48 |
| 2 | +4.938 (tension) | 1481.48 |
| 3 | -5.185 (compression) | -518.52 |

| Mode | Frequency (Hz) |
|---|---|
| 1 | 3182.03 |
| 2 | 7682.83 |

## Verification

| Check | Result |
|---|---|
| Static equilibrium (reactions + applied load) | Zero within rounding |
| Cantilever test: 2 elements of 1 m, 100 mm2, E = 200 GPa, 1 kN at the free end | u = 0.100 mm, matches hand calculation (FL/EA) |

## Note on units

The lab reference sheet for this problem lists stiffness values 100 times larger than EA/L in SI
(for element 1: 207e9 x 5e-4 / 0.2 = 5.175 x 10^8 N/m, not 5.175 x 10^10). The solver converts all
inputs to SI internally, so this mistake cannot occur.

## Tools

Python, NumPy, SciPy, VS Code

## Possible extensions

- Distributed and thermal loads
- Lumped mass matrix comparison
- Validation against ANSYS or HyperMesh / OptiStruct rod models
- Plots of displacement, stress and mode shapes
- 2D truss elements
