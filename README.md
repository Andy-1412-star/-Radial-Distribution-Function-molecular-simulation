# Radial Distribution Function From Scratch

A teaching-oriented computational physics project for implementing
Radial Distribution Functions (RDF) and static structure factors
from first principles in Python.

This project is a GitHub-ready Python implementation of the **Radial Distribution Function (RDF)**. The goal is not just to make a plot, but to help students understand:

- what RDF means physically,
- why the normalization matters,
- how periodic boundary conditions enter the calculation,
- how to implement the core algorithm by hand,
- and how the same workflow extends to `XYZ`, `LAMMPS`, and `GROMACS` trajectories.

It also includes a teaching-oriented bridge from real-space structure
$g(r)$ to reciprocal-space structure $S(k)$.

---

## Overview

The project is designed as an entry point for students in:

- computational physics
- molecular dynamics
- statistical mechanics
- atomistic simulation
- computational chemistry and materials science

The code **does not call a prebuilt RDF routine**. Instead, it implements the RDF workflow explicitly:

1. pair distance calculation
2. histogram binning
3. shell-volume normalization
4. number-density normalization
5. periodic boundary condition handling
6. minimum image convention
7. coordination-number integration
8. numerical structure-factor integration from $g(r)$

---

## What This Project Demonstrates

This project demonstrates:

- implementing a physics algorithm from first principles,
- handling atomistic trajectory data,
- applying periodic boundary conditions,
- connecting real-space RDFs to reciprocal-space structure factors,
- building reproducible scientific Python workflows,
- and developing computational physics tooling for molecular simulation.

The repository is designed to bridge:

```text
physics theory -> numerical implementation -> trajectory analysis -> structural interpretation
```

---

## Physics Background

The **Radial Distribution Function** $g(r)$ describes how particle density varies as a function of distance from a reference particle.

Intuitively:

> If I stand on one atom, how likely am I to find another atom at distance $r$, compared with an ideal uniformly random system?

That is why RDF is often summarized as:

$$
g(r)=\frac{\text{actual pair density}}{\text{ideal gas pair density}}
$$

### What does $g(r)$ mean?

- $g(r)=1$: the local particle distribution looks like a random uniform system
- $g(r)>1$: particles are more likely to be found at that distance
- $g(r)<1$: particles are less likely to be found at that distance

### Why is $g(r)$ small at short distance?

At very small $r$, atoms cannot overlap. Because of repulsion and excluded volume, the probability of finding another atom extremely close to the reference atom is near zero, so $g(r)$ usually starts near zero.

### Why do liquids show peaks and damped oscillations?

Liquids have **short-range order**:

- a strong **first peak** from nearest neighbors,
- a weaker **second peak** from the next shell,
- oscillations that gradually decay.

At long distance, the system forgets the original reference atom, so:

$$
g(r)\rightarrow 1
$$

### Ideal gas, liquid, crystal

- **Ideal gas**: $g(r)$ is nearly flat and close to 1
- **Liquid**: clear first peak, second peak, and damped oscillations
- **Crystal**: sharp peaks at specific shell distances, often extending to long range

---

## Mathematical Formulation

For an isotropic system, the RDF depends only on distance $r$.

If we consider a spherical shell between $r$ and $r+dr$, the shell volume is:

$$
dV_{\text{shell}} = 4\pi r^2dr
$$

If the average number density is:

$$
\rho = \frac{N}{V}
$$

then the expected number of neighbors in an ideal uniform system is:

$$
dN_{\text{ideal}}=\rho\,4\pi r^2dr
$$

The real system has:

$$
dN_{\text{real}}=\rho\,g(r)\,4\pi r^2dr
$$

Therefore:

$$
g(r)=\frac{dN_{\text{real}}}{dN_{\text{ideal}}}
$$

### Why shell volume is required

If we only count pairs, bins at larger radius contain more volume and therefore more pairs even in a completely random system. Dividing by the shell volume removes this geometric effect.

In the code, the histogram bins have finite width, so the exact shell volume is used:

$$
V_{\text{shell}}=\frac{4}{3}\pi\left(r_{\text{outer}}^3-r_{\text{inner}}^3\right)
$$

### Why number density is required

If two systems have different densities, raw pair counts are not directly comparable. Dividing by the number density $\rho$ converts the result into a dimensionless structural measure.

### Coordination number

The cumulative coordination number is:

$$
N(r)=4\pi\rho\int_0^r g(r')r'^2\,dr'
$$

This tells us the average number of neighbors within radius $r$.

In practice:

- the **first peak** marks the most probable nearest-neighbor distance,
- the **first minimum** after that peak is often used as the boundary of the first coordination shell,
- integrating up to that first minimum gives the **first coordination number**.

---

## From RDF to Structure Factor

The RDF describes structure in real space. The static structure factor $S(k)$
describes structure in reciprocal space. The two are linked by a spherical
Fourier transform:

$$
S(k)=1+4\pi\rho\int_0^\infty [g(r)-1]\frac{\sin(kr)}{kr}r^2\,dr
$$

This formula is one of the most useful bridges in molecular simulation:

- $g(r)-1$ measures the deviation from an ideal gas,
- $\frac{\sin(kr)}{kr}$ is the isotropic Fourier kernel,
- $\rho$ sets the density scale,
- integrating over $r$ converts local structure into wavevector-dependent structure.

### Physical meaning of $S(k)$

- $S(k)=1$: ideal-gas-like behavior at that wavevector
- large peaks in $S(k)$: strong structural correlations at the corresponding length scale
- small-$k$ behavior: connected to long-wavelength density fluctuations

In experiment, $S(k)$ is closely related to X-ray and neutron scattering.

---

## PBC and Minimum Image Convention

### Why periodic boundary conditions matter

Finite simulation boxes create edge effects. Atoms near the boundary have missing neighbors if we treat the box as isolated.

Periodic boundary conditions (PBC) solve this by tiling space with copies of the simulation box, so the system behaves more like bulk matter.

### Minimum image convention

With PBC, each particle has infinitely many periodic images. The physically relevant distance is the **shortest periodic-image distance**. This is the minimum image convention.

In a cubic box of length $L$, each displacement component is wrapped to:

$$
[-L/2,L/2)
$$

This logic is implemented manually in [`src/utils.py`](src/utils.py).

---

## Bin Width and Statistical Noise

The RDF is a statistical quantity, so binning matters:

- small `bin_width`: high resolution, but noisier
- large `bin_width`: smoother, but fine structure is blurred

This project keeps the core RDF unsmoothed, but uses a small moving average during analysis to stabilize automatic peak and minimum detection.

---

## Features

- Handwritten RDF core algorithm
- `A-A`, `A-B`, and `B-B` RDF
- `XYZ` single-frame and multi-frame support
- `LAMMPS dump` support
- `GROMACS .gro` support
- Optional `MDAnalysis` support for `.xtc` and `.trr`
- Automatic first-peak detection
- Automatic first-minimum detection
- Automatic first coordination number
- Automatic species detection and all-pairs RDF batch export
- Pair-specific RDF parameter configuration via JSON
- Static structure factor $S(k)$ derived from RDF
- Species-resolved $S(k)$ batch export for multi-component systems
- Lennard-Jones fluid trajectory generation from a simple Monte Carlo simulation
- CSV output
- Publication-style RDF plots with annotations
- Water RDF example: `O-O`, `O-H`, `H-H`

---

## Installation

Use Python `3.11`.

### Core installation

```bash
pip install -r requirements.txt
```

### Development tools

```bash
pip install -r requirements-dev.txt
```

### Optional GROMACS trajectory support

This project does **not** require `MDAnalysis` for the core RDF algorithm.

If you want to read `GROMACS .xtc` or `.trr`, install:

```bash
pip install MDAnalysis
```

If `MDAnalysis` is missing, the code will raise:

```text
GROMACS trajectory support requires MDAnalysis. Please install it using pip install MDAnalysis.
```

---

## Quick Start

From the project root:

```bash
python examples/run_rdf_example.py
python examples/run_water_rdf_example.py
python examples/run_lammps_rdf_example.py
python examples/run_auto_pairs_rdf_example.py
python examples/run_structure_factor_example.py
python examples/run_water_structure_factor_example.py
python examples/run_lj_fluid_example.py
python examples/run_lj_acceptance_scan.py
```

If `MDAnalysis` is installed, you can also run:

```bash
python examples/run_gromacs_rdf_example.py
```

---

## Examples

### 1. Argon-like RDF

```bash
python examples/run_rdf_example.py
```

This reads [`data/example_argon.xyz`](data/example_argon.xyz) and writes:

- [`results/rdf_Ar_Ar.csv`](results/rdf_Ar_Ar.csv)
- [`results/rdf_Ar_Ar.png`](results/rdf_Ar_Ar.png)

### 2. Water RDF Analysis

```bash
python examples/run_water_rdf_example.py
```

This reads [`data/example_water.xyz`](data/example_water.xyz) and computes:

- $g_{OO}(r)$
- $g_{OH}(r)$
- $g_{HH}(r)$

Outputs:

- [`results/rdf_O_O.csv`](results/rdf_O_O.csv)
- [`results/rdf_O_H.csv`](results/rdf_O_H.csv)
- [`results/rdf_H_H.csv`](results/rdf_H_H.csv)
- [`results/rdf_O_O.png`](results/rdf_O_O.png)
- [`results/rdf_O_H.png`](results/rdf_O_H.png)
- [`results/rdf_H_H.png`](results/rdf_H_H.png)

#### Physical interpretation of water RDFs

- **O-O RDF**: the first peak corresponds to the most probable intermolecular oxygen-oxygen separation.
- **O-H RDF**: highlights the local hydrogen environment and can reflect hydrogen-bond-related structure.
- **H-H RDF**: shows spatial correlation between hydrogen atoms.

### 3. Automatic All-Pairs RDF

```bash
python examples/run_auto_pairs_rdf_example.py
```

This workflow automatically:

1. scans the trajectory to detect every atomic species present,
2. builds all unique unordered pairs such as `O-O`, `O-H`, and `H-H`,
3. computes one RDF for each pair,
4. writes one CSV and one PNG per pair,
5. writes a summary table with peak positions and coordination numbers.

For the bundled water example, the script writes into:

- `results/auto_pairs/rdf_O_O.csv`
- `results/auto_pairs/rdf_O_H.csv`
- `results/auto_pairs/rdf_H_H.csv`
- `results/auto_pairs/rdf_summary.csv`

This is useful when you want to explore a new trajectory without manually listing every pair first.

#### Pair-specific parameter configuration

Real systems often benefit from using different RDF settings for different pairs. For example:

- `O-O` may need a larger `r_max`
- `O-H` may benefit from a smaller `bin_width`
- `H-H` may use a different short-range window

The automatic script accepts a JSON configuration file. A bundled example lives at:

- [`data/example_pair_parameters_water.json`](data/example_pair_parameters_water.json)

Example structure:

```json
{
  "default": {
    "r_max": 6.0,
    "bin_width": 0.05,
    "smoothing_window": 5,
    "title_prefix": "Automatic Pair RDF",
    "output_prefix": "rdf"
  },
  "pairs": {
    "O-O": {
      "r_max": 6.0,
      "bin_width": 0.05,
      "smoothing_window": 5,
      "title": "Water O-O RDF",
      "output_stem": "rdf_water_O_O"
    },
    "O-H": {
      "r_max": 4.0,
      "bin_width": 0.025,
      "smoothing_window": 7,
      "title": "Water O-H RDF",
      "output_stem": "rdf_water_O_H"
    },
    "H-H": {
      "r_max": 4.0,
      "bin_width": 0.025,
      "smoothing_window": 7,
      "title": "Water H-H RDF",
      "output_stem": "rdf_water_H_H"
    }
  }
}
```

Usage:

```bash
python examples/run_auto_pairs_rdf_example.py --pair-config data/example_pair_parameters_water.json
```

If a pair is not listed in the config file, the script falls back to the `default` values.

The config can control, per pair:

- `r_max`
- `bin_width`
- `smoothing_window`
- plot `title`
- output filename stem via `output_stem`

### 4. LAMMPS Trajectory Analysis

```bash
python examples/run_lammps_rdf_example.py
```

This reads [`data/example_lammps.dump`](data/example_lammps.dump).

The current parser supports a simple orthogonal-box format with:

```text
ITEM: TIMESTEP
ITEM: NUMBER OF ATOMS
ITEM: BOX BOUNDS
ITEM: ATOMS id type x y z
```

It also supports multiple frames and trajectory averaging.

If the atom type is numeric, you can pass a mapping such as:

```python
type_mapping = {
    1: "O",
    2: "H",
}
```

The example script exposes that through `--type-map`.

Example:

```bash
python examples/run_lammps_rdf_example.py --atom-a O --atom-b H --type-map 1:O 2:H
```

If the file format is unsupported, the parser raises a clear error message explaining which section or coordinate columns are missing.

### 5. GROMACS Trajectory Analysis

#### GRO example without MDAnalysis

```bash
python examples/run_gromacs_rdf_example.py
```

This reads [`data/example_water.gro`](data/example_water.gro) directly using pure Python.

#### User-defined pair RDF

```bash
python examples/run_gromacs_rdf_example.py --atom-a O --atom-b H
```

#### XTC/TRR example with MDAnalysis installed

```bash
python examples/run_gromacs_rdf_example.py --topology-file system.gro --trajectory-file traj.xtc --atom-a O --atom-b O
```

The core RDF is still computed by **our own code**. `MDAnalysis` is used only as an optional trajectory reader.

### 6. Structure Factor Example

```bash
python examples/run_structure_factor_example.py
```

This script first computes the Argon-like `Ar-Ar` RDF, then evaluates:

$$
S(k)=1+4\pi\rho\int_0^\infty [g(r)-1]\frac{\sin(kr)}{kr}r^2\,dr
$$

using the project's own numerical integration code.

Outputs:

- `results/structure_factor_Ar_Ar.csv`
- `results/structure_factor_Ar_Ar.png`

This keeps the project teaching-oriented: you can see how one structural signal
appears both in real space through $g(r)$ and in reciprocal space through $S(k)$.

### 7. Water Structure Factor Example

```bash
python examples/run_water_structure_factor_example.py
```

This example computes species-resolved structure factors for:

- `O-O`
- `O-H`
- `H-H`

using the corresponding water RDFs as input.

Outputs are written to:

- `results/structure_factor_water/structure_factor_water_O_O.csv`
- `results/structure_factor_water/structure_factor_water_O_H.csv`
- `results/structure_factor_water/structure_factor_water_H_H.csv`
- `results/structure_factor_water/structure_factor_summary.csv`

This is a useful next step for multi-component systems because it shows how
different pair correlations contribute to reciprocal-space structure.

### 8. Lennard-Jones Fluid Example

```bash
python examples/run_lj_fluid_example.py
```

This example takes the project one step further: instead of only reading a
trajectory, it first **generates** one using a small Metropolis Monte Carlo
simulation of a Lennard-Jones fluid in reduced units.

The script then:

1. saves the sampled trajectory as `lj_fluid.xyz`,
2. computes $g(r)$ for the `LJ-LJ` fluid,
3. computes $S(k)$ from that RDF.

Outputs are written to:

- `results/lj_fluid/lj_fluid.xyz`
- `results/lj_fluid/rdf_LJ_LJ.csv`
- `results/lj_fluid/rdf_LJ_LJ.png`
- `results/lj_fluid/structure_factor_LJ_LJ.csv`
- `results/lj_fluid/structure_factor_LJ_LJ.png`

This turns the repository into a more complete teaching pipeline:

```text
simulation -> trajectory -> RDF -> structure factor
```

### 9. Lennard-Jones Acceptance Scan

```bash
python examples/run_lj_acceptance_scan.py
```

This example scans the Monte Carlo trial displacement and measures the
acceptance ratio. It is useful because it teaches a practical simulation lesson:

- too small a displacement gives very high acceptance but slow exploration,
- too large a displacement gives low acceptance and inefficient sampling.

Outputs are written to:

- `results/lj_acceptance_scan/lj_acceptance_scan.csv`
- `results/lj_acceptance_scan/lj_acceptance_scan.png`

---

## Example Figures

Generated figures include:

- RDF curve $g(r)$ vs $r$
- structure factor curve $S(k)$ vs $k$
- ideal-gas reference line $g(r)=1$
- first peak marker
- first minimum marker
- shaded coordination-number region up to the first minimum

---

### O-O RDF

![](results/rdf_O_O.png)

---

### O-H RDF

![](results/rdf_O_H.png)

---

### H-H RDF

![](results/rdf_H_H.png)

---

### Lennard-Jones Fluid RDF

![](results/lj_fluid/rdf_LJ_LJ.png)

---

### Structure Factor Example

![](results/structure_factor_Ar_Ar.png)

Example output figures are written into `results/`.

---

## Docs

The repository also includes a short teaching note and notebook:

- [`docs/lj_rdf_sk_tutorial.md`](docs/lj_rdf_sk_tutorial.md)
- [`notebooks/01_rdf_and_structure_factor_walkthrough.ipynb`](notebooks/01_rdf_and_structure_factor_walkthrough.ipynb)

This document connects the whole learning path:

```text
Monte Carlo simulation -> trajectory -> RDF -> structure factor
```

---

## Simulation Workflow

```text
Monte Carlo / Molecular Dynamics
                ↓
          Atomic Trajectory
                ↓
      Radial Distribution Function g(r)
                ↓
       Structure Factor S(k)
                ↓
        Structural Interpretation
```

This repository is intended to help students understand how simulation data becomes physically meaningful structural observables.

---

## Project Structure

```text
rdf_project/
├── README.md
├── requirements.txt
├── requirements-dev.txt
├── .gitignore
├── LICENSE
├── docs/
│   └── lj_rdf_sk_tutorial.md
├── notebooks/
│   └── 01_rdf_and_structure_factor_walkthrough.ipynb
├── data/
│   ├── example_argon.xyz
│   ├── example_water.xyz
│   ├── example_lammps.dump
│   ├── example_pair_parameters_water.json
│   └── example_water.gro
├── src/
│   ├── rdf.py
│   ├── io_utils.py
│   ├── plotting.py
│   ├── analysis.py
│   ├── lj_simulation.py
│   ├── structure_factor.py
│   ├── trajectory.py
│   └── utils.py
├── examples/
│   ├── run_rdf_example.py
│   ├── run_water_rdf_example.py
│   ├── run_lammps_rdf_example.py
│   ├── run_gromacs_rdf_example.py
│   ├── run_auto_pairs_rdf_example.py
│   ├── run_structure_factor_example.py
│   ├── run_water_structure_factor_example.py
│   ├── run_lj_fluid_example.py
│   └── run_lj_acceptance_scan.py
└── results/
```

### Module roles

- [`src/rdf.py`](src/rdf.py): hand-written RDF algorithm
- [`src/io_utils.py`](src/io_utils.py): file parsing for `XYZ`, `LAMMPS`, and `GROMACS`
- [`src/analysis.py`](src/analysis.py): first peak, first minimum, coordination number
- [`src/lj_simulation.py`](src/lj_simulation.py): simple Lennard-Jones Metropolis Monte Carlo sampling
- [`src/structure_factor.py`](src/structure_factor.py): numerical $S(k)$ integration from RDF
- [`src/plotting.py`](src/plotting.py): figure generation and annotations
- [`src/trajectory.py`](src/trajectory.py): high-level workflow helpers for scripts
- [`src/utils.py`](src/utils.py): PBC, minimum image convention, shell volumes, smoothing

### Automatic batch workflow

The batch interface lives in [`src/trajectory.py`](src/trajectory.py):

- `detect_species(...)` finds every unique atomic label across the trajectory
- `build_auto_pair_settings(...)` constructs all unique unordered pairs
- `run_multiple_pair_rdfs(...)` computes the batch
- `export_batch_rdf_results(...)` writes per-pair files and one summary table

---

## Future Roadmap

- Partial and species-resolved structure factors
- Lennard-Jones fluid simulation
- Monte Carlo simulation
- Diffusion coefficient analysis
- Velocity autocorrelation function
- Hydrogen bond network analysis
- Time-dependent RDF
- Quantum molecular dynamics observables
- Integration with OpenMM / LAMMPS / GROMACS

---

## References

1. D. Frenkel and B. Smit, *Understanding Molecular Simulation*
2. M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*
3. J.-P. Hansen and I. R. McDonald, *Theory of Simple Liquids*
4. MDAnalysis documentation
5. LAMMPS documentation
6. GROMACS documentation

---

## Project Positioning

This repository is best thought of as:

> A computational physics and molecular simulation project focused on implementing structural observables from first principles in Python, with emphasis on RDFs, structure factors, trajectory analysis, and statistical mechanics interpretation.
