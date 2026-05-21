# Lennard-Jones Fluid Tutorial

This note explains the teaching path behind the Lennard-Jones part of the project:

```text
Monte Carlo simulation -> sampled trajectory -> RDF g(r) -> structure factor S(k)
```

## 1. Why add a Lennard-Jones fluid?

Earlier parts of the repository focused on **reading trajectories** and then
computing structural observables. The Lennard-Jones example adds a missing
piece: it shows how a simple simulation can generate the trajectory in the
first place.

That makes the repository more complete for learning:

- you generate atom positions yourself,
- analyze them in real space with `g(r)`,
- and then analyze the same structure in reciprocal space with `S(k)`.

## 2. The Lennard-Jones potential

The pair potential is

```text
U(r) = 4 epsilon [ (sigma / r)^12 - (sigma / r)^6 ]
```

Interpretation:

- the `r^-12` term gives strong short-range repulsion
- the `r^-6` term gives intermediate-range attraction
- `sigma` sets the length scale
- `epsilon` sets the energy scale

In this project we work in **reduced Lennard-Jones units**, where:

- distance is measured in units of `sigma`
- energy is measured in units of `epsilon`
- temperature is measured in units of `epsilon / k_B`

## 3. Why Monte Carlo?

The repository uses a small **Metropolis Monte Carlo** sampler instead of full
molecular dynamics because the algorithm is easier to explain line by line.

One Monte Carlo move is:

1. pick one particle
2. propose a random displacement
3. compute the energy change `ΔU`
4. accept the move if `ΔU <= 0`
5. otherwise accept with probability `exp(-ΔU / T*)`

Here `T*` is the reduced temperature.

## 4. Acceptance ratio and trial displacement

The trial displacement size matters a lot:

- too small: moves are accepted often, but the system explores configuration space slowly
- too large: moves are rejected often, so the simulation wastes steps

That is why the repository includes the acceptance-scan example:

```bash
python examples/run_lj_acceptance_scan.py
```

This produces:

- a CSV table of acceptance ratios
- a plot of acceptance ratio versus trial step size

This is a classic simulation diagnostic.

## 5. From the simulated fluid to RDF

After sampling a trajectory, the project computes:

```text
g(r) = actual pair density / ideal gas pair density
```

For a Lennard-Jones fluid, the RDF should typically show:

- `g(r) ≈ 0` at short range because particles cannot overlap
- a first peak near the preferred nearest-neighbor distance
- damped oscillations as liquid-like order decays

## 6. From RDF to S(k)

Once `g(r)` is available, the code evaluates:

```text
S(k) = 1 + 4πρ ∫ [g(r) - 1] sin(kr)/(kr) r² dr
```

So the same simulated trajectory can be described in two complementary ways:

- `g(r)` in real space
- `S(k)` in reciprocal space

## 7. Suggested learning order

If you want to study this part of the repository, read the files in this order:

1. `src/lj_simulation.py`
2. `examples/run_lj_fluid_example.py`
3. `src/rdf.py`
4. `src/structure_factor.py`
5. `examples/run_lj_acceptance_scan.py`

## 8. What to try next

- change the reduced temperature and see how the RDF peak changes
- change the density and compare liquid-like structure
- compare low-density and high-density `S(k)`
- test how acceptance ratio changes when you change `max_displacement`
