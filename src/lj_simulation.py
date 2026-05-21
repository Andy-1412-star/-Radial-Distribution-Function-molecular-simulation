"""Teaching-oriented Lennard-Jones fluid simulation using Metropolis Monte Carlo."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List

import numpy as np

from utils import AtomFrame, minimum_image_displacement, parse_box_length


@dataclass
class LJSimulationResult:
    """Store sampled frames and summary statistics from an LJ Monte Carlo run."""

    frames: List[AtomFrame]
    acceptance_ratio: float
    box_length: float
    number_density: float
    temperature: float
    sigma: float
    epsilon: float


@dataclass
class LJAcceptanceScanPoint:
    """Store one Monte Carlo acceptance-ratio scan result."""

    max_displacement: float
    acceptance_ratio: float
    sampled_frames: int


def simulate_lj_fluid(
    n_atoms: int = 32,
    reduced_density: float = 0.75,
    reduced_temperature: float = 1.2,
    n_equilibration_sweeps: int = 100,
    n_production_sweeps: int = 200,
    max_displacement: float = 0.18,
    sample_interval: int = 10,
    sigma: float = 1.0,
    epsilon: float = 1.0,
    cutoff: float = 2.5,
    random_seed: int = 7,
) -> LJSimulationResult:
    """Generate a simple Lennard-Jones fluid trajectory.

    Notes
    -----
    This simulation uses reduced Lennard-Jones units:

    - distance in units of ``sigma``
    - energy in units of ``epsilon``
    - temperature in units of ``epsilon / k_B``

    The goal is educational clarity, not production-level performance.
    """

    if n_atoms <= 0:
        raise ValueError("n_atoms must be positive.")
    if reduced_density <= 0.0:
        raise ValueError("reduced_density must be positive.")
    if reduced_temperature <= 0.0:
        raise ValueError("reduced_temperature must be positive.")
    if n_equilibration_sweeps < 0 or n_production_sweeps <= 0:
        raise ValueError("Sweep counts must be non-negative, with production > 0.")
    if max_displacement <= 0.0:
        raise ValueError("max_displacement must be positive.")
    if sample_interval <= 0:
        raise ValueError("sample_interval must be positive.")

    rng = np.random.default_rng(random_seed)
    box_length = (n_atoms / reduced_density) ** (1.0 / 3.0)
    positions = initialize_cubic_positions(n_atoms=n_atoms, box_length=box_length)
    symbols = np.array(["LJ"] * n_atoms, dtype=str)

    total_attempts = 0
    total_accepts = 0
    sampled_frames: List[AtomFrame] = []

    for sweep in range(n_equilibration_sweeps + n_production_sweeps):
        for particle_index in range(n_atoms):
            total_attempts += 1
            old_position = positions[particle_index].copy()
            old_energy = particle_site_energy(
                particle_index=particle_index,
                positions=positions,
                box_length=box_length,
                sigma=sigma,
                epsilon=epsilon,
                cutoff=cutoff,
            )

            trial_displacement = rng.uniform(-max_displacement, max_displacement, size=3)
            positions[particle_index] = (positions[particle_index] + trial_displacement) % box_length

            new_energy = particle_site_energy(
                particle_index=particle_index,
                positions=positions,
                box_length=box_length,
                sigma=sigma,
                epsilon=epsilon,
                cutoff=cutoff,
            )
            delta_energy = new_energy - old_energy
            if delta_energy <= 0.0 or rng.random() < np.exp(-delta_energy / reduced_temperature):
                total_accepts += 1
            else:
                positions[particle_index] = old_position

        production_sweep = sweep - n_equilibration_sweeps
        if production_sweep >= 0 and production_sweep % sample_interval == 0:
            sampled_frames.append(
                AtomFrame(
                    symbols=symbols.copy(),
                    positions=positions.copy(),
                    metadata={
                        "box_length": f"{box_length:.6f}",
                        "sweep": str(sweep),
                        "reduced_density": f"{reduced_density:.6f}",
                        "reduced_temperature": f"{reduced_temperature:.6f}",
                    },
                )
            )

    return LJSimulationResult(
        frames=sampled_frames,
        acceptance_ratio=total_accepts / total_attempts if total_attempts > 0 else 0.0,
        box_length=box_length,
        number_density=reduced_density,
        temperature=reduced_temperature,
        sigma=sigma,
        epsilon=epsilon,
    )


def scan_acceptance_vs_displacement(
    displacement_values: Iterable[float],
    n_atoms: int = 32,
    reduced_density: float = 0.75,
    reduced_temperature: float = 1.2,
    n_equilibration_sweeps: int = 60,
    n_production_sweeps: int = 120,
    sample_interval: int = 12,
    sigma: float = 1.0,
    epsilon: float = 1.0,
    cutoff: float = 2.5,
    random_seed: int = 7,
) -> List[LJAcceptanceScanPoint]:
    """Evaluate how Monte Carlo acceptance changes with trial step size.

    This is a useful teaching diagnostic because it shows the tradeoff:

    - very small displacements are accepted often but explore configuration space slowly
    - very large displacements are rejected often and waste Monte Carlo steps
    """

    scan_results: List[LJAcceptanceScanPoint] = []
    for index, max_displacement in enumerate(displacement_values):
        simulation = simulate_lj_fluid(
            n_atoms=n_atoms,
            reduced_density=reduced_density,
            reduced_temperature=reduced_temperature,
            n_equilibration_sweeps=n_equilibration_sweeps,
            n_production_sweeps=n_production_sweeps,
            max_displacement=max_displacement,
            sample_interval=sample_interval,
            sigma=sigma,
            epsilon=epsilon,
            cutoff=cutoff,
            random_seed=random_seed + index,
        )
        scan_results.append(
            LJAcceptanceScanPoint(
                max_displacement=max_displacement,
                acceptance_ratio=simulation.acceptance_ratio,
                sampled_frames=len(simulation.frames),
            )
        )
    return scan_results


def initialize_cubic_positions(n_atoms: int, box_length: float) -> np.ndarray:
    """Place particles on a simple cubic grid inside the box."""

    n_side = int(np.ceil(n_atoms ** (1.0 / 3.0)))
    grid = np.linspace(0.0, box_length, n_side, endpoint=False, dtype=float)
    spacing = box_length / n_side
    offset = 0.5 * spacing

    positions = []
    for x in grid:
        for y in grid:
            for z in grid:
                positions.append([x + offset, y + offset, z + offset])
                if len(positions) == n_atoms:
                    return np.asarray(positions, dtype=float)

    raise RuntimeError("Failed to initialize the requested number of particles.")


def particle_site_energy(
    particle_index: int,
    positions: np.ndarray,
    box_length: float,
    sigma: float,
    epsilon: float,
    cutoff: float,
) -> float:
    """Return the interaction energy involving one selected particle."""

    origin = positions[particle_index]
    others = np.delete(positions, particle_index, axis=0)
    if len(others) == 0:
        return 0.0

    displacements = others - origin
    displacements = minimum_image_displacement(displacements, parse_box_length(box_length))
    distances = np.linalg.norm(displacements, axis=1)
    return float(np.sum(lennard_jones_pair_energy(distances, sigma=sigma, epsilon=epsilon, cutoff=cutoff)))


def lennard_jones_pair_energy(
    distances: np.ndarray,
    sigma: float,
    epsilon: float,
    cutoff: float,
) -> np.ndarray:
    """Evaluate the shifted Lennard-Jones pair potential on many distances."""

    energies = np.zeros_like(distances, dtype=float)
    valid = distances < cutoff
    if not np.any(valid):
        return energies

    safe_distances = distances[valid]
    inverse_r = sigma / safe_distances
    inverse_r6 = inverse_r**6
    inverse_r12 = inverse_r6**2
    raw_energy = 4.0 * epsilon * (inverse_r12 - inverse_r6)

    cutoff_inverse = sigma / cutoff
    cutoff_inverse6 = cutoff_inverse**6
    cutoff_inverse12 = cutoff_inverse6**2
    cutoff_shift = 4.0 * epsilon * (cutoff_inverse12 - cutoff_inverse6)
    energies[valid] = raw_energy - cutoff_shift
    return energies
