"""Manual RDF implementation for teaching statistical mechanics concepts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional

import numpy as np

from io_utils import frame_box_lengths
from utils import AtomFrame, box_volume, pair_distances, shell_volumes


@dataclass
class RDFResult:
    """Store the radial distribution function and related metadata."""

    r: np.ndarray
    g_r: np.ndarray
    bin_edges: np.ndarray
    counts: np.ndarray
    shell_volume: np.ndarray
    pair_label: str
    atom_type_a: str
    atom_type_b: str
    frames_used: int
    box_length: np.ndarray
    neighbor_density: float
    center_count_average: float
    neighbor_count_average: float


def _expected_shell_counts(
    n_a: int,
    n_b: int,
    volume: float,
    shell_volume: np.ndarray,
    same_species: bool,
) -> tuple[np.ndarray, float]:
    """Return the ideal-gas reference counts used in RDF normalization.

    The RDF is

    ``g(r) = observed pair density / ideal uniform pair density``.

    For unlike species ``A-B``:

    ``expected = N_A * rho_B * V_shell``

    because every A center sees a uniform density ``rho_B = N_B / V`` of
    possible B neighbors.

    For like species ``A-A`` we count each pair only once, so we divide by 2:

    ``expected = 1/2 * N_A * rho_A(neighbor) * V_shell``

    with ``rho_A(neighbor) = (N_A - 1)/V``.
    """

    if same_species:
        if n_a < 2:
            return np.zeros_like(shell_volume), 0.0
        neighbor_density = (n_a - 1) / volume
        expected = 0.5 * n_a * neighbor_density * shell_volume
    else:
        if n_a == 0 or n_b == 0:
            return np.zeros_like(shell_volume), 0.0
        neighbor_density = n_b / volume
        expected = n_a * neighbor_density * shell_volume

    return expected, float(neighbor_density)


def compute_rdf(
    frames: List[AtomFrame],
    atom_type_a: str,
    atom_type_b: str,
    r_max: float,
    bin_width: float,
    box_length: Optional[float | Iterable[float]] = None,
    use_pbc: bool = True,
) -> RDFResult:
    """Compute an RDF by manually counting pairs and normalizing shell by shell.

    Notes
    -----
    The algorithm is intentionally explicit:

    1. Select atoms of type A and type B.
    2. Compute all pair distances.
    3. Bin those distances into a histogram.
    4. Divide by the ideal number of pairs in each spherical shell.

    The final normalization is what turns a raw pair histogram into a
    dimensionless structural quantity.
    """

    if r_max <= 0.0:
        raise ValueError("r_max must be positive.")
    if bin_width <= 0.0:
        raise ValueError("bin_width must be positive.")

    same_species = atom_type_a == atom_type_b
    bin_edges = np.arange(0.0, r_max + bin_width, bin_width, dtype=float)
    if bin_edges[-1] < r_max:
        bin_edges = np.append(bin_edges, r_max)
    r = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    shell_volume = shell_volumes(bin_edges)

    total_counts = np.zeros_like(r)
    total_expected = np.zeros_like(r)
    center_counts = []
    neighbor_counts = []
    neighbor_densities = []
    resolved_box = None

    for frame in frames:
        frame_box = box_length
        if frame_box is None:
            metadata_box = frame_box_lengths(frame)
            if metadata_box is not None:
                frame_box = metadata_box
        if use_pbc and frame_box is None:
            raise ValueError(
                "box_length is required for PBC. Provide it explicitly or store "
                "box_length=... in the XYZ comment line."
            )

        indices_a = frame.subset_indices(atom_type_a)
        indices_b = frame.subset_indices(atom_type_b)
        positions_a = frame.positions[indices_a]
        positions_b = frame.positions[indices_b]

        distances = pair_distances(
            positions_a=positions_a,
            positions_b=positions_b,
            box_length=frame_box,
            use_pbc=use_pbc,
            same_species=same_species,
        )
        counts, _ = np.histogram(distances, bins=bin_edges)

        volume = box_volume(frame_box) if frame_box is not None else _bounding_box_volume(frame.positions)
        expected, neighbor_density = _expected_shell_counts(
            n_a=len(indices_a),
            n_b=len(indices_b),
            volume=volume,
            shell_volume=shell_volume,
            same_species=same_species,
        )

        total_counts += counts
        total_expected += expected
        center_counts.append(len(indices_a))
        neighbor_counts.append(len(indices_b))
        neighbor_densities.append(neighbor_density)
        if frame_box is not None:
            resolved_box = np.asarray(frame_box if not np.isscalar(frame_box) else [frame_box] * 3, dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        g_r = np.divide(
            total_counts,
            total_expected,
            out=np.zeros_like(total_counts, dtype=float),
            where=total_expected > 0.0,
        )

    box_array = resolved_box if resolved_box is not None else np.asarray([np.nan, np.nan, np.nan], dtype=float)
    return RDFResult(
        r=r,
        g_r=g_r,
        bin_edges=bin_edges,
        counts=total_counts,
        shell_volume=shell_volume,
        pair_label=f"{atom_type_a}-{atom_type_b}",
        atom_type_a=atom_type_a,
        atom_type_b=atom_type_b,
        frames_used=len(frames),
        box_length=box_array,
        neighbor_density=float(np.mean(neighbor_densities)) if neighbor_densities else 0.0,
        center_count_average=float(np.mean(center_counts)) if center_counts else 0.0,
        neighbor_count_average=float(np.mean(neighbor_counts)) if neighbor_counts else 0.0,
    )


def _bounding_box_volume(positions: np.ndarray) -> float:
    """Fallback non-PBC volume estimate from the coordinate bounding box."""

    span = positions.max(axis=0) - positions.min(axis=0)
    if np.any(span <= 0.0):
        raise ValueError("Cannot infer a positive volume from the coordinate bounding box.")
    return float(np.prod(span))
