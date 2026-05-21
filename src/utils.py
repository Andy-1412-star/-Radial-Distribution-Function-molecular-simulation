"""Utility helpers shared across the RDF teaching project."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np


@dataclass
class AtomFrame:
    """Container for one XYZ frame.

    Parameters
    ----------
    symbols
        Atomic species labels such as ``Ar`` or ``O``.
    positions
        Cartesian coordinates with shape ``(n_atoms, 3)``.
    metadata
        Optional key-value data parsed from the XYZ comment line.
    """

    symbols: np.ndarray
    positions: np.ndarray
    metadata: Dict[str, str] = field(default_factory=dict)

    def subset_indices(self, atom_type: str) -> np.ndarray:
        """Return atom indices for one element/species."""
        return np.where(self.symbols == atom_type)[0]


def ensure_directory(path: str | Path) -> Path:
    """Create a directory if it does not exist and return it as a ``Path``."""
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    return target


def parse_box_length(box_length: float | Iterable[float]) -> np.ndarray:
    """Convert a scalar or iterable box length into a 3D vector.

    A cubic box may be passed as a single float. A rectangular box may be
    passed as an iterable of length 3.
    """

    if np.isscalar(box_length):
        return np.full(3, float(box_length), dtype=float)

    box = np.asarray(list(box_length), dtype=float)
    if box.shape != (3,):
        raise ValueError("box_length must be a float or an iterable of length 3.")
    return box


def box_volume(box_length: float | Iterable[float]) -> float:
    """Return the simulation cell volume."""
    return float(np.prod(parse_box_length(box_length)))


def minimum_image_displacement(
    delta: np.ndarray,
    box_length: float | Iterable[float],
) -> np.ndarray:
    """Apply the minimum image convention to displacement vectors.

    For each Cartesian component, the vector is wrapped into the interval
    ``[-L/2, L/2)``. This chooses the nearest periodic image of the atom.
    """

    box = parse_box_length(box_length)
    return delta - box * np.round(delta / box)


def pair_distances(
    positions_a: np.ndarray,
    positions_b: np.ndarray,
    box_length: Optional[float | Iterable[float]] = None,
    use_pbc: bool = True,
    same_species: bool = False,
) -> np.ndarray:
    """Compute pair distances between two sets of atoms.

    Parameters
    ----------
    positions_a, positions_b
        Arrays of Cartesian coordinates.
    box_length
        Simulation box length. Required if ``use_pbc=True``.
    use_pbc
        Whether periodic boundary conditions should be used.
    same_species
        If ``True``, only unique pairs ``i < j`` are kept.
    """

    if positions_a.size == 0 or positions_b.size == 0:
        return np.empty(0, dtype=float)

    displacements = positions_b[None, :, :] - positions_a[:, None, :]
    if use_pbc:
        if box_length is None:
            raise ValueError("box_length must be provided when use_pbc=True.")
        displacements = minimum_image_displacement(displacements, box_length)

    distances = np.linalg.norm(displacements, axis=-1)
    if same_species:
        upper_triangle = np.triu_indices(len(positions_a), k=1)
        return distances[upper_triangle]

    return distances.ravel()


def shell_volumes(bin_edges: np.ndarray) -> np.ndarray:
    """Compute the exact volume of each spherical shell.

    The shell between ``r_inner`` and ``r_outer`` has volume

    ``V_shell = 4/3 pi (r_outer^3 - r_inner^3)``.
    """

    r_outer = bin_edges[1:]
    r_inner = bin_edges[:-1]
    return (4.0 / 3.0) * np.pi * (r_outer**3 - r_inner**3)


def moving_average(values: np.ndarray, window: int = 5) -> np.ndarray:
    """Return a simple centered moving-average smoothed signal."""

    if window <= 1:
        return values.copy()
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(values, kernel, mode="same")
