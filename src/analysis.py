"""Analysis helpers for RDF peak finding and coordination numbers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.signal import find_peaks

from rdf import RDFResult
from utils import moving_average


@dataclass
class RDFAnalysis:
    """Derived structural information from an RDF curve."""

    first_peak_r: Optional[float]
    first_peak_g: Optional[float]
    first_minimum_r: Optional[float]
    first_minimum_g: Optional[float]
    coordination_number: Optional[float]
    cumulative_coordination: np.ndarray


def cumulative_coordination_number(result: RDFResult) -> np.ndarray:
    r"""Compute cumulative coordination number by shell-by-shell integration.

    The continuous formula is

    ``N(r) = 4 pi rho \int_0^r g(r') r'^2 dr'``.

    Numerically, each histogram bin already represents a shell of finite volume,
    so the discrete analogue is simply

    ``delta N_i = rho * g_i * V_shell,i``.
    """

    shell_population = result.neighbor_density * result.g_r * result.shell_volume
    return np.cumsum(shell_population)


def analyze_rdf(result: RDFResult, smoothing_window: int = 5) -> RDFAnalysis:
    """Find the first peak, first minimum, and first-shell coordination number."""

    smoothed = moving_average(result.g_r, window=smoothing_window)
    peak_indices, _ = find_peaks(smoothed, height=1.0)
    peak_indices = peak_indices[result.g_r[peak_indices] > 0.0]

    first_peak_r = None
    first_peak_g = None
    first_minimum_r = None
    first_minimum_g = None
    coordination_number = None
    cumulative_cn = cumulative_coordination_number(result)

    if peak_indices.size > 0:
        peak_index = int(peak_indices[0])
    else:
        peak_index = int(np.argmax(result.g_r))
        if result.g_r[peak_index] <= 0.0:
            peak_index = -1

    if peak_index >= 0:
        first_peak_r = float(result.r[peak_index])
        first_peak_g = float(result.g_r[peak_index])

        if peak_index < len(result.g_r) - 2:
            search_region = smoothed[peak_index + 1 :]
            minimum_offset = int(np.argmin(search_region))
            minimum_index = peak_index + 1 + minimum_offset
            first_minimum_r = float(result.r[minimum_index])
            first_minimum_g = float(result.g_r[minimum_index])
            coordination_number = float(cumulative_cn[minimum_index])

    return RDFAnalysis(
        first_peak_r=first_peak_r,
        first_peak_g=first_peak_g,
        first_minimum_r=first_minimum_r,
        first_minimum_g=first_minimum_g,
        coordination_number=coordination_number,
        cumulative_coordination=cumulative_cn,
    )
