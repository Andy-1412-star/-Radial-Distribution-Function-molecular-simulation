"""Static structure factor S(k) computed from a previously calculated RDF."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import simpson

from rdf import RDFResult


@dataclass
class StructureFactorResult:
    """Store a static structure factor curve and its metadata."""

    k: np.ndarray
    s_k: np.ndarray
    pair_label: str
    number_density: float
    rdf_cutoff: float


def compute_structure_factor_from_rdf(
    rdf_result: RDFResult,
    k_min: float,
    k_max: float,
    k_step: float,
) -> StructureFactorResult:
    r"""Compute an isotropic static structure factor from RDF data.

    The relation between real-space and reciprocal-space structure is

    ``S(k) = 1 + 4 \pi \rho \int_0^\infty [g(r) - 1] \frac{\sin(kr)}{kr} r^2 dr``.

    In this project, the integral is evaluated numerically on the discrete RDF
    grid. The RDF is available only up to a finite cutoff, so the integral is
    truncated at the largest sampled radius.
    """

    if k_min < 0.0:
        raise ValueError("k_min must be non-negative.")
    if k_max <= k_min:
        raise ValueError("k_max must be larger than k_min.")
    if k_step <= 0.0:
        raise ValueError("k_step must be positive.")
    if rdf_result.neighbor_density <= 0.0:
        raise ValueError("A positive number density is required to compute S(k).")

    k_values = np.arange(k_min, k_max + 0.5 * k_step, k_step, dtype=float)
    r_values = rdf_result.r
    h_r = rdf_result.g_r - 1.0
    integrand_prefactor = h_r * r_values**2

    s_k_values = np.empty_like(k_values)
    for index, k_value in enumerate(k_values):
        kernel = _sinc_like_kernel(k_value, r_values)
        integral = simpson(integrand_prefactor * kernel, x=r_values)
        s_k_values[index] = 1.0 + 4.0 * np.pi * rdf_result.neighbor_density * integral

    return StructureFactorResult(
        k=k_values,
        s_k=s_k_values,
        pair_label=rdf_result.pair_label,
        number_density=rdf_result.neighbor_density,
        rdf_cutoff=float(r_values[-1]) if len(r_values) > 0 else 0.0,
    )


def structure_factor_dataframe(result: StructureFactorResult) -> pd.DataFrame:
    """Convert a structure factor result into a CSV-ready dataframe."""

    return pd.DataFrame(
        {
            "k": result.k,
            "s_k": result.s_k,
        }
    )


def save_structure_factor_csv(
    result: StructureFactorResult,
    output_path: str | Path,
) -> None:
    """Write an S(k) curve to CSV."""

    dataframe = structure_factor_dataframe(result)
    dataframe.to_csv(output_path, index=False)


def compute_structure_factor_batch(
    rdf_results: Dict[Tuple[str, str], RDFResult],
    k_min: float,
    k_max: float,
    k_step: float,
) -> Dict[Tuple[str, str], StructureFactorResult]:
    """Compute S(k) for a batch of previously calculated RDF results."""

    batch_results: Dict[Tuple[str, str], StructureFactorResult] = {}
    for pair_key, rdf_result in rdf_results.items():
        batch_results[pair_key] = compute_structure_factor_from_rdf(
            rdf_result=rdf_result,
            k_min=k_min,
            k_max=k_max,
            k_step=k_step,
        )
    return batch_results


def summarize_structure_factor_batch(
    batch_results: Dict[Tuple[str, str], StructureFactorResult],
) -> pd.DataFrame:
    """Create a compact summary table for a batch of S(k) calculations."""

    rows = []
    for (atom_type_a, atom_type_b), result in batch_results.items():
        peak_index = int(np.argmax(result.s_k)) if len(result.s_k) > 0 else -1
        peak_k = float(result.k[peak_index]) if peak_index >= 0 else np.nan
        peak_s_k = float(result.s_k[peak_index]) if peak_index >= 0 else np.nan
        rows.append(
            {
                "pair": result.pair_label,
                "atom_type_a": atom_type_a,
                "atom_type_b": atom_type_b,
                "peak_k": peak_k,
                "peak_s_k": peak_s_k,
                "rdf_cutoff": result.rdf_cutoff,
                "number_density": result.number_density,
            }
        )
    return pd.DataFrame(rows).reset_index(drop=True)


def export_structure_factor_batch(
    batch_results: Dict[Tuple[str, str], StructureFactorResult],
    output_directory: str | Path,
    file_prefix: str = "structure_factor",
    summary_filename: str = "structure_factor_summary.csv",
    title_prefix: str = "Static Structure Factor",
) -> pd.DataFrame:
    """Write CSV, PNG, and summary output for a batch of S(k) calculations."""

    output_dir = Path(output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    for (atom_type_a, atom_type_b), result in batch_results.items():
        stem = f"{file_prefix}_{atom_type_a}_{atom_type_b}"
        save_structure_factor_csv(result, output_dir / f"{stem}.csv")
        plot_structure_factor(
            result,
            output_path=output_dir / f"{stem}.png",
            title=f"{title_prefix}: {atom_type_a}-{atom_type_b}",
        )

    summary = summarize_structure_factor_batch(batch_results)
    summary.to_csv(output_dir / summary_filename, index=False)
    return summary


def plot_structure_factor(
    result: StructureFactorResult,
    output_path: Optional[str | Path] = None,
    title: Optional[str] = None,
) -> None:
    """Plot the static structure factor S(k)."""

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(8, 5), dpi=160)
    ax.plot(result.k, result.s_k, lw=2.2, color="tab:purple", label=f"S(k) {result.pair_label}")
    ax.axhline(1.0, lw=1.2, ls="--", color="gray", alpha=0.8, label="Ideal gas limit S(k)=1")
    ax.set_xlabel("k", fontsize=12)
    ax.set_ylabel("S(k)", fontsize=12)
    ax.set_title(title or f"Static Structure Factor: {result.pair_label}", fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)
    fig.tight_layout()

    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def _sinc_like_kernel(k_value: float, r_values: np.ndarray) -> np.ndarray:
    """Return sin(kr)/(kr), with the correct k->0 limit."""

    kr = k_value * r_values
    kernel = np.ones_like(r_values, dtype=float)
    nonzero = np.abs(kr) > 1.0e-12
    kernel[nonzero] = np.sin(kr[nonzero]) / kr[nonzero]
    return kernel
