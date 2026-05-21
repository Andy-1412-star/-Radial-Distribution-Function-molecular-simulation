"""Plotting utilities for RDF results."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt

from analysis import RDFAnalysis
from rdf import RDFResult


def plot_rdf(
    result: RDFResult,
    analysis: Optional[RDFAnalysis] = None,
    output_path: str | Path | None = None,
    title: str | None = None,
) -> None:
    """Plot ``g(r)`` and optionally annotate the first shell features."""

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(8, 5), dpi=160)
    ax.plot(result.r, result.g_r, lw=2.2, color="tab:blue", label=f"RDF {result.pair_label}")
    ax.axhline(1.0, lw=1.2, ls="--", color="gray", alpha=0.8, label="Ideal gas limit g(r)=1")

    if analysis is not None and analysis.first_minimum_r is not None:
        coordination_mask = result.r <= analysis.first_minimum_r
        ax.fill_between(
            result.r[coordination_mask],
            result.g_r[coordination_mask],
            color="tab:orange",
            alpha=0.18,
            label="First-shell coordination region",
        )

    if analysis is not None and analysis.first_peak_r is not None:
        ax.scatter(
            [analysis.first_peak_r],
            [analysis.first_peak_g],
            color="tab:red",
            s=45,
            zorder=3,
            label=f"First peak: r={analysis.first_peak_r:.3f}",
        )
        ax.annotate(
            f"Peak\nr={analysis.first_peak_r:.3f}",
            xy=(analysis.first_peak_r, analysis.first_peak_g),
            xytext=(10, 12),
            textcoords="offset points",
            fontsize=10,
            color="tab:red",
        )

    if analysis is not None and analysis.first_minimum_r is not None:
        ax.scatter(
            [analysis.first_minimum_r],
            [analysis.first_minimum_g],
            color="tab:green",
            s=45,
            zorder=3,
            label=f"First minimum: r={analysis.first_minimum_r:.3f}",
        )
        text = f"Minimum\nr={analysis.first_minimum_r:.3f}"
        if analysis.coordination_number is not None:
            text += f"\nCN={analysis.coordination_number:.2f}"
        ax.annotate(
            text,
            xy=(analysis.first_minimum_r, analysis.first_minimum_g),
            xytext=(10, -30),
            textcoords="offset points",
            fontsize=10,
            color="tab:green",
        )

    ax.set_xlabel("r", fontsize=12)
    ax.set_ylabel("g(r)", fontsize=12)
    ax.set_title(title or f"Radial Distribution Function: {result.pair_label}", fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)
    fig.tight_layout()

    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
