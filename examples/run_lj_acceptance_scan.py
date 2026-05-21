"""Scan Lennard-Jones Monte Carlo acceptance ratio versus trial displacement."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from lj_simulation import scan_acceptance_vs_displacement
from utils import ensure_directory


def main() -> None:
    """Run an acceptance-ratio parameter scan for the LJ Monte Carlo sampler."""

    results_dir = ensure_directory(PROJECT_ROOT / "results" / "lj_acceptance_scan")
    displacement_values = np.linspace(0.05, 0.45, 9)
    scan_results = scan_acceptance_vs_displacement(
        displacement_values=displacement_values,
        n_atoms=32,
        reduced_density=0.75,
        reduced_temperature=1.2,
        n_equilibration_sweeps=60,
        n_production_sweeps=120,
        sample_interval=12,
        random_seed=21,
    )

    dataframe = pd.DataFrame(
        {
            "max_displacement": [item.max_displacement for item in scan_results],
            "acceptance_ratio": [item.acceptance_ratio for item in scan_results],
            "sampled_frames": [item.sampled_frames for item in scan_results],
        }
    )
    csv_path = results_dir / "lj_acceptance_scan.csv"
    dataframe.to_csv(csv_path, index=False)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(8, 5), dpi=160)
    ax.plot(
        dataframe["max_displacement"],
        dataframe["acceptance_ratio"],
        marker="o",
        lw=2.0,
        color="tab:orange",
        label="Acceptance ratio",
    )
    ax.axhline(0.5, color="gray", ls="--", lw=1.2, label="Reference 0.5")
    ax.set_xlabel("Maximum trial displacement", fontsize=12)
    ax.set_ylabel("Acceptance ratio", fontsize=12)
    ax.set_title("Lennard-Jones Monte Carlo: Acceptance vs Displacement", fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)
    fig.tight_layout()
    png_path = results_dir / "lj_acceptance_scan.png"
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)

    best_index = int(np.argmin(np.abs(dataframe["acceptance_ratio"] - 0.5)))
    best_row = dataframe.iloc[best_index]

    print("Lennard-Jones acceptance scan completed.")
    print(f"CSV saved to: {csv_path}")
    print(f"Figure saved to: {png_path}")
    print(
        "Closest-to-0.5 acceptance point: "
        f"max_displacement = {best_row['max_displacement']:.3f}, "
        f"acceptance_ratio = {best_row['acceptance_ratio']:.4f}"
    )


if __name__ == "__main__":
    main()
