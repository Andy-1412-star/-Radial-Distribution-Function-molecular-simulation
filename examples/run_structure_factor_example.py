"""Compute S(k) from the Argon-like RDF as a teaching example."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from structure_factor import (
    compute_structure_factor_from_rdf,
    plot_structure_factor,
    save_structure_factor_csv,
)
from trajectory import load_frames, run_rdf_analysis
from utils import ensure_directory


def main() -> None:
    """Compute an Ar-Ar RDF and transform it into S(k)."""

    results_dir = ensure_directory(PROJECT_ROOT / "results")
    frames = load_frames(PROJECT_ROOT / "data" / "example_argon.xyz", file_format="xyz")
    rdf_result, rdf_analysis = run_rdf_analysis(
        frames=frames,
        atom_type_a="Ar",
        atom_type_b="Ar",
        r_max=5.0,
        bin_width=0.05,
        box_length=10.52,
        use_pbc=True,
        smoothing_window=5,
    )
    structure_result = compute_structure_factor_from_rdf(
        rdf_result=rdf_result,
        k_min=0.1,
        k_max=20.0,
        k_step=0.1,
    )

    csv_path = results_dir / "structure_factor_Ar_Ar.csv"
    png_path = results_dir / "structure_factor_Ar_Ar.png"
    save_structure_factor_csv(structure_result, csv_path)
    plot_structure_factor(
        structure_result,
        output_path=png_path,
        title="Argon-like Static Structure Factor from RDF",
    )

    print("Structure factor calculation completed.")
    print(f"Pair type: {rdf_result.pair_label}")
    print(f"RDF first peak: r = {rdf_analysis.first_peak_r:.4f}" if rdf_analysis.first_peak_r is not None else "RDF first peak: not found")
    print(f"CSV saved to: {csv_path}")
    print(f"Figure saved to: {png_path}")
    print(f"Density used in S(k): {structure_result.number_density:.6f}")
    print(f"RDF cutoff used in the integral: r_max = {structure_result.rdf_cutoff:.4f}")


if __name__ == "__main__":
    main()
