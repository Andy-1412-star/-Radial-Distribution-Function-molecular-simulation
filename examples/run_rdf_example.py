"""End-to-end RDF example for an Argon-like cubic system."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from trajectory import export_rdf_results, load_frames, run_rdf_analysis
from utils import ensure_directory


def main() -> None:
    """Compute an Ar-Ar RDF and save the outputs."""

    data_file = PROJECT_ROOT / "data" / "example_argon.xyz"
    results_dir = ensure_directory(PROJECT_ROOT / "results")

    frames = load_frames(data_file, file_format="xyz")
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

    csv_path = results_dir / "rdf_Ar_Ar.csv"
    png_path = results_dir / "rdf_Ar_Ar.png"
    export_rdf_results(
        result=rdf_result,
        analysis=rdf_analysis,
        csv_path=csv_path,
        figure_path=png_path,
        title="Argon-like Ar-Ar Radial Distribution Function",
    )

    print("RDF calculation completed.")
    print(f"Frames used: {rdf_result.frames_used}")
    print(f"Pair type: {rdf_result.pair_label}")
    print(f"CSV saved to: {csv_path}")
    print(f"Figure saved to: {png_path}")
    if rdf_analysis.first_peak_r is not None:
        print(f"First peak position: r = {rdf_analysis.first_peak_r:.4f}")
        print(f"First peak height: g(r) = {rdf_analysis.first_peak_g:.4f}")
    if rdf_analysis.first_minimum_r is not None:
        print(f"First minimum position: r = {rdf_analysis.first_minimum_r:.4f}")
        print(f"First minimum value: g(r) = {rdf_analysis.first_minimum_g:.4f}")
    if rdf_analysis.coordination_number is not None:
        print(f"First-shell coordination number: {rdf_analysis.coordination_number:.4f}")


if __name__ == "__main__":
    main()
