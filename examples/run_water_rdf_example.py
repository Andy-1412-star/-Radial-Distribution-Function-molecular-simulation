"""Compute O-O, O-H, and H-H RDFs for a small water-like trajectory."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from trajectory import export_rdf_results, load_frames, run_multiple_pair_rdfs
from utils import ensure_directory


def main() -> None:
    """Run three RDF calculations for a water-like XYZ trajectory."""

    frames = load_frames(PROJECT_ROOT / "data" / "example_water.xyz", file_format="xyz")
    results_dir = ensure_directory(PROJECT_ROOT / "results")

    pair_settings = [
        ("O", "O", 6.0, 0.05),
        ("O", "H", 4.0, 0.05),
        ("H", "H", 4.0, 0.05),
    ]
    pair_results = run_multiple_pair_rdfs(
        frames=frames,
        pair_settings=pair_settings,
        box_length=12.0,
        use_pbc=True,
        smoothing_window=5,
    )

    print("Water RDF calculation completed.")
    print(f"Frames used: {len(frames)}")

    for atom_a, atom_b, _, _ in pair_settings:
        result, analysis = pair_results[(atom_a, atom_b)]
        csv_path = results_dir / f"rdf_{atom_a}_{atom_b}.csv"
        png_path = results_dir / f"rdf_{atom_a}_{atom_b}.png"
        export_rdf_results(
            result=result,
            analysis=analysis,
            csv_path=csv_path,
            figure_path=png_path,
            title=f"Water-like RDF: {atom_a}-{atom_b}",
        )

        print(f"{atom_a}-{atom_b}: CSV -> {csv_path.name}, PNG -> {png_path.name}")
        if analysis.first_peak_r is not None:
            print(f"  first peak r = {analysis.first_peak_r:.4f}, g(r) = {analysis.first_peak_g:.4f}")
        if analysis.first_minimum_r is not None:
            print(f"  first minimum r = {analysis.first_minimum_r:.4f}, g(r) = {analysis.first_minimum_g:.4f}")
        if analysis.coordination_number is not None:
            print(f"  first-shell coordination number = {analysis.coordination_number:.4f}")


if __name__ == "__main__":
    main()
