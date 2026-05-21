"""Compute species-resolved S(k) curves for a small water-like trajectory."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from structure_factor import compute_structure_factor_batch, export_structure_factor_batch
from trajectory import load_frames, run_multiple_pair_rdfs
from utils import ensure_directory


def main() -> None:
    """Compute O-O, O-H, and H-H structure factors from their RDFs."""

    results_dir = ensure_directory(PROJECT_ROOT / "results" / "structure_factor_water")
    frames = load_frames(PROJECT_ROOT / "data" / "example_water.xyz", file_format="xyz")
    pair_settings = [
        ("O", "O", 6.0, 0.05),
        ("O", "H", 4.0, 0.025),
        ("H", "H", 4.0, 0.025),
    ]
    rdf_batch = run_multiple_pair_rdfs(
        frames=frames,
        pair_settings=pair_settings,
        box_length=12.0,
        use_pbc=True,
        smoothing_window=5,
    )
    rdf_results_only = {pair_key: payload[0] for pair_key, payload in rdf_batch.items()}
    structure_batch = compute_structure_factor_batch(
        rdf_results=rdf_results_only,
        k_min=0.1,
        k_max=25.0,
        k_step=0.1,
    )
    summary = export_structure_factor_batch(
        batch_results=structure_batch,
        output_directory=results_dir,
        file_prefix="structure_factor_water",
        summary_filename="structure_factor_summary.csv",
        title_prefix="Water Static Structure Factor",
    )

    print("Water structure factor batch completed.")
    print(f"Results directory: {results_dir}")
    print("Summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
