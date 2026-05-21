"""Demonstrate RDF calculation directly from a GROMACS GRO trajectory."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from trajectory import export_rdf_results, load_frames, run_rdf_analysis
from utils import ensure_directory


def parse_arguments() -> argparse.Namespace:
    """Parse CLI arguments for GRO or MDAnalysis-backed GROMACS input."""

    parser = argparse.ArgumentParser(description="Compute RDF from GROMACS-style files.")
    parser.add_argument(
        "--topology-file",
        type=Path,
        default=PROJECT_ROOT / "data" / "example_water.gro",
        help="Topology file such as .gro or .tpr.",
    )
    parser.add_argument(
        "--trajectory-file",
        type=Path,
        default=None,
        help="Optional trajectory file such as .xtc or .trr. If omitted, the topology file is treated as a GRO trajectory.",
    )
    parser.add_argument("--atom-a", default="O", help="Center atom type.")
    parser.add_argument("--atom-b", default="H", help="Neighbor atom type.")
    parser.add_argument("--r-max", type=float, default=0.4, help="Maximum RDF radius in nm.")
    parser.add_argument("--bin-width", type=float, default=0.005, help="Histogram bin width in nm.")
    return parser.parse_args()


def main() -> None:
    """Read a GROMACS-compatible trajectory and compute one RDF."""

    arguments = parse_arguments()
    results_dir = ensure_directory(PROJECT_ROOT / "results")

    if arguments.trajectory_file is None:
        frames = load_frames(arguments.topology_file, file_format="gro")
    else:
        frames = load_frames(
            arguments.trajectory_file,
            file_format=arguments.trajectory_file.suffix.lstrip("."),
            topology_file=arguments.topology_file,
        )

    result, analysis = run_rdf_analysis(
        frames=frames,
        atom_type_a=arguments.atom_a,
        atom_type_b=arguments.atom_b,
        r_max=arguments.r_max,
        bin_width=arguments.bin_width,
        use_pbc=True,
        smoothing_window=5,
    )
    csv_path = results_dir / f"rdf_gromacs_{arguments.atom_a}_{arguments.atom_b}.csv"
    png_path = results_dir / f"rdf_gromacs_{arguments.atom_a}_{arguments.atom_b}.png"
    export_rdf_results(
        result=result,
        analysis=analysis,
        csv_path=csv_path,
        figure_path=png_path,
        title=f"GROMACS RDF: {arguments.atom_a}-{arguments.atom_b}",
    )

    print("GROMACS GRO RDF example completed.")
    print(f"Frames used: {result.frames_used}")
    print(f"Pair type: {result.pair_label}")
    print("Units: nm")
    print(f"CSV saved to: {csv_path}")
    print(f"Figure saved to: {png_path}")
    if analysis.first_peak_r is not None:
        print(f"First peak position: r = {analysis.first_peak_r:.4f} nm")
    if analysis.coordination_number is not None:
        print(f"First-shell coordination number: {analysis.coordination_number:.4f}")


if __name__ == "__main__":
    main()
