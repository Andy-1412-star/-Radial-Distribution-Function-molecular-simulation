"""Demonstrate RDF calculation directly from a LAMMPS dump trajectory."""

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
    """Parse command-line arguments for the LAMMPS RDF example."""

    parser = argparse.ArgumentParser(description="Compute RDF from a simple LAMMPS dump.")
    parser.add_argument(
        "--dump-file",
        type=Path,
        default=PROJECT_ROOT / "data" / "example_lammps.dump",
        help="Path to the LAMMPS dump file.",
    )
    parser.add_argument("--atom-a", default="O", help="Center atom type.")
    parser.add_argument("--atom-b", default="O", help="Neighbor atom type.")
    parser.add_argument("--type-map", nargs="*", default=["1:O", "2:H"], help="Mappings such as 1:O 2:H.")
    parser.add_argument("--r-max", type=float, default=6.0, help="Maximum RDF radius.")
    parser.add_argument("--bin-width", type=float, default=0.05, help="Histogram bin width.")
    return parser.parse_args()


def parse_type_mapping(items: list[str]) -> dict[int, str]:
    """Convert CLI mapping strings into a LAMMPS type dictionary."""

    mapping: dict[int, str] = {}
    for item in items:
        if ":" not in item:
            raise ValueError(f"Invalid type mapping '{item}'. Use the form integer:symbol.")
        raw_key, symbol = item.split(":", maxsplit=1)
        mapping[int(raw_key)] = symbol
    return mapping


def main() -> None:
    """Read a LAMMPS dump and compute one RDF."""

    arguments = parse_arguments()
    results_dir = ensure_directory(PROJECT_ROOT / "results")
    frames = load_frames(
        arguments.dump_file,
        file_format="lammps-dump",
        atom_type_map=parse_type_mapping(arguments.type_map),
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
    csv_path = results_dir / f"rdf_lammps_{arguments.atom_a}_{arguments.atom_b}.csv"
    png_path = results_dir / f"rdf_lammps_{arguments.atom_a}_{arguments.atom_b}.png"
    export_rdf_results(
        result=result,
        analysis=analysis,
        csv_path=csv_path,
        figure_path=png_path,
        title=f"LAMMPS RDF: {arguments.atom_a}-{arguments.atom_b}",
    )

    print("LAMMPS dump RDF example completed.")
    print(f"Frames used: {result.frames_used}")
    print(f"Pair type: {result.pair_label}")
    print(f"CSV saved to: {csv_path}")
    print(f"Figure saved to: {png_path}")
    if analysis.first_peak_r is not None:
        print(f"First peak position: r = {analysis.first_peak_r:.4f}")
    if analysis.coordination_number is not None:
        print(f"First-shell coordination number: {analysis.coordination_number:.4f}")


if __name__ == "__main__":
    main()
