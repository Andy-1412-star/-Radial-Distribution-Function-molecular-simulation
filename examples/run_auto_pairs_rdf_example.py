"""Automatically detect atomic species and compute all unique pair RDFs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from trajectory import (
    build_auto_pair_jobs,
    detect_species,
    export_batch_rdf_job_results,
    load_frames,
    load_pair_parameter_config,
    run_pair_jobs,
)
from utils import ensure_directory


def parse_arguments() -> argparse.Namespace:
    """Parse CLI arguments for automatic all-pairs RDF generation."""

    parser = argparse.ArgumentParser(
        description="Automatically detect species and compute RDFs for all unique pairs."
    )
    parser.add_argument(
        "--input-file",
        type=Path,
        default=PROJECT_ROOT / "data" / "example_water.xyz",
        help="Input structure or trajectory file.",
    )
    parser.add_argument(
        "--file-format",
        default="xyz",
        help="Input format: xyz, lammps-dump, gro, xtc, or trr.",
    )
    parser.add_argument(
        "--default-r-max",
        type=float,
        default=6.0,
        help="Default RDF cutoff for automatically generated pairs.",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.05,
        help="Default RDF bin width for automatically generated pairs.",
    )
    parser.add_argument(
        "--pair-config",
        type=Path,
        default=PROJECT_ROOT / "data" / "example_pair_parameters_water.json",
        help=(
            "Optional JSON configuration file that provides default RDF parameters "
            "and pair-specific overrides."
        ),
    )
    parser.add_argument(
        "--box-length",
        type=float,
        default=12.0,
        help="Optional cubic box length override.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "auto_pairs",
        help="Directory where all RDF CSVs, PNGs, and the summary table are written.",
    )
    return parser.parse_args()


def main() -> None:
    """Run all unique pair RDFs for an automatically detected species list."""

    arguments = parse_arguments()
    frames = load_frames(arguments.input_file, file_format=arguments.file_format)
    species = detect_species(frames)
    default_r_max = arguments.default_r_max
    default_bin_width = arguments.bin_width
    default_smoothing_window = 5
    default_title_prefix = "Automatic Pair RDF"
    default_output_prefix = "rdf"
    pair_overrides = None
    pair_metadata_overrides = None

    if arguments.pair_config is not None:
        config = load_pair_parameter_config(arguments.pair_config)
        if config["default_r_max"] is not None:
            default_r_max = config["default_r_max"]
        if config["default_bin_width"] is not None:
            default_bin_width = config["default_bin_width"]
        if config["default_smoothing_window"] is not None:
            default_smoothing_window = config["default_smoothing_window"]
        if config["default_title_prefix"] is not None:
            default_title_prefix = config["default_title_prefix"]
        if config["default_output_prefix"] is not None:
            default_output_prefix = config["default_output_prefix"]
        pair_overrides = config["pair_parameter_overrides"]
        pair_metadata_overrides = config["pair_metadata_overrides"]

    pair_jobs = build_auto_pair_jobs(
        frames=frames,
        default_r_max=default_r_max,
        default_bin_width=default_bin_width,
        pair_parameter_overrides=pair_overrides,
        default_smoothing_window=default_smoothing_window,
        pair_metadata_overrides=pair_metadata_overrides,
        default_title_prefix=default_title_prefix,
        default_output_prefix=default_output_prefix,
    )
    batch_results = run_pair_jobs(
        frames=frames,
        pair_jobs=pair_jobs,
        box_length=arguments.box_length,
        use_pbc=True,
        default_smoothing_window=default_smoothing_window,
    )

    ensure_directory(arguments.output_dir)
    summary = export_batch_rdf_job_results(
        batch_results=batch_results,
        output_directory=arguments.output_dir,
        default_file_prefix=default_output_prefix,
        summary_filename="rdf_summary.csv",
        default_title_prefix=default_title_prefix,
    )

    print("Automatic all-pairs RDF calculation completed.")
    print(f"Detected species: {', '.join(species)}")
    print(f"Unique pairs computed: {len(pair_jobs)}")
    if arguments.pair_config is not None:
        print(f"Pair-parameter config: {arguments.pair_config}")
    print(f"Results directory: {arguments.output_dir}")
    print("Summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
