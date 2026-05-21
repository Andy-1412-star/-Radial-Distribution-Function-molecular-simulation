"""Simulate a small Lennard-Jones fluid and analyze RDF plus S(k)."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from io_utils import write_xyz
from lj_simulation import simulate_lj_fluid
from structure_factor import (
    compute_structure_factor_from_rdf,
    plot_structure_factor,
    save_structure_factor_csv,
)
from trajectory import export_rdf_results, run_rdf_analysis
from utils import ensure_directory


def main() -> None:
    """Run a small LJ Monte Carlo simulation and analyze its structure."""

    results_dir = ensure_directory(PROJECT_ROOT / "results" / "lj_fluid")
    simulation = simulate_lj_fluid(
        n_atoms=32,
        reduced_density=0.75,
        reduced_temperature=1.2,
        n_equilibration_sweeps=120,
        n_production_sweeps=240,
        max_displacement=0.18,
        sample_interval=12,
        random_seed=7,
    )

    trajectory_path = results_dir / "lj_fluid.xyz"
    write_xyz(simulation.frames, trajectory_path)

    rdf_result, rdf_analysis = run_rdf_analysis(
        frames=simulation.frames,
        atom_type_a="LJ",
        atom_type_b="LJ",
        r_max=3.5,
        bin_width=0.05,
        box_length=simulation.box_length,
        use_pbc=True,
        smoothing_window=5,
    )
    export_rdf_results(
        result=rdf_result,
        analysis=rdf_analysis,
        csv_path=results_dir / "rdf_LJ_LJ.csv",
        figure_path=results_dir / "rdf_LJ_LJ.png",
        title="Lennard-Jones Fluid RDF",
    )

    structure_result = compute_structure_factor_from_rdf(
        rdf_result=rdf_result,
        k_min=0.1,
        k_max=20.0,
        k_step=0.1,
    )
    save_structure_factor_csv(structure_result, results_dir / "structure_factor_LJ_LJ.csv")
    plot_structure_factor(
        structure_result,
        output_path=results_dir / "structure_factor_LJ_LJ.png",
        title="Lennard-Jones Fluid Structure Factor",
    )

    print("Lennard-Jones fluid example completed.")
    print(f"Sampled frames: {len(simulation.frames)}")
    print(f"Acceptance ratio: {simulation.acceptance_ratio:.4f}")
    print(f"Trajectory saved to: {trajectory_path}")
    print(f"RDF first peak: r = {rdf_analysis.first_peak_r:.4f}" if rdf_analysis.first_peak_r is not None else "RDF first peak: not found")
    print(f"S(k) CSV saved to: {results_dir / 'structure_factor_LJ_LJ.csv'}")
    print(f"S(k) figure saved to: {results_dir / 'structure_factor_LJ_LJ.png'}")


if __name__ == "__main__":
    main()
