"""High-level trajectory and RDF workflows for example scripts and notebooks."""

from __future__ import annotations

from itertools import combinations_with_replacement
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from analysis import RDFAnalysis, analyze_rdf
from io_utils import load_trajectory
from plotting import plot_rdf
from rdf import RDFResult, compute_rdf
from utils import AtomFrame, ensure_directory


PairParameterMap = Dict[Tuple[str, str], Tuple[float, float]]
PairMetadataMap = Dict[Tuple[str, str], Dict[str, Any]]


def load_frames(
    file_path: str | Path,
    file_format: str,
    atom_type_map: Optional[Dict[str | int, str]] = None,
    topology_file: Optional[str | Path] = None,
) -> List[AtomFrame]:
    """Load trajectory frames from one of the supported formats."""

    return load_trajectory(
        file_path=file_path,
        file_format=file_format,
        atom_type_map=atom_type_map,
        topology_file=topology_file,
    )


def run_rdf_analysis(
    frames: List[AtomFrame],
    atom_type_a: str,
    atom_type_b: str,
    r_max: float,
    bin_width: float,
    box_length: Optional[float | Iterable[float]] = None,
    use_pbc: bool = True,
    smoothing_window: int = 5,
) -> Tuple[RDFResult, RDFAnalysis]:
    """Compute an RDF and immediately derive peak/minimum/CN analysis."""

    rdf_result = compute_rdf(
        frames=frames,
        atom_type_a=atom_type_a,
        atom_type_b=atom_type_b,
        r_max=r_max,
        bin_width=bin_width,
        box_length=box_length,
        use_pbc=use_pbc,
    )
    rdf_analysis = analyze_rdf(rdf_result, smoothing_window=smoothing_window)
    return rdf_result, rdf_analysis


def run_multiple_pair_rdfs(
    frames: List[AtomFrame],
    pair_settings: Sequence[Tuple[str, str, float, float]],
    box_length: Optional[float | Iterable[float]] = None,
    use_pbc: bool = True,
    smoothing_window: int = 5,
) -> Dict[Tuple[str, str], Tuple[RDFResult, RDFAnalysis]]:
    """Compute a batch of pair RDFs for one trajectory."""

    results: Dict[Tuple[str, str], Tuple[RDFResult, RDFAnalysis]] = {}
    for atom_type_a, atom_type_b, r_max, bin_width in pair_settings:
        results[(atom_type_a, atom_type_b)] = run_rdf_analysis(
            frames=frames,
            atom_type_a=atom_type_a,
            atom_type_b=atom_type_b,
            r_max=r_max,
            bin_width=bin_width,
            box_length=box_length,
            use_pbc=use_pbc,
            smoothing_window=smoothing_window,
        )
    return results


def run_pair_jobs(
    frames: List[AtomFrame],
    pair_jobs: Sequence[Dict[str, Any]],
    box_length: Optional[float | Iterable[float]] = None,
    use_pbc: bool = True,
    default_smoothing_window: int = 5,
) -> Dict[Tuple[str, str], Dict[str, Any]]:
    """Compute a batch of RDF jobs that may carry pair-specific metadata."""

    results: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for job in pair_jobs:
        atom_type_a = str(job["atom_type_a"])
        atom_type_b = str(job["atom_type_b"])
        smoothing_window = int(job.get("smoothing_window", default_smoothing_window))
        rdf_result, rdf_analysis = run_rdf_analysis(
            frames=frames,
            atom_type_a=atom_type_a,
            atom_type_b=atom_type_b,
            r_max=float(job["r_max"]),
            bin_width=float(job["bin_width"]),
            box_length=box_length,
            use_pbc=use_pbc,
            smoothing_window=smoothing_window,
        )
        results[(atom_type_a, atom_type_b)] = {
            "result": rdf_result,
            "analysis": rdf_analysis,
            "title": job.get("title"),
            "output_stem": job.get("output_stem"),
            "smoothing_window": smoothing_window,
        }
    return results


def detect_species(frames: List[AtomFrame]) -> List[str]:
    """Return unique species labels in order of first appearance."""

    if not frames:
        raise ValueError("Cannot detect species from an empty trajectory.")

    seen = set()
    species: List[str] = []
    for frame in frames:
        for symbol in frame.symbols.tolist():
            if symbol not in seen:
                seen.add(symbol)
                species.append(symbol)
    if not species:
        raise ValueError("No atomic species were found in the provided frames.")
    return species


def build_auto_pair_settings(
    frames: List[AtomFrame],
    default_r_max: float,
    default_bin_width: float,
    pair_parameter_overrides: Optional[PairParameterMap] = None,
) -> List[Tuple[str, str, float, float]]:
    """Create all unique species-pair RDF settings automatically.

    The returned list includes same-species pairs and mixed pairs once each:
    ``A-A, A-B, B-B`` rather than duplicating ``B-A``.
    """

    if default_r_max <= 0.0:
        raise ValueError("default_r_max must be positive.")
    if default_bin_width <= 0.0:
        raise ValueError("default_bin_width must be positive.")

    species = detect_species(frames)
    overrides = pair_parameter_overrides or {}
    pair_settings: List[Tuple[str, str, float, float]] = []

    for atom_type_a, atom_type_b in combinations_with_replacement(species, 2):
        parameters = overrides.get((atom_type_a, atom_type_b))
        if parameters is None:
            parameters = overrides.get((atom_type_b, atom_type_a))
        r_max, bin_width = parameters if parameters is not None else (default_r_max, default_bin_width)
        pair_settings.append((atom_type_a, atom_type_b, r_max, bin_width))

    return pair_settings


def build_auto_pair_jobs(
    frames: List[AtomFrame],
    default_r_max: float,
    default_bin_width: float,
    pair_parameter_overrides: Optional[PairParameterMap] = None,
    default_smoothing_window: int = 5,
    pair_metadata_overrides: Optional[PairMetadataMap] = None,
    default_title_prefix: Optional[str] = None,
    default_output_prefix: str = "rdf",
) -> List[Dict[str, Any]]:
    """Create all unique species-pair RDF jobs with per-pair metadata support."""

    pair_settings = build_auto_pair_settings(
        frames=frames,
        default_r_max=default_r_max,
        default_bin_width=default_bin_width,
        pair_parameter_overrides=pair_parameter_overrides,
    )
    metadata_overrides = pair_metadata_overrides or {}
    jobs: List[Dict[str, Any]] = []

    for atom_type_a, atom_type_b, r_max, bin_width in pair_settings:
        metadata = metadata_overrides.get((atom_type_a, atom_type_b))
        if metadata is None:
            metadata = metadata_overrides.get((atom_type_b, atom_type_a), {})
        output_stem = metadata.get("output_stem", f"{default_output_prefix}_{atom_type_a}_{atom_type_b}")
        title = metadata.get("title")
        if title is None and default_title_prefix:
            title = f"{default_title_prefix}: {atom_type_a}-{atom_type_b}"

        jobs.append(
            {
                "atom_type_a": atom_type_a,
                "atom_type_b": atom_type_b,
                "r_max": r_max,
                "bin_width": bin_width,
                "smoothing_window": int(metadata.get("smoothing_window", default_smoothing_window)),
                "title": title,
                "output_stem": output_stem,
            }
        )

    return jobs


def load_pair_parameter_config(
    config_path: str | Path,
) -> Dict[str, Any]:
    """Read per-pair RDF parameters from a JSON configuration file.

    Expected structure
    ------------------
    {
      "default": {
        "r_max": 6.0,
        "bin_width": 0.05,
        "smoothing_window": 5,
        "title_prefix": "Automatic Pair RDF",
        "output_prefix": "rdf"
      },
      "pairs": {
        "O-O": {"r_max": 6.0, "bin_width": 0.05},
        "O-H": {
          "r_max": 4.0,
          "bin_width": 0.02,
          "smoothing_window": 7,
          "title": "Water O-H RDF",
          "output_stem": "rdf_water_O_H"
        }
      }
    }
    """

    path = Path(config_path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Pair-parameter config file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in pair-parameter config file: {path}") from exc

    if not isinstance(payload, dict):
        raise ValueError("Pair-parameter config must be a JSON object.")

    default_section = payload.get("default", {})
    if default_section is None:
        default_section = {}
    if not isinstance(default_section, dict):
        raise ValueError("The 'default' section of the pair config must be an object.")

    default_r_max = _optional_positive_float(default_section.get("r_max"), "default.r_max")
    default_bin_width = _optional_positive_float(default_section.get("bin_width"), "default.bin_width")
    default_smoothing_window = _optional_positive_int(
        default_section.get("smoothing_window"),
        "default.smoothing_window",
    )
    default_title_prefix = _optional_string(default_section.get("title_prefix"), "default.title_prefix")
    default_output_prefix = _optional_string(default_section.get("output_prefix"), "default.output_prefix")

    pairs_section = payload.get("pairs", {})
    if not isinstance(pairs_section, dict):
        raise ValueError("The 'pairs' section of the pair config must be an object.")

    overrides: PairParameterMap = {}
    metadata_overrides: PairMetadataMap = {}
    for pair_label, parameters in pairs_section.items():
        if not isinstance(pair_label, str) or "-" not in pair_label:
            raise ValueError(
                "Each pair override key must look like 'A-B', for example 'O-H' or 'Ar-Ar'."
            )
        if not isinstance(parameters, dict):
            raise ValueError(f"Pair override '{pair_label}' must map to an object.")

        atom_type_a, atom_type_b = [item.strip() for item in pair_label.split("-", maxsplit=1)]
        if not atom_type_a or not atom_type_b:
            raise ValueError(f"Pair override '{pair_label}' is invalid. Use non-empty labels.")

        r_max = _required_positive_float(parameters.get("r_max"), f"pairs.{pair_label}.r_max")
        bin_width = _required_positive_float(
            parameters.get("bin_width"),
            f"pairs.{pair_label}.bin_width",
        )
        overrides[(atom_type_a, atom_type_b)] = (r_max, bin_width)
        metadata_overrides[(atom_type_a, atom_type_b)] = {
            "smoothing_window": _optional_positive_int(
                parameters.get("smoothing_window"),
                f"pairs.{pair_label}.smoothing_window",
            ),
            "title": _optional_string(parameters.get("title"), f"pairs.{pair_label}.title"),
            "output_stem": _optional_string(
                parameters.get("output_stem"),
                f"pairs.{pair_label}.output_stem",
            ),
        }

    return {
        "default_r_max": default_r_max,
        "default_bin_width": default_bin_width,
        "default_smoothing_window": default_smoothing_window,
        "default_title_prefix": default_title_prefix,
        "default_output_prefix": default_output_prefix,
        "pair_parameter_overrides": overrides,
        "pair_metadata_overrides": metadata_overrides,
    }


def summarize_batch_results(
    batch_results: Dict[Tuple[str, str], Tuple[RDFResult, RDFAnalysis]],
) -> pd.DataFrame:
    """Create a compact summary table for a set of RDF calculations."""

    rows = []
    for (atom_type_a, atom_type_b), (result, analysis) in batch_results.items():
        rows.append(
            {
                "pair": result.pair_label,
                "atom_type_a": atom_type_a,
                "atom_type_b": atom_type_b,
                "frames_used": result.frames_used,
                "first_peak_r": analysis.first_peak_r,
                "first_peak_g": analysis.first_peak_g,
                "first_minimum_r": analysis.first_minimum_r,
                "first_minimum_g": analysis.first_minimum_g,
                "coordination_number": analysis.coordination_number,
            }
        )
    return pd.DataFrame(rows).reset_index(drop=True)


def rdf_dataframe(result: RDFResult, analysis: RDFAnalysis) -> pd.DataFrame:
    """Return a tidy dataframe ready to be written to CSV."""

    return pd.DataFrame(
        {
            "r": result.r,
            "g_r": result.g_r,
            "pair_counts": result.counts,
            "shell_volume": result.shell_volume,
            "coordination_number": analysis.cumulative_coordination,
        }
    )


def export_rdf_results(
    result: RDFResult,
    analysis: RDFAnalysis,
    csv_path: str | Path,
    figure_path: str | Path,
    title: Optional[str] = None,
) -> None:
    """Write CSV and figure outputs for one RDF result."""

    dataframe = rdf_dataframe(result, analysis)
    dataframe.to_csv(csv_path, index=False)
    plot_rdf(
        result=result,
        analysis=analysis,
        output_path=figure_path,
        title=title,
    )


def export_batch_rdf_results(
    batch_results: Dict[Tuple[str, str], Tuple[RDFResult, RDFAnalysis]],
    output_directory: str | Path,
    file_prefix: str = "rdf",
    summary_filename: str = "rdf_summary.csv",
    title_prefix: str = "Radial Distribution Function",
) -> pd.DataFrame:
    """Write CSV, PNG, and one summary CSV for a batch of pair RDFs."""

    output_dir = ensure_directory(output_directory)

    for (atom_type_a, atom_type_b), (result, analysis) in batch_results.items():
        stem = f"{file_prefix}_{atom_type_a}_{atom_type_b}"
        export_rdf_results(
            result=result,
            analysis=analysis,
            csv_path=output_dir / f"{stem}.csv",
            figure_path=output_dir / f"{stem}.png",
            title=f"{title_prefix}: {atom_type_a}-{atom_type_b}",
        )

    summary = summarize_batch_results(batch_results)
    summary.to_csv(output_dir / summary_filename, index=False)
    return summary


def export_batch_rdf_job_results(
    batch_results: Dict[Tuple[str, str], Dict[str, Any]],
    output_directory: str | Path,
    default_file_prefix: str = "rdf",
    summary_filename: str = "rdf_summary.csv",
    default_title_prefix: str = "Radial Distribution Function",
) -> pd.DataFrame:
    """Write CSV, PNG, and summary output for pair jobs with per-pair metadata."""

    output_dir = ensure_directory(output_directory)
    classic_batch: Dict[Tuple[str, str], Tuple[RDFResult, RDFAnalysis]] = {}

    for pair_key, payload in batch_results.items():
        result: RDFResult = payload["result"]
        analysis: RDFAnalysis = payload["analysis"]
        output_stem = payload.get("output_stem") or f"{default_file_prefix}_{pair_key[0]}_{pair_key[1]}"
        title = payload.get("title") or f"{default_title_prefix}: {pair_key[0]}-{pair_key[1]}"
        export_rdf_results(
            result=result,
            analysis=analysis,
            csv_path=output_dir / f"{output_stem}.csv",
            figure_path=output_dir / f"{output_stem}.png",
            title=title,
        )
        classic_batch[pair_key] = (result, analysis)

    summary = summarize_batch_results(classic_batch)
    summary.to_csv(output_dir / summary_filename, index=False)
    return summary


def _optional_positive_float(value: object, field_name: str) -> Optional[float]:
    """Validate an optional positive floating-point field."""

    if value is None:
        return None
    return _required_positive_float(value, field_name)


def _required_positive_float(value: object, field_name: str) -> float:
    """Validate a required positive floating-point field."""

    if value is None:
        raise ValueError(f"Missing required field '{field_name}' in pair-parameter config.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Field '{field_name}' must be a positive number.") from exc
    if number <= 0.0:
        raise ValueError(f"Field '{field_name}' must be positive.")
    return number


def _optional_positive_int(value: object, field_name: str) -> Optional[int]:
    """Validate an optional positive integer field."""

    if value is None:
        return None
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Field '{field_name}' must be a positive integer.") from exc
    if number <= 0:
        raise ValueError(f"Field '{field_name}' must be positive.")
    return number


def _optional_string(value: object, field_name: str) -> Optional[str]:
    """Validate an optional string field."""

    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"Field '{field_name}' must be a string.")
    cleaned = value.strip()
    if not cleaned:
        raise ValueError(f"Field '{field_name}' cannot be empty.")
    return cleaned
