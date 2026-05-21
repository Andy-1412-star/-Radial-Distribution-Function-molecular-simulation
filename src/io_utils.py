"""Input/output helpers for XYZ, LAMMPS, and GROMACS trajectories."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np

from utils import AtomFrame


def parse_xyz_comment(comment: str) -> Dict[str, str]:
    """Parse key-value metadata from an XYZ comment line."""

    metadata: Dict[str, str] = {}
    for part in comment.split(";"):
        item = part.strip()
        if "=" not in item:
            continue
        key, value = item.split("=", maxsplit=1)
        metadata[key.strip()] = value.strip()
    return metadata


def read_xyz(file_path: str | Path) -> List[AtomFrame]:
    """Read an XYZ file that may contain one frame or many frames."""

    path = Path(file_path)
    lines = path.read_text(encoding="utf-8").splitlines()

    frames: List[AtomFrame] = []
    line_index = 0
    while line_index < len(lines):
        if not lines[line_index].strip():
            line_index += 1
            continue

        n_atoms = int(lines[line_index].strip())
        if line_index + 1 >= len(lines):
            raise ValueError(f"Missing comment line near line {line_index + 1}.")

        comment_line = lines[line_index + 1].strip()
        metadata = parse_xyz_comment(comment_line)

        start = line_index + 2
        end = start + n_atoms
        if end > len(lines):
            raise ValueError("XYZ file ended before all coordinates were read.")

        symbols = []
        positions = []
        for raw_line in lines[start:end]:
            parts = raw_line.split()
            if len(parts) < 4:
                raise ValueError(f"Malformed XYZ coordinate line: {raw_line}")
            symbols.append(parts[0])
            positions.append([float(parts[1]), float(parts[2]), float(parts[3])])

        frames.append(
            AtomFrame(
                symbols=np.asarray(symbols, dtype=str),
                positions=np.asarray(positions, dtype=float),
                metadata=metadata,
            )
        )
        line_index = end

    if not frames:
        raise ValueError(f"No frames found in XYZ file: {path}")

    return frames


def write_xyz(frames: Iterable[AtomFrame], file_path: str | Path) -> None:
    """Write one or more frames to a multi-frame XYZ file."""

    path = Path(file_path)
    lines: List[str] = []
    for frame_index, frame in enumerate(frames, start=1):
        n_atoms = len(frame.symbols)
        if n_atoms != len(frame.positions):
            raise ValueError("Frame symbols and positions must have the same length.")

        if frame.metadata:
            comment = "; ".join(f"{key}={value}" for key, value in frame.metadata.items())
        else:
            comment = f"Frame {frame_index}"

        lines.append(str(n_atoms))
        lines.append(comment)
        for symbol, position in zip(frame.symbols, frame.positions):
            lines.append(
                f"{symbol} {float(position[0]):.6f} {float(position[1]):.6f} {float(position[2]):.6f}"
            )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_lammps_dump(
    file_path: str | Path,
    atom_type_map: Optional[Dict[str | int, str]] = None,
) -> List[AtomFrame]:
    """Read a simple orthogonal LAMMPS dump trajectory.

    Supported atom columns:

    - ``element x y z``
    - ``type x y z``
    - ``type xs ys zs`` with scaled coordinates

    Parameters
    ----------
    atom_type_map
        Optional mapping from LAMMPS type IDs to human-readable symbols,
        for example ``{1: "O", 2: "H"}``.
    """

    path = Path(file_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    frames: List[AtomFrame] = []
    i = 0

    while i < len(lines):
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue

        timestep = lines[i + 1].strip()
        if lines[i + 2].strip() != "ITEM: NUMBER OF ATOMS":
            raise ValueError("Unexpected LAMMPS dump format: missing NUMBER OF ATOMS.")
        n_atoms = int(lines[i + 3].strip())

        if not lines[i + 4].startswith("ITEM: BOX BOUNDS"):
            raise ValueError("Unexpected LAMMPS dump format: missing BOX BOUNDS.")

        bounds = []
        for offset in range(3):
            low, high, *_ = lines[i + 5 + offset].split()
            bounds.append((float(low), float(high)))
        box_lengths = np.array([high - low for low, high in bounds], dtype=float)

        atom_header = lines[i + 8].strip()
        if not atom_header.startswith("ITEM: ATOMS"):
            raise ValueError("Unexpected LAMMPS dump format: missing ATOMS section.")
        columns = atom_header.split()[2:]
        column_index = {name: idx for idx, name in enumerate(columns)}

        symbols = []
        positions = []
        atom_start = i + 9
        atom_end = atom_start + n_atoms
        if atom_end > len(lines):
            raise ValueError("LAMMPS dump ended before all atom lines were read.")

        for raw_line in lines[atom_start:atom_end]:
            parts = raw_line.split()
            symbol = _symbol_from_lammps_parts(parts, column_index, atom_type_map)
            position = _position_from_lammps_parts(parts, column_index, bounds, box_lengths)
            symbols.append(symbol)
            positions.append(position)

        metadata = {
            "box_length_x": f"{box_lengths[0]:.8f}",
            "box_length_y": f"{box_lengths[1]:.8f}",
            "box_length_z": f"{box_lengths[2]:.8f}",
            "timestep": timestep,
        }
        frames.append(
            AtomFrame(
                symbols=np.asarray(symbols, dtype=str),
                positions=np.asarray(positions, dtype=float),
                metadata=metadata,
            )
        )
        i = atom_end

    if not frames:
        raise ValueError(f"No frames found in LAMMPS dump file: {path}")

    return frames


def read_gro(file_path: str | Path) -> List[AtomFrame]:
    """Read a GROMACS ``.gro`` file with one or multiple frames.

    Notes
    -----
    GRO coordinates are stored in nanometers. This reader keeps those units.
    Therefore, when you compute RDF from GRO data, choose ``r_max`` and
    ``bin_width`` in nanometers as well.
    """

    path = Path(file_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    frames: List[AtomFrame] = []
    i = 0

    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue

        title = lines[i].strip()
        n_atoms = int(lines[i + 1].strip())
        atom_start = i + 2
        atom_end = atom_start + n_atoms
        box_line_index = atom_end
        if box_line_index >= len(lines):
            raise ValueError("GRO file ended before the box line was read.")

        symbols = []
        positions = []
        for raw_line in lines[atom_start:atom_end]:
            atom_name = raw_line[10:15].strip()
            x = float(raw_line[20:28].strip())
            y = float(raw_line[28:36].strip())
            z = float(raw_line[36:44].strip())
            symbols.append(_infer_symbol_from_atom_name(atom_name))
            positions.append([x, y, z])

        box_values = [float(value) for value in lines[box_line_index].split()]
        if len(box_values) < 3:
            raise ValueError("GRO box line must contain at least three numbers.")

        metadata = {
            "box_length_x": f"{box_values[0]:.8f}",
            "box_length_y": f"{box_values[1]:.8f}",
            "box_length_z": f"{box_values[2]:.8f}",
            "title": title,
            "units": "nm",
        }
        frames.append(
            AtomFrame(
                symbols=np.asarray(symbols, dtype=str),
                positions=np.asarray(positions, dtype=float),
                metadata=metadata,
            )
        )
        i = box_line_index + 1

    if not frames:
        raise ValueError(f"No frames found in GRO file: {path}")

    return frames


def read_gromacs_with_mdanalysis(
    topology_file: str | Path,
    trajectory_file: Optional[str | Path] = None,
) -> List[AtomFrame]:
    """Read GROMACS-compatible trajectories through MDAnalysis.

    This optional route is useful for formats such as ``.xtc`` and ``.trr``.
    It is only used if the user has installed ``MDAnalysis`` separately.
    """

    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise ImportError(
            "GROMACS trajectory support requires MDAnalysis. Please install it using "
            "pip install MDAnalysis."
        ) from exc

    universe = mda.Universe(str(topology_file), str(trajectory_file) if trajectory_file else None)
    frames: List[AtomFrame] = []

    for ts in universe.trajectory:
        dimensions = ts.dimensions[:3]
        metadata = {
            "box_length_x": f"{dimensions[0]:.8f}",
            "box_length_y": f"{dimensions[1]:.8f}",
            "box_length_z": f"{dimensions[2]:.8f}",
            "frame": str(ts.frame),
            "source": "MDAnalysis",
        }
        symbols = _symbols_from_mdanalysis(universe)
        positions = universe.atoms.positions.astype(float).copy()
        frames.append(AtomFrame(symbols=symbols, positions=positions, metadata=metadata))

    if not frames:
        raise ValueError("No frames were read by MDAnalysis.")

    return frames


def load_trajectory(
    file_path: str | Path,
    file_format: str,
    atom_type_map: Optional[Dict[str | int, str]] = None,
    topology_file: Optional[str | Path] = None,
) -> List[AtomFrame]:
    """Generic dispatcher for supported trajectory formats."""

    normalized = file_format.lower()
    if normalized == "xyz":
        return read_xyz(file_path)
    if normalized in {"lammps", "dump", "lammps-dump"}:
        return read_lammps_dump(file_path, atom_type_map=atom_type_map)
    if normalized in {"gro", "gromacs-gro"}:
        return read_gro(file_path)
    if normalized in {"xtc", "trr", "gromacs"}:
        if topology_file is None:
            raise ValueError(
                "Reading GROMACS .xtc or .trr files requires a topology file "
                "such as .gro or .tpr."
            )
        return read_gromacs_with_mdanalysis(topology_file, file_path)
    raise ValueError(f"Unsupported file format: {file_format}")


def frame_box_lengths(frame: AtomFrame) -> Optional[np.ndarray]:
    """Extract box lengths from metadata if available."""

    if "box_length" in frame.metadata:
        box = float(frame.metadata["box_length"])
        return np.array([box, box, box], dtype=float)

    keys = ("box_length_x", "box_length_y", "box_length_z")
    if all(key in frame.metadata for key in keys):
        return np.array([float(frame.metadata[key]) for key in keys], dtype=float)

    return None


def _symbol_from_lammps_parts(
    parts: List[str],
    column_index: Dict[str, int],
    atom_type_map: Optional[Dict[str | int, str]],
) -> str:
    """Infer the species label for one LAMMPS atom record."""

    if "element" in column_index:
        return parts[column_index["element"]]

    if "type" not in column_index:
        raise ValueError("LAMMPS ATOMS section must contain either 'element' or 'type'.")

    raw_type = parts[column_index["type"]]
    if atom_type_map is None:
        return raw_type

    if raw_type in atom_type_map:
        return atom_type_map[raw_type]

    try:
        integer_type = int(raw_type)
    except ValueError:
        integer_type = None
    if integer_type is not None and integer_type in atom_type_map:
        return atom_type_map[integer_type]

    return raw_type


def _position_from_lammps_parts(
    parts: List[str],
    column_index: Dict[str, int],
    bounds: Iterable[tuple[float, float]],
    box_lengths: np.ndarray,
) -> np.ndarray:
    """Extract Cartesian position from one LAMMPS atom record."""

    if all(key in column_index for key in ("x", "y", "z")):
        return np.array([float(parts[column_index["x"]]), float(parts[column_index["y"]]), float(parts[column_index["z"]])], dtype=float)

    if all(key in column_index for key in ("xu", "yu", "zu")):
        return np.array([float(parts[column_index["xu"]]), float(parts[column_index["yu"]]), float(parts[column_index["zu"]])], dtype=float)

    if all(key in column_index for key in ("xs", "ys", "zs")):
        lows = np.array([low for low, _ in bounds], dtype=float)
        scaled = np.array(
            [float(parts[column_index["xs"]]), float(parts[column_index["ys"]]), float(parts[column_index["zs"]])],
            dtype=float,
        )
        return lows + scaled * box_lengths

    raise ValueError("LAMMPS ATOMS section must contain x/y/z, xu/yu/zu, or xs/ys/zs.")


def _infer_symbol_from_atom_name(atom_name: str) -> str:
    """Infer a chemical symbol from a GRO atom name."""

    letters = "".join(char for char in atom_name if char.isalpha())
    if not letters:
        return atom_name.strip() or "X"

    upper = letters.upper()
    if upper.startswith("OW") or upper == "O":
        return "O"
    if upper.startswith("HW") or upper == "H":
        return "H"

    if len(letters) == 1:
        return letters.upper()
    return letters[0].upper() + letters[1:].lower()


def _symbols_from_mdanalysis(universe) -> np.ndarray:
    """Extract per-atom symbols from an MDAnalysis universe."""

    if hasattr(universe.atoms, "elements") and len(universe.atoms.elements) == len(universe.atoms):
        return np.asarray(universe.atoms.elements, dtype=str)
    return np.asarray([_infer_symbol_from_atom_name(name) for name in universe.atoms.names], dtype=str)
