#!/usr/bin/env python3
"""
extract_tddft_cubes.py

Recorre los confórmeros TDDFT del proyecto PAPER1 y genera, para cada uno:
  - Cubos HOMO/LUMO del estado fundamental (desde tddft.gbw)
  - Cubos hole/particle (NTO) de cada estado excitado (desde tddft.sN.nto)

Los pares NTO adicionales al dominante se incluyen solo si superan un umbral
relativo al coeficiente dominante (--nto-threshold), evitando generar cubos
redundantes en estados de carácter limpio y capturando estados con carácter
mixto (ej. 0.52/0.45) cuando corresponde. El "participation ratio" (PR) de
cada estado se reporta en el CSV final como diagnóstico de cuántos pares
NTO son realmente necesarios para describirlo.

Uso:
    python extract_tddft_cubes.py --root PAPER1 --dry-run
    python extract_tddft_cubes.py --root PAPER1 --grid 80 --force
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger("extract_tddft_cubes")

# --------------------------------------------------------------------------- #
# Estructuras de datos
# --------------------------------------------------------------------------- #


@dataclass
class NTOPair:
    hole: int
    particle: int
    coef: float


@dataclass
class ExcitedState:
    state_idx: int
    ntos: list[NTOPair] = field(default_factory=list)

    @property
    def participation_ratio(self) -> Optional[float]:
        """PR = 1 / sum(lambda_i^2). lambda_i normalizado a partir de coef."""
        if not self.ntos:
            return None
        total = sum(p.coef for p in self.ntos)
        if total <= 0:
            return None
        lambdas = [p.coef / total for p in self.ntos]
        denom = sum(l ** 2 for l in lambdas)
        if denom <= 0:
            return None
        return 1.0 / denom


@dataclass
class ConformerJob:
    conf_dir: Path
    molecule: str
    run_label: str
    conf_id: str
    homo_idx: Optional[int] = None
    lumo_idx: Optional[int] = None
    is_open_shell: bool = False
    states: dict[int, ExcitedState] = field(default_factory=dict)


@dataclass
class ReportRow:
    molecule: str
    run_label: str
    conf_id: str
    kind: str  # "ground_state" | "state1", "state2", ...
    homo_or_hole: Optional[int]
    lumo_or_particle: Optional[int]
    coef: Optional[float]
    participation_ratio: Optional[float]
    fosc: Optional[float]
    cube_files: str
    status: str
    error: str = ""


# --------------------------------------------------------------------------- #
# Descubrimiento de carpetas
# --------------------------------------------------------------------------- #


def discover_conformers(root: Path) -> list[Path]:
    """Encuentra carpetas conf_* bajo */runs/*/05_tddft/ con tddft.out + tddft.gbw."""
    pattern = "*/runs/*/05_tddft/conf_*"
    candidates = sorted(root.glob(pattern))
    valid = []
    for conf_dir in candidates:
        if not conf_dir.is_dir():
            continue
        out_file = conf_dir / "tddft.out"
        gbw_file = conf_dir / "tddft.gbw"
        if out_file.exists() and gbw_file.exists():
            valid.append(conf_dir)
        else:
            missing = [
                name
                for name, f in (("tddft.out", out_file), ("tddft.gbw", gbw_file))
                if not f.exists()
            ]
            logger.warning(
                "Saltando %s: faltan %s", conf_dir, ", ".join(missing)
            )
    return valid


def parse_conformer_path(conf_dir: Path, root: Path) -> tuple[str, str, str]:
    """Extrae (molecula, run_label, conf_id) de PAPER1/MOL/runs/B/05_tddft/conf_6."""
    rel = conf_dir.relative_to(root)
    parts = rel.parts
    # parts: (MOL, "runs", RUN_LABEL, "05_tddft", "conf_X")
    molecule = parts[0]
    run_label = parts[2] if len(parts) > 2 else "unknown"
    conf_id = conf_dir.name
    return molecule, run_label, conf_id


# --------------------------------------------------------------------------- #
# Parseo de tddft.out
# --------------------------------------------------------------------------- #

ORBITAL_ENERGIES_HEADER = re.compile(r"ORBITAL ENERGIES")
# Filas tipo:   6   2.0000    -0.522130      -14.2079
ORBITAL_ROW = re.compile(
    r"^\s*(\d+)\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)\s*$"
)
NTO_STATE_HEADER = re.compile(r"NATURAL TRANSITION ORBITALS FOR STATE\s+(\d+)")
# Filas tipo:   54a -> 55a  :  n= 0.523456
NTO_ROW = re.compile(
    r"^\s*(\d+)([ab])\s*->\s*(\d+)([ab])\s*:\s*n\s*=\s*([\-\d.]+)"
)

ABSORPTION_HEADER = re.compile(
    r"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
)
# Fila tipo:  0-1A  ->  1-1A    4.203068   33900.0   295.0   0.752016216 ...
ABSORPTION_ROW = re.compile(
    r"^\s*\d+-\S+\s*->\s*(\d+)-\S+\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)"
)


def parse_oscillator_strengths(out_text: str) -> dict[int, float]:
    """Extrae fosc por estado desde la primera tabla 'ABSORPTION SPECTRUM VIA
    TRANSITION ELECTRIC DIPOLE MOMENTS' (dipolo eléctrico, no velocidad ni CD).
    Devuelve {estado: fosc}."""
    lines = out_text.splitlines()
    header_positions = [i for i, l in enumerate(lines) if ABSORPTION_HEADER.search(l)]
    if not header_positions:
        return {}

    start = header_positions[0]
    fosc_by_state: dict[int, float] = {}
    dash_count = 0
    for line in lines[start:]:
        if "-----" in line:
            dash_count += 1
            if dash_count >= 3:
                break  # tercera línea de guiones = fin de la tabla de datos
            continue
        if dash_count < 2:
            continue  # todavía en encabezado (nombre, unidades)
        m = ABSORPTION_ROW.match(line)
        if m:
            state_idx = int(m.group(1))
            fosc = float(m.group(5))
            fosc_by_state[state_idx] = fosc
    return fosc_by_state


def select_dominant_state(
    fosc_by_state: dict[int, float],
    fosc_min: float,
) -> Optional[int]:
    """Selecciona el único estado 'dominante' (el de deltamax): el de mayor
    fosc, solo si ese máximo supera fosc_min. Si no hay tabla de fosc,
    devuelve None con significado especial de 'sin filtrar' (fallback inseguro
    ya no aplica aquí: sin tabla no hay forma de saber cuál es el dominante,
    así que no se extrae ningún NTO y se advierte)."""
    if not fosc_by_state:
        return None
    brightest = max(fosc_by_state, key=fosc_by_state.get)
    if fosc_by_state[brightest] <= fosc_min:
        return None
    return brightest


def parse_ground_state_frontier_orbitals(
    out_text: str,
) -> tuple[Optional[int], Optional[int], bool]:
    """Devuelve (homo_idx, lumo_idx, is_open_shell) desde el último bloque
    ORBITAL ENERGIES en tddft.out."""
    lines = out_text.splitlines()
    header_positions = [
        i for i, l in enumerate(lines) if ORBITAL_ENERGIES_HEADER.search(l)
    ]
    if not header_positions:
        return None, None, False

    # Tomar el bloque a partir de la última ocurrencia
    start = header_positions[-1]
    block = lines[start:]

    # Heurística open-shell: ORCA imprime "SPIN UP ORBITALS" / "SPIN DOWN ORBITALS"
    is_open_shell = any(
        "SPIN UP ORBITALS" in l or "SPIN DOWN ORBITALS" in l for l in block
    )

    homo_idx: Optional[int] = None
    lumo_idx: Optional[int] = None
    prev_occ_idx: Optional[int] = None

    for line in block:
        m = ORBITAL_ROW.match(line)
        if not m:
            # Si ya encontramos HOMO/LUMO y llegamos a una sección nueva, cortar
            if homo_idx is not None and lumo_idx is not None:
                break
            continue
        idx = int(m.group(1))
        occ = float(m.group(2))
        if occ > 0:
            prev_occ_idx = idx
        elif occ == 0 and prev_occ_idx is not None and lumo_idx is None:
            homo_idx = prev_occ_idx
            lumo_idx = idx
            if not is_open_shell:
                break

    return homo_idx, lumo_idx, is_open_shell


def parse_nto_states(out_text: str) -> dict[int, ExcitedState]:
    """Extrae, para cada estado excitado, TODOS los pares NTO reportados
    (con su coeficiente), ordenados por coeficiente descendente."""
    lines = out_text.splitlines()
    states: dict[int, ExcitedState] = {}

    current_state: Optional[int] = None
    for line in lines:
        header_m = NTO_STATE_HEADER.search(line)
        if header_m:
            current_state = int(header_m.group(1))
            states[current_state] = ExcitedState(state_idx=current_state)
            continue
        if current_state is None:
            continue
        row_m = NTO_ROW.match(line)
        if row_m:
            hole = int(row_m.group(1))
            particle = int(row_m.group(3))
            coef = float(row_m.group(5))
            states[current_state].ntos.append(
                NTOPair(hole=hole, particle=particle, coef=coef)
            )

    # Ordenar cada lista por coeficiente descendente (por si el output no
    # viniera ya ordenado)
    for state in states.values():
        state.ntos.sort(key=lambda p: p.coef, reverse=True)

    return states


def select_ntos_for_state(
    state: ExcitedState,
    threshold: float,
    min_coef: float,
    max_pairs: Optional[int],
) -> list[NTOPair]:
    """Selecciona los pares NTO a extraer: el dominante siempre, más los que
    cumplan coef >= threshold * coef_dominante Y coef >= min_coef."""
    if not state.ntos:
        return []
    dominant = state.ntos[0]
    selected = [dominant]
    for pair in state.ntos[1:]:
        if pair.coef >= min_coef and pair.coef >= threshold * dominant.coef:
            selected.append(pair)
    if max_pairs is not None:
        selected = selected[:max_pairs]
    return selected


def match_ntos_to_existing_files(
    conf_dir: Path, states: dict[int, ExcitedState]
) -> dict[int, ExcitedState]:
    """Cruza los estados parseados con los archivos tddft.sN.nto que
    realmente existen, para no inventar estados no calculados."""
    existing_state_ids = set()
    for f in conf_dir.glob("tddft.s*.nto"):
        m = re.match(r"tddft\.s(\d+)\.nto", f.name)
        if m:
            existing_state_ids.add(int(m.group(1)))

    filtered = {
        idx: state for idx, state in states.items() if idx in existing_state_ids
    }
    missing = set(states.keys()) - existing_state_ids
    if missing:
        logger.warning(
            "%s: estados en tddft.out sin archivo .nto correspondiente: %s",
            conf_dir,
            sorted(missing),
        )
    return filtered


# --------------------------------------------------------------------------- #
# Construcción y ejecución de orca_plot
# --------------------------------------------------------------------------- #


def build_orca_plot_stdin(orbital_indices: list[int], grid: int) -> str:
    """Arma el bloque de comandos interactivo de orca_plot para una lista de
    índices de orbitales, reutilizando la sesión."""
    commands = ["5", "7", "4", str(grid)]
    for idx in orbital_indices:
        commands += ["1", "1", "2", str(idx), "3", "0", "11"]
    commands.append("12")
    return "\n".join(commands) + "\n"


def run_orca_plot(
    orca_plot_bin: str,
    target_file: Path,
    orbital_indices: list[int],
    grid: int,
    timeout: int,
) -> tuple[bool, str]:
    """Ejecuta orca_plot en modo interactivo vía stdin sobre target_file.
    Devuelve (exito, mensaje)."""
    stdin_text = build_orca_plot_stdin(orbital_indices, grid)
    try:
        result = subprocess.run(
            [orca_plot_bin, target_file.name, "-i"],
            input=stdin_text,
            cwd=target_file.parent,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return False, f"timeout tras {timeout}s"
    except FileNotFoundError:
        return False, f"binario orca_plot no encontrado: {orca_plot_bin}"

    if result.returncode != 0:
        return False, f"returncode={result.returncode}: {result.stderr[-500:]}"
    if "PLOTTING FINISHED" not in result.stdout:
        return False, "no se encontró 'PLOTTING FINISHED' en stdout"
    return True, "ok"


def find_generated_cubes(target_file: Path, orbital_indices: list[int]) -> list[Path]:
    """orca_plot nombra los cubos como basename.mo{N}a.cube (o .b. si beta).
    Los busca en la carpeta del target_file."""
    found = []
    basename = target_file.stem  # ej. "tddft" o "tddft.s6"
    for idx in orbital_indices:
        for spin in ("a", "b"):
            candidate = target_file.parent / f"{basename}.mo{idx}{spin}.cube"
            if candidate.exists():
                found.append(candidate)
    return found


# --------------------------------------------------------------------------- #
# Procesamiento de un confórmero
# --------------------------------------------------------------------------- #


def process_conformer(
    conf_dir: Path,
    root: Path,
    args: argparse.Namespace,
) -> list[ReportRow]:
    rows: list[ReportRow] = []
    molecule, run_label, conf_id = parse_conformer_path(conf_dir, root)

    out_text = (conf_dir / "tddft.out").read_text(errors="replace")

    homo_idx, lumo_idx, is_open_shell = parse_ground_state_frontier_orbitals(out_text)
    if is_open_shell:
        logger.warning(
            "%s: sistema open-shell detectado, asumiendo alpha por defecto "
            "(revisar manualmente si esto no es correcto)",
            conf_dir,
        )

    cubes_dir = conf_dir / "cubes"
    if args.dry_run:
        cubes_dir_display = str(cubes_dir)
    else:
        cubes_dir.mkdir(exist_ok=True)
        cubes_dir_display = str(cubes_dir)

    # --- Estado fundamental: HOMO/LUMO ------------------------------------
    if homo_idx is None or lumo_idx is None:
        rows.append(
            ReportRow(
                molecule=molecule,
                run_label=run_label,
                conf_id=conf_id,
                kind="ground_state",
                homo_or_hole=homo_idx,
                lumo_or_particle=lumo_idx,
                coef=None,
                participation_ratio=None,
                fosc=None,
                cube_files="",
                status="error",
                error="no se pudo identificar HOMO/LUMO en ORBITAL ENERGIES",
            )
        )
    else:
        gbw_target = conf_dir / "tddft.gbw"
        cube_targets = [
            conf_dir / f"tddft.mo{homo_idx}a.cube",
            conf_dir / f"tddft.mo{lumo_idx}a.cube",
        ]
        already_exist = all(c.exists() for c in cube_targets)
        final_names = [cubes_dir / "HOMO_mo{}a.cube".format(homo_idx),
                        cubes_dir / "LUMO_mo{}a.cube".format(lumo_idx)]
        final_exist = all(c.exists() for c in final_names)

        if final_exist and not args.force:
            status, error = "skipped_exists", ""
            cube_list = final_names
        elif args.dry_run:
            status, error = "dry_run", ""
            cube_list = final_names
        else:
            ok, msg = run_orca_plot(
                args.orca_plot, gbw_target, [homo_idx, lumo_idx], args.grid, args.timeout
            )
            if not ok:
                status, error = "error", msg
                cube_list = []
            else:
                generated = find_generated_cubes(gbw_target, [homo_idx, lumo_idx])
                cube_list = []
                labels = {homo_idx: "HOMO", lumo_idx: "LUMO"}
                for g in generated:
                    m = re.search(r"\.mo(\d+)([ab])\.cube$", g.name)
                    if not m:
                        continue
                    idx = int(m.group(1))
                    label = labels.get(idx, f"mo{idx}")
                    dest = cubes_dir / f"{label}_mo{idx}{m.group(2)}.cube"
                    shutil.move(str(g), str(dest))
                    cube_list.append(dest)
                status, error = ("ok", "") if cube_list else ("error", "cubos no encontrados tras orca_plot")

        rows.append(
            ReportRow(
                molecule=molecule,
                run_label=run_label,
                conf_id=conf_id,
                kind="ground_state",
                homo_or_hole=homo_idx,
                lumo_or_particle=lumo_idx,
                coef=None,
                participation_ratio=None,
                fosc=None,
                cube_files=";".join(str(c.name) for c in cube_list),
                status=status,
                error=error,
            )
        )

    # --- Estados excitados: pares NTO --------------------------------------
    all_states = parse_nto_states(out_text)
    states = match_ntos_to_existing_files(conf_dir, all_states)

    fosc_by_state = parse_oscillator_strengths(out_text)
    dominant_state = select_dominant_state(fosc_by_state, args.fosc_min)

    if not fosc_by_state:
        logger.warning(
            "%s: no se encontró tabla de fosc, no se puede identificar el "
            "estado dominante (deltamax) — se omite la extracción de NTO "
            "para este confórmero", conf_dir,
        )
    elif dominant_state is None:
        brightest = max(fosc_by_state, key=fosc_by_state.get)
        logger.warning(
            "%s: ningún estado supera fosc_min=%.3f (máximo: estado %d, "
            "fosc=%.4f) — se omite la extracción de NTO para este confórmero",
            conf_dir, args.fosc_min, brightest, fosc_by_state[brightest],
        )

    for state_idx in sorted(states):
        fosc = fosc_by_state.get(state_idx)

        if state_idx != dominant_state:
            rows.append(
                ReportRow(
                    molecule=molecule,
                    run_label=run_label,
                    conf_id=conf_id,
                    kind=f"state{state_idx}",
                    homo_or_hole=None,
                    lumo_or_particle=None,
                    coef=None,
                    participation_ratio=states[state_idx].participation_ratio,
                    fosc=fosc,
                    cube_files="",
                    status="skipped_not_dominant",
                    error="",
                )
            )
            continue

        state = states[state_idx]
        selected = select_ntos_for_state(
            state, args.nto_threshold, args.nto_min_coef, args.nto_max_pairs
        )
        pr = state.participation_ratio

        if not selected:
            rows.append(
                ReportRow(
                    molecule=molecule,
                    run_label=run_label,
                    conf_id=conf_id,
                    kind=f"state{state_idx}",
                    homo_or_hole=None,
                    lumo_or_particle=None,
                    coef=None,
                    participation_ratio=pr,
                    fosc=fosc,
                    cube_files="",
                    status="error",
                    error="sin pares NTO parseados",
                )
            )
            continue

        nto_target = conf_dir / f"tddft.s{state_idx}.nto"
        orbital_indices = []
        for pair in selected:
            orbital_indices += [pair.hole, pair.particle]

        final_names = []
        for i, pair in enumerate(selected, start=1):
            suffix = "" if i == 1 else f"_alt{i}"
            final_names.append(
                cubes_dir / f"state{state_idx}_hole{suffix}_mo{pair.hole}a.cube"
            )
            final_names.append(
                cubes_dir / f"state{state_idx}_particle{suffix}_mo{pair.particle}a.cube"
            )
        final_exist = all(c.exists() for c in final_names)

        if final_exist and not args.force:
            status, error = "skipped_exists", ""
            cube_list = final_names
        elif args.dry_run:
            status, error = "dry_run", ""
            cube_list = final_names
        else:
            ok, msg = run_orca_plot(
                args.orca_plot, nto_target, orbital_indices, args.grid, args.timeout
            )
            if not ok:
                status, error = "error", msg
                cube_list = []
            else:
                generated = find_generated_cubes(nto_target, orbital_indices)
                cube_list = []
                idx_to_role = {}
                for i, pair in enumerate(selected, start=1):
                    suffix = "" if i == 1 else f"_alt{i}"
                    idx_to_role[pair.hole] = f"hole{suffix}"
                    idx_to_role[pair.particle] = f"particle{suffix}"
                for g in generated:
                    m = re.search(r"\.mo(\d+)([ab])\.cube$", g.name)
                    if not m:
                        continue
                    idx = int(m.group(1))
                    role = idx_to_role.get(idx, f"mo{idx}")
                    dest = cubes_dir / f"state{state_idx}_{role}_mo{idx}{m.group(2)}.cube"
                    shutil.move(str(g), str(dest))
                    cube_list.append(dest)
                status, error = ("ok", "") if cube_list else ("error", "cubos no encontrados tras orca_plot")

        rows.append(
            ReportRow(
                molecule=molecule,
                run_label=run_label,
                conf_id=conf_id,
                kind=f"state{state_idx}",
                homo_or_hole=selected[0].hole,
                lumo_or_particle=selected[0].particle,
                coef=selected[0].coef,
                participation_ratio=pr,
                fosc=fosc,
                cube_files=";".join(str(c.name) for c in cube_list),
                status=status,
                error=error,
            )
        )

        caracter = "mixto" if len(selected) > 1 else "puro"
        logger.info(
            "%s estado %d (dominante, fosc=%.4f): carácter %s (PR=%.2f), "
            "%d par(es) NTO extraídos",
            conf_dir, state_idx, fosc or -1, caracter, pr or -1, len(selected),
        )

    return rows


# --------------------------------------------------------------------------- #
# Reporte final
# --------------------------------------------------------------------------- #


def write_report(rows: list[ReportRow], out_path: Path) -> None:
    fieldnames = [
        "molecule", "run_label", "conf_id", "kind",
        "homo_or_hole", "lumo_or_particle", "coef", "participation_ratio",
        "fosc", "cube_files", "status", "error",
    ]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "molecule": row.molecule,
                    "run_label": row.run_label,
                    "conf_id": row.conf_id,
                    "kind": row.kind,
                    "homo_or_hole": row.homo_or_hole,
                    "lumo_or_particle": row.lumo_or_particle,
                    "coef": row.coef,
                    "participation_ratio": (
                        f"{row.participation_ratio:.3f}"
                        if row.participation_ratio is not None
                        else ""
                    ),
                    "fosc": f"{row.fosc:.6f}" if row.fosc is not None else "",
                    "cube_files": row.cube_files,
                    "status": row.status,
                    "error": row.error,
                }
            )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extrae cubos HOMO/LUMO y NTO hole/particle desde runs TDDFT."
    )
    parser.add_argument("--root", type=Path, default=Path("PAPER1"),
                         help="raíz de búsqueda (default: PAPER1)")
    parser.add_argument("--grid", type=int, default=80,
                         help="resolución de grilla para orca_plot (default: 80)")
    parser.add_argument(
        "--orca-plot",
        default="/home/kmilozeta/orca_6_1_1_avx2/orca_plot",
        help="ruta al binario orca_plot",
    )
    parser.add_argument(
        "--nto-threshold", type=float, default=0.15,
        help="fracción relativa al coeficiente dominante para incluir un par "
             "NTO adicional (default: 0.15, calibrado contra conf_5/B: separa "
             "limpiamente PR<=1.37 (puro) de PR>=1.39 (mixto) en ese conformero; "
             "ej. dominante=0.64 -> se incluyen pares con coef >= 0.096)",
    )
    parser.add_argument(
        "--nto-min-coef", type=float, default=0.05,
        help="piso absoluto de coeficiente para considerar un par NTO real, "
             "evita ruido numérico (default: 0.05)",
    )
    parser.add_argument(
        "--nto-max-pairs", type=int, default=2,
        help="tope duro de pares NTO por estado, aunque el umbral sugiera "
             "más (default: 2, evita nomenclatura de carpeta compleja)",
    )
    parser.add_argument(
        "--fosc-min", type=float, default=0.1,
        help="fuerza de oscilador mínima para considerar un estado como el "
             "'dominante' (deltamax) de un confórmero (default: 0.1). Solo se "
             "extraen NTO del único estado con mayor fosc, y solo si ese "
             "máximo supera este umbral; el resto queda como "
             "status=skipped_not_dominant en el reporte. Si el máximo no "
             "supera el umbral, o si no se encuentra la tabla ABSORPTION "
             "SPECTRUM en el .out, no se extrae ningún NTO para ese confórmero.",
    )
    parser.add_argument("--force", action="store_true",
                         help="regenerar aunque ya existan cubos")
    parser.add_argument("--dry-run", action="store_true",
                         help="solo lista qué haría, sin correr orca_plot")
    parser.add_argument("--timeout", type=int, default=600,
                         help="timeout por llamada a orca_plot en segundos (default: 600)")
    parser.add_argument(
        "--report", type=Path, default=None,
        help="ruta del CSV de salida (default: <root>/tddft_cubes_report.csv)",
    )
    parser.add_argument("-v", "--verbose", action="store_true",
                         help="logging DEBUG")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    root = args.root.resolve()
    if not root.exists():
        logger.error("La raíz de búsqueda no existe: %s", root)
        return 1

    report_path = args.report or (root / "tddft_cubes_report.csv")

    conformers = discover_conformers(root)
    if not conformers:
        logger.warning("No se encontraron confórmeros TDDFT bajo %s", root)
        return 0

    logger.info("Encontrados %d confórmeros TDDFT válidos", len(conformers))

    all_rows: list[ReportRow] = []
    n_errors = 0
    for conf_dir in conformers:
        try:
            rows = process_conformer(conf_dir, root, args)
            all_rows.extend(rows)
            n_errors += sum(1 for r in rows if r.status == "error")
        except Exception as exc:  # noqa: BLE001 - no queremos matar el batch
            logger.exception("Fallo procesando %s", conf_dir)
            all_rows.append(
                ReportRow(
                    molecule="", run_label="", conf_id=conf_dir.name,
                    kind="conformer", homo_or_hole=None, lumo_or_particle=None,
                    coef=None, participation_ratio=None, fosc=None, cube_files="",
                    status="error", error=str(exc),
                )
            )
            n_errors += 1

    if not args.dry_run:
        write_report(all_rows, report_path)
        logger.info("Reporte escrito en %s", report_path)
    else:
        logger.info("[dry-run] no se escribió reporte; %d filas planificadas", len(all_rows))

    logger.info(
        "Terminado: %d confórmeros, %d filas, %d errores",
        len(conformers), len(all_rows), n_errors,
    )
    return 1 if n_errors else 0


if __name__ == "__main__":
    sys.exit(main())
