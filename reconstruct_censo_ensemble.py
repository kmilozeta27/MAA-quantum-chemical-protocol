#!/usr/bin/env python3
"""
reconstruct_censo_ensemble.py

Reconstruye la salida de la etapa SCREENING de CENSO a partir de
screening_recovered.csv (generado por recover_screening.py), sin
1_SCREENING.out real, porque CENSO abortó la etapa completa al
encontrar 2 fallos (CONF341, CONF346) con ignore_failed=False.

Genera DOS archivos, porque pipeline.py necesita AMBOS para no
re-lanzar CENSO desde cero al resumir:

  1) censo_ensemble_sorted.xyz
     Las geometrías filtradas, en orden de estabilidad. Se construye
     reutilizando extract_conformers_by_id() de tu propio pipeline.py
     (import directo del módulo), así el formato de salida es
     BYTE-IDÉNTICO al que el pipeline generaría por su camino normal.

  2) 1_SCREENING.out
     Archivo sintético con la tabla de resultados en el mismo formato
     de columnas que parsean parse_censo_screening_out() y
     parse_boltzmann_weights(): CONF# E(DFT) ΔGsolv GmRRHO Gtot ΔGtot
     Boltzmann%. Sin esto, stage_censo_prescreen() NO reconoce el
     screening como completo (revisa 1_SCREENING.out explícitamente)
     y select_nconf_boltzmann() no tiene de dónde leer los pesos.

APROXIMACIÓN A TENER EN CUENTA:
    La columna GmRRHO (corrección de energía libre vibracional) se
    fija en 0.00 como placeholder — no la recuperamos porque el
    screening real de CENSO no vuelve a computar RRHO (lo reutiliza
    de prescreening) y esos archivos NumFreq de screening nunca se
    generaron. Gtot = E_solv (SP con CPCM) tal como quedó en el CSV,
    que es la misma cantidad que CENSO habría usado como Gtot salvo
    por ese término RRHO. El ORDEN de estabilidad y los pesos
    Boltzmann relativos deberían ser muy similares a los reales, pero
    no son bit-idénticos a los que CENSO habría producido. Si tu
    threshold de Boltzmann (--boltzmann-threshold) es sensible a
    diferencias de <0.5 kcal/mol, ten esto presente al revisar qué
    confórmeros terminan pasando a OPT.

USO (correr desde el directorio PAPER1, donde vive pipeline.py):

    python3 reconstruct_censo_ensemble.py \\
        --csv       MGLYVAL/runs/R/01b_censo_prescreen/1_SCREENING/screening_recovered.csv \\
        --source-xyz MGLYVAL/runs/R/01_crest/crest_conformers.xyz \\
        --censo-dir MGLYVAL/runs/R/01b_censo_prescreen \\
        --temperature 298.15

Esto escribe:
    MGLYVAL/runs/R/01b_censo_prescreen/censo_ensemble_sorted.xyz
    MGLYVAL/runs/R/01b_censo_prescreen/1_SCREENING.out
"""
import argparse
import csv
import importlib.util
import logging
import math
import re
import sys
from pathlib import Path

R_KCAL = 0.0019872041  # kcal/(mol*K)
HARTREE_TO_KCAL = 627.5094740631


def load_pipeline_module(pipeline_path: Path):
    spec = importlib.util.spec_from_file_location("pipeline_mod", pipeline_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ensure_conf_labels(source_xyz: Path, log: logging.Logger) -> Path:
    """
    crest_conformers.xyz "plano" (el formato normal de CREST) no trae
    ningún token CONF<n> en el comentario -- solo la energía
    (ej. "        -78.73620006"). extract_conformers_by_id() de
    pipeline.py busca literalmente ese token y, si no lo encuentra,
    no escribe nada (esto es justo lo que pasó).

    CREST siempre ordena crest_conformers.xyz por energía ascendente, y
    tanto CREST como CENSO usan esa POSICIÓN como numeración implícita:
    frame 1 (1-indexed) = CONF1, frame 2 = CONF2, etc. Esta es la misma
    convención de fallback que tu propio pipeline.py ya usa en
    stage_opt(): `_xyz_id_map.get(i, f"CONF{i+1}")`.

    Si el archivo YA trae tokens CONF# reales en los comentarios (p.ej.
    porque le pasaste un ensemble ya filtrado/etiquetado), no se toca
    nada y se devuelve el path original intacto.

    Devuelve el path a usar (el original, o una copia temporal con los
    labels inyectados por posición).
    """
    raw = source_xyz.read_text(errors="replace").splitlines()

    if any(re.search(r"\bCONF\d+", line) for line in raw):
        log.info(f"[RECON] {source_xyz.name} ya trae tokens CONF# -- no se modifica")
        return source_xyz

    log.info(f"[RECON] {source_xyz.name} no trae tokens CONF# (formato CREST "
             f"plano) -- inyectando CONF<n> por posición (orden = energía "
             f"ascendente, convención estándar CREST/CENSO)")

    out_lines = []
    i = 0
    frame_n = 0
    n = len(raw)
    while i < n:
        stripped = raw[i].strip()
        if stripped.isdigit():
            nat = int(stripped)
            frame_n += 1
            out_lines.append(raw[i])
            comment = raw[i + 1] if i + 1 < n else ""
            out_lines.append(f"CONF{frame_n} {comment.strip()}")
            i += 2
            for _ in range(nat):
                if i < n:
                    out_lines.append(raw[i])
                    i += 1
        else:
            out_lines.append(raw[i])
            i += 1

    labeled_path = source_xyz.with_name(source_xyz.stem + "_labeled.xyz")
    labeled_path.write_text("\n".join(out_lines) + "\n")
    log.info(f"[RECON] {frame_n} frames etiquetados -> {labeled_path.name}")
    return labeled_path


def read_recovered_csv(csv_path: Path):
    """
    Lee screening_recovered.csv (ya viene ordenado ascendente por
    E_total_Eh, es decir, por estabilidad, igual que
    parse_censo_screening_out() ordena por Gtot).
    """
    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(dict(
                conf=row["CONF"],
                e_gas=float(row["E_gas_Eh"]) if row["E_gas_Eh"] else None,
                e_solv=float(row["E_solv_Eh"]),
                gsolv_kcal=float(row["dGsolv_kcal"]) if row["dGsolv_kcal"] else None,
                e_total=float(row["E_total_Eh"]),
            ))
    if not rows:
        sys.exit(f"No se encontraron filas en {csv_path}")
    return rows


def write_synthetic_screening_out(rows, out_path: Path, temperature: float):
    """
    Escribe 1_SCREENING.out sintético, con las mismas 7 columnas que
    parse_censo_screening_out()/parse_boltzmann_weights() esperan.
    """
    gtot_min = min(r["e_total"] for r in rows)
    rt = R_KCAL * temperature

    # Boltzmann weight (%) usando Gtot=E_total como proxy (ver docstring:
    # falta el término GmRRHO real).
    dgtot_list = [(r["e_total"] - gtot_min) * HARTREE_TO_KCAL for r in rows]
    boltz_factors = [math.exp(-dg / rt) for dg in dgtot_list]
    z = sum(boltz_factors)
    weights = [100.0 * bf / z for bf in boltz_factors]

    lines = []
    lines.append("PARTIAL RECOVERY — reconstruido manualmente desde "
                  "screening_recovered.csv (ver reconstruct_censo_ensemble.py).")
    lines.append("GmRRHO fijado en 0.00 (placeholder, no recuperado). "
                  "Ver docstring del script para el detalle de la aproximación.")
    lines.append("")
    lines.append(f"{'CONF#':<10}{'E(DFT)':>14}{'dGsolv':>12}{'GmRRHO':>10}"
                  f"{'Gtot':>14}{'dGtot':>10}{'Boltzmann%':>12}")
    for r, dgtot, w in zip(rows, dgtot_list, weights):
        e_dft = r["e_gas"] if r["e_gas"] is not None else r["e_total"]
        gsolv = r["gsolv_kcal"] if r["gsolv_kcal"] is not None else 0.0
        lines.append(
            f"{r['conf']:<10}{e_dft:>14.6f}{gsolv:>12.2f}{0.00:>10.2f}"
            f"{r['e_total']:>14.6f}{dgtot:>10.2f}{w:>12.2f}"
        )
    out_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, type=Path,
                     help="screening_recovered.csv (de recover_screening.py)")
    ap.add_argument("--source-xyz", required=True, type=Path,
                     help="multi-XYZ de CREST usado como -i de CENSO "
                          "(crest_conformers.xyz, comentarios 'CONF<n> ...')")
    ap.add_argument("--censo-dir", required=True, type=Path,
                     help="01b_censo_prescreen/ (o el subdir T<temp>K/ en modo "
                          "--multitemp) — ahí se escriben ambos archivos de salida")
    ap.add_argument("--temperature", type=float, default=298.15,
                     help="Temperatura usada en el run de CENSO, para el peso "
                          "Boltzmann sintético (default: 298.15, revisa tu "
                          "censo.log -T)")
    ap.add_argument("--pipeline", type=Path, default=Path("pipeline.py"),
                     help="Ruta a pipeline.py (default: ./pipeline.py)")
    args = ap.parse_args()

    for p, label in [(args.csv, "--csv"), (args.source_xyz, "--source-xyz"),
                      (args.pipeline, "--pipeline")]:
        if not p.exists():
            sys.exit(f"No existe {label}: {p}")
    if not args.censo_dir.exists():
        sys.exit(f"No existe --censo-dir: {args.censo_dir}")

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    log = logging.getLogger("reconstruct")

    pipeline = load_pipeline_module(args.pipeline)

    rows = read_recovered_csv(args.csv)
    conf_ids = [r["conf"] for r in rows]
    log.info(f"[RECON] {len(conf_ids)} confórmeros leídos de {args.csv.name}, "
              f"orden: {conf_ids[0]} (más estable) -> {conf_ids[-1]} (menos estable)")

    # 1) censo_ensemble_sorted.xyz — vía el propio pipeline.py
    labeled_source = ensure_conf_labels(args.source_xyz, log)
    output_xyz = args.censo_dir / "censo_ensemble_sorted.xyz"
    ok = pipeline.extract_conformers_by_id(
        source_xyz=labeled_source,
        conf_ids=conf_ids,
        output_xyz=output_xyz,
        log=log,
    )
    if not ok:
        sys.exit("[RECON] extract_conformers_by_id() no escribio nada -- revisa "
                  "que --source-xyz tenga comentarios 'CONF<n>' y que esos IDs "
                  "existan ahi.")

    from ase.io import read
    n_written = len(read(str(output_xyz), index=":"))
    log.info(f"[RECON] {output_xyz} -> {n_written}/{len(conf_ids)} confórmeros escritos")
    if n_written != len(conf_ids):
        log.warning(f"[RECON]  ATENCION: faltan {len(conf_ids) - n_written} "
                     f"-- revisa el warning de extract_conformers_by_id arriba "
                     f"(probablemente CONF#s del CSV no encontrados en --source-xyz)")

    # 2) 1_SCREENING.out sintético
    screening_out = args.censo_dir / "1_SCREENING.out"
    write_synthetic_screening_out(rows, screening_out, args.temperature)
    log.info(f"[RECON] {screening_out} escrito ({len(rows)} filas, "
              f"T={args.temperature} K)")

    # Verificación: que el propio pipeline.py pueda re-parsear lo que acabamos
    # de escribir, igual que lo haria al resumir.
    reparsed = pipeline.parse_censo_screening_out(screening_out)
    if reparsed == conf_ids:
        log.info("[RECON] OK: parse_censo_screening_out() reproduce el mismo "
                  "orden de CONF# -- listo para que el pipeline resuma en OPT.")
    else:
        log.warning("[RECON] parse_censo_screening_out() devolvió un orden "
                     "distinto al esperado -- revisa manualmente antes de resumir.")
