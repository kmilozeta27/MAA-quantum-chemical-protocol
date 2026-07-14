#!/usr/bin/env python3
"""
recover_screening.py

Reconstruye manualmente la tabla de resultados de la etapa SCREENING de
CENSO a partir de los outputs de ORCA que ya quedaron en disco, sin
necesidad de re-correr los ~209 confórmeros que sí terminaron antes de
que el run muriera por el colgado/OOM de CONF346.

CENSO nunca escribió estos resultados a censo_ensemble_sorted.xyz ni al
JSON porque la etapa se aborta como bloque completo (ignore_failed=False)
en cuanto detecta los 2 fallos — pero los .out de ORCA de los jobs que sí
terminaron siguen físicamente en:

    01b_censo_prescreen/1_SCREENING/CONF<n>/gsolv/*.out

USO (correr DESDE el directorio 01b_censo_prescreen, o pasar --root):

    # 1) Modo inspección — mira la estructura real de archivos antes de
    #    confiar en el parseo automático. Corre esto primero.
    python3 recover_screening.py --root 1_SCREENING --inspect

    # 2) Modo extracción — una vez confirmada la estructura en el paso 1.
    python3 recover_screening.py --root 1_SCREENING --out screening_recovered.csv

Si el modo --inspect muestra 2 archivos .out por CONF y el script no los
puede clasificar automáticamente como gas/solv (te lo va a decir), ajusta
la función classify_files() con los nombres reales que use tu instalación
de CENSO.
"""
import argparse
import re
import sys
from pathlib import Path

HARTREE_TO_KCAL = 627.5094740631

FINAL_SP_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
TERM_OK_RE = re.compile(r"ORCA TERMINATED NORMALLY")


def find_conf_dirs(root: Path):
    dirs = [d for d in root.glob("CONF*/gsolv") if d.is_dir()]
    return sorted(dirs, key=lambda d: int(re.search(r"CONF(\d+)", d.parent.name).group(1)))


def parse_out_file(path: Path):
    text = path.read_text(errors="replace")
    energies = [float(e) for e in FINAL_SP_RE.findall(text)]
    terminated_ok = bool(TERM_OK_RE.search(text))
    return energies, terminated_ok


def classify_files(outs):
    """Intenta separar gas-phase vs CPCM/solvente por nombre de archivo.
    Ajusta esto si tu CENSO usa otra convención de nombres (revisa con --inspect)."""
    gas, solv, other = None, None, []
    for o in outs:
        name = o.name.lower()
        if "gas" in name or "vac" in name:
            gas = o
        elif "cpcm" in name or "solv" in name or "water" in name:
            solv = o
        else:
            other.append(o)
    return gas, solv, other


def inspect(root: Path, n: int):
    dirs = find_conf_dirs(root)
    print(f"Encontradas {len(dirs)} carpetas CONF*/gsolv bajo {root}\n")
    if not dirs:
        print("No se encontró nada. ¿--root apunta a la carpeta 1_SCREENING correcta?")
        return
    for d in dirs[:n]:
        outs = sorted(d.glob("*.out"))
        print(f"--- {d} ---")
        if not outs:
            print("  (sin archivos .out)")
        for o in outs:
            energies, ok = parse_out_file(o)
            print(f"  {o.name}: {len(energies)} x 'FINAL SINGLE POINT ENERGY' -> {energies}  "
                  f"| terminó OK: {ok}")
        print()
    print("Revisa el patrón de arriba:")
    print("- Si cada CONF tiene 1 solo .out con 1 energía -> esa ya es la SP")
    print("  con CPCM (solvatada); se usa directo como energía de ranking.")
    print("- Si hay 2 .out (gas + cpcm) -> el modo extracción intentará")
    print("  distinguirlos por nombre de archivo (classify_files()). Si no puede,")
    print("  te lo marcará como AMBIGUO en vez de adivinar.")


def extract(root: Path, out_csv: Path, failed_confs: set):
    dirs = find_conf_dirs(root)
    rows = []
    problems = []

    for d in dirs:
        conf = d.parent.name
        if conf in failed_confs:
            problems.append((conf, "conocido como fallido según censo.log — excluido"))
            continue

        outs = sorted(d.glob("*.out"))
        if not outs:
            problems.append((conf, "sin archivos .out en gsolv/"))
            continue

        if len(outs) == 1:
            energies, ok = parse_out_file(outs[0])
            if not energies:
                size = outs[0].stat().st_size
                last_line = outs[0].read_text(errors="replace").strip().splitlines()[-1:] or [""]
                problems.append((conf, f"{outs[0].name}: sin 'FINAL SINGLE POINT ENERGY' "
                                        f"({size} bytes, última línea: {last_line[0][:60]!r})"))
                continue
            if not ok:
                problems.append((conf, f"{outs[0].name}: no dice 'ORCA TERMINATED NORMALLY' "
                                        f"(tiene energía, pero revísalo a mano)"))
            e_final = energies[-1]  # última ocurrencia, por si hubo restart de SCF
            rows.append(dict(conf=conf, n_files=1, e_gas=None, e_solv=e_final,
                              gsolv_kcal=None, e_total=e_final, ok=ok))

        elif len(outs) == 2:
            gas, solv, other = classify_files(outs)
            if gas and solv:
                e_gas_list, ok_gas = parse_out_file(gas)
                e_solv_list, ok_solv = parse_out_file(solv)

                if e_gas_list and e_solv_list:
                    e_gas, e_solv = e_gas_list[-1], e_solv_list[-1]
                    gsolv_kcal = (e_solv - e_gas) * HARTREE_TO_KCAL
                    rows.append(dict(conf=conf, n_files=2, e_gas=e_gas, e_solv=e_solv,
                                      gsolv_kcal=gsolv_kcal, e_total=e_solv,
                                      ok=(ok_gas and ok_solv)))
                elif e_solv_list and not e_gas_list:
                    # gas job murió/truncado (típico del OOM), pero solv está OK ->
                    # usamos solo la energía solvatada como proxy de ranking.
                    e_solv = e_solv_list[-1]
                    rows.append(dict(conf=conf, n_files=1, e_gas=None, e_solv=e_solv,
                                      gsolv_kcal=None, e_total=e_solv, ok=ok_solv))
                    problems.append((conf, f"gas ({gas.name}) truncado/sin energía, "
                                            f"usada solo energía solv como fallback"))
                elif e_gas_list and not e_solv_list:
                    problems.append((conf, f"solv ({solv.name}) truncado/sin energía "
                                            f"y gas no sirve solo para ranking -> excluido"))
                else:
                    problems.append((conf, "2 .out pero ninguno tiene energía legible"))
            else:
                problems.append((conf, f"2 .out pero no los pude clasificar gas/solv: "
                                        f"{[o.name for o in outs]} -> ajusta classify_files()"))
        else:
            problems.append((conf, f"{len(outs)} archivos .out (caso inesperado): "
                                    f"{[o.name for o in outs]}"))

    rows.sort(key=lambda r: r["e_total"])

    with open(out_csv, "w") as f:
        f.write("CONF,n_out_files,E_gas_Eh,E_solv_Eh,dGsolv_kcal,E_total_Eh,terminated_ok\n")
        for r in rows:
            f.write(
                f"{r['conf']},{r['n_files']},"
                f"{r['e_gas'] if r['e_gas'] is not None else ''},"
                f"{r['e_solv']},"
                f"{r['gsolv_kcal'] if r['gsolv_kcal'] is not None else ''},"
                f"{r['e_total']},{r['ok']}\n"
            )

    print(f"OK: {len(rows)} confórmeros recuperados -> {out_csv}")
    if problems:
        print(f"\n{len(problems)} confórmeros con problemas (NO incluidos en el CSV):")
        for conf, msg in problems:
            print(f"  {conf}: {msg}")
    total_expected = len(dirs)
    print(f"\nResumen: {len(rows)}/{total_expected} recuperados, "
          f"{len(problems)}/{total_expected} con problemas o excluidos.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path,
                     help="Carpeta 1_SCREENING (contiene CONF<n>/gsolv/)")
    ap.add_argument("--inspect", action="store_true",
                     help="Solo inspecciona la estructura de archivos, no extrae nada")
    ap.add_argument("--inspect-n", type=int, default=5,
                     help="Cuántos CONF mostrar en modo --inspect (default: 5)")
    ap.add_argument("--out", type=Path, default=Path("screening_recovered.csv"),
                     help="CSV de salida (default: screening_recovered.csv)")
    ap.add_argument("--failed", nargs="*", default=["CONF341", "CONF346"],
                     help="Confórmeros conocidos como fallidos, se excluyen "
                          "(default: CONF341 CONF346)")
    args = ap.parse_args()

    if not args.root.exists():
        sys.exit(f"No existe la carpeta: {args.root}")

    if args.inspect:
        inspect(args.root, args.inspect_n)
    else:
        extract(args.root, args.out, set(args.failed))
