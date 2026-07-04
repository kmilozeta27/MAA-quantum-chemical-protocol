#!/usr/bin/env python3
"""
Recovery script para CENSO 3.0.x — reconstruye 0_PRESCREENING.json
a partir de los .out ya calculados en disco, sin volver a correr nada.

Por qué esto es necesario:
CENSO (en esta arquitectura basada en Dask) solo escribe el JSON consolidado
de resultados (0_PRESCREENING.json) AL FINAL de la parte, cuando todos los
confórmeros terminaron en la misma ejecución del proceso. Si el proceso se
mata antes (KeyboardInterrupt, crash, etc.), ese JSON nunca se genera, y
CENSO no tiene ningún mecanismo de "saltar cálculo si el .out ya existe":
siempre vuelve a escribir el .inp y a correr desde cero. El único punto de
entrada para reutilizar resultados antiguos es el flag `--reload archivo.json`,
que exige que TODOS los conformeros de la ensemble actual estén presentes
en ese JSON.

Este script reconstruye ese JSON a mano, parseando exactamente lo mismo
que parsea CENSO internamente:
  - sp/sp.out              -> "FINAL SINGLE POINT ENERGY" (último token) -> energy
  - xtb_gsolv/gas.out      -> línea "| TOTAL ENERGY", token[3]            -> energy_gas
  - xtb_gsolv/solv.out     -> línea "| TOTAL ENERGY", token[3]            -> energy_solv
  - gsolv = energy_solv - energy_gas
  - grrho = 0.0  (prescreening no calcula grrho)
  - gtot  = energy + gsolv + grrho

Uso:
    python recover_prescreening.py /home/kmilozeta/MGLYVAL/runs/SB/01b_censo_prescreen

Esto generará:
    .../01b_censo_prescreen/0_PRESCREENING_recovered.json

que luego se usa así:
    cd /home/kmilozeta/MGLYVAL/runs/SB/01b_censo_prescreen
    censo --prescreening -i crest_conformers.xyz -c -1 -u 0 --solvent water \
          --maxcores 12 -T 298.15 \
          --reload 0_PRESCREENING_recovered.json

El script reporta al final qué confórmeros NO pudo recuperar (faltantes o
con error de parseo), para que sepas si falta relanzar algo puntual antes
de usar --reload (recuerda: read_output exige el set COMPLETO).
"""

import json
import re
import sys
from pathlib import Path

PART_DIR_NAME = "0_PRESCREENING"
OUT_JSON_NAME = "0_PRESCREENING_recovered.json"


def parse_orca_sp_energy(outfile: Path) -> float | None:
    """Replica el parseo de OrcaProc._sp: busca FINAL SINGLE POINT ENERGY."""
    if not outfile.exists():
        return None
    try:
        text = outfile.read_text(errors="ignore")
    except OSError:
        return None

    if "ORCA TERMINATED NORMALLY" not in text:
        return None

    for line in text.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                return float(line.split()[-1])
            except (ValueError, IndexError):
                return None
    return None


def parse_xtb_total_energy(outfile: Path) -> float | None:
    """Replica el parseo de XtbProc._sp: busca '| TOTAL ENERGY', token[3]."""
    if not outfile.exists():
        return None
    try:
        text = outfile.read_text(errors="ignore")
    except OSError:
        return None

    for line in text.splitlines():
        if "| TOTAL ENERGY" in line:
            try:
                return float(line.split()[3])
            except (ValueError, IndexError):
                return None
    return None


def natural_conf_sort_key(name: str):
    """Ordena CONF1, CONF2, ..., CONF579 numéricamente (no lexicográficamente)."""
    m = re.search(r"\d+", name)
    return int(m.group()) if m else 0


def recover(prescreening_dir: Path) -> None:
    if not prescreening_dir.is_dir():
        sys.exit(f"ERROR: no existe el directorio {prescreening_dir}")

    conf_dirs = sorted(
        [d for d in prescreening_dir.iterdir() if d.is_dir() and d.name.startswith("CONF")],
        key=lambda d: natural_conf_sort_key(d.name),
    )

    if not conf_dirs:
        sys.exit(f"ERROR: no se encontraron carpetas CONF* en {prescreening_dir}")

    print(f"Encontradas {len(conf_dirs)} carpetas de confórmeros en {prescreening_dir}\n")

    data: dict[str, dict[str, float]] = {}
    missing: list[str] = []
    errors: list[tuple[str, str]] = []

    for conf_dir in conf_dirs:
        name = conf_dir.name

        sp_out = conf_dir / "sp" / "sp.out"
        gas_out = conf_dir / "xtb_gsolv" / "gas.out"
        solv_out = conf_dir / "xtb_gsolv" / "solv.out"

        energy = parse_orca_sp_energy(sp_out)
        if energy is None:
            missing.append(name)
            errors.append((name, f"sp.out ausente, sin TERMINATED NORMALLY, o sin FINAL SINGLE POINT ENERGY ({sp_out})"))
            continue

        energy_gas = parse_xtb_total_energy(gas_out)
        energy_solv = parse_xtb_total_energy(solv_out)

        if energy_gas is None or energy_solv is None:
            missing.append(name)
            faltante = "gas.out" if energy_gas is None else "solv.out"
            errors.append((name, f"no se pudo parsear {faltante} en xtb_gsolv/"))
            continue

        gsolv = energy_solv - energy_gas
        grrho = 0.0
        gtot = energy + gsolv + grrho

        data[name] = {
            "energy": energy,
            "gsolv": gsolv,
            "grrho": grrho,
            "gtot": gtot,
        }

    print(f"Recuperados correctamente: {len(data)} / {len(conf_dirs)}\n")

    if errors:
        print("--- Confórmeros NO recuperados ---")
        for name, reason in errors:
            print(f"  {name}: {reason}")
        print()

    if not data:
        sys.exit("No se pudo recuperar ningún confórmero. Revisa las rutas.")

    # Estructura exacta que espera EnsembleData.read_output:
    # un dict con clave "data" (formato "nested" que el parser acepta),
    # más metadata informativa (no es requerida por el parser pero no molesta).
    dump = {
        "part_name": "prescreening",
        "data": data,
    }

    out_path = prescreening_dir / OUT_JSON_NAME
    out_path.write_text(json.dumps(dump, indent=4))

    print(f"JSON de recovery escrito en: {out_path}")
    print(f"Total confórmeros incluidos: {len(data)}")

    if missing:
        print(
            f"\n⚠️  ATENCIÓN: {len(missing)} confórmero(s) quedaron fuera del JSON: "
            f"{', '.join(missing)}"
        )
        print(
            "   read_output() de CENSO exige que TODOS los confórmeros de la ensemble\n"
            "   actual (los del crest_conformers.xyz que le pases con -i) estén en el\n"
            "   JSON, o lanzará RuntimeError. Antes de usar --reload, debes:\n"
            "     (a) completar manualmente esos cálculos y volver a correr este script, o\n"
            "     (b) quitar esos confórmeros del crest_conformers.xyz de entrada, o\n"
            "     (c) generar para ellos un .json separado de un run que sí los complete\n"
            "         y pasar ambos archivos a --reload (se cargan en orden)."
        )
    else:
        print("\n✅ Todos los confórmeros recuperados. Listo para usar --reload.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Uso: python {Path(__file__).name} <ruta_a_01b_censo_prescreen>")

    recover(Path(sys.argv[1]).resolve())
