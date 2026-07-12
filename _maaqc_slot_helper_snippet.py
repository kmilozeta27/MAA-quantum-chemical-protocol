# ─────────────────────────────────────────────────────────────────────────
# MAAQC-PATCH: coordinación de slots de CPU entre SP/OPT paralelos de CENSO
# ─────────────────────────────────────────────────────────────────────────
# Ver README.md junto a este parche para el contexto completo. Resumen:
# CENSO congela ENVIRON una sola vez (params.py) y todas las llamadas
# paralelas a ORCA lo heredan idéntico -- sin esto, 2+ SP simultáneos
# reciben el mismo MAAQC_PE_LIST y colisionan en los mismos cores físicos.
#
# pipeline.py fija MAAQC_CENSO_SLOTS antes de lanzar CENSO, con el pool
# particionado en N grupos de PUs separados por '|', ej:
#   MAAQC_CENSO_SLOTS="0,1,4,5|2,3,6,7"   (2 slots de 4 PUs cada uno)
#
# Este helper deriva, de forma determinística y sin coordinación en tiempo
# real, a qué slot pertenece cada conformero (hash estable de
# job.conf.name módulo N), y construye el MAAQC_PE_LIST correspondiente
# (traducido de PU a índice de Core, mismo mecanismo que build_omp_env()
# en pipeline.py -- PRTE interpreta pe-list como índice de CORE, no PU).
import hashlib

# Misma tabla que CORE_TO_PUS/PU_TO_CORE en pipeline.py (ThinkPad T14,
# Intel Ultra 7 165U, ver historial de CPU pinning). Duplicada aquí porque
# este archivo vive en el paquete CENSO, proceso independiente de
# pipeline.py -- no se puede importar directamente.
_MAAQC_CORE_TO_PUS: dict[int, list[int]] = {
    0: [0, 1], 1: [2, 3],                        # P-cores, SMT
    2: [4], 3: [5], 4: [6], 5: [7],
    6: [8], 7: [9], 8: [10], 9: [11],            # E-cores
    10: [12], 11: [13],                          # LP E-cores
}
_MAAQC_PU_TO_CORE: dict[int, int] = {
    pu: core for core, pus in _MAAQC_CORE_TO_PUS.items() for pu in pus
}


def _maaqc_build_slot_env(conf_name: str) -> dict[str, str] | None:
    """
    Devuelve un env (copia de ENVIRON + MAAQC_PE_LIST del slot que le
    corresponde a `conf_name`), o None si MAAQC_CENSO_SLOTS no está
    fijada (fail-open: sin esta variable, el comportamiento es idéntico
    al de antes del parche -- _make_call cae a ENVIRON global tal cual).

    El slot se deriva con un hash estable (md5 del nombre, no hash()
    built-in de Python -- ese varía entre procesos por PYTHONHASHSEED
    aleatorio, lo cual rompería la reproducibilidad entre corridas).
    """
    slots_raw = os.environ.get("MAAQC_CENSO_SLOTS")
    if not slots_raw:
        return None

    groups = [g for g in slots_raw.split("|") if g]
    if not groups:
        return None

    n_slots = len(groups)
    digest = hashlib.md5(conf_name.encode("utf-8")).hexdigest()
    slot_idx = int(digest, 16) % n_slots

    try:
        pu_list = [int(x) for x in groups[slot_idx].split(",") if x]
    except ValueError:
        # MAAQC_CENSO_SLOTS mal formado -- fail-open, no romper el cálculo
        # por un problema de pinning.
        return None

    if not pu_list:
        return None

    pe_list = [_MAAQC_PU_TO_CORE.get(pu) for pu in pu_list]
    if any(c is None for c in pe_list):
        # PU fuera de la tabla conocida -- fail-open.
        return None

    env = ENVIRON.copy()
    env["MAAQC_PE_LIST"] = ",".join(str(c) for c in pe_list)
    return env
