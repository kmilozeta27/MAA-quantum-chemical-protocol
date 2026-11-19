# Sesión 2026-11-19: CENSO (protocolo, integración en maaqc y corrida de C′)

## 1. Diagnóstico del uso de CENSO

Revisé el código de CENSO 3.0.9.dev1 (`~/CENSO`), `~/.censo2rc` y las corridas previas de DG y AP.

**Problemas en `~/.censo2rc`:**
- `gsolv_included = False` con ORCA corre dos SP por conformero (gas y solvente) para obtener el mismo Gtot (E_gas + (E_solv − E_gas) = E_solv). Con `True` el screening cuesta la mitad, y el screening fue ~80 % del tiempo en DG.
- `sm = alpb` en `[prescreening]` no tiene efecto: CENSO ignora sin avisar las claves desconocidas. El Gsolv del prescreening sale de `[general] sm_rrho`, que por defecto es GBSA.
- Falta `[general]`, así que `ignore_failed` queda en False.
- Sin `[optimization]` ni `[refinement]`, el programa por defecto es TURBOMOLE y `-O`/`-R` fallan en la validación.

**Lo que mostraron los datos de DG:**
- El filtro descartó poco: 88 → 82 → 71 conformeros en 4DG_A.
- Sobre geometrías GFN2 el paisaje es plano: 61 de 70 conformeros quedan a menos de 2 kcal/mol.
- Hubo conformeros a 3.5–3.9 kcal/mol en el prescreening que terminaron a menos de 2.5 kcal/mol en el screening, así que el umbral de 4.0 no debe bajarse.
- maaqc elegía los 3 conformeros a optimizar con ωB97X-D3 a partir de esos SP sobre geometrías GFN2, una base poco fiable.

## 2. Protocolo adoptado (`censo2rc_refine`, en la raíz del repo)

| Etapa | Nivel | Umbral |
|---|---|---|
| Prescreening | B97-D3/def2-SV(P) en gas + ALPB(GFN2) | 4.0 kcal/mol |
| Screening | r²SCAN-3c CPCM, `gsolv_included = True` | 3.5 kcal/mol + ajuste automático |
| Optimización | r²SCAN-3c CPCM, ANCOPT con macrociclos | 3.0 kcal/mol, gradthr 0.01 |
| Refinamiento | ωB97X-D3/def2-TZVP SMD (mismo funcional y solvente que el SP final de maaqc) | 95 % de población Boltzmann |

**Recursos** (Core Ultra 7 165U, 12 núcleos físicos, 30 GB): `--maxcores 12 --omp-min 4`, es decir, 3 trabajos de 4 núcleos. Los microestados se corren uno tras otro, no en paralelo.

## 3. Cambios en `maaqc.py`

- **`--censo-refine`:** CENSO corre P+S+O+R usando `censo2rc_refine`.
  - Si ya existe un 1_SCREENING o 2_OPTIMIZATION terminado, retoma desde ahí (`-O -R` o `-R` sobre el `.xyz` de esa parte, con salida en `censo_refine.log`).
  - Si falta 3_REFINEMENT, la etapa falla en vez de volver al screening sin avisar.
  - Siempre pasa `--ignore-failed`, porque la opción de línea de comandos de CENSO sobrescribe el valor del rc.
  - El límite de tiempo es 120 h.
- **`--censo-rc RCFILE`:** se pasa a CENSO como `--inprc`.
- **`find_censo_result()`:** toma el resultado más avanzado (3_REFINEMENT, luego 2_OPTIMIZATION, luego 1_SCREENING). Para O y R usa solo las geometrías r²SCAN-3c de CENSO, nunca las de CREST. `_find_screening_out` y `--boltzmann-threshold` leen esos mismos pesos.
- **Arreglo de Open MPI en `_run_censo_locked`:**
  - *Síntoma:* CPU casi en cero.
  - *Causa 1:* CENSO lanza ORCA con `stderr=PIPE` y no lee la tubería hasta el final. Open MPI 5.0.10 escribe ~260 KB de avisos UCX/TCP por SP, la tubería (64 KB) se llena y `prterun` queda bloqueado para siempre.
  - *Causa 2:* cada `mpirun` ata sus procesos a los núcleos 0–3, así que los 3 trabajos se apilan sobre los mismos 4 núcleos.
  - *Solución:* `OMPI_MCA_pml=ob1`, `OMPI_MCA_btl=self,sm`, `OMPI_MCA_osc=^ucx`, `PRTE_MCA_hwloc_default_binding_policy=none`.
  - *Resultado:* un SP de Cp pasa de 85 s a 15 s con 3 en paralelo, sin cuelgues. **Si se llama `censo` a mano, hay que exportar estas variables antes.**
- **Tests:** `test_maaqc.py` suma 3 tests nuevos; pytest = 68 passed.

## 4. Corrida de C′ (Cp)

```bash
python maaqc.py run --censo-standalone MGLYVAL/microstates_new/merged/Cp/crest_ensemble.xyz --out MGLYVAL/microstates_new/censo/Cp --charge -2 --temperature 298.15 --censo-refine
```

- **Entrada:** 197 conformeros, carga −2 (según `.CHRG`).
- **Prescreening:** terminado en 37 min (2227 s); 197 → 40 conformeros.
- **Estado a las 09:49:** screening en curso; optimización y refinamiento pendientes.

## 5. Pendiente

- Terminar C′ y revisar `3_REFINEMENT.out`: cuántos conformeros quedan y si hay duplicados tras optimizar, ya que CENSO 3 no los elimina.
- Correr los otros 6 microestados (D′, R′, 010_N8H/N9H, 100_N8H/N9H) cuando estén fusionados; Dp tiene 336 conformeros.
- Decidir si `~/.censo2rc` se reemplaza por `censo2rc_refine`.
- Commit de los cambios (no se ha hecho).
