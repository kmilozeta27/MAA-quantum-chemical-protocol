# Archivo en disco externo TOSHIBA — inventario y plan (2026-11-17)

Disco: `/run/media/kmilozeta/TOSHIBA`, Toshiba MQ04UBF100 (HDD USB de 2.5"), NTFS vía ntfs-3g, 932 GB (160 GB usados, 772 GB libres). Inventario hecho en solo lectura; **no se movió ni borró nada del disco**. Listado completo (≈134 000 archivos) en el scratchpad de la sesión.

## Qué hay (proyecto, 93 GB en total)

| Carpeta | Tamaño | Contenido |
|---|---|---|
| `PAPER1/` | 60.3 GB, 90 634 archivos | Árbol "canónico": `MGLYVAL/runs/{A,B,C,D,R}` (33.8 GB), `SHINORINE/runs/{SA,SB,SD}` (14 GB), `CALIBRATION/{ASPARTIC,GLUTAMIC,GLYCINE,VALINE}` (9.5 GB), `PALYTHINE/runs/PR`, `QCG/`, scripts, `calib_v2.csv`, `.git` (1.6 GB) |
| `PAPER1.2/` | 32.9 GB, 43 320 archivos | Instantánea más nueva (hasta 10-oct): `ASPARTIC/` (sin el prefijo `CALIBRATION/`), `GADUSOL/` (3.6 GB, 4-deoxigadusol con nombres `GAD_dia{1,2}_{A,B,T}`), `MGLYVAL/runs/{B,R}` parcial, `calibration_v3.json`, `calib_v3.csv`, `maaqc.py` más nuevo, `.git` (0.6 GB) |
| `04c_sp_smd`, `(2)`, `(3)`, `(4)`, `QCG-2/` | casi vacíos | Restos sueltos de copias: solo un `sp_boltzmann_weights.json` cada uno; sus subcarpetas `conf_N` dan **error de E/S** |
| Otros (Música, GarageBand, Logic, respaldos del iPhone, etc.) | — | **Ajenos al proyecto; no se tocan** |

Por extensión (todo el proyecto): `.gbw` 17.0 GB, `.densities` 16.6 GB, `.nto` 7.0 GB, `.tmp` 5.6 GB, `.out` 5.4 GB, `.ges` 3.9 GB.

## Hallazgos

1. **Datos que existen solo en el portátil**: `AP/` (calibrador 4-aminopiridina), `DG/` (4-deoxigadusol, nombres `4DG_*`), y el `MGLYVAL/runs/B` local con `07_adc2` y `gold_01` (gold standard). El disco no tiene `AP` ni `DG`. **Ya copiados** (ver abajo).
2. **`PAPER1.2` es casi todo redundante con `PAPER1`**: 28 880 archivos (25 GB) tienen misma ruta (tras mapear `ASPARTIC` ↔ `CALIBRATION/ASPARTIC`), mismo tamaño y misma fecha. Falta confirmar con hash antes de tratarlos como duplicados.
3. **Lo único que `PAPER1.2` aporta y `PAPER1` no tiene** (7 979 archivos, 4.64 GB): `GADUSOL/` completo (GAD_dia1_A/B/T, GAD_dia2_A/B/T), 105 archivos más en `MGLYVAL/runs/B`, `calibration_v3.json`, `calib_v3.csv`, `calibration.csv/json` y un `maaqc.py` y `pka_mglyval_raw.log` con tamaño distinto.
4. **No hay duplicados grandes dentro de `PAPER1`** (0 grupos de archivos >50 MB con mismo nombre y tamaño). Los duplicados reales son entre las dos carpetas y entre los `.git`.
5. **Salud del disco**: error de E/S al listar `04c_sp_smd/conf_2` y otros `conf_N` sueltos, y un nombre de archivo con salto de línea. Es un síntoma de entradas NTFS corruptas (probablemente copias interrumpidas) y, en un HDD USB, puede ser también desgaste. `smartctl` necesita permisos de root y no pude leerlo.
6. **Los 134 GB de datos del portátil no están en `PAPER1`**: el repo local pesa 1.6 GB (casi todo `.git`). Lo grande del portátil (Música 28 GB, `orca_6_1_1_avx2` 17 GB, miniforge 13 GB, Descargas 12 GB) no es dato del proyecto.

## Hecho (aditivo, verificado)

`/run/media/kmilozeta/TOSHIBA/MAAqc_ARCHIVO/00_snapshot_portatil_2026-11-17/`: copia de `AP`, `DG`, `MGLYVAL`, `docs` (sin `docs/papers`), `maaqc.py`, `test_maaqc.py`, `README.md`, `TODO.md`, etc. 514 archivos, 33 MB; verificación por checksum: 0 diferencias. Incluye `GIT_HEAD.txt` y `GIT_STATUS.txt`.

## Plan propuesto (nada de esto está ejecutado)

**Orden estricto; cada paso solo si el anterior terminó bien.**

0. **Salud del disco primero.** Con permisos de root: `sudo smartctl -a -d sat /dev/sda`, y reparar las entradas NTFS con `chkdsk` desde Windows o `sudo ntfsfix` con el disco desmontado. No reorganizar en un disco que da errores de E/S sin tener antes la segunda copia.
1. **Segundo disco, copia 1:1 y verificada** de `PAPER1`, `PAPER1.2` y `MAAqc_ARCHIVO` (`rsync -a` + pasada `-c` por checksum). Hasta que no exista, el disco actual es la única copia de SHINORINE, PALYTHINE, CALIBRATION y MGLYVAL A/C/D/R. Los `.gbw` y `.densities` no se pueden regenerar sin repetir el cálculo.
2. **Hash completo** de la intersección `PAPER1` ∩ `PAPER1.2` (~25 GB, ~15–25 min por USB) para confirmar duplicados.
3. **Fusión sin renombrar** hacia `MAAqc_ARCHIVO/01_canonico/` conservando los nombres actuales (`CALIBRATION/ASPARTIC`, `MGLYVAL/runs/A…`), para no romper `maaqc.py calibrate --base`: `PAPER1` como base, más lo único de `PAPER1.2` (4.64 GB), más `AP`, `DG` y el MGLYVAL/B más nuevo del portátil. Si un archivo difiere entre versiones se conservan ambos (sufijo `.v2`), no se pisa.
4. **Cuarentena, no borrado**: los duplicados confirmados (≈25 GB) y los restos sueltos (`04c_sp_smd*`, `QCG-2`) se mueven a `MAAqc_ARCHIVO/02_cuarentena/` con un `MANIFIESTO.csv` (ruta original, hash, motivo). Se borra solo cuando el usuario lo diga, tras abrir las corridas importantes.
5. **Candidatos a borrar (requieren tu confirmación explícita)**: `.tmp` (5.6 GB, scratch de ORCA). **No** se proponen para borrar `.gbw`, `.densities`, `.nto` ni `.cis`; si hace falta espacio, se comprimen con `zstd` (el algoritmo no pierde información) en lugar de borrarlos.
6. **Portátil**: tras los pasos 1–4 y verificación, el `PAPER1` local puede quedarse solo con el repo (código, docs, `.out/.inp/.xyz` de texto, que están en git). No se elimina nada local hasta confirmar que está en el disco *y* en la copia.
7. **Papers**: `docs/papers/` (404 MB, 164 PDFs) se copia también a `MAAqc_ARCHIVO/03_papers_zotero/`.

## Decisiones pendientes del usuario

- ¿Segundo disco ya disponible y montado? Condiciona el paso 1.
- ¿Se confirma la fusión **sin renombrar** (recomendada)?
- ¿Se autoriza borrar los `.tmp` (5.6 GB) tras la cuarentena?
- ¿Qué hacer con lo ajeno al proyecto en el portátil que ocupa espacio (Música 28 GB, Descargas 12 GB)? No se tocó.

## Diagnóstico y reparación del TOSHIBA (2026-11-17)

**Hardware: sano.** SMART leído sin root vía `udisks` (`gdbus ... org.freedesktop.UDisks2.Drive.Ata`): `SmartFailing=false`, 0 atributos fallando, 0 sectores defectuosos, último autotest `success`, 32 °C. No hay síntoma de falla física.

**El problema es del sistema de archivos NTFS.** El error de E/S aparece al listar subcarpetas `conf_N` de los restos sueltos `04c_sp_smd*`, más un nombre de archivo con salto de línea. Es el patrón típico de entradas de índice de directorio corruptas, normalmente por una copia interrumpida o por desmontar el disco sin expulsarlo. Está acotado a esos restos: el árbol `PAPER1/` y `PAPER1.2/` se lista y se copia sin errores (ver el log de consolidación).

**Qué se puede y qué no se puede hacer desde Linux.** `ntfsfix` no es un `chkdsk`: arregla algunos problemas, limpia la marca de "sucio" y deja el volumen listo para que Windows lo revise, pero no reconstruye índices de directorio dañados. La herramienta que sí lo hace es `chkdsk` de Windows.

**Orden recomendado (el paso 1 ya está hecho):**
1. Respaldo completo y verificado en `EC IEII/KMILO/MAAqc_BACKUP_FINAL` (hecho). Hasta tener además la segunda copia, no se repara nada en el TOSHIBA.
2. Revisión sin escribir, requiere root (no la pude ejecutar yo): 
   ```bash
   udisksctl unmount -b /dev/sda1
   sudo ntfsfix -n /dev/sda1
   ```
   (`-n` = no actuar; solo informa).
3. Reparación real: conectar el disco a un Windows y ejecutar `chkdsk X: /f` (más rápido) y, si hay dudas sobre sectores, `chkdsk X: /f /r` (varias horas en 1 TB).
4. Como las carpetas dañadas son solo restos sin datos legibles, la alternativa más simple, **después de tener las dos copias verificadas**, es borrarlas desde Windows tras el `chkdsk`.
5. Prevención: expulsar siempre con "Expulsar/Desmontar de forma segura" (`udisksctl power-off -b /dev/sda`) antes de desconectar el disco.
6. Si se quiere dejar de depender de NTFS en Linux: reformatear a ext4 (solo Linux) o exFAT (como el disco nuevo), tras confirmar las dos copias.

## Respaldo final (2026-11-17) — resumen

Destino: `/run/media/kmilozeta/EC IEII/KMILO/MAAqc_BACKUP_FINAL/` (disco XS1000, exFAT, 954 GB). Estructura, reglas y pendientes en su `README_BACKUP.md`. Resultado de la consolidación: el respaldo de octubre de ese disco era un subconjunto exacto de `TOSHIBA/PAPER1` (63 975 archivos idénticos); se reubicó con `mv` y se añadieron SHINORINE, PALYTHINE, QCG (24 909 + 6 451 + 92 archivos), GADUSOL y lo nuevo de `PAPER1.2`, AP, DG, el gold standard (07_adc2, gold_01) y el código del portátil. Los `.git` de TOSHIBA no se copiaron (todos sus commits están en el repo local, 184); en su lugar hay un `git bundle` verificado. Con el disco nuevo en uso: ~89 GB de 954 GB.
