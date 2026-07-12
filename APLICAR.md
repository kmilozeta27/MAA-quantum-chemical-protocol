# Cómo aplicar el parche (paso a paso, sin necesitar `patch`)

## 1. Backup del archivo original

```bash
cp /home/kmilozeta/CENSO/src/censo/processing/orca_processor.py \
   /home/kmilozeta/CENSO/src/censo/processing/orca_processor.py.bak-pre-slots-patch
```

## 2. Agregar el import de ENVIRON y hashlib (arriba del archivo)

Buscar este bloque (líneas ~31-36 del archivo original):

```python
from ..params import (
    WARNLEN,
    USER_ASSETS_PATH,
    OrcaSolvMod,
    Prog,
)
```

Reemplazarlo por:

```python
import hashlib

from ..params import (
    WARNLEN,
    USER_ASSETS_PATH,
    OrcaSolvMod,
    Prog,
    ENVIRON,
)
```

## 3. Insertar el helper de slots

Pegar el contenido completo de `_maaqc_slot_helper_snippet.py` (adjunto)
inmediatamente después de la línea `logger = setup_logger(__name__)`
(línea ~40 del original), antes de `@final class OrcaProc`.

## 4. Parchear `_sp()` — reemplazar la llamada a ORCA

Buscar (dentro de `_sp()`, línea ~506-508 del original):

```python
        # call orca
        call = [config.paths.orca, f"{filename}.inp"]
        returncode, _ = self._make_call(call, outputpath, jobdir)
```

Reemplazar por:

```python
        # call orca
        call = [config.paths.orca, f"{filename}.inp"]
        # MAAQC-PATCH: env explícito por-conformero para evitar que SP
        # paralelos de CENSO compitan por los mismos cores (ver README.md
        # junto a este parche). Fail-open: si MAAQC_CENSO_SLOTS no está
        # fijada, _maaqc_build_slot_env devuelve None y _make_call cae a
        # ENVIRON como antes del parche.
        _maaqc_env = _maaqc_build_slot_env(job.conf.name)
        returncode, _ = self._make_call(call, outputpath, jobdir, env=_maaqc_env)
```

## 5. Parchear `opt()` — reemplazar la llamada a ORCA

Buscar (dentro de `opt()`, línea ~653-655 del original):

```python
        # call orca
        call = [config.paths.orca, f"{filename}.inp"]
        returncode, _ = self._make_call(call, outputpath, jobdir)
```

⚠️ Este bloque es **textualmente idéntico** al de `_sp()` (paso 4) — tu
editor puede confundir cuál es cuál si usás buscar/reemplazar automático.
Confirmá que estás editando el que está dentro de `def opt(` (unas 150
líneas más abajo que el de `_sp()`) y no lo apliques dos veces al mismo.

Reemplazar por:

```python
        # call orca
        call = [config.paths.orca, f"{filename}.inp"]
        # MAAQC-PATCH: ver nota idéntica en _sp() más arriba.
        _maaqc_env = _maaqc_build_slot_env(job.conf.name)
        returncode, _ = self._make_call(call, outputpath, jobdir, env=_maaqc_env)
```

## 6. Verificar sintaxis antes de correr nada real

```bash
conda activate sciencebitch
python3 -c "import ast; ast.parse(open('/home/kmilozeta/CENSO/src/censo/processing/orca_processor.py').read())" && echo "SINTAXIS OK"
python3 -c "from censo.processing.orca_processor import OrcaProc; print('import OK')"
```

Si el segundo comando falla con ImportError, típicamente es el import de
ENVIRON mal ubicado (paso 2) — revisar que quedó dentro del bloque
`from ..params import (...)` y no duplicado en otro lado.
