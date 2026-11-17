# docs/papers

Copia ordenada de los PDFs de la biblioteca Zotero (carpeta `storage/` del snap de Zotero).

- Generado con un script Python a partir de `zotero.sqlite` (copia de solo lectura en un directorio temporal; integrity_check ok). No se uso el `.bak`.
- Origen (`~/snap/zotero-snap/common/Zotero/`) NO se toco: los PDFs se copiaron (nunca movidos) y cada copia se verifico por tamano y sha256.
- Un paper por carpeta (tema principal); los temas secundarios y las colecciones de Zotero figuran en `INDEX.md` / `index.csv`.
- Archivos: `INDEX.md`, `index.csv` (con clave Zotero, ruta original y sha256), `references.bib` (BibTeX basico), `duplicados_y_huerfanos.md`.
- Los PDFs estan en `.gitignore` (`docs/papers/**/*.pdf`); solo se versionan los indices.
- Los papers de filosofía (Kant, qualia, estética) se retiraron de este repo el 2026-11-17 y viven en el disco de respaldo, `MAAqc_BACKUP_FINAL/05_otros_no_proyecto/`; los originales siguen en Zotero.
