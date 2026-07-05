# Fixtures de test

Usar una molécula chica (glicina) para correr el pipeline completo en
menos de 1 minuto sin ORCA real (mockear la llamada) o con GFN2-xTB real
(rápido). **NO usar MGLYVAL** — es demasiado grande para CI y es además
el sistema real del paper, mezclar fixtures de test con datos de
investigación activa es una fuente fácil de confusión.

Pendiente: agregar `glycine.xyz` + outputs `.out` de referencia
(generados una vez, a mano, y congelados) para los tests de parsing en
`io/parsing.py` y `orchestrator/state.py`.
