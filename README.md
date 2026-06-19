# MAA-quantum-chemical-protocol
Pipeline computacional sistemático para la caracterización estructural y electrónica de mycosporine-like amino acids (MAAs), diseñado para reemplazar inputs geométricos no validados (PubChem 2D, minimización por campo de fuerza) en estudios de docking molecular y bioactividad por estructuras quantum-mecánicamente optimizadas, con estado de protonación termodinámicamente asignado y muestreo conformacional en fase solución.
 
Proyecto desarrollado para FAE 2026 (Feria Antártica Escolar / INACH), con manuscrito en preparación para *Marine Drugs*.
 
## Motivación
 
Los estudios de docking molecular publicados que postulan a las MAAs como inhibidores de ACE2 y otros blancos farmacológicos humanos parten de estructuras 2D sin optimización cuántica, sin estado de protonación validado termodinámicamente y sin análisis conformacional en solución. Esto es relevante porque trabajo experimental previo (titulaciones UV-Vis combinadas con cálculos RI-CC2/TZVP) estableció que la especie absorbente a pH fisiológico es la forma protonada y estabilizada por zwitterión del cromóforo, no el núcleo neutro. Este proyecto aborda ese vacío metodológico aplicando un protocolo validado de química cuántica a M-Gly:Val (mycosporine-glycine:valine), con foco adicional en el contexto antártico (274 K, alta exposición UV-B).
 
## Protocolo
 
Pipeline de cinco etapas, orquestado en Python:
 
1. **CREST** (iMTD-GC / GFN2-xTB / ALPB-agua) — muestreo conformacional
2. **OPT** — optimización de geometría DFT (ωB97X-D3/def2-SVP, CPCM agua)
3. **NumFreq** — verificación por frecuencias vibracionales numéricas
4. **SP** — energía de punto simple (ωB97X-D3/def2-TZVPD, CPCM agua)
5. **TDDFT** — propiedades electrónicas en estado excitado (transiciones, NTO, ESP, HOMO-LUMO)
Adicionalmente: ciclos termodinámicos de pKa para asignación de estado de protonación a pH fisiológico, y análisis de población de Boltzmann a 298.15 K y 274 K (este último relevante para el contexto antártico).
 
### Stack de software
 
| Componente | Versión |
|---|---|
| ORCA | 6.1.1 |
| CREST | 3.0.2 |
| xTB | 6.7.1 |
| Multiwfn | 2026.6.2 |
 
**Nivel de teoría:** ωB97X-D3 · def2-SVP (optimización/frecuencias) · def2-TZVPD (SP/TD-DFT) · CPCM(agua) implícito · RIJCOSX · TightSCF
 
## Resultados destacados (M-Gly:Val)
 
## Estructura del repositorio
 
```
.
├── pipeline.py              # Orquestador principal (CREST → OPT → NumFreq → SP → TDDFT)
├── results_collector.py     # Recolección y consolidación de resultados
├── MGLYVAL/                 # Directorio de trabajo del compuesto de referencia
└── ...
```

## Estado del proyecto
 
En desarrollo activo. El manuscrito asociado se encuentra en revisión interna (estructura IMRAD, validación de hipótesis y nomenclatura en curso).
 
## Autoría
 
Camilo Zúñiga Jara, en colaboración con Dr. Ricardo Veloso (UFRO) y Dr. Desmond Macleod (UAutónoma).
