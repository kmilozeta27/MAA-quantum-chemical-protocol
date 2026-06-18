# LITERATURA CONSOLIDADA — Proyecto FAE 2026 / Marine Drugs
## Reporte generado: 2026-06-14 | Papers procesados: 9 | Orden: relevancia descendente

---
---

# REPORTE BIBLIOGRÁFICO — Pracht et al. 2024

## 1. Datos de Identificación y Citación

**Título completo**: CREST—A program for the exploration of low-energy molecular chemical space  
**Autores**: Pracht, P., Grimme, S., Bannwarth, C., Bohle, F., Ehlert, S., Feldmann, G., Gorges, J., Müller, M., Neudecker, T., Plett, C., Spicher, S., Steinbach, P., Wesołowski, P. A., Zeller, F.  
**Año**: 2024  
**Revista**: The Journal of Chemical Physics  
**Volumen / Número / Páginas**: 160, 114110  
**DOI**: https://doi.org/10.1063/5.0197592  
**Tipo de estudio**: métodos / computacional  
**Formato de cita (APA)**: Pracht, P., Grimme, S., Bannwarth, C., Bohle, F., Ehlert, S., Feldmann, G., Gorges, J., Müller, M., Neudecker, T., Plett, C., Spicher, S., Steinbach, P., Wesołowski, P. A., & Zeller, F. (2024). CREST—A program for the exploration of low-energy molecular chemical space. *The Journal of Chemical Physics*, *160*, 114110. https://doi.org/10.1063/5.0197592  
**Formato de cita (Vancouver)**: Pracht P, Grimme S, Bannwarth C, Bohle F, Ehlert S, Feldmann G, et al. CREST—A program for the exploration of low-energy molecular chemical space. J Chem Phys. 2024;160:114110. doi:10.1063/5.0197592

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★★★  
**Tipo de relevancia**: [marcar todos los que apliquen]
- [x] Metodología computacional transferible (benchmark funcionales / ciclos termodinámicos / CREST)

**Justificación de relevancia**: Este es el paper de referencia primario para el software CREST 3.0.2, herramienta que constituye el primer paso del protocolo del proyecto (muestreo conformacional). Debe citarse obligatoriamente en §3.x Métodos como la referencia canónica del programa. Describe el algoritmo de iMTD-GC y las capacidades de GFN2-xTB que son exactamente las utilizadas en el pipeline. La versión 2024 refactorizó el backend de cálculo con mejoras significativas de velocidad para moléculas de tamaño farmacéutico (como los MAAs), lo que hace esta cita más actualizada y precisa que la referencia original de 2020.

---

## 3. Puntos Clave para el Proyecto

- **Software CREST**: El programa implementa muestreo conformacional automatizado, optimización de geometría, cálculo de entropía molecular absoluta e identificación de sitios de protonación/desprotonación, todos los cuales son directamente utilizados en el protocolo del proyecto.
- **Base teórica GFNn-xTB**: CREST opera sobre los Hamiltonianos GFNn-xTB y el campo de fuerza GFN-FF, lo que permite un balance costo-precisión apropiado para el muestreo de espacio conformacional de moléculas de tamaño de fármaco antes de refinamiento DFT.
- **Refactorización 2024**: El backend de cálculo refactorizado (publicado en esta versión 2024) provee aceleraciones significativas para moléculas de tamaño medio como los MAAs, haciendo el pipeline más eficiente que con versiones anteriores.
- **Interoperabilidad**: El programa acepta interfaces con cualquier software de química cuántica y de campo de fuerzas mediante creación de archivos de entrada estándar, lo que soporta directamente la integración con ORCA 6.1.1 usada en el proyecto.
- **Paralelización single-node**: Los cálculos se configuran para ejecución concurrente con paralelización eficiente en un único nodo, consistente con el hardware ThinkPad T14 usado en el proyecto.
- **pKa y protonación automatizada**: CREST incluye procedimientos automatizados para la identificación de sitios de protonación y desprotonación, relevante para los ciclos termodinámicos de pKa del protocolo.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
GFN2-xTB (Hamiltoniano semiempírico), GFN-FF (campo de fuerza), metadinámica iMTD-GC, modelo de solvente ALPB (para cálculos conformacionales en agua).

**Moléculas estudiadas**:  
Aplicaciones generales incluyendo moléculas farmacéuticas de tamaño medio; no se estudian MAAs específicamente en este paper.

**Tipo de datos experimentales**: No aplica (paper de métodos/software).

**Limitaciones metodológicas documentadas por los autores**:  
La calidad de la exploración conformacional está limitada por la calidad del potencial utilizado; paisajes de energía de alta fidelidad (DFT) son costosos para muestreos a gran escala. Las interfaces a software externo deben crearse manualmente. El paper reconoce que el GFN2-xTB puede perder features del paisaje de energía que solo un nivel DFT capturaría.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §3.x Methods — CREST | Cita principal del software; describir el algoritmo iMTD-GC y el nivel GFN2-xTB | "Conformational sampling was performed using CREST 3.0.2 [Pracht et al., 2024], employing the iMTD-GC algorithm at the GFN2-xTB level with implicit solvation." |
| §2.x Background teórico | Justificar la elección de CREST frente a otros métodos de muestreo | "CREST provides an efficient, automated exploration of molecular chemical space at the semiempirical level, balancing computational cost and accuracy for drug-sized molecules [Pracht et al., 2024]." |

---

## 6. Fragmentos Parafraseables

1. CREST implements automated conformational sampling using the iMTD-GC algorithm driven by the GFN2-xTB Hamiltonian, providing an efficient route to the low-energy molecular chemical space of drug-sized molecules with minimal user input. [Pracht et al., 2024]

2. The 2024 refactored calculation backend introduced significant speed improvements for sampling of small to medium-sized molecules, enabling more sophisticated multi-layer setups including QM/MM and minimum energy crossing point calculations. [Pracht et al., 2024]

3. CREST includes automated procedures for determining molecular protonation and deprotonation sites, making it directly applicable to studies requiring a systematic evaluation of possible ionic states prior to higher-level quantum chemical treatment. [Pracht et al., 2024]

---

## 7. Notas Adicionales

Este paper es la versión 2024 que reemplaza al paper original de CREST (Pracht et al., *Phys. Chem. Chem. Phys.* 22, 7169, 2020). La versión 2024 es la que corresponde a CREST 3.x utilizado en el proyecto. Siempre citar esta versión más reciente en el manuscrito. Open access bajo CC BY.

---
---

# REPORTE BIBLIOGRÁFICO — Matsuyama et al. 2015

## 1. Datos de Identificación y Citación

**Título completo**: pH-Independent Charge Resonance Mechanism for UV Protective Functions of Shinorine and Related Mycosporine-like Amino Acids  
**Autores**: Matsuyama, K., Matsumoto, J., Yamamoto, S., Nagasaki, K., Inoue, Y., Nishijima, M., Mori, T.  
**Año**: 2015  
**Revista**: The Journal of Physical Chemistry A  
**Volumen / Número / Páginas**: 119(49), 12722–12729  
**DOI**: https://doi.org/10.1021/acs.jpca.5b09988  
**Tipo de estudio**: mixto (experimental + computacional)  
**Formato de cita (APA)**: Matsuyama, K., Matsumoto, J., Yamamoto, S., Nagasaki, K., Inoue, Y., Nishijima, M., & Mori, T. (2015). pH-Independent Charge Resonance Mechanism for UV Protective Functions of Shinorine and Related Mycosporine-like Amino Acids. *The Journal of Physical Chemistry A*, *119*(49), 12722–12729. https://doi.org/10.1021/acs.jpca.5b09988  
**Formato de cita (Vancouver)**: Matsuyama K, Matsumoto J, Yamamoto S, Nagasaki K, Inoue Y, Nishijima M, et al. pH-Independent Charge Resonance Mechanism for UV Protective Functions of Shinorine and Related Mycosporine-like Amino Acids. J Phys Chem A. 2015;119(49):12722–9. doi:10.1021/acs.jpca.5b09988

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★★★  
**Tipo de relevancia**:
- [x] Fotoquímica y espectroscopía de MAAs (experimental)
- [x] Química computacional de MAAs (TD-DFT / conformacional / pKa)

**Justificación de relevancia**: Este es uno de los papers más importantes para el proyecto porque establece experimentalmente y computacionalmente el mecanismo de absorción UV de los MAAs en función del estado de protonación. Demuestra que los grupos amino de los MAAs adoptan la forma zwitteriónica en solución acuosa a pH fisiológico (pH 4–10), y que esta forma zwitteriónica activa el mecanismo de resonancia de carga que explica la absorción en la región UVA. Este resultado es directamente relevante para la elección del estado de protonación que justifica el proyecto (Claim 2 del argumento central) y provee contexto experimental para la validación del protocolo computacional.

---

## 3. Puntos Clave para el Proyecto

- **Independencia de pH**: Los máximos de absorción de mycosporine-glycine y shinorine permanecen prácticamente invariantes en el rango pH 4–10, con solo un pequeño desplazamiento hipsocrómico a pH 1–2, demostrando que el estado de protonación predominante en condiciones fisiológicas es el responsable de la absorción UVA.
- **Mecanismo de resonancia de carga**: La naturaleza zwitteriónica del residuo amino facilita la protonación de los cromóforos 3-aminociclohexenona y 1-amino-3-iminocyclohexeno, y la resonancia de carga en las especies protonadas explica las transiciones de baja energía en la región UVA.
- **Validación computacional**: Cálculos RI-CC2/TZVP sobre sistemas modelo en su forma protonada reproducen bien las energías de transición y las fuerzas de oscilador observadas, con sobreestimaciones sistemáticas menores.
- **El sistema neutro no absorbe en UV-A**: Los cromóforos simples de 3-hidroxi- y 3-aminociclohexenona, así como 1-amino-3-iminocyclohexeno, no poseen la capacidad de absorber luz UV por encima de ~300 nm por sí solos; es el estado protonado-zwitteriónico el que activa esta capacidad.
- **Desactivación no radiativa eficiente**: Las investigaciones de fluorescencia de shinorine en agua confirmaron la efectiva desactivación no radiativa del estado singlete excitado, consistente con la función fotoprotectora.
- **Efecto del pH extremo**: A pH 1–2, la protonación parcial del anión carboxilato en el residuo amino produce un pequeño desplazamiento al azul, explicado por la teoría.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
RI-CC2/TZVP sobre sistemas modelo (sin el residuo amino completo); espectroscopía UV en rango pH 1–10; espectroscopía de fluorescencia de shinorine en agua; HPLC; MS.

**Moléculas estudiadas**:  
Mycosporine-glycine (2b), shinorine (3b), porphyra-334 (3c), 4-deoxygadusol (1a), y sus formas modelo sin residuo amino (2a, 3a, 1a).

**Tipo de datos experimentales**:  
λmax UV en función de pH para mycosporine-glycine y shinorine; espectros de fluorescencia de shinorine; espectros de excitación; datos de HPLC y MS para caracterización.

**Limitaciones metodológicas documentadas por los autores**:  
Los cálculos se realizaron sobre sistemas modelo sin el residuo amino completo (se truncaron por razones computacionales), lo que limita la capacidad de capturar efectos electrónicos del sustituyente amino completo. Los autores reconocen sobreestimaciones sistemáticas menores en las energías de transición y fuerzas de oscilador con RI-CC2.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.3 State of quantum chemical... | Precedente experimental del mecanismo de resonancia de carga y del rol del estado de protonación | "Prior experimental and theoretical work established that the zwitterionic form of MAAs in aqueous solution at physiological pH is responsible for their UV-A absorption via a charge resonance mechanism [Matsuyama et al., 2015]." |
| §3.5 pKa cycles / protonation | Justificar por qué la determinación del estado de protonación es no trivial y requiere metodología rigurosa | "This pH-dependence underscores the importance of accurately determining the protonation state of MAAs in solution, as the absorbing species differs fundamentally from the neutral form [Matsuyama et al., 2015]." |
| §4.3 Results TD-DFT | Comparar λmax simulados con valores experimentales reportados | "Experimental UV absorption maxima for shinorine and mycosporine-glycine in aqueous solution at neutral pH have been reported at 334 and 310 nm, respectively [Matsuyama et al., 2015]." |
| §5.3 Discussion | Conectar resultados computacionales del proyecto con el mecanismo de resonancia de carga | "Our electronic structure results, including NTO and ESP analyses, are consistent with the charge resonance mechanism proposed for related MAAs [Matsuyama et al., 2015]." |

---

## 6. Fragmentos Parafraseables

1. The UV absorption maxima of both mycosporine-glycine and shinorine in aqueous solution are essentially invariant over a wide pH range (pH 4–10), with only a minor blue-shift at strongly acidic conditions (pH 1–2), indicating that the predominant absorbing species at physiological pH has a well-defined protonation state. [Matsuyama et al., 2015]

2. The zwitterionic character of the amino acid residue in MAAs facilitates protonation of the chromophoric core moieties, enabling a charge resonance mechanism in the protonated species that accounts for the allowed low-energy UV-A transition; the unsubstituted chromophores alone are incapable of absorbing above approximately 300 nm. [Matsuyama et al., 2015]

3. Ab initio calculations at the RI-CC2/TZVP level on model MAA systems reproduced the observed transition energies and oscillator strengths with only minor systematic overestimations, confirming the charge resonance interpretation of the UV-protective absorption. [Matsuyama et al., 2015]

4. Fluorescence spectral studies of shinorine in water confirmed effective nonradiative deactivation of the singlet excited state following photoexcitation, consistent with the efficient heat-dissipation mechanism central to MAA photoprotection. [Matsuyama et al., 2015]

---

## 7. Notas Adicionales

Este paper es citado por Hatakeyama & Nakamura 2025 (ref. 14 en ese trabajo) como la referencia canónica del mecanismo de resonancia de carga y del rol del estado zwitteriónico. Tiene datos de λmax experimentales precisos que pueden usarse para la tabla de validación en §4.3. Acceso verificado (ACS). No es open access; verificar disponibilidad institucional.

---
---

# REPORTE BIBLIOGRÁFICO — Woolley et al. 2018

## 1. Datos de Identificación y Citación

**Título completo**: Unravelling the Photoprotection Properties of Mycosporine Amino Acid Motifs  
**Autores**: Woolley, J. M., Staniforth, M., Horbury, M. D., Richings, G. W., Wills, M., Stavros, V. G.  
**Año**: 2018  
**Revista**: The Journal of Physical Chemistry Letters  
**Volumen / Número / Páginas**: 9(11), 3043–3048  
**DOI**: https://doi.org/10.1021/acs.jpclett.8b00921  
**Tipo de estudio**: mixto (experimental + computacional)  
**Formato de cita (APA)**: Woolley, J. M., Staniforth, M., Horbury, M. D., Richings, G. W., Wills, M., & Stavros, V. G. (2018). Unravelling the Photoprotection Properties of Mycosporine Amino Acid Motifs. *The Journal of Physical Chemistry Letters*, *9*(11), 3043–3048. https://doi.org/10.1021/acs.jpclett.8b00921  
**Formato de cita (Vancouver)**: Woolley JM, Staniforth M, Horbury MD, Richings GW, Wills M, Stavros VG. Unravelling the Photoprotection Properties of Mycosporine Amino Acid Motifs. J Phys Chem Lett. 2018;9(11):3043–8. doi:10.1021/acs.jpclett.8b00921

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★★  
**Tipo de relevancia**:
- [x] Fotoquímica y espectroscopía de MAAs (experimental)
- [x] Química computacional de MAAs (TD-DFT / conformacional / pKa)

**Justificación de relevancia**: Paper fundamental del grupo Stavros que por primera vez aplica espectroscopía de absorción transitoria en femtosegundos (TEAS) a subestructuras cromofóricas de MAAs (motivos ACyO y NN), combinado con cálculos ab initio del paisaje de energía de estado excitado. Demuestra experimentalmente que el núcleo ciclohexenimina promueve relajación ultrarrápida vía intersección cónica accesible, mientras que el núcleo ciclohexenona queda atrapado en el mínimo de estado excitado. Estos datos establecen la base mecanicista que cualquier estudio computacional de MAAs debe reproducir o ser consistente con, y son directamente relevantes para la discusión de los resultados TD-DFT del proyecto.

---

## 3. Puntos Clave para el Proyecto

- **Dos mecanismos distintos según el núcleo**: El núcleo ciclohexenona (ACyO) queda atrapado en el mínimo del estado excitado, mientras que el núcleo ciclohexenimine (NN) experimenta relajación ultrarrápida hacia S₀ a través de una intersección cónica accesible mediada por una flexión fuera del plano.
- **Escala temporal de relajación**: Para el sistema NN (núcleo ciclohexenimine), la vida media del estado excitado S₁ es del orden de pocos cientos de femtosegundos, consistente con una intersección cónica de baja barrera.
- **Validación de estudios teóricos previos**: Los resultados experimentales son consistentes con cálculos teóricos previos sobre porphyra-334 y palythine que predicen la flexión fuera del plano como coordenada de relajación hacia la CI.
- **Enfoque bottom-up**: El estudio adopta una estrategia molecular bottom-up, estudiando motivos cromofóricos simplificados para aislar el efecto del núcleo, independientemente de los residuos amino laterales.
- **Absorción en UVA/UVB**: Ambos motivos muestran absorción fuerte en la región UVA/UVB (400–280 nm), confirmando que el núcleo por sí mismo es el responsable de la función fotoprotectora en primer orden.
- **Implicación para diseño**: La identificación del motivo ciclohexenimine como responsable de la relajación ultrarrápida eficiente provee un criterio de diseño racional para nuevos fotoprotectores basados en MAAs.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
TEAS (femtosecond transient electronic absorption spectroscopy); cálculos ab initio del paisaje de energía de estado excitado sobre sistemas modelo (nivel de teoría Pendiente de verificación en el texto completo); análisis GLOTARAN para descomposición de espectros transitorios.

**Moléculas estudiadas**:  
ACyO (3-aminociclohex-2-en-1-ona, motivo núcleo ciclohexenona), NN ((Z)-N-(3-(butilamino)-2-metilciclohex-2-en-1-ilideno)-butan-1-aminio 4-metilbencenosulfonato, motivo núcleo ciclohexenimine).

**Tipo de datos experimentales**:  
Espectros de absorción UV-vis en estado estacionario de ACyO y NN; mapas 2D de absorción transitoria; espectros de absorción de especie de estado excitado en función del tiempo de retardo; constantes de tiempo de decaimiento.

**Limitaciones metodológicas documentadas por los autores**:  
Los resultados se obtienen sobre motivos simplificados, no sobre MAAs completos; los residuos amino laterales pueden modificar cuantitativamente las vidas medias y los paisajes de energía. Los autores reconocen que la extensión a MAAs completos en solución acuosa requerirá estudios adicionales.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.1 MAAs as multifunctional... | Establecer función fotoprotectora y mecanismo de disipación de energía del núcleo | "The cyclohexenimine core of MAAs promotes ultrafast nonradiative decay to the electronic ground state, establishing the molecular basis for their photoprotective function [Woolley et al., 2018]." |
| §1.3 State of quantum chemical... | Citar como precedente experimental para el carácter π→π* de S₁ | "Experimental femtosecond spectroscopy on MAA chromophore motifs confirmed that S₁ relaxation proceeds through a conical intersection accessible via out-of-plane ring flexing [Woolley et al., 2018]." |
| §5.x Discussion | Conectar caracteres NTO del proyecto con el mecanismo experimental establecido | "The π→π* character identified in NTO analysis of S₁ for M-Gly:Val is consistent with the photoprotective relaxation mechanism experimentally elucidated for the cyclohexenimine core [Woolley et al., 2018]." |

---

## 6. Fragmentos Parafraseables

1. Femtosecond transient absorption spectroscopy on isolated MAA chromophore motifs revealed that the cyclohexenimine core undergoes rapid internal conversion to the electronic ground state through an energetically accessible conical intersection reached via an out-of-plane ring flexing coordinate, while the cyclohexenone core becomes trapped in the excited-state minimum. [Woolley et al., 2018]

2. These experimental findings, obtained using a systematic bottom-up approach on simplified MAA chromophore analogues, are consistent with prior theoretical predictions of the excited-state potential energy landscape for porphyra-334 and palythine. [Woolley et al., 2018]

3. The strong UV absorption of both cyclohexenone and cyclohexenimine core motifs in the UVA/UVB spectral region (400–280 nm) establishes the chromophoric core as the primary determinant of the UV-filtering function, independent of the amino acid substituents. [Woolley et al., 2018]

---

## 7. Notas Adicionales

Este paper es una Letter de JPCL, de alto impacto dentro del campo. Es open access (CC-BY). Es el paper fundacional del enfoque bottom-up del grupo Stavros para MAAs, y establece la base experimental para los trabajos posteriores de Whittock et al. 2022 y Hatakeyama & Nakamura 2025. Debe citarse en conjunto con estos trabajos cuando se discuta el mecanismo de disipación de energía.

---
---

# REPORTE BIBLIOGRÁFICO — Hatakeyama & Nakamura 2025

## 1. Datos de Identificación y Citación

**Título completo**: "Button-on-a-String" Mechanism in Water, the Ultrafast UV-to-Heat Conversion by Mycosporine-like Amino Acid Porphyra-334 of Natural Sunscreen Compound  
**Autores**: Hatakeyama, M., Nakamura, S.  
**Año**: 2025  
**Revista**: ACS Physical Chemistry Au  
**Volumen / Número / Páginas**: 5, 274–282  
**DOI**: https://doi.org/10.1021/acsphyschemau.4c00107  
**Tipo de estudio**: computacional (dinámica no adiabática / MD cuántica)  
**Formato de cita (APA)**: Hatakeyama, M., & Nakamura, S. (2025). "Button-on-a-String" Mechanism in Water, the Ultrafast UV-to-Heat Conversion by Mycosporine-like Amino Acid Porphyra-334 of Natural Sunscreen Compound. *ACS Physical Chemistry Au*, *5*, 274–282. https://doi.org/10.1021/acsphyschemau.4c00107  
**Formato de cita (Vancouver)**: Hatakeyama M, Nakamura S. "Button-on-a-String" Mechanism in Water, the Ultrafast UV-to-Heat Conversion by Mycosporine-like Amino Acid Porphyra-334 of Natural Sunscreen Compound. ACS Phys Chem Au. 2025;5:274–82. doi:10.1021/acsphyschemau.4c00107

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★★  
**Tipo de relevancia**:
- [x] Fotoquímica y espectroscopía de MAAs (experimental)
- [x] Química computacional de MAAs (TD-DFT / conformacional / pKa)

**Justificación de relevancia**: Este es el paper más reciente del grupo Nakamura sobre dinámica no adiabática de MAAs en agua, y el primero en caracterizar el mecanismo de conversión UV→calor de porphyra-334 en entorno acuoso explícito mediante dinámica cuántica en primeros principios. El mecanismo propuesto (rotación del anillo ciclohexenimine mientras se mantienen los puentes de hidrógeno con el agua) es directamente relevante para la discusión del estado excitado S₁ del proyecto y para justificar por qué la geometría del estado excitado y las coordenadas de relajación dependen fuertemente del entorno acuoso. Este paper también establece explícitamente que los aminoácidos adoptan la forma zwitteriónica en agua, lo cual valida la elección de estados de protonación del proyecto.

---

## 3. Puntos Clave para el Proyecto

- **Mecanismo "botón-en-cuerda"**: La conversión ultrarrápida UV→calor de porphyra-334 en agua procede mediante rotación del anillo ciclohexenimine intramolecular dentro de varios cientos de femtosegundos post-excitación UV, sin ruptura de los puentes de hidrógeno con el agua circundante.
- **Rol de los aminoácidos laterales**: Los aminoácidos laterales de porphyra-334 mantienen sus posiciones de centro de masa durante la rotación no radiativa del anillo, actuando como el "hilo" que preserva los puentes de hidrógeno; su flexibilidad estructural es esencial para el mecanismo.
- **Dependencia del agua**: La capacidad fotoprotectora de los MAAs depende de la presencia de agua tanto para la absorción UV (modificando la π-conjugación mediante la forma zwitteriónica) como para la disipación de energía post-absorción (mantenimiento de puentes de hidrógeno durante la rotación del anillo).
- **Forma zwitteriónica en agua**: Los aminoácidos de los MAAs adoptan la forma zwitteriónica en agua (citando Matsuyama et al. 2015), lo que modifica la π-conjugación hacia la estructura de resonancia de carga y aumenta la longitud de onda de absorción.
- **Mecanismo intrínseco**: El mecanismo surge de la propia estructura de porphyra-334, no de las moléculas de agua circundantes, aunque el agua es necesaria para la estabilización y preservación de los puentes de hidrógeno durante la rotación.
- **Más de 70 tipos de MAAs**: El paper reafirma que existen más de 70 tipos de MAAs en la naturaleza, variando en los aminoácidos laterales, lo que diversifica sus propiedades fotoprotectoras.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
Dinámica no adiabática en primeros principios (nonadiabatic molecular dynamics) con moléculas de agua explícitas; simulaciones de dinámica molecular; mecánica cuántica en primeros principios (nivel de teoría específico Pendiente de verificación del texto completo). Software: TINKER, potencial TINKER.

**Moléculas estudiadas**:  
Porphyra-334 en agua explícita (incluye box de solvente).

**Tipo de datos experimentales**: No aplica (paper computacional puro).

**Limitaciones metodológicas documentadas por los autores**:  
El estudio se centra en porphyra-334 específicamente; la generalización a otros MAAs con diferentes aminoácidos laterales requeriría estudios adicionales. El rol del agua explícita es esencial, lo que implica que modelos de solvente implícito solos podrían no capturar la dinámica completa.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.1 MAAs as multifunctional... | Describir el mecanismo molecular de la función fotoprotectora | "The photoprotective mechanism of MAAs has been attributed to ultrafast nonradiative decay following UV absorption, dissipating photon energy as heat on a sub-picosecond timescale [Hatakeyama & Nakamura, 2025]." |
| §1.3 State of quantum chemical... | Justificar la necesidad de estructuras precisas en agua para estudios de dinámica | "Accurate simulation of the MAA UV-to-heat conversion mechanism requires the correct aqueous protonation state and conformational geometry, as both directly control the ring rotation coordinate and hydrogen bond network [Hatakeyama & Nakamura, 2025]." |
| §5.x Discussion | Contextualizar los resultados TD-DFT del proyecto en relación con el mecanismo de relajación | "The π→π* character of S₁ identified in this work is consistent with the charge resonance form required for the ring rotation mechanism proposed for porphyra-334 [Hatakeyama & Nakamura, 2025]." |

---

## 6. Fragmentos Parafraseables

1. Porphyra-334 undergoes nonradiative decay within several hundred femtoseconds after UV absorption in water by rotating its intramolecular cyclohexenimine ring while maintaining hydrogen bonds with surrounding water molecules, a process the authors describe as a "button-on-a-string" mechanism. [Hatakeyama & Nakamura, 2025]

2. The amino acid substituents of porphyra-334 preserve their center-of-mass positions during ring rotation, providing the structural flexibility needed to maintain the hydrogen bond network with water throughout the nonradiative decay process. [Hatakeyama & Nakamura, 2025]

3. In aqueous solution, the amino acid residues of MAAs adopt a zwitterionic form that modifies the intramolecular π-conjugation toward the charge resonance structure, thereby increasing the UV absorption wavelength and enabling the UV-A protective function. [Hatakeyama & Nakamura, 2025]

---

## 7. Notas Adicionales

Este paper está bajo CC-BY-NC-ND 4.0. Es la publicación más reciente del grupo Nakamura (Hokkaido University) sobre dinámica de MAAs, grupo que ha producido una serie sobre porphyra-334 (2019, 2022, 2025). El paper de 2025 es el que debe citarse como la referencia de dinámica más actualizada. Ya estaba identificado en el SKILL.md de fae-maa-project como paper "en mano".

---
---

# REPORTE BIBLIOGRÁFICO — Whittock et al. 2022

## 1. Datos de Identificación y Citación

**Título completo**: Investigating the Ultrafast Dynamics and Long-Term Photostability of an Isomer Pair, Usujirene and Palythene, from the Mycosporine-like Amino Acid Family  
**Autores**: Whittock, A. L., Woolley, J. M., Auckloo, N., Corre, C., Stavros, V. G.  
**Año**: 2022  
**Revista**: Molecules  
**Volumen / Número / Páginas**: 27(7), 2272  
**DOI**: https://doi.org/10.3390/molecules27072272  
**Tipo de estudio**: mixto (experimental + computacional)  
**Formato de cita (APA)**: Whittock, A. L., Woolley, J. M., Auckloo, N., Corre, C., & Stavros, V. G. (2022). Investigating the Ultrafast Dynamics and Long-Term Photostability of an Isomer Pair, Usujirene and Palythene, from the Mycosporine-like Amino Acid Family. *Molecules*, *27*(7), 2272. https://doi.org/10.3390/molecules27072272  
**Formato de cita (Vancouver)**: Whittock AL, Woolley JM, Auckloo N, Corre C, Stavros VG. Investigating the Ultrafast Dynamics and Long-Term Photostability of an Isomer Pair, Usujirene and Palythene, from the Mycosporine-like Amino Acid Family. Molecules. 2022;27(7):2272. doi:10.3390/molecules27072272

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★★  
**Tipo de relevancia**:
- [x] Fotoquímica y espectroscopía de MAAs (experimental)
- [x] Química computacional de MAAs (TD-DFT / conformacional / pKa)

**Justificación de relevancia**: Este paper del grupo Stavros aplica espectroscopía de absorción transitoria en femtosegundos y métodos computacionales al par isomérico usujirene/palythene, dos MAAs con núcleo ciclohexenimine y un doble enlace conjugado adicional. Es relevante porque (a) aporta datos experimentales de fotodinámica ultrarrápida sobre MAAs completos (no solo motivos), (b) incluye cálculos DFT/PBE0/6-311++G** de geometrías neutras y zwitteriónicas en agua implícita, mostrando que el confórmero zwitteriónico de palythene es el más estable, y (c) demuestra alta fotostabilidad a largo plazo. Estos datos proporcionan contexto comparativo para la discusión de fotodinámica de M-Gly:Val.

---

## 3. Puntos Clave para el Proyecto

- **Vidas medias ultrarrápidas de S₁**: Los estados excitados electrónicos de usujirene y palythene persisten solo unos pocos cientos de femtosegundos antes de que la población se canalice a través de una intersección cónica energéticamente accesible, con posterior transferencia de energía vibracional al solvente.
- **Estabilidad conformacional zwitteriónica**: A nivel DFT/PBE0/6-311++G** con solvente implícito (COSMO), el confórmero zwitteriónico de palythene es el más estable energéticamente, mientras que el confórmero neutro es menos estable; resultado coherente con Matsuyama et al. 2015.
- **Recuperación incompleta**: Una fracción menor del par isomérico no recuperó su estado original dentro de 3 ns tras la fotoexcitación, consistente con una pequeña probabilidad de fotoisomerización.
- **Alta fotostabilidad a largo plazo**: Experimentos de irradiación continua a una sola longitud de onda y con simulador solar mostraron altos niveles de fotostabilidad en ambas moléculas, validando su función como filtros UV duraderos.
- **Extinción molar alta**: Usujirene y palythene presentan ε ≈ 50,000 dm³ mol⁻¹ cm⁻¹ en UVA, entre los valores más altos del conjunto de MAAs conocidos.
- **Extracción natural**: Las moléculas fueron obtenidas por extracción de producto natural y purificación, no por síntesis química, lo que indica disponibilidad a partir de fuentes biológicas.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
DFT/PBE0/6-311++G** con COSMO (solvente implícito agua) para optimizaciones; CC2/aug-cc-pVDZ para energías verticales de excitación (Pendiente de verificación); TEAS (femtosecond transient electronic absorption spectroscopy) en solución acuosa; irradiación continua con láser monocromático y simulador solar.

**Moléculas estudiadas**:  
Usujirene y palythene (par isomérico geométrico, ambos MAAs con núcleo ciclohexenimine), formas neutra y zwitteriónica.

**Tipo de datos experimentales**:  
Espectros UV-vis en estado estacionario de usujirene y palythene; espectros de absorción transitoria (2D tiempo-longitud de onda); constantes de tiempo de decaimiento; cuantificación de fotostabilidad a largo plazo.

**Limitaciones metodológicas documentadas por los autores**:  
El nivel de teoría DFT/PBE0 puede no capturar completamente la dinámica del estado excitado; los cálculos de geometría se realizaron en la aproximación de Franck-Condon sin relajación del estado excitado. El solvente implícito no reproduce el efecto de los puentes de hidrógeno explícitos.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.1 MAAs as multifunctional... | Citar como evidencia de fotostabilidad alta y capacidad antioxidante | "MAAs such as usujirene and palythene exhibit both ultrafast excited-state deactivation and high long-term photostability, properties that underpin their potential as natural UV filters [Whittock et al., 2022]." |
| §1.3 State of quantum chemical... | Citar como ejemplo de enfoque computacional previo sobre MAAs (DFT sin CREST) | "Previous computational studies on MAAs employed geometry optimization of specific protonation states without systematic conformational sampling [Whittock et al., 2022], a limitation addressed in the present work." |
| §5.x Discussion | Comparar estrategia conformacional del proyecto con la de este paper | "In contrast to approaches that optimize a single MAA conformer [Whittock et al., 2022], the present protocol employs Boltzmann-weighted ensemble averaging over the CREST conformational landscape." |

---

## 6. Fragmentos Parafraseables

1. The electronic excited state lifetimes of usujirene and palythene in aqueous solution are on the order of a few hundred femtoseconds, after which the excited state population is channeled through a conical intersection with subsequent vibrational energy transfer to the solvent. [Whittock et al., 2022]

2. Geometry optimization at the DFT/PBE0 level with implicit water solvation showed that the zwitterionic conformer of palythene is the lowest-energy form, consistent with the proposed role of zwitterionic amino acid residues in enabling the UV-A absorption of MAAs. [Whittock et al., 2022]

3. Both usujirene and palythene exhibited high levels of photostability under continuous UV irradiation at a single wavelength and under solar-simulated conditions, supporting their efficacy as durable UV filters in biological systems. [Whittock et al., 2022]

---

## 7. Notas Adicionales

Open access (CC BY), MDPI Molecules. El paper forma parte de la serie del grupo Stavros (Warwick) sobre MAAs. Whittock et al. 2022 es la secuela directa de Woolley et al. 2018, extendiendo el estudio desde motivos simplificados a MAAs completos. También figura en Hatakeyama & Nakamura 2025 como ref. 16, confirmando reconocimiento mutuo entre grupos de trabajo.

---
---

# REPORTE BIBLIOGRÁFICO — Losantos et al. 2015

## 1. Datos de Identificación y Citación

**Título completo**: Photochemistry and photophysics of mycosporine-like amino acids and gadusols, nature's ultraviolet screens  
**Autores**: Losantos, R., Sampedro, D., Churio, M. S.  
**Año**: 2015  
**Revista**: Pure and Applied Chemistry  
**Volumen / Número / Páginas**: 87(9–10), 979–996  
**DOI**: https://doi.org/10.1515/pac-2015-0304  
**Tipo de estudio**: revisión  
**Formato de cita (APA)**: Losantos, R., Sampedro, D., & Churio, M. S. (2015). Photochemistry and photophysics of mycosporine-like amino acids and gadusols, nature's ultraviolet screens. *Pure and Applied Chemistry*, *87*(9–10), 979–996. https://doi.org/10.1515/pac-2015-0304  
**Formato de cita (Vancouver)**: Losantos R, Sampedro D, Churio MS. Photochemistry and photophysics of mycosporine-like amino acids and gadusols, nature's ultraviolet screens. Pure Appl Chem. 2015;87(9-10):979–96. doi:10.1515/pac-2015-0304

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★★  
**Tipo de relevancia**:
- [x] Fondo de MAAs / revisión general
- [x] Fotoquímica y espectroscopía de MAAs (experimental)
- [x] Química computacional de MAAs (TD-DFT / conformacional / pKa)

**Justificación de relevancia**: Revisión comprensiva (basada en ponencia ICP-16) sobre fotoquímica y fotofísica de MAAs y gadusoles. Cubre estudios experimentales y teóricos sobre el mecanismo fotoprotector molecular, síntesis de MAAs, y diseño racional de derivados. Es relevante como revisión de referencia para §1.1 y §1.3: provee datos de λmax, coeficientes de extinción, discusión de multifuncionalidad, y análisis de estudios computacionales previos (incluyendo trabajo del propio grupo Sampedro). La cobertura de gadusoles como intermediarios en la biosíntesis también es útil para el contexto biosintético.

---

## 3. Puntos Clave para el Proyecto

- **Rango espectral de MAAs**: Las bandas de absorción UV de los MAAs van desde 320 nm hasta ca. 360 nm, con coeficientes de absorción molar entre 28,000 y 60,000 M⁻¹ cm⁻¹, siendo considerados los absorbentes UVA más potentes en la naturaleza.
- **Estructura química diversificada**: Los MAAs consisten en un núcleo ciclohexenona o ciclohexenimine; los MAAs con núcleo ciclohexenimina con glicina en posición 1 y un segundo sustituyente muestran λmax entre 320–360 nm, mientras que las micosporinas con núcleo ciclohexenona absorben entre 310–320 nm.
- **Gadusoles y origen evolutivo**: Los gadusoles y deoxygadusoles (ciclohexenonas naturales) se proponen como intermediarios en la biosíntesis de MAAs y como proto-fotoprotectores en organismos primitivos anteriores a la atmósfera oxigenada.
- **Multifuncionalidad**: Además de la función fotoprotectora, se han postulado funciones osmorreguladoras, de transducción UV en fotosíntesis, de reservorio de nitrógeno intrace lular y de conectividad ecológica para los MAAs.
- **Síntesis química**: La revisión cubre esfuerzos sintéticos para obtener MAAs y estructuras relacionadas simplificadas, útil para contextualizar la disponibilidad de moléculas para estudios experimentales.
- **Estudios computacionales de referencia**: La revisión cita estudios computacionales de Sampedro (2011) sobre exploración del paisaje de energía de potencial de MAAs como punto de partida para la literatura computacional de este campo.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
Revisión de estudios experimentales y computacionales previos; sin metodología computacional propia original (es una revisión).

**Moléculas estudiadas**:  
Shinorine, porphyra-334, palythine, mycosporine-glycine, mycosporine-taurine, palythene, usujirene, palythinol, gadusol, 6-deoxygadusol.

**Tipo de datos experimentales**:  
λmax para múltiples MAAs (tabla de valores); coeficientes de extinción molar; datos de fotodegradación y fotoisomerización; estudios de fluorescencia (rendimiento cuántico, barrido espectral); estudios de fotoestabilidad.

**Limitaciones metodológicas documentadas por los autores**:  
Como revisión, los autores señalan que la química y bases moleculares del potencial fotoprotector de MAAs y gadusoles han sido examinadas en menor extensión que su distribución y presencia biológica.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.1 MAAs as multifunctional... | Datos de rango espectral, ε molar, y multifuncionalidad | "MAAs exhibit UV absorption maxima spanning 320 to 360 nm, with molar absorption coefficients reaching up to 60,000 M⁻¹ cm⁻¹, placing them among the most effective UV-absorbing compounds in nature [Losantos et al., 2015]." |
| §1.3 State of quantum chemical... | Proveer contexto histórico de los estudios computacionales | "Early computational explorations of MAA photoprotection established the relevance of the excited-state potential energy landscape [Losantos et al., 2015], motivating the need for more rigorous structural inputs." |
| §2.x Background | Describir estructura química general de MAAs y gadusoles | "MAAs are low-molecular-weight, water-soluble compounds sharing a cyclohexenone or cyclohexenimine core; those bearing a cyclohexenimine ring with amino acid substituents absorb at longer wavelengths (320–360 nm) than the cyclohexenone-bearing mycosporines (310–320 nm) [Losantos et al., 2015]." |

---

## 6. Fragmentos Parafraseables

1. Mycosporine-like amino acids are water-soluble, low-molecular-weight compounds that absorb UV radiation with maxima between 320 and 360 nm and molar absorption coefficients in the range 28,000–60,000 M⁻¹ cm⁻¹, placing them among the strongest UVA-absorbing natural compounds known. [Losantos et al., 2015]

2. Beyond their well-documented UV-screening role, MAAs have been proposed to serve additional biological functions including osmoregulation, nitrogen storage, and participation in UV transduction during photosynthesis, though their primary ecological significance is as photoprotective agents in marine ecosystems. [Losantos et al., 2015]

3. The molecular basis of MAA photoprotection has been examined through a combination of laser spectroscopy and computational methods, revealing excited-state relaxation pathways that account for the remarkable photostability of these compounds under environmental UV conditions. [Losantos et al., 2015]

---

## 7. Notas Adicionales

Paper de acceso a través de De Gruyter (no OA por defecto). Es una conference paper de ICP-16, lo que significa que puede tener revisión menos rigurosa que un article original; sin embargo, los autores (Sampedro, Churio) son referentes establecidos en el campo. Ya figuraba en el SKILL.md de fae-maa-project como "Losantos et al. 2015 — ecological adaptations." El título exacto del SKILL.md no coincide exactamente ("ecological adaptations") con el contenido del paper (fotoquímica/fotofísica), lo que sugiere que en el SKILL.md se usó una descripción imprecisa del contenido real.

---
---

# REPORTE BIBLIOGRÁFICO — Harada et al. 2022

## 1. Datos de Identificación y Citación

**Título completo**: Molecular Design Learned from the Natural Product Porphyra-334: Molecular Generation via Chemical Variational Autoencoder versus Database Mining via Similarity Search, A Comparative Study  
**Autores**: Harada, Y., Hatakeyama, M., Maeda, S., Gao, Q., Koizumi, K., Sakamoto, Y., Ono, Y., Nakamura, S.  
**Año**: 2022  
**Revista**: ACS Omega  
**Volumen / Número / Páginas**: 7(10), 8581–8590  
**DOI**: https://doi.org/10.1021/acsomega.1c06453  
**Tipo de estudio**: computacional (ML + química cuántica)  
**Formato de cita (APA)**: Harada, Y., Hatakeyama, M., Maeda, S., Gao, Q., Koizumi, K., Sakamoto, Y., Ono, Y., & Nakamura, S. (2022). Molecular Design Learned from the Natural Product Porphyra-334: Molecular Generation via Chemical Variational Autoencoder versus Database Mining via Similarity Search, A Comparative Study. *ACS Omega*, *7*(10), 8581–8590. https://doi.org/10.1021/acsomega.1c06453  
**Formato de cita (Vancouver)**: Harada Y, Hatakeyama M, Maeda S, Gao Q, Koizumi K, Sakamoto Y, et al. Molecular Design Learned from the Natural Product Porphyra-334. ACS Omega. 2022;7(10):8581–90. doi:10.1021/acsomega.1c06453

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★  
**Tipo de relevancia**:
- [x] Fotoquímica y espectroscopía de MAAs (experimental)
- [x] Metodología computacional transferible (benchmark funcionales / ciclos termodinámicos / CREST)
- [x] Otro: Diseño molecular basado en MAAs + ML

**Justificación de relevancia**: Este paper del grupo Nakamura explora el diseño molecular de nuevos compuestos inspirados en porphyra-334 mediante un autoencoder variacional químico (VAE) y búsqueda por similitud en bases de datos. Es relevante en dos sentidos: (a) utiliza propiedades de función de onda cuántica como descriptores para el VAE, lo que muestra una integración de química cuántica con ML aplicada a MAAs, y (b) demuestra que moléculas con capacidad de absorción de luz comparable a porphyra-334 pueden diseñarse computacionalmente. La relevancia para el proyecto es media-alta porque refuerza el argumento de que las propiedades electrónicas cuánticamente calculadas (como las del proyecto) son la base del diseño funcional. No es una referencia directa de metodología del protocolo, pero es útil en §1.1 y §5.x como ejemplo de cómo la comprensión computacional de MAAs habilita aplicaciones.

---

## 3. Puntos Clave para el Proyecto

- **Porphyra-334 como plantilla de diseño**: El paper usa porphyra-334 como molécula semilla para generar nuevos compuestos con capacidad de absorción UV similar, demostrando su relevancia como prototipo molecular para el campo del fotoprotector.
- **VAE química supera búsqueda por similitud**: El método de VAE química generó 138 moléculas que satisfacen los criterios de absorción de luz frente a 6 encontradas por búsqueda de similitud en base de datos, mostrando la superioridad del enfoque generativo para exploración de espacio químico.
- **Propiedades de función de onda como descriptores**: El uso de propiedades de función de onda cuántica como descriptores del VAE permitió encontrar moléculas con geometrías inesperadas pero propiedades ópticas comparables a porphyra-334.
- **Reducción del espesor de la capa de ozono**: El paper contextualiza la necesidad de nuevos fotoprotectores en el contexto del adelgazamiento de la capa de ozono, relevante para la justificación ambiental del proyecto.
- **MAAs como compuestos anti-UV**: El paper reafirma que los MAAs incluyen porphyra-334 y previenen el daño inducido por radiación UV, contribuyendo al fondo de §1.1.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
Chemical Variational Autoencoder (VAE) usando SELFIES como representación molecular; descriptores de función de onda cuántica (no especificados en el texto truncado); DFTB (Density Functional Tight Binding) para evaluación de propiedades de absorción; búsqueda por similitud en bases de datos químicas.

**Moléculas estudiadas**:  
Porphyra-334 (semilla principal), MAAs del conjunto GSEEDS; moléculas generadas por VAE (GVAE) y por búsqueda de similitud (GSIM).

**Tipo de datos experimentales**: No aplica (paper computacional + ML).

**Limitaciones metodológicas documentadas por los autores**:  
Los resultados del VAE dependen del conjunto de entrenamiento y de los descriptores de función de onda elegidos. Las geometrías generadas no fueron validadas experimentalmente; el papel del estado de protonación y de la conformación molecular en los descriptores del VAE no se aborda explícitamente.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.1 MAAs as multifunctional... | Mencionar como evidencia del potencial de MAAs para inspirar diseño de nuevos fotoprotectores | "Porphyra-334 has been used as a molecular template for the design of novel UV-absorbing compounds through quantum-chemistry-guided machine learning approaches [Harada et al., 2022]." |
| §5.x Discussion | Implicaciones del estudio para diseño de fotoprotectores | "The electronic structure characterization developed in the present work provides physically rigorous descriptors that could enhance quantum-chemistry-guided molecular design approaches for MAA-inspired photoprotectors [Harada et al., 2022]." |

---

## 6. Fragmentos Parafraseables

1. Using porphyra-334 as a seed molecule and quantum chemical wave function properties as descriptors, a chemical variational autoencoder generated over 130 candidate molecules satisfying the UV absorption criteria of the natural product, substantially outperforming similarity-search-based database mining approaches. [Harada et al., 2022]

2. The design of environmentally benign photoprotective molecules that replicate the UV-filtering properties of natural MAAs requires computational support integrating quantum chemical properties with generative molecular design strategies. [Harada et al., 2022]

---

## 7. Notas Adicionales

Paper del grupo Nakamura (Hokkaido), el mismo que publicó Hatakeyama & Nakamura 2025. El autor Makoto Hatakeyama aparece también en este paper como coautor, mostrando continuidad del grupo en la intersección ML + fotoquímica de MAAs. Relevancia media para el proyecto: útil como "implicación aplicada" en §5.x pero no es referencia metodológica primaria. ACS Omega es OA.

---
---

# REPORTE BIBLIOGRÁFICO — Baker et al. 2017

## 1. Datos de Identificación y Citación

**Título completo**: Photoprotection: extending lessons learned from studying natural sunscreens to the design of artificial sunscreen constituents  
**Autores**: Baker, L. A., Marchetti, B., Karsili, T. N. V., Stavros, V. G., Ashfold, M. N. R.  
**Año**: 2017  
**Revista**: Chemical Society Reviews  
**Volumen / Número / Páginas**: 46(12), 3770–3791  
**DOI**: https://doi.org/10.1039/c7cs00102a  
**Tipo de estudio**: revisión  
**Formato de cita (APA)**: Baker, L. A., Marchetti, B., Karsili, T. N. V., Stavros, V. G., & Ashfold, M. N. R. (2017). Photoprotection: extending lessons learned from studying natural sunscreens to the design of artificial sunscreen constituents. *Chemical Society Reviews*, *46*(12), 3770–3791. https://doi.org/10.1039/c7cs00102a  
**Formato de cita (Vancouver)**: Baker LA, Marchetti B, Karsili TNV, Stavros VG, Ashfold MNR. Photoprotection: extending lessons learned from studying natural sunscreens to the design of artificial sunscreen constituents. Chem Soc Rev. 2017;46(12):3770–91. doi:10.1039/c7cs00102a

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★★  
**Tipo de relevancia**:
- [x] Fondo de MAAs / revisión general
- [x] Fotoquímica y espectroscopía de MAAs (experimental)

**Justificación de relevancia**: Revisión amplia de alto impacto (Chem. Soc. Rev.) sobre mecanismos de fotoprotección en moléculas naturales y protectores solares comerciales, con énfasis en espectroscopía ultrarrápida de bomba-sonda y química cuántica. Provee el marco conceptual general sobre (a) qué hace una buena molécula filtro UV (alta sección eficaz, mecanismo de disipación de energía sin pérdida de integridad molecular), (b) los mecanismos de desactivación del estado excitado (IC ultrarrápida, IC a estados oscuros), y (c) la "controversia del protector solar" como motor de investigación en el campo. Útil en §1.1 como referencia de amplio contexto para justificar la relevancia del estudio de MAAs como fotoprotectores.

---

## 3. Puntos Clave para el Proyecto

- **Criterios de filtro UV eficiente**: Los requisitos clave para una molécula filtro UV incluyen alta sección eficaz de absorción en UVA y UVB, y disponibilidad de mecanismos de disipación de la energía fotónica sin pérdida de integridad molecular (fotoestabilidad).
- **Historia de los protectores solares**: La revisión ubica los MAAs en la historia más amplia de la fotoprotección química, desde las formulaciones antiguas (aceite de arroz, jazmín, aceite de oliva) hasta los protectores comerciales modernos.
- **Controversia del protector solar**: La controversia sobre la seguridad de algunos componentes de protectores solares ha generado interés en alternativas naturales, contextualizando la relevancia de los MAAs como candidatos.
- **Espectroscopía ultrarrápida como herramienta central**: La revisión describe las técnicas de espectroscopía de bomba-sonda en femtosegundos como las herramientas experimentales clave para estudiar la dinámica de desactivación de moléculas fotoprotectoras.
- **Extensión natural→artificial**: El argumento central de la revisión —que entender la fotofísica de moléculas naturales puede guiar el diseño de protectores artificiales mejorados— es paralelo al argumento del proyecto (estructuras cuánticamente validadas como base para diseño futuro).

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
Revisión de métodos experimentales (espectroscopía de absorción transitoria, espectroscopía de velocidades de iones, espectroscopía de resolución temporal) y computacionales (TDDFT, CASPT2, MRCI) usados en la literatura del campo.

**Moléculas estudiadas**:  
Revisión amplia: ácidos fenólicos, cinamatos, benzofenoles, MAAs (mencionados en contexto de fotoprotección natural), melanina, nucleobases de ADN/ARN, aminoácidos aromáticos.

**Tipo de datos experimentales**: No aplica (revisión).

**Limitaciones metodológicas documentadas por los autores**:  
Como revisión, los autores señalan que el campo de la fotodinámica molecular es relativamente pequeño pero de crecimiento rápido; muchos mecanismos permanecen bajo debate.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.1 MAAs as multifunctional... | Contexto amplio sobre fotoprotección natural y MAAs como protectores solares naturales | "Natural photoprotective compounds, including MAAs, satisfy the key requirements of effective UV filters: high absorption cross-section and efficient non-destructive energy dissipation [Baker et al., 2017]." |
| §1.3 State of quantum chemical... | Mencionar la convergencia de métodos experimentales y computacionales en el campo | "The integration of ultrafast spectroscopy and quantum chemical calculations has become the standard framework for characterizing photoprotective mechanisms at the molecular level [Baker et al., 2017]." |

---

## 6. Fragmentos Parafraseables

1. Effective UV filter molecules for photoprotective applications require both a large absorption cross-section across the UVA and UVB spectral regions and the availability of one or more mechanisms for dissipating the absorbed photon energy without compromising the molecular integrity of the filter compound. [Baker et al., 2017]

2. Growing concerns about the safety of some commercial sunscreen components have stimulated research into natural photoprotective molecules, including MAAs, as alternatives that combine high photostability with biological compatibility. [Baker et al., 2017]

---

## 7. Notas Adicionales

Open access (Creative Commons Attribution 3.0). Chem. Soc. Rev. tiene alto impacto (IF ~40+); esta revisión es una referencia canónica del campo de la fotofísica de protectores solares. Autores incluyen V. G. Stavros (Warwick), quien también es coautor de Woolley et al. 2018 y Whittock et al. 2022; hay coherencia en el grupo de trabajo. Citar con precaución para evitar sobre-citar revisiones amplias cuando existen papers originales más específicos.

---
---

# REPORTE BIBLIOGRÁFICO — Roberts & Stavros 2014

## 1. Datos de Identificación y Citación

**Título completo**: The role of πσ* states in the photochemistry of heteroaromatic biomolecules and their subunits: insights from gas-phase femtosecond spectroscopy  
**Autores**: Roberts, G. M., Stavros, V. G.  
**Año**: 2014  
**Revista**: Chemical Science  
**Volumen / Número / Páginas**: 5(5), 1698–1722  
**DOI**: https://doi.org/10.1039/c3sc53175a  
**Tipo de estudio**: revisión (Perspective)  
**Formato de cita (APA)**: Roberts, G. M., & Stavros, V. G. (2014). The role of πσ* states in the photochemistry of heteroaromatic biomolecules and their subunits: insights from gas-phase femtosecond spectroscopy. *Chemical Science*, *5*(5), 1698–1722. https://doi.org/10.1039/c3sc53175a  
**Formato de cita (Vancouver)**: Roberts GM, Stavros VG. The role of πσ* states in the photochemistry of heteroaromatic biomolecules and their subunits: insights from gas-phase femtosecond spectroscopy. Chem Sci. 2014;5(5):1698–722. doi:10.1039/c3sc53175a

---

## 2. Relevancia para el Proyecto MAA

**Puntuación de relevancia**: ★★  
**Tipo de relevancia**:
- [x] Fotoquímica y espectroscopía de MAAs (experimental)
- [x] Otro: Contexto mecanicista sobre estados πσ* en biomoléculas heteroaromáticas

**Justificación de relevancia**: Artículo de perspectiva de Chemical Science sobre el papel de los estados πσ* disociativos en la fotoquímica de biomoléculas heteroaromáticas (bases de ADN/ARN, aminoácidos aromáticos), estudiados mediante espectroscopía de femtosegundos en fase gaseosa. No aborda MAAs directamente, pero provee el marco conceptual general sobre fotodinámica de biomoléculas que es el antecedente intelectual de los estudios de Stavros sobre MAAs (Woolley et al. 2018, Whittock et al. 2022). Relevancia secundaria para el proyecto: podría citarse en §1.3 como referencia metodológica amplia si se discute el uso de espectroscopía ultrarrápida como herramienta de validación de mecanismos fotoprotectores.

---

## 3. Puntos Clave para el Proyecto

- **Estados πσ* como mediadores de fotodegradación**: En biomoléculas heteroaromáticas, los estados πσ* electrónicamente excitados y disociativos contribuyen a la fotodegradación al promover transferencia de hidrógeno foto-inducida o disociación N-H/O-H.
- **Resistencia fotoquímica como criterio selectivo evolutivo**: La review propone que la fotorresistencia de las bases de ADN/ARN y aminoácidos aromáticos puede reflejar selección evolutiva para moléculas que disipan eficientemente la energía UV sin daño fotoquímico.
- **Espectroscopía de bomba-sonda en fase gaseosa**: La revisión describe los fundamentos de la espectroscopía de velocidades de iones con resolución temporal como herramienta para estudiar la dinámica de desactivación fotoinducida en sistemas modelo aislados.
- **Estrategia bottom-up**: El enfoque de estudiar primero biomoléculas simples aisladas y luego escalar a sistemas más complejos es el mismo adoptado por Woolley et al. 2018 para MAAs.

---

## 4. Metodología Relevante

**Nivel de teoría / métodos usados**:  
Revisión de métodos: espectroscopía de velocidades de iones con resolución temporal (TR-VMI), espectroscopía de fotoelectrones con resolución temporal; perspectiva computacional sobre estados πσ*.

**Moléculas estudiadas**:  
Bases de ADN/ARN (adenina, uracilo, etc.), aminoácidos aromáticos (fenilalanina, tirosina, triptófano), fenol, imidazol, indol, pirimidinas. No MAAs.

**Tipo de datos experimentales**: Revisión de datos de la literatura.

**Limitaciones metodológicas documentadas por los autores**:  
Los experimentos en fase gaseosa aíslan la molécula del entorno, limitando la transferibilidad directa a condiciones solución. La relevancia de los mecanismos en fase gaseosa para sistemas en solución acuosa (como los MAAs) requiere estudios puente adicionales.

---

## 5. Ubicación en el Manuscript

| Sección del paper MAA | Uso específico | Frase de transición sugerida |
|----------------------|----------------|------------------------------|
| §1.3 State of quantum chemical... | Citar como antecedente metodológico del campo de la fotodinámica de biomoléculas | "The role of electronically excited states in UV photoprotection has been systematically studied in model biological chromophores using gas-phase femtosecond spectroscopy [Roberts & Stavros, 2014], an approach subsequently extended to MAAs in aqueous solution." |

---

## 6. Fragmentos Parafraseables

1. Gas-phase femtosecond spectroscopy studies on heteroaromatic biomolecules and their chromophore subunits have revealed that photostability often arises from rapid internal conversion pathways involving electronically excited states, providing a conceptual framework applicable to the study of natural photoprotective compounds such as MAAs. [Roberts & Stavros, 2014]

---

## 7. Notas Adicionales

Artículo de perspectiva en Chemical Science (RSC), open access. Útil como referencia metodológica de contexto amplio pero no es una cita obligatoria para el proyecto MAA. Si se incluye, debe ser en la introducción como referencia amplia junto con Baker et al. 2017, no como cita central. El paper es un precedente directo del programa de investigación del grupo Stavros sobre MAAs que culmina en Woolley et al. 2018.

---
---

# ÍNDICE CONSOLIDADO DE RELEVANCIA

| Rank | Primer autor | Año | Revista | Estrellas | Uso principal en manuscript |
|------|-------------|-----|---------|-----------|----------------------------|
| 1 | Pracht | 2024 | J Chem Phys | ★★★★★ | §3.x Methods — cita canónica CREST |
| 2 | Matsuyama | 2015 | J Phys Chem A | ★★★★★ | §1.3, §3.5, §4.3, §5.3 — mecanismo de protonación / λmax experimental |
| 3 | Woolley | 2018 | J Phys Chem Lett | ★★★★ | §1.1, §1.3, §5.x — mecanismo fotoprotector experimental (TEAS + ab initio) |
| 4 | Hatakeyama & Nakamura | 2025 | ACS Phys Chem Au | ★★★★ | §1.1, §1.3, §5.x — dinámica no adiabática en agua; mecanismo "botón-en-cuerda" |
| 5 | Whittock | 2022 | Molecules | ★★★★ | §1.1, §1.3, §5.x — fotodinámica ultrarrápida + fotostabilidad de MAAs completos |
| 6 | Losantos | 2015 | Pure Appl Chem | ★★★★ | §1.1, §1.3, §2.x — revisión de fotoquímica/fotofísica de MAAs; datos de λmax y ε |
| 7 | Harada | 2022 | ACS Omega | ★★★ | §1.1, §5.x — diseño molecular de MAAs + ML cuántico |
| 8 | Baker | 2017 | Chem Soc Rev | ★★★ | §1.1, §1.3 — contexto amplio de fotoprotección natural vs. artificial |
| 9 | Roberts & Stavros | 2014 | Chem Sci | ★★ | §1.3 — antecedente metodológico espectroscopía ultrarrápida en biomoléculas |

---

## Gaps de Literatura Identificados

A partir de los papers procesados, se confirman los siguientes gaps de citas que el proyecto aún necesita antes de finalizar §1.2 y §5.3:

1. **Bioactividad de MAAs (antiviral, anti-inflamatorio, antioxidante, docking SARS-CoV-2)**: Ninguno de los 9 papers cubre este tema directamente. Es el gap más urgente para el argumento central del paper (Claim 1). Buscar: "mycosporine docking SARS-CoV-2", "MAA antiviral", "mycosporine anti-inflammatory 2020–2025".

2. **M-Gly:Val específicamente**: Ninguno de los 9 papers cita o estudia Mycosporine-glycine:valine directamente. Se necesitan papers de identificación, distribución, o propiedades de M-Gly:Val.

3. **Benchmark ωB97X-D3 para cromóforos orgánicos**: No hay papers de benchmark del funcional usado en el proyecto. Buscar: grupo Jacquemin, benchmark TD-DFT range-separated hybrids, organic chromophores.

4. **Contexto antártico / Southern Ocean UV-B**: Ninguno de los 9 papers cubre el contexto antártico. Necesario para Claim 4. Buscar: Karentz 1991, ozone depletion Southern Ocean, MAAs Antarctic organisms.

5. **Ciclos termodinámicos de pKa con DFT**: No hay papers metodológicos de pKa computacional. Necesario para §3.5. Buscar: Born-Haber pKa cycle DFT, Pracht & Grimme 2021 pKa.

---

*Fin del archivo consolidado — 9 reportes individuales + índice + gaps identificados*
