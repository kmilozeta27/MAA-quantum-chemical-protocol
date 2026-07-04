# MAA-quantum-chemical-protocol
A systematic computational pipeline for the structural and electronic characterization of mycosporine-like amino acids (MAAs), designed to replace unvalidated geometry inputs (PubChem 2D, force-field minimization)—typically used in molecular docking and bioactivity studies—with quantum-mechanically optimized structures featuring thermodynamically assigned protonation states and solution-phase conformational sampling.
 
Project developed for FAE 2026 (Antarctic School Fair / INACH), with a manuscript currently in preparation for *Marine Drugs*.
 
## Motivation
 
Published molecular docking studies proposing MAAs as inhibitors of ACE2 and other human pharmacological targets rely on 2D structures lacking quantum optimization, thermodynamic validation of protonation states, and solution-phase conformational analysis. This is significant because prior experimental work (UV-Vis titrations combined with RI-CC2/TZVP calculations) established that the species absorbing at physiological pH is the protonated, zwitterion-stabilized form of the chromophore, rather than the neutral core (Matsuyama et al., 2015). This project addresses this methodological gap by applying a validated quantum chemistry protocol to M-Gly:Val (mycosporine-glycine:valine), with an additional focus on the Antarctic context (274 K, high UV-B exposure) and to other relevant MAAs such as shinorine, palythine, mycosporine-glycine and prophyra-334.
 
## Protocol
 
Five-stage pipeline, orchestrated in Python:
 
1. **CREST** (iMTD-GC / GFN2-xTB / ALPB-water) — conformational sampling
2. **CENSO** (pre-screening, screening and optimization) — conformational energetic sorting
3. **OPT** — DFT geometry optimization (ωB97X-D3/def2-SVP, CPCM water)
4. **NumFreq** — verification via numerical vibrational frequencies
5. **SP** — single-point energy (ωB97X-D3/def2-TZVPD, CPCM water)
6. **TDDFT** — excited-state electronic properties (transitions, NTO, ESP, HOMO-LUMO)
7. **ADC2** (def2-SVP, gas-phase) — calibration of the TD-DFT results (systematic overstimations and blue-shift)
Additionally: pKa thermodynamic cycles for protonation state assignment at physiological pH, and Boltzmann population analysis at 298.15 K and 274 K (the latter relevant to the Antarctic context).
 
**Nivel de teoría:** ωB97X-D3 · ACD2 · def2-SVP (optimización/frecuencias) · def2-TZVPD (SP/TD-DFT/ACD2) · CPCM(agua) implícito · Gas-Phase (ACD2) · RIJCOSX · TightSCF
 
## Repo structure
 
```
.
├── pipeline.py              # Orchestrator script (CREST → CENSO → OPT → NumFreq → SP → TDDFT → ACD2)
├── protonation_cycle.py     # Theoretical pKa calculations and determination of the relevant protonation state
├── results_collector.py     # Collection and consolidation of results across all explored states
├── MGLYVAL/                 # Work directory for the subject of study (M-Gly:Val)
└── runs/                    # Relevant Files (.inp, .out, .xyz)
```

## Project state
 
In active development.
 
## Authorship
 
Camilo Zúñiga Jara, in collaboration with Dr. Ricardo Veloso (UFRO) and Dr. Desmond Macleod (UAutónoma).
