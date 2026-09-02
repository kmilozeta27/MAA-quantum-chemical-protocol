# MAA-quantum-chemical-protocol
A systematic computational pipeline for the structural and electronic characterization of mycosporine-like amino acids (MAAs), designed to replace unvalidated geometry inputs (PubChem 2D, force-field minimization)—typically used in molecular docking and bioactivity studies—with quantum-mechanically optimized structures featuring thermodynamically assigned protonation states and solution-phase conformational sampling.

Manuscript currently in preparation for *Marine Drugs*.
 
## Motivation
 
Several published molecular docking studies propose MAAs as inhibitors of human pharmacological targets (e.g., ACE2, Sars-CoV2), building their structural inputs from representations that carry no quantum-mechanical validation, no thermodynamic treatment of protonation state, and no solution-phase conformational analysis. This matters because it isn't a merely technical omission. Prior experimental work combining UV-Vis titrations with RI-CC2/TZVP calculations (Matsuyama et al., 2015) showed that the species responsible for absorption at physiological pH is the protonated, zwitterion-stabilized form of the chromophore, not the neutral core that a 2D-to-3D conversion would generate. A docking study built on the wrong protonation state is, in a strict sense, modeling the wrong molecule. This project applies a validated quantum-chemistry protocol to close that gap for M-Gly:Val (mycosporine-glycine:valine), with an explicit extension to Antarctic environmental conditions (274 K, high UV-B exposure), and a design intended to generalize to other MAAs of interest such as shinorine, palythine, mycosporine-glycine, porphyra-334, etc.
 
## Protocol
 
Five-stage pipeline, orchestrated in Python:
 
1. **CREST** (iMTD-GC / GFN2-xTB / ALPB-water) — conformational sampling
2. **CENSO** (pre-screening, screening and optimization) — conformational energetic sorting
3. **OPT** — DFT geometry optimization (ωB97X-D3/def2-SVP, CPCM water)
4. **NumFreq** — verification via numerical vibrational frequencies
5. **SP** — single-point energy (ωB97X-D3/def2-TZVPD, SMD water)
6. **TDDFT** — excited-state electronic properties (transitions, NTO, ESP, HOMO-LUMO)
7. **ADC2** (def2-SVP, gas-phase) — calibration of the TD-DFT results (systematic overstimations and blue-shift)
- pKa thermodynamic cycles for protonation state assignment at physiological pH, and Boltzmann population analysis at 298.15 K and 274 K (the latter relevant to the Antarctic context).
 
**Theory Level:** ωB97X-D3 · ADC2 · def2-SVP · def2-TZVPD (SP/TD-DFT/ACD2) · SMD(water) implícito · CPCM(water) implícito · Gas-Phase · RIJCOSX · TightSCF
 
## Repo structure
 
```
.
├── pipeline.py              # Orchestrator script (CREST → CENSO → OPT → NumFreq → SP → TDDFT → ACD2)
├── protonation_cycle.py     # Theoretical pKa calculations and determination of the relevant protonation state
├── results_collector.py     # Collection and consolidation of results across all explored states
├── CALIBRATIONS/            # Molecules used for pKa cycle calibration and ADC2 cross-check
├── MGLYVAL/                 # Work directory for the subject of study (M-Gly:Val)
└── runs/                    # Relevant Files (.inp, .out, .xyz)
```
## Why this design
The pipeline is built the way production quantum-chemistry software is built, not the way a one-off calculation is run: each stage is a discrete, resumable unit with a defined input/output contract, so a failure at NumFreq or a stalled TD-DFT job doesn't require re-running CREST from scratch. This mirrors the orchestration logic of tools like CREST/CENSO/ORCA themselves

## Expected impact
 
- **Methodological:** a reusable, generalizable protocol for assigning validated 3D structures to MAAs. Directly transferable to shinorine, palythine, mycosporine-glycine, and porphyra-334, and in principle to other small, protonation-sensitive photoactive natural products.
- **For bioactivity research:** a corrected structural basis for any future docking or QSAR study involving MAAs, addressing a validation gap in the existing published literature rather than adding another isolated data point to it.
- **For astro/photobiology:** explicit characterization of MAA photochemical behavior under Antarctic-relevant conditions (274 K, high UV-B), relevant to understanding UV-photoprotection in polar marine organisms and, more broadly, chromophore stability under extreme-environment radiation regimes.
- **Output:** manuscript in preparation for *Marine Drugs*.

## Project state
 
In active development.
 
## Authorship
 
Camilo Zúñiga Jara, in collaboration with Mg. Tania Tapia (UFRO), Dr. Ricardo Veloso (UFRO) and Dr. Desmond Macleod (UAutónoma de Santiago).