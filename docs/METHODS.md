# Computational methods

This document summarises Section 4 of the manuscript in preparation. Numerical results are not
reported here.

## 1. Workflow overview

The whole protocol runs from one program, `maaqc.py`. The subcommands are `run` (computation
stages), `status`, `collect` (analysis), `pka` (protonation cycle and speciation), `calibrate`
(calibration against reference compounds) and `cubes` (orbital and NTO cube extraction). Before
running, each stage checks whether its output is already complete, and resumes instead of
recomputing.

## 2. Structures and protonation-state convention

Starting geometries come from PubChem and are converted to 3D automatically. Candidate
protonation states are labelled with a single-letter suffix, and the pipeline reads the formal
charge from that suffix:

| Label | Species | Charge |
|---|---|---|
| A | Cation: ammonium/amidinium with all carboxyl groups protonated | +1 |
| R | Neutral reference tautomer | 0 |
| D | Zwitterion (tautomer of R) | 0 |
| B | Monoanion | −1 |
| C | Dianion | −2 |

The general chain is **A⁺ ⇌ (R ⇌ D) ⇌ B⁻ ⇌ C²⁻**. R and D have the same charge and the same
atom count, so they are merged into one Boltzmann-weighted neutral macrostate. Their tautomeric
fractions do not depend on pH. Compounds with a single carboxyl group have a shorter chain.
4-Deoxygadusol ionises through its enol. Before any calculation, the pipeline checks the charge
and atom-count relationship of each pair of states.

## 3. Conformational sampling and ensemble refinement

```bash
crest crest_input.xyz --T <N> --chrg <q> --uhf 0 --alpb water
censo --prescreening --screening -i <ensemble> -c <q> -u 0 \
      --solvent water --maxcores <N> --omp-min 4 -T <temperature>
```

Conformers go on to DFT according to their CENSO screening Boltzmann weight. The threshold is
5 %, with at least three conformers kept. If screening weights are not available, an
energy-window fallback is used. A multi-temperature mode runs a separate CENSO ensemble for
each temperature.

## 4. DFT optimisation, frequencies and single points (ORCA 6.1.1)

```text
! wB97X-D3 def2-SVP   CPCM(water) RIJCOSX TightSCF TightOpt           (optimisation)
! wB97X-D3 def2-SVP   CPCM(water) RIJCOSX TightSCF NumFreq            (frequencies)
! wB97X-D3 def2-TZVPD def2/J CPCM(water) RIJCOSX TightSCF defgrid3    (aqueous SP)
! wB97X-D3 def2-TZVPD def2/J RIJCOSX TightSCF defgrid3                (gas-phase SP)
! wB97X-D3 def2-TZVPD def2/J SMD(water) RIJCOSX TightSCF defgrid3     (SMD SP)
```

Thermal corrections come from one GoodVibes run per frequency output. The harmonic (G) and
quasi-harmonic (qh-G, Grimme interpolation, 100 cm⁻¹ cutoff) values are read from the same
table, so the two differ only in the damping. Free energies are evaluated at **250.0 K,
274.15 K (Antarctic surface ocean), 298.15 K and 310.0 K**.

## 5. Thermodynamic pKa cycle

The ensemble free energy is used in place of the lowest conformer:

$$G_\mathrm{ens} = G_\mathrm{min} - RT\ln\sum_i \exp\!\left(-\frac{G_i - G_\mathrm{min}}{RT}\right)$$

The aqueous free energy of the proton is assembled from its components:

$$G(\mathrm{H^+,aq}) = \underbrace{-6.28}_{\text{gas (Tissandier 1998)}} + \underbrace{1.894}_{RT\ln 24.5} + \underbrace{(-265.9)}_{\Delta G_\mathrm{solv}\ \text{(Kelly 2006)}} = -270.286\ \mathrm{kcal\,mol^{-1}}$$

Its gas-phase and standard-state terms are made temperature-dependent in closed form for the
ideal monatomic gas. They reduce exactly to the values above at 298.15 K.

$$\Delta G_\mathrm{aq} = G^\mathrm{aq}_\mathrm{total}(\mathrm{A^-}) - G^\mathrm{aq}_\mathrm{total}(\mathrm{AH}) + G(\mathrm{H^+,aq}), \qquad \mathrm{p}K_a = \frac{\Delta G_\mathrm{aq}}{RT\ln 10}$$

Once the proton reference is unified, the Born cycle and the direct continuum route are
algebraically identical. Their agreement is therefore reported as an *implementation check*,
not as independent validation. pKa values are computed with both CPCM and SMD single points.

**About the proton constant.** Some implementations of the direct cycle use −264.61 kcal/mol as
the total proton free energy. That value is a solvation free energy only. Using it as a total
shifts every pKa by about +4.2 units, which happened in an early round of this project.

## 6. Class-resolved empirical calibration

| Class | Charge change | Group | Calibrators |
|---|---|---|---|
| `carboxyl` | +1 → 0 | first COOH | Gly, Val, Asp(α), Glu(α) |
| `carboxyl_2` | 0 → −1 | second COOH | Asp(β), Glu(γ) |
| `amino` | 0 → −1 | NH₃⁺, no neighbouring charge | Gly, Val |
| `amino_dianion` | −1 → −2 | NH₃⁺ flanked by two COO⁻ | Glu |
| `pyridinium` | +1 → 0 | aromatic N–H⁺ | 4-aminopyridine |

- Each class gets its own least-squares fit. A free slope is used only if the class has at
  least three calibrators spanning at least one pKa unit. Otherwise the fit is offset-only.
- RMSE and leave-one-out RMSE are reported per class. A class with a single calibrator is
  marked as having undefined RMSE and is not applied without an explicit override.
- **The pipeline refuses to apply a calibrator across charge classes.** For glutamic acid's
  −1 → −2 step, the raw cycle deviates from experiment by about 8 pKa units. This class cannot be
  transferred from the `amino` class.
- The primary specification (solvation model, quasi-harmonic corrections, automatic fit mode) is
  fixed *before* any MAA pKa is examined. All other combinations are reported as sensitivity
  analyses.

## 7. Uncertainty propagation

A Monte Carlo over the full protonation chain (10 000 samples by default) samples three error
sources:

- the electronic energy, σ = 0.5 kcal/mol per state;
- the thermal correction, σ = 0.4 kcal/mol per state;
- G(H⁺,aq), σ = 2.5 kcal/mol.

The proton perturbation is drawn **once per sample and shared across all deprotonation steps**,
because it is the same physical quantity in each step. Sampling it independently per step
would inflate the uncertainty of derived quantities such as pKa differences and species
fractions. For calibrated values, the proton term is replaced by the class residual error so
that the same error is not counted twice. The intervals are conditional on the solvation model.

## 8. Speciation

Species fractions come from the chained Henderson–Hasselbalch expression over pH 0–14. They
are reported at pH 7.0 (physiological) and pH 8.1 (seawater). The R/D tautomeric split is
applied as a constant factor.

## 9. Excited states and electronic structure

```text
! wB97X-D3 def2-TZVPD def2-TZVPD/C CPCM(water) RIJCOSX TightSCF
%tddft NRoots 15  MaxDim 150  ETol 1e-6  RTol 1e-6
       TDA false  Triplets false  DoCD true  DoNTO true end
```

- Ensemble spectra are built by Gaussian broadening of each conformer's stick spectrum, weighted
  by Boltzmann population.
- NTO pairs for the brightest state are selected with a dual criterion: the dominant pair is
  always kept, and any other pair is kept if its coefficient exceeds 15 % of the dominant one
  and an absolute floor of 0.05 (at most two pairs). The participation ratio is reported.
- The electrostatic potential is analysed in Multiwfn on the 0.001 a.u. isodensity surface,
  using the wavefunction converted with `orca_2mkl`. Reported quantities: surface extrema,
  positive and negative areas, average potential, variance and the molecular polarity index.

## 10. Software and hardware

ORCA 6.1.1 · CREST 3.0.2 · xTB 6.7.1 · CENSO 3.0.9 · GoodVibes · Multiwfn · ASE · NumPy ·
pandas · matplotlib. All calculations ran on a Lenovo ThinkPad T14 (Intel Core Ultra 7,
32 GB RAM, 512 GB SSD).

## Key references

- Pracht *et al.*, *J. Chem. Phys.* **160**, 114110 (2024): CREST.
- Neese *et al.*, *J. Chem. Phys.* **152**, 224108 (2020); Neese, *WIREs CMS* **15**, e70019 (2025): ORCA.
- Luchini *et al.*, *F1000Research* **9**, 291 (2020): GoodVibes.
- Lu & Chen, *J. Comput. Chem.* **33**, 580 (2012): Multiwfn.
- Ho & Coote, *Theor. Chem. Acc.* **125**, 3 (2010): continuum pKa calculations.
- Matsuyama *et al.*, *J. Phys. Chem. A* **119**, 12722 (2015): protonation and UV absorption of shinorine and related MAAs.
- Hatakeyama & Nakamura, *ACS Phys. Chem. Au* **5**, 274 (2025): explicit-solvent photodynamics of porphyra-334.
