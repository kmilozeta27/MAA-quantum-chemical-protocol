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
censo --prescreening --screening --optimization --refinement --ignore-failed \
      --inprc censo2rc_refine -i <ensemble> -c <q> -u 0 \
      --solvent water --maxcores <N> --omp-min 4 -T <temperature>
```

CENSO 3 runs four stages (configuration in [`censo2rc_refine`](../censo2rc_refine);
`maaqc.py run --censo-refine`):

| Stage | Level | Cut-off |
|---|---|---|
| Pre-screening | B97-D3/def2-SV(P), gas phase, + ALPB(GFN2) Gsolv | 4.0 kcal/mol |
| Screening | r²SCAN-3c, CPCM (`gsolv_included = True`) | 3.5 kcal/mol, adjusted automatically |
| Optimization | r²SCAN-3c, CPCM, ANCOPT macrocycles | 3.0 kcal/mol |
| Refinement | ωB97X-D3/def2-TZVP, SMD | 95 % cumulative Boltzmann population |

Conformers go on to DFT from the r²SCAN-3c geometries, according to their refined Boltzmann
weight (threshold 5 %, at least three conformers), never from the CREST geometries. If the
refinement is missing, the pipeline stops instead of falling back to screening energies. A run
that already finished screening or optimization is resumed from that stage. The first version
of the protocol (pre-screening and screening only, selection from single points on GFN2-xTB
geometries) is still the default without `--censo-refine`. CENSO does not remove duplicates
after optimization, so refined ensembles must be checked for them. A multi-temperature mode
runs a separate CENSO ensemble for each temperature. If CENSO is missing or fails, the pipeline
stops with an error instead of silently continuing with the unfiltered CREST ensemble; the old
behaviour requires the explicit `--allow-no-censo` flag.

ORCA inside CENSO needs `OMPI_MCA_pml=ob1`, `OMPI_MCA_btl=self,sm`, `OMPI_MCA_osc=^ucx` and
`PRTE_MCA_hwloc_default_binding_policy=none` (set by `maaqc.py`). Without them, Open MPI 5.0.10
fills the stderr pipe that CENSO does not read, and ORCA hangs at 0 % CPU; it also binds
every parallel job to the same four cores.

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

**Small imaginary modes.** GoodVibes by default removes imaginary modes from the vibrational
sums. Here GoodVibes is called with its mode-inversion option (`--invert` in v4,
`--invertifreq` in v3) at 20 cm⁻¹: imaginary modes with |ν| < 20 cm⁻¹ are treated as real
modes of equal |ν| and then go through the quasi-RRHO treatment (Grimme interpolation,
100 cm⁻¹ cutoff). A mode more negative than −20 cm⁻¹ is not inverted (GoodVibes removes it)
and flags a structure to re-optimise. If the installed GoodVibes has no inversion option, a
warning is logged and the default (removal) applies. Structures are not re-optimised for small
imaginary modes.

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

| Class | Charge change | Group | Calibrators (n) | Raw error, pKa units (v4, SMD/direct, GoodVibes) |
|---|---|---|---|---|
| `carboxyl` | +1 → 0 | first COOH | Gly, Val, Asp(α), Glu(α) (4) | +1.03 to +1.51 |
| `carboxyl_2` | 0 → −1 | second COOH | Asp(β), Glu(γ) (2) | +2.06, +1.40 |
| `amino` | 0 → −1 | NH₃⁺, no neighbouring charge | Gly, Val (2) | +1.81, +3.15 |
| `amino_dianion` | −1 → −2 | NH₃⁺ flanked by two COO⁻ | Glu (1) | +6.01 |
| `pyridinium` | +1 → 0 | aromatic N–H⁺ | 4-aminopyridine (1) | +1.64 |

Values are the errors of `calibration/calibration_v4.json` (pKa experimental values from IUPAC where a Reliable entry exists, §11). The earlier Born/CPCM
calibration put the `amino_dianion` error near 8 units. A class with n = 1 (`amino_dianion`,
`pyridinium`) has undefined RMSE and no cross-validation.

**Known mismatch.** The B → C step of mycosporine-glycine:valine deprotonates a vinylogous
amidinium, not an sp³ NH₃⁺, so `amino_dianion` (Glu) is an imperfect calibrator for it.

**Calibration source.** `--empirical-correction` requires a current calibration file: pass
`--calibration-file` or provide `calibration/calibration_v4.json`. If neither is found the
program stops with an error. The table embedded in `maaqc.py` is obsolete and is used only with
the explicit `--allow-builtin-calibration` flag.

- Each class gets its own least-squares fit. A free slope is used only if the class has at
  least five calibrators (n ≥ 5), a calculated-pKa range of at least 3 units, and
  AICc(linear) < AICc(offset) − 2 (ΔAICc > 2 in favour of the linear model). Otherwise the fit
  is an offset per class. The frozen rule is in §11.
- RMSE and leave-one-out RMSE are reported per class. A class with a single calibrator is
  marked as having undefined RMSE and is not applied without an explicit override.
- **The pipeline refuses to apply a calibrator across charge classes.** For glutamic acid's
  −1 → −2 step, the raw cycle deviates from experiment by about 6 pKa units with SMD and the
  direct cycle (about 8 with the earlier Born/CPCM cycle). This class cannot be transferred from
  the `amino` class.
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

**Wave-function reference (status).** STEOM-DLPNO-CCSD has been run on the full
mycosporine-glycine:valine B molecule with def2-SVP, def2-SVPD, ma-def2-SVP and CPCM(water),
giving 3.41–3.62 eV for S1 (ωB97X-D3: 4.20 eV; experiment: 3.70 eV, 335 nm). The def2-TZVP
increment and the bt-PNO-EOM-CCSD reference are not available (both attempts failed), and the
ADC(2) run was cut by a timeout and is not used. Extractors accept an ORCA output only after
checking for normal termination, and for STEOM they read the STEOM spectrum rather than the
preceding CIS one. Ensemble conformer folders are identified by their numeric id
(`conf_5` ≡ `conf_005`); folders with a suffix (e.g. `conf_6_PBE0`) are auxiliary calculations
and are excluded from the ensemble.

## 10. Software and hardware

ORCA 6.1.1 · CREST 3.0.2 · xTB 6.7.1 · CENSO 3.0.9 · GoodVibes · Multiwfn · ASE · NumPy ·
pandas · matplotlib. All calculations ran on a single laptop, 32 GB RAM.

## 11. Frozen specification pka-spec-v1.1

**tag: pka-spec-v1.1** (supersedes `pka-spec-v1`). This section fixes the primary protocol
before any MGLYVAL pKa is examined. Rationale and source numbers:
the calibration plan (hereafter CP, an internal document), §3 and §4. Calibrators and MAA use
exactly the same specification.

**Primary protocol.**

- Solvation: **SMD(water)** is the primary model (ωB97X-D3/def2-TZVPD SP on the geometry
  optimised with CPCM). **CPCM(water)** is reported only as a sensitivity analysis (CP §4).
- Thermochemistry: **qh-GoodVibes** corrections (Grimme qRRHO, 100 cm⁻¹ cutoff, harmonic H),
  298.15 K. Calibration is done at 298.15 K only. Pure harmonic is a sensitivity analysis.
- Cycle: direct method, G(H⁺,aq) = −270.286 kcal/mol (§5).
- Sampling and DFT: as in §3 and §4 (CREST → CENSO, ωB97X-D3/def2-SVP CPCM, no imaginary
  frequencies with |ν| > 20 cm⁻¹).
- Macrostates: G_macro = −RT ln Σᵢ Σ_conf exp(−G/RT) over all microstates of the same charge,
  with degeneracy gᵢ (CP §4).

**Calibration model: offset per class.** The default mode is offset only. A free slope is
allowed only with n ≥ 5 calibrators, a calculated-pKa range ≥ 3 units, and
AICc(linear) < AICc(offset) − 2 (exact criterion as stated in CP §3.1). RMSE_LOO is reported
in all modes. A calibrator is never applied to another charge class (cross-class protection,
§6). Guanidinium groups are excluded from the fit (CP §5, D7).

**Macro-pKa** (CP §1; each macrostate groups all microstates of the same charge):

- **pKaM1** (+1 → 0): A → {D, D′, R, R′}. First COOH with the amidinium present.
- **pKaM2** (0 → −1): {D, D′, …} → {B, 010, 100}. Second COOH.
- **pKaM3** (−1 → −2): {B, …} → {C, C′}. Amidinium N–H with 2×COO⁻.

Micro-pKa values (A→D, A→D′, D→B, D′→B, A→R) are diagnostic only. The primary pKaM3 value
(double difference or dianion class alone) follows decision D1′ of CP §5.

**Acceptance criteria** (CP §3.4, unmodified):

- RMSE_LOO ≤ 1.0 per class (≤ 0.5 for `carboxyl`);
- |mean external error| ≤ 1.0 (external validation H1–H4 and porphyra-334, outside the fit);
- micro-pKa cycle closure ≤ 0.5;
- additivity: |[δ(pyridinium, 0→−1) − δ(pyridinium, +1→0)] − [δ(ammonium, 0→−1) − δ(ammonium, +1→0)]| ≤ 1.0;
  if it fails, the double difference is kept as a sensitivity analysis only.

**Change rule.** Any change after this specification (solvation model, thermal corrections,
fit rule, definition of sites or macro-pKa, acceptance criteria) requires a new tag
(`pka-spec-v2`, …) and an explicit record. Other combinations are reported as sensitivity
analyses, never as a silent replacement of the primary one.

## Key references

- Pracht *et al.*, *J. Chem. Phys.* **160**, 114110 (2024): CREST.
- Neese *et al.*, *J. Chem. Phys.* **152**, 224108 (2020); Neese, *WIREs CMS* **15**, e70019 (2025): ORCA.
- Luchini *et al.*, *F1000Research* **9**, 291 (2020): GoodVibes.
- Lu & Chen, *J. Comput. Chem.* **33**, 580 (2012): Multiwfn.
- Ho & Coote, *Theor. Chem. Acc.* **125**, 3 (2010): continuum pKa calculations.
- Matsuyama *et al.*, *J. Phys. Chem. A* **119**, 12722 (2015): protonation and UV absorption of shinorine and related MAAs.
- Hatakeyama & Nakamura, *ACS Phys. Chem. Au* **5**, 274 (2025): explicit-solvent photodynamics of porphyra-334.
