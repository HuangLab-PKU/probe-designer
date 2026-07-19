# Thermodynamic parameters — values, locations, and sources

Provenance reference for every Tm/ΔG parameter in the designer. When reviewing
or improving the thermodynamics, start here: each value below lists **what it
is**, **where it is set in code**, and **the literature/protocol source** it
comes from. Introduced 2026-07-17 alongside the reaction-buffer refactor
(`probe_designer/chemistry.py`); see the audits
`experiments/20260715_tm_deltag_methods_audit/` (analysis) and
`experiments/20260717_tm_buffer_config_impl/` (implementation + impact).

## 1. Nearest-neighbor (NN) parameter tables — per duplex chemistry

| Chemistry | Table (Biopython) | Set at | Source |
|---|---|---|---|
| DNA probe : mRNA target (direct-mRNA, iLock) | `R_DNA_NN1` | `chemistry.dna_revcomp_to_rna` + `mt.Tm_NN` in `filtering/__init__.py`, `ext/tcr/probe.py`, `rt_primer.py` | **Sugimoto et al. 1995**, *Biochemistry* 34:11211 (RNA/DNA hybrid NN). |
| DNA : DNA (cDNA chemistry) | `DNA_NN4` | same call sites | **SantaLucia & Hicks 2004**, *Annu Rev Biophys* 33:415 (= SantaLucia "unified" 1998, *PNAS* 95:1460; WC values Allawi & SantaLucia 1997, *Biochemistry* 36:10581). |
| RNA : RNA (RNA-probe branch) | `RNA_NN1` | `filtering/__init__.py` | **Freier et al. 1986**, *PNAS* 83:9373. (Newer: Xia et al. 1998 = `RNA_NN2`; candidate upgrade.) |

**Strand convention (critical).** `R_DNA_NN1` is asymmetric; Biopython requires
the **RNA strand** ("Note that `seq` must be the RNA sequence" — installed
`Bio/SeqUtils/MeltingTemp.py:855,885`, v1.85). We therefore feed it the
reverse complement of the DNA arm (`dna_revcomp_to_rna`). Passing the DNA strand
(the pre-2026-07-17 bug in TCR/RT) biases Tm by ~4 °C. Verified against the
audit: `experiments/20260715_tm_deltag_methods_audit/`.

**Modern hybrid alternative (not yet adopted — candidate R6).** Sugimoto 1995
is 1 M-NaCl-derived; at physiological salt it errs ~4.9 °C in Tm.
**Banerjee et al. 2020**, *NAR* 48:12042 (physiological-salt hybrid NN, ~1.1 °C
error) is the modern set; see the audit for the upgrade plan.

## 2. Salt / Mg²⁺ correction

| Parameter | Value | Set at | Source |
|---|---|---|---|
| Salt-correction method | `saltcorr = 5` | `chemistry.ReactionConditions.saltcorr` → `tm_nn_kwargs()` | **SantaLucia 1998** entropic correction (Biopython method 5), with the **von Ahsen et al. 2001** (*Clin Chem* 47:1956) Mg²⁺→Na⁺-equivalent conversion applied automatically when Mg≠0 (Biopython `MeltingTemp.salt_correction:534`). Recommended combination per **von Ahsen, Wittwer & Schütz 2011**, *Brief Bioinform* 12:514. |
| Alternative (high Mg) | `saltcorr = 7` | same | **Owczarzy et al. 2008**, *Biochemistry* 47:5336 (divalent decision tree). Slightly more conservative at ≥8 mM Mg²⁺; user-selectable. |

## 3. Formamide correction

| Parameter | Value | Set at | Source |
|---|---|---|---|
| Model | linear (`chem_correction`, `fmdmethod=1`): `Tm -= factor × %formamide` | `ReactionConditions.apply_formamide` | **McConaughy et al. 1969**, *Biochemistry* 8:3289 (via Biopython `chem_correction`). |
| Coefficient | `formamide_deg_per_pct = 0.5` °C/% | `ReactionConditions` | Chosen to **match the existing cross-lig / NUPACK screen model** (`qc/cross_ligation.py`, `ext/nupack/config.py` both used 0.5 °C/%). Literature range is 0.6–0.72 (Biopython default 0.65); 0.5 is on the conservative end. **Tunable** — revisit if calibration data warrants. |

The same coefficient drives the ΔG-screen effective temperature
(`effective_celsius = lab_temp_c + formamide_pct × formamide_deg_per_pct`), so
the Tm-domain and ΔG-domain treatments stay numerically consistent.

## 4. Reaction buffers

### 4a. Hybridization (padlock arm anneal) — default `ReactionConditions()`

| Parameter | Value | Source |
|---|---|---|
| Monovalent (K⁺) | 75 mM | Protocol `spatial_workshop/docs/protocols/rca.md` v5.3 §2: Ampligase 1× (25 mM KCl, vendor) + 50 mM added KCl. |
| Mg²⁺ | 10 mM | Ampligase 1× buffer (Lucigen/Biosearch), vendor literature. |
| Formamide | 20 % v/v | rca.md v5.3 §2. |
| Probe conc. | 100 nM (0.1 µM) | rca.md v5.3 §2 (range 0.05–0.2 µM, median). |
| Anneal temp | 45 °C | rca.md v5.3 §2 (after 55 °C denaturation). |

*Note:* Ampligase Tris (~20 mM) contributes a minor (~+10 mM Na-equiv) term not
modelled (kept as `monovalent_mM = 75`, K⁺ only). SplintR ligation itself runs at
37 °C in ~10 mM Mg²⁺, no formamide; the 45 °C anneal is the binding-limiting step.

### 4b. Reverse transcription — `rt_primer.MAXIMA_H_MINUS_RT`

| Parameter | Value | Source |
|---|---|---|
| Enzyme / temp | Maxima H Minus RT, **50 °C** | User protocol; Thermo Scientific EP0751/2/3, MAN0012047. |
| Monovalent (K⁺) | 75 mM | 1× of Thermo 5× RT buffer (375 mM KCl → 75 mM). |
| Mg²⁺ | 3 mM | 1× of 5× buffer (15 mM MgCl₂ → 3 mM). |
| Formamide | 0 % | No formamide in reverse transcription. |

Thermo 5× RT buffer: 250 mM Tris-HCl (pH 8.3), 375 mM KCl, 15 mM MgCl₂, 50 mM DTT
(<https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012047_TS_Maxima_H_Minus_Reverse_Transcriptase_2000U_UG.pdf>).
Tris (~50 mM at 1×, ~+25 mM Na-equiv) omitted as a minor term.

### 4c. Cross-ligation / NUPACK screens (unchanged; to be unified — Phase 2)

`qc/cross_ligation.py` and `ext/nupack/config.py` independently encode the same
2026-05-26 lab buffer (75 mM monovalent, 10 mM Mg²⁺, 0.1 µM, 55 °C
formamide-effective). Phase 2 will point these at `ReactionConditions`.

## 5. Tm window anchoring

| Parameter | Value | Set at | Source |
|---|---|---|---|
| Min arm Tm | `lab_temp_c + tm_margin_c` = 45 + 5 = **50 °C** | `ReactionConditions.min_arm_tm`; `FilterConfig.min_tm`; `configs/*.yaml` | **Krzywkowski & Nilsson 2017**, *NAR* 45:e161: padlock arm Tm must exceed the reaction temperature; SNP-padlock lineage uses +5–10 °C. Anchored to the 45 °C hyb anneal. |
| Max arm Tm | `lab_temp_c + 25` = **70 °C** | `FilterConfig.max_tm` | Upper bound for specificity (avoid over-stable arms). |
| RT primer window | **[55, 75]** °C | `rt_primer` defaults; `ext/*/config.py rt_primer_tm_range` | RT temp 50 °C + [5, 25] margin. |
| Concentration term | Tm carries `R·ln(C_T/x)`, x=4 (non-self-comp.) | Biopython `Tm_NN` (`dnac1=dnac2=strand_nM`) | **SantaLucia & Hicks 2004**, Eq. 3. |

**Impact of re-anchoring:** modelling Mg²⁺ (+~13 °C) and 20 % formamide (−10 °C)
nets a **+3 °C** Tm shift; the anchored [50, 70] window fits 99.6 % of shipped
arms (`experiments/20260717_tm_buffer_config_impl/output/impact_summary.txt`).

**Soft Tm scoring (2026-07-19).** The hard [min_tm, max_tm] gate is **off by
default** (`FilterConfig.enforce_tm_gate = False`); absolute arm Tm is instead a
**soft, two-sided** scoring term (`scoring.compute_target_score` `tm_proximity`)
that peaks at the reaction-anchored target (`min_tm` = `lab_temp_c + tm_margin_c`,
≈ 50 °C — a few °C above the 45 °C hyb temp so arms stay bound, per Krzywkowski
2017) and falls off in both directions. Candidates rank by score rather than
being dropped on Tm. `min_tm`/`max_tm` are the scoring target/reference; set
`enforce_tm_gate = True` to restore the hard cutoff. TCR's per-clone selection
(`ext/tcr/probe.filter_by_chemistry_tm`) still hard-gates — separate follow-up.

## 7. Tm as a reference annotation

`filtering/thermo_profile.py` (`compute_tm_profile` / `cached_tm_profile`)
precomputes the per-position arm Tm along a genome/RNA reference under a given
`ReactionConditions` + arm length + chemistry, cached to `.npy` keyed by
`ReactionConditions.signature()`. Design / filtering / visualization can read
this annotation instead of recomputing Tm per candidate. Mirrors the RNAplfold
accessibility cache (`filtering/accessibility.cached_profile`).

**Genome-annotation format (IGV-compatible).** `write_bedgraph()` exports a
profile as **bedGraph** (`chrom  start  end  value`), which IGV/UCSC read
directly — store on NAS beside the reference (one track per condition + arm
length + chemistry, named by the buffer signature). For genome-scale tracks
convert to indexed **BigWig** (UCSC `bedGraphToBigWig` with a chrom.sizes, or
`pyBigWig`). BigWig is the standard for dense per-base numeric tracks like Tm.

**Generating the batch.** `probe_designer/annotate.build_reference_annotations`
and the `probe-design annotate --fasta <ref> --out-dir <NAS>` CLI write the
first-batch tracks for each reference at the current hyb conditions (buffer
flags default to rca.md v5.3): **arm Tm** and **accessibility** (RNAplfold
p_unpaired from full-mRNA folding, at the formamide-effective temperature
`effective_celsius`). Extensible — add more track builders to `annotate.py`.

**Auto-emit from design runs.** `probe-design design --annotations-dir <NAS>`
writes these per-transcript tracks as a byproduct of a run (the annotation
library grows with the genes you design), plus a **canonical-transcript genome
projection** under `<dir>/genome/` for IGV overlay on the assembly.
`genome_projection.py` maps a per-transcript profile to genomic coordinates via
the exon structure + strand and picks the representative transcript (canonical
flag, else longest). Only the canonical transcript is projected — isoforms
disagree at shared exon positions and accessibility is per-isoform; transcript
space stays the lossless source of truth. Windowed metrics are placed at the
window anchor (a junction-spanning arm is annotated at its anchor, not split).

## 6. Tool versions (env `probe-design`)

Biopython **1.85** · primer3-py **2.3.0** · ViennaRNA **2.7.0** · NUPACK **4.0.1.8**.
Biopython changed its default `saltcorr` (5→0) across versions, so it is **set
explicitly** everywhere (never relied on as a default).
