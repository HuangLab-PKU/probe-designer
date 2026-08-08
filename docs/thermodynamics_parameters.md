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

## 3. Co-solvent correction (flexible)

Co-solvents depress Tm ~linearly (`Tm -= Σ coeff × %`). Formamide has a dedicated
field; any other solvent goes in `ReactionConditions.solvents` keyed to the
`SOLVENT_TM_COEFF` registry (`register_solvent(name, coeff)` adds one).

| Solvent | Coefficient (°C/%) | Source |
|---|---|---|
| Formamide | `formamide_deg_per_pct = 0.5` | **McConaughy 1969** (via Biopython `chem_correction`). Matches the cross-lig/NUPACK screens (both 0.5); literature 0.6–0.72, Biopython default 0.65 — 0.5 is conservative, **tunable**. |
| DMSO | `SOLVENT_TM_COEFF["dmso"] = 0.6` | **von Ahsen 2001** (typical 0.5–0.75). |
| other (betaine, glycol, …) | user-supplied via `register_solvent` | No built-in NN term — supply a **lab-calibrated** coefficient (approximate). |

`ReactionConditions.solvent_tm_depression` sums these; the same value drives the
ΔG-screen effective temperature (`effective_celsius = lab_temp_c + depression`),
keeping the Tm-domain and ΔG-domain treatments consistent.

**Flexibility & its limit.** Tm-from-buffer is only as flexible as the correction
*physics*: monovalent cations (Na/K/Tris) are **equivalent** for Tm (ionic
strength) — set the total `monovalent_mM` or use `ReactionConditions.from_buffer(
na_mM=, k_mM=, tris_mM=)`; **Mg²⁺ is the only modelled divalent** (approximate a
novel divalent as Mg²⁺-equivalent); solvents beyond formamide/DMSO need an
empirical coefficient. A component with no physical model **fails fast** rather
than silently returning a wrong Tm.

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

### 4c. Cross-ligation / NUPACK / orthogonality screens (unified 2026-07-20)

`qc/cross_ligation.py` and `ext/nupack/config.py` previously each hard-coded the
lab buffer; both now **derive** their constants from the shared
`ReactionConditions` defaults, so the screens and the Tm path cannot drift apart
(values unchanged: 75 mM monovalent, 10 mM Mg²⁺, 0.1 µM, 55 °C formamide-effective).

`filtering/pairwise_duplex.py` (panel orthogonality) previously scored
padlock–padlock **DNA:DNA** dimers with ViennaRNA's **RNA** Turner parameters;
it now uses **primer3** (SantaLucia unified DNA NN) at the same buffer. ΔG
magnitudes are therefore not comparable to pre-2026-07-20 values; the −12 kcal/mol
default is a permissive log-only threshold, re-tune before using it to drop probes.

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
**soft, two-sided** scoring term that peaks at the reaction-anchored target
(`min_tm` = `lab_temp_c + tm_margin_c`, ≈ 50 °C — a few °C above the 45 °C hyb
temp so arms stay bound, per Krzywkowski 2017) and falls off in both directions.
Candidates rank by score rather than being dropped on Tm. `min_tm`/`max_tm` are
the scoring target/reference; set `enforce_tm_gate = True` to restore the hard
cutoff. TCR's per-clone selection also went soft in the same change.

### 5a. One rule, one definition — `policy.ThermoPolicy` (2026-07-21)

The arm-Tm rules used to exist **twice**: `filtering/` held the thresholds it
rejects on, `scoring/` held its own copies to rank with — and they had drifted
(`compute_target_score` defaulted to `min_arm_tm=50 / max_tm_diff=10` regardless
of `FilterConfig`, whose own `min_tm` default of 45 contradicted every shipped
YAML's 50). That is audit **R7**, and its cause is duplication rather than any
one wrong number.

`ThermoPolicy` now owns both the thresholds **and** their normalized 0–1 scoring
shapes; `scoring` keeps only the composite *weights* (how much Tm proximity is
worth against isoform coverage is a ranking decision, not a thermodynamic one),
and `filtering` keeps its per-call override surface. Build one with
`ThermoPolicy.resolve(filter_config)` or `ThermoPolicy.from_reaction(reaction)`
— the latter re-anchors the target when you change `lab_temp_c`, so raising the
hyb temperature moves what counts as a good arm without editing a threshold.

| Parameter | Value | Set at | Source |
|---|---|---|---|
| Two-arm Tm balance | **5 °C** | `ThermoPolicy.max_tm_diff`; `FilterConfig.max_tm_diff`; `MutationConfig.max_tm_diff` | **Larsson et al. 2010** *Nat Methods* 7:395; **Krzywkowski & Nilsson 2017** — the two arms must be comparable; field practice ≤ 3–5 °C. |
| Tm proximity falloff | **20 °C** | `ThermoPolicy.tm_falloff_c` | Deliberately wide: hybrid-Tm model error is a few °C (Banerjee 2020 reports 1.1 °C at best), so a narrow falloff would rank on noise. |

**Impact of R4 (balance 10 → 5 °C).** Measured before changing anything, on 232
shipped probes at the corrected buffer: as a **hard gate**, 5 °C rejects **17.7 %**
(19.7 % of cDNA, 3.4 % of dRNA); 10 °C rejects 0.9 %. So the tolerance tightened
as a **scoring reference** and the gate became soft
(`enforce_tm_diff_gate = False`, matching `enforce_tm_gate`) — ranking gains the
discrimination field practice asks for, and no probe is discarded. Set
`enforce_tm_diff_gate = True` for a hard cutoff.
(`experiments/20260717_tm_buffer_config_impl/output/phase3_summary.txt`)

### 5b. Target accessibility — `chemistry.FoldingConditions` (2026-07-21)

Accessibility asks whether a candidate window is *open* in the folded transcript
(local RNAplfold `p_unpaired`, Lange 2012), which beats a per-window self-fold
MFE and puts us ahead of MERFISH/Xenium (they filter only *probe* self-structure).
The geometry that question is asked over used to be three loose numbers repeated
across `FilterConfig`, `search_strategies`, `annotate`, `accessibility` and two
YAMLs — the same duplication the buffer had before `ReactionConditions`, with the
same consequence: `build_canonical_genome_annotations` silently ignored the
caller's geometry, and the filter folded at a nominal 37 °C while the annotation
writer folded at the solvent-effective temperature, so the two disagreed about
the same transcript.

| Parameter | Value | Set at | Source |
|---|---|---|---|
| Span `L` | **100 nt** | `FoldingConditions.span`; `FilterConfig.plfold_span` | **Lange et al. 2012**, *NAR* 40:5215 — the ViennaRNA `W = L = 70` default is artifact-prone; recommended `L ≈ 100`. |
| Window `W` | **150 nt** (`= L + 50`) | `FoldingConditions.window`; `FilterConfig.plfold_window` | Lange 2012, `W = L + 50`. |
| Folding temperature | **tracks the reaction** (`effective_celsius` ≈ 55 °C) | `FoldingConditions.temperature_c = None`; `FilterConfig.plfold_temperature` | Folding must happen in the same buffer the arm anneals in; `None` means "derive", so the accessibility that *filters* is the accessibility that gets *annotated*. Pin a number to decouple. |
| Aggregation | mean `p_unpaired` over the arm footprint, expression-weighted across isoforms | `accessibility.mean_open_probability` / `aggregate_open_probability` | Lange 2012 (average over the footprint). |

**Impact of R9 (`W=70,L=40` → `W=150,L=100`).** On 5 real transcripts the old
geometry **over-stated** mean `p_unpaired` by **0.10** (0.64 → 0.54), with
|Δ| > 0.1 at **half** of all bases — a span shorter than the real base-pair
distance cannot represent a long-range stem, so a site locked behind one reads
as open. Gate impact is negligible (**0.1 %** of 40-nt windows cross a 0.30
accessibility gate), and the gate is off by default (`min_accessibility = 0`).
`FoldingConditions.signature()` keys the geometry as well as the temperature, so
profiles computed under different `W`/`L` cannot collide in a cache or an
annotation directory.
(`experiments/20260717_tm_buffer_config_impl/output/phase3_summary.txt`)

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

**Declaration (R10, 2026-07-21).** `primer3-py` backs the cross-ligation screen
and the panel orthogonality screen — both production paths that import it
unguarded — but was declared in neither `requirements.txt` nor `pyproject.toml`,
so a fresh install only worked if the environment happened to provide it. It is
now a hard dependency (`primer3-py>=2.0.0`) in both. `tests/unit/test_packaging_deps.py`
lints the *class* of bug (imported-but-undeclared) across the thermodynamic
libraries rather than re-asserting this one case. NUPACK stays undeclared on
purpose: it needs academic registration and a manual install, and every call
site is gated behind `has_nupack()`.
