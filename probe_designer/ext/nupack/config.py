"""Default NUPACK 4 model + tube parameters for HuangLab lab buffer.

These reflect the actual padlock probe hybridization conditions
(2026-05-26, user-confirmed):

* 1× Ampligase buffer = 20 mM Tris-HCl pH 8.3 + 25 mM KCl + 10 mM MgCl₂
* + 50 mM extra KCl → 75 mM total monovalent (K⁺ ≡ Na⁺ in NN models)
* + 20% formamide → destabilizes duplex; modeled as +10 °C effective
  temperature shift (0.5 °C / %FA × 20%, midpoint of literature range)
* Reaction at 45 °C lab temperature → 55 °C effective in no-formamide
  NUPACK DNA model
* Probe concentration 0.05–0.2 μM, default 0.1 μM (midpoint)

NUPACK 4 has no native formamide parameter; we encode formamide
destabilization as an effective ``celsius``.
"""
from __future__ import annotations

from probe_designer.chemistry import LIGASE_CLAMP_NT, ReactionConditions

# Single source of truth for the lab buffer: the shared ReactionConditions
# defaults (protocol rca.md v5.3). These screen constants are DERIVED from it so
# the cross-lig / NUPACK screens and the Tm path never drift apart (audit P5).
_LAB = ReactionConditions()
_NK = _LAB.nupack_kwargs()

# Monovalent salt (K+) — passed as 'sodium' to NUPACK (NN treats Na+/K+ alike).
DEFAULT_SODIUM_M: float = _NK["sodium"]

# Divalent (Mg2+) from Ampligase buffer.
DEFAULT_MAGNESIUM_M: float = _NK["magnesium"]

# Effective temperature in the no-formamide NN model (lab temp + formamide shift).
DEFAULT_CELSIUS: float = _NK["celsius"]

# Informational only — not passed to NUPACK.
DEFAULT_LAB_TEMP_C: float = _LAB.lab_temp_c
DEFAULT_FORMAMIDE_PCT: float = _LAB.formamide_pct
DEFAULT_FORMAMIDE_DEG_PER_PCT: float = _LAB.formamide_deg_per_pct

# NUPACK ensemble — 'stacking' = coaxial stacking on (DEFAULT, the whole point
# of using NUPACK over primer3 for nick-junction energetics).
# 'nostacking' = ablation test; produces less-stable ternary ΔG.
DEFAULT_ENSEMBLE: str = "stacking"

# Probe strand concentration (each of arm3, arm5, splint), from the shared buffer.
# 0.1 μM is the midpoint of the panel-dependent range 0.05–0.2 μM.
DEFAULT_STRAND_CONC_M: float = _LAB.strand_nM * 1e-9

# Productive ternary complex fraction threshold:
# [arm3·splint·arm5] / [splint]_total > threshold flags as cross-lig risk.
# 1e-4 ≈ "1 in 10000 splint molecules in productive ternary at equilibrium"
# — empirically picked to match v2.2 sensitivity on validated cross-lig pairs.
DEFAULT_FRACTION_THRESHOLD: float = 1e-4

# Ligase clamp: n bases each side of the nick, contiguously paired in the MFE
# ternary. Imported, not restated — the register scan, this check and the BLAST
# verdict have to mean the same thing by "ligation-competent", and three private
# copies would let them disagree without anything failing.
DEFAULT_VICINITY_N: int = LIGASE_CLAMP_NT


__all__ = [
    "DEFAULT_SODIUM_M", "DEFAULT_MAGNESIUM_M", "DEFAULT_CELSIUS",
    "DEFAULT_LAB_TEMP_C", "DEFAULT_FORMAMIDE_PCT", "DEFAULT_FORMAMIDE_DEG_PER_PCT",
    "DEFAULT_ENSEMBLE", "DEFAULT_STRAND_CONC_M",
    "DEFAULT_FRACTION_THRESHOLD", "DEFAULT_VICINITY_N",
]
