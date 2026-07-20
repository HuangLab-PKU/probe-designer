"""Optional NUPACK 4 extension for cross-ligation ternary-complex confirmation.

This package adds a Tier 3 ensemble-aware ligation-competence check on top
of qc/cross_ligation.py's v2.2 split-fragment screen. NUPACK 4 brings:

* **Coaxial stacking** (Model ensemble='stacking' default) — the ~−1 to
  −2 kcal/mol energy that stabilizes a nick between adjacent oligos on
  the splint. primer3 binary heterodimer scoring misses this entirely.
* **Multi-strand tube_analysis** — the equilibrium concentration of the
  productive (arm3·splint·arm5) ternary complex at lab conditions.
* **Ensemble sampling** — secondary alignments contribute proportional
  to their Boltzmann weight, not just the MFE.

NUPACK 4 is an OPTIONAL dependency (academic registration + download at
nupack.org). All ``import nupack`` calls are lazy (inside function bodies);
``from probe_designer.ext.nupack.check import screen_ternary`` succeeds
without nupack installed — only the actual call raises a clear ImportError.

See ``[[project-cross-lig-v2-3]]`` for the full v2.3 design.
"""

__all__ = ["check", "config", "cli"]
