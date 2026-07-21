"""Reaction-condition (buffer) model shared by every thermodynamic call.

Historically the designer computed all melting temperatures with Biopython's
``Tm_NN`` defaults (Na=50 mM, Mg=0, 50 nM oligo), which does not match the lab
reaction buffer and — because Mg2+ was ignored — understated arm Tm by ~10 C
(audit: ``experiments/20260715_tm_deltag_methods_audit/``). This module makes
the buffer a first-class, user-passable object so every Tm/ΔG site computes at
the real hybridization conditions, and so the (previously duplicated) cross-lig
and NUPACK buffers share a single source of truth.

Defaults follow the current canonical protocol (spatial-workshop
``docs/protocols/rca.md`` v5.3): Ampligase 1x (25 mM K+, 10 mM Mg2+) + 50 mM
added KCl + 20% formamide, annealed at 45 C.

Flexibility & its limits: Tm-from-buffer is only as flexible as the correction
*physics*. Monovalent cations (Na/K/Tris) are equivalent for Tm (ionic strength)
— pass the total via ``monovalent_mM`` or the granular ``from_buffer`` helper.
Mg2+ is the only modelled divalent; approximate a novel divalent as Mg2+-equiv.
Co-solvents depress Tm linearly: formamide has a dedicated field, any other
solvent goes in ``solvents`` keyed to ``SOLVENT_TM_COEFF`` (``register_solvent``
adds an empirical coefficient for a denaturant with no built-in NN term). The
same solvent depression raises the effective screen temperature, keeping the
Tm-domain and ΔG-domain treatments consistent.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

from Bio.SeqUtils import MeltingTemp as mt

# Valid Biopython MeltingTemp salt-correction methods (see salt_correction).
_VALID_SALTCORR = frozenset({1, 2, 3, 4, 5, 6, 7})

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
               "U": "A",   # accept RNA input
               "-": "-"}   # preserve gaps in aligned sequences

# Tm depression coefficient (degC per % v/v) for co-solvents / denaturants.
# Formamide keeps a dedicated field (historical + default); ANY other solvent
# goes in ReactionConditions.solvents keyed here. This is the extension point
# for a denaturant Biopython doesn't build in (Tm-from-buffer is only as flexible
# as the correction physics — beyond salt/formamide/DMSO there is no NN model, so
# a novel solvent needs an empirical coefficient supplied here). Values are
# midpoints of the literature range; override per lab as needed.
SOLVENT_TM_COEFF: Dict[str, float] = {
    "formamide": 0.5,   # McConaughy 1969 (repo-consistent; literature 0.6-0.72)
    "dmso": 0.6,        # von Ahsen 2001 (typical 0.5-0.75)
}


def register_solvent(name: str, deg_per_pct: float) -> None:
    """Register/override a co-solvent's Tm depression coefficient (degC per %).

    Model a denaturant not built into the NN corrections (e.g. betaine,
    ethylene glycol) by supplying an empirical coefficient, then list it in a
    ReactionConditions ``solvents`` mapping. This is an approximation — there is
    no first-principles NN term for it — so use a lab-calibrated value.
    """
    SOLVENT_TM_COEFF[name.strip().lower()] = float(deg_per_pct)


def reverse_complement(seq: str) -> str:
    """Reverse complement, uppercase. The canonical implementation for this package.

    ``U`` reads as ``T`` (RNA input is accepted), a gap ``-`` is preserved so
    aligned sequences survive the round trip, and anything else maps to ``N``
    rather than raising — callers here feed validated ACGT, and a permissive
    mapping keeps a batch of thousands from dying on one stray character.

    This replaced seven near-identical private copies scattered across the
    package, which had quietly diverged on three axes: whether the input was
    uppercased first, whether an unknown base raised ``KeyError`` or became
    ``N``, and whether ``-`` and ``U`` were handled at all. Divergence on the
    case axis in particular meant the same call was safe in one module and a
    crash in another.

    Args:
        seq: nucleotide sequence, 5'->3' (case-insensitive).

    Returns:
        The reverse complement, uppercase.
    """
    return "".join(_COMPLEMENT.get(base, "N") for base in reversed(seq.upper()))


def dna_revcomp_to_rna(seq: str) -> str:
    """The reverse complement of a DNA arm, named for its RNA-sense intent.

    For a DNA probe arm hybridizing to mRNA, the reverse complement *is* the
    target (mRNA) sequence. Biopython's ``R_DNA_NN1`` table requires the RNA
    strand ("Note that ``seq`` must be the RNA sequence"), so this is the string
    that must be fed to ``Tm_NN`` for a DNA:RNA hybrid. Thymine is intentionally
    kept (the table is keyed on T), so "RNA-sense" here does not mean
    U-containing.

    Mechanically identical to :func:`reverse_complement`; kept as a separate
    name because at the call sites the intent — "make the Tm table's required
    strand" — is what needs to be readable.
    """
    return reverse_complement(seq)


@dataclass
class ReactionConditions:
    """Ionic + additive composition of a hybridization/ligation reaction.

    All concentrations are per the buffer the probe arm actually anneals in.
    Instances are cheap, immutable-in-spirit value objects; construct one per
    chemistry (dRNA hyb differs from cDNA) and pass it into the Tm/ΔG engines
    via the ``*_kwargs`` helpers.

    Attributes:
        monovalent_mM: Monovalent cation (K+) concentration (mM).
        mg_mM: Mg2+ concentration (mM).
        dntp_mM: dNTP concentration (mM); chelates Mg2+ when present.
        strand_nM: Oligo/strand concentration (nM).
        formamide_pct: Formamide (% v/v); destabilizes duplexes.
        formamide_deg_per_pct: Tm depression per % formamide (C/%). Also the
            effective-temperature offset coefficient for the ΔG screens.
        solvents: Extra co-solvents beyond formamide, ``{name: % v/v}`` — each
            depresses Tm by ``SOLVENT_TM_COEFF[name] * %`` (register_solvent to
            add a coefficient). Novel divalent ions have no NN model; approximate
            as Mg2+-equivalent via ``mg_mM``. Monovalents (Na/K/Tris) are
            equivalent for Tm — pass the total via ``monovalent_mM`` or the
            granular ``from_buffer(na_mM=, k_mM=, tris_mM=)`` helper.
        lab_temp_c: Hybridization anneal temperature (C); the arm-Tm anchor.
        tm_margin_c: Required arm-Tm margin above ``lab_temp_c`` (Krzywkowski 2017).
        saltcorr: Biopython MeltingTemp salt-correction method (1-7). 5 =
            SantaLucia 1998 entropic (+ von Ahsen Mg-equivalent); 7 = Owczarzy
            2008 divalent decision tree (more accurate at high Mg2+).
        hybrid_nn_model: Which DNA:RNA nearest-neighbor set to use —
            ``"sugimoto1995"`` (default, the software default everywhere, so our
            numbers stay comparable with colleagues') or ``"banerjee2020"`` (more
            accurate at physiological salt: 1.1 C Tm error vs 4.9 C). Affects
            only the hybrid chemistries; cDNA (DNA:DNA) and RNA:RNA are
            unchanged. See nn_tables.
    """

    monovalent_mM: float = 75.0
    mg_mM: float = 10.0
    dntp_mM: float = 0.0
    strand_nM: float = 100.0
    formamide_pct: float = 20.0
    formamide_deg_per_pct: float = 0.5
    solvents: Dict[str, float] = field(default_factory=dict)
    lab_temp_c: float = 45.0
    tm_margin_c: float = 5.0
    saltcorr: int = 5
    hybrid_nn_model: str = "sugimoto1995"

    @classmethod
    def from_buffer(
        cls, *, na_mM: float = 0.0, k_mM: float = 0.0, tris_mM: float = 0.0, **kwargs
    ) -> "ReactionConditions":
        """Build from granular monovalent cations (Na/K/Tris).

        Monovalents are equivalent for Tm (the salt correction sees only the
        total ionic strength), so they are folded into ``monovalent_mM`` as
        ``Na + K + Tris/2`` (von Ahsen convention). All other fields pass through.
        """
        return cls(monovalent_mM=na_mM + k_mM + tris_mM / 2.0, **kwargs)

    def __post_init__(self) -> None:
        """Validate the composition, failing fast with an actionable message."""
        for name in ("monovalent_mM", "mg_mM", "dntp_mM", "formamide_pct"):
            value = getattr(self, name)
            if value < 0:
                raise ValueError(f"{name} must be >= 0, got {value}")
        if self.strand_nM <= 0:
            raise ValueError(f"strand_nM must be > 0, got {self.strand_nM}")
        if self.monovalent_mM <= 0:
            raise ValueError(
                "monovalent_mM must be > 0 (salt correction takes ln[Mon+]); "
                f"got {self.monovalent_mM}"
            )
        if self.saltcorr not in _VALID_SALTCORR:
            raise ValueError(
                f"saltcorr must be one of {sorted(_VALID_SALTCORR)}, got {self.saltcorr}"
            )
        # Imported here: nn_tables imports nothing from this module, but keeping
        # the dependency one-directional at module scope would still read oddly.
        from probe_designer.nn_tables import HYBRID_MODELS

        if self.hybrid_nn_model not in HYBRID_MODELS:
            raise ValueError(
                f"hybrid_nn_model must be one of {sorted(HYBRID_MODELS)}, "
                f"got {self.hybrid_nn_model!r}"
            )
        for name, pct in self.solvents.items():
            if pct < 0:
                raise ValueError(f"solvent {name!r} percent must be >= 0, got {pct}")

    @property
    def min_arm_tm(self) -> float:
        """Minimum acceptable arm Tm: anneal temperature plus the margin (C)."""
        return self.lab_temp_c + self.tm_margin_c

    @property
    def solvent_tm_depression(self) -> float:
        """Total Tm depression from all co-solvents (C): formamide + solvents.

        Raises if a listed solvent has no coefficient (register_solvent first) —
        Tm-from-buffer is only as flexible as the correction physics.
        """
        dep = self.formamide_pct * self.formamide_deg_per_pct
        for name, pct in self.solvents.items():
            key = name.strip().lower()
            if key not in SOLVENT_TM_COEFF:
                raise ValueError(
                    f"unknown solvent {name!r}; call register_solvent(name, "
                    "deg_per_pct) with an empirical coefficient first"
                )
            dep += SOLVENT_TM_COEFF[key] * pct
        return dep

    @property
    def effective_celsius(self) -> float:
        """Solvent-shifted temperature for the ΔG/complex screens (C).

        Co-solvents are modelled as a temperature increase for the two-strand
        screens (raising stringency), equal to the Tm depression applied to arm
        Tm — both use ``solvent_tm_depression``.
        """
        return self.lab_temp_c + self.solvent_tm_depression

    def tm_nn_kwargs(self) -> dict:
        """Keyword arguments for ``Bio.SeqUtils.MeltingTemp.Tm_NN``.

        Monovalent salt is passed as K+ (the protocol is KCl-based); Biopython
        sums Na+K+Tris and folds Mg2+ in via the von Ahsen equivalent for
        ``saltcorr`` 5/6, or the Owczarzy 2008 tree for 7.
        """
        return {
            "Na": 0.0,
            "K": self.monovalent_mM,
            "Tris": 0.0,
            "Mg": self.mg_mM,
            "dNTPs": self.dntp_mM,
            "dnac1": self.strand_nM,
            "dnac2": self.strand_nM,
            "saltcorr": self.saltcorr,
        }

    def apply_solvents(self, tm: float) -> float:
        """Depress a computed Tm for all co-solvents (formamide + ``solvents``).

        Linear per solvent (McConaughy ``fmdmethod=1``): ``Tm -= sum(coeff * %)``,
        matching Biopython ``chem_correction`` for the formamide term.
        """
        return tm - self.solvent_tm_depression

    # Backward-compat alias (pre-2026-07-20 name; formamide was the only solvent).
    apply_formamide = apply_solvents

    def primer3_kwargs(self) -> dict:
        """Buffer kwargs for primer3 ``calc_heterodimer`` / ``calc_tm``.

        primer3 expects monovalent (mM), divalent (mM), dNTP (mM), DNA conc (nM),
        and a simulation temperature (C). The formamide-effective temperature is
        used so the ΔG screen stringency matches the reaction.
        """
        return {
            "mv_conc": self.monovalent_mM,
            "dv_conc": self.mg_mM,
            "dntp_conc": self.dntp_mM,
            "dna_conc": self.strand_nM,
            "temp_c": self.effective_celsius,
        }

    def nupack_kwargs(self) -> dict:
        """Buffer kwargs for a NUPACK 4 ``Model`` (molar salts, Celsius)."""
        return {
            "sodium": self.monovalent_mM / 1000.0,
            "magnesium": self.mg_mM / 1000.0,
            "celsius": self.effective_celsius,
        }

    def signature(self) -> str:
        """Compact, filesystem-safe key of the parameters that affect Tm.

        Used to key precomputed per-position Tm profiles / reference annotations
        (see filtering.thermo_profile). Excludes ``lab_temp_c`` / ``tm_margin_c``
        (they set the gate/target, not the Tm value itself).
        """
        sv = "".join(
            f"-{n}{p:g}" for n, p in sorted(self.solvents.items())
        )
        return (
            f"mv{self.monovalent_mM:g}_mg{self.mg_mM:g}_dn{self.dntp_mM:g}"
            f"_st{self.strand_nM:g}_fa{self.formamide_pct:g}"
            f"_fd{self.formamide_deg_per_pct:g}_sc{self.saltcorr}{sv}"
        )


# ViennaRNA's own convention when no reaction context is available. Real runs
# pass a ReactionConditions and fold at its solvent-effective temperature.
DEFAULT_FOLD_TEMPERATURE_C = 37.0


@dataclass(frozen=True)
class FoldingConditions:
    """RNAplfold window geometry + temperature for target-accessibility profiles.

    The companion of :class:`ReactionConditions` for the structure side: where
    that object answers "what buffer does the arm anneal in", this one answers
    "over how much sequence context, and at what temperature, do we ask whether
    the site is open".

    Geometry defaults follow Lange 2012 (NAR 40:5215), which shows ViennaRNA's
    W = L = 70 default is artifact-prone for accessibility: a span shorter than
    the real base-pair distance cannot represent a long-range stem, so a site
    locked behind one reads as open. The recommendation is L ~ 100 with
    W = L + 50. (Pre-2026-07-21 this package used W = 70, L = 40 — shorter still.)

    Attributes:
        window: RNAplfold ``W``, the sliding window size (nt).
        span: RNAplfold ``L``, the maximum base-pair span within the window (nt).
            Must not exceed ``window``.
        temperature_c: Folding temperature (C). ``None`` (default) means *track
            the reaction* — resolved to :attr:`ReactionConditions.effective_celsius`
            so the accessibility used to filter candidates is the same
            accessibility written into the reference annotation. Set a number to
            pin it independently of the buffer.
    """

    window: int = 150
    span: int = 100
    temperature_c: Optional[float] = None

    def __post_init__(self) -> None:
        """Validate the geometry, failing fast with an actionable message."""
        if self.window <= 0:
            raise ValueError(f"window (W) must be > 0, got {self.window}")
        if self.span <= 0:
            raise ValueError(f"span (L) must be > 0, got {self.span}")
        if self.span > self.window:
            raise ValueError(
                f"span (L={self.span}) must not exceed window (W={self.window}); "
                "RNAplfold cannot pair further apart than the window it sees"
            )

    def temperature_for(self, reaction: Optional["ReactionConditions"] = None) -> float:
        """Resolve the folding temperature against a reaction (C).

        An explicit ``temperature_c`` always wins; otherwise the reaction's
        solvent-effective temperature is used, falling back to the bare
        ViennaRNA convention when no reaction is supplied.
        """
        if self.temperature_c is not None:
            return float(self.temperature_c)
        if reaction is not None:
            return reaction.effective_celsius
        return DEFAULT_FOLD_TEMPERATURE_C

    def plfold_kwargs(self, reaction: Optional["ReactionConditions"] = None) -> dict:
        """Keyword arguments for ``filtering.accessibility.compute_plfold_profile``.

        Mirrors :meth:`ReactionConditions.tm_nn_kwargs` — splat it into the
        primitive so callers never hand-carry the geometry triple again.
        """
        return {
            "window": self.window,
            "span": self.span,
            "temperature": self.temperature_for(reaction),
        }

    def signature(self, reaction: Optional["ReactionConditions"] = None) -> str:
        """Compact, filesystem-safe key of everything that changes the profile.

        The geometry belongs in the key, not just the temperature: widening the
        span shifts mean p_unpaired by ~0.10 on real transcripts, so two runs
        with different W/L must not collide in a cache or an annotation
        directory. Mirrors :meth:`ReactionConditions.signature`.
        """
        return f"W{self.window}_L{self.span}_T{self.temperature_for(reaction):g}"
