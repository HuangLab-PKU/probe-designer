"""Nearest-neighbor parameter sets, and the one place that picks between them.

Which NN table applies is decided by the chemistry — DNA:RNA hybrid for
direct-mRNA and iLock, DNA:DNA for cDNA, RNA:RNA for an RNA probe — and that
decision was written out separately in four modules (``filtering``,
``filtering.thermo_profile``, ``ext.tcr.probe``, ``rt_primer``). Adding a second
hybrid model to only some of them would have meant a config saying
``banerjee2020`` while the TCR path quietly stayed on Sugimoto, so the mapping
lives here once and every Tm call goes through :func:`melting_temperature`.

**Salt reference — the subtle part.** Biopython's ``Tm_NN`` salt corrections
assume the parameters were derived at 1 M Na+ and correct *downward* to the
buffer you pass. That is right for Sugimoto 1995, SantaLucia 1998 and Freier
1986. It is **wrong for Banerjee 2020**, whose parameters were measured *at*
100 mM NaCl: applying the 1 M correction on top of them corrects twice. Measured
here, that mistake costs **−14.5 °C** on real 20-mer arms, against the ≈ −4.1 °C
the two models genuinely differ by. So a model carries the salt it was measured
at, and :func:`melting_temperature` uses the salt model only for the *relative*
shift from that reference to the actual buffer.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from Bio.SeqUtils import MeltingTemp as mt

# Banerjee et al. 2020, NAR 48:12042, Table 2 — RNA/DNA hybrid NN parameters
# measured at 100 mM NaCl, AS CORRECTED by NAR 2021 49:10796 (the correction
# fixed ΔS values that had been derived from rounded ΔG°37; it changes the
# initiation entropies materially, −4.9 → −6.4 and −7.0 → −8.4).
#
# Notation: the paper writes both strands 5'->3' (rAC/dGT = RNA 5'-AC-3' paired
# with DNA 5'-GT-3'). Biopython keys are RNA 5'->3' (U written as T) over the
# *complement* read 3'->5', so the DNA dimer is reversed: rAC/dGT -> "AC/TG".
# The mapped key set is identical to Biopython's R_DNA_NN1, and every entry
# satisfies ΔG°37 = ΔH° − 310.15·ΔS° to within 0.015 kcal/mol (rounding only) —
# both checked in tests/unit/test_nn_tables.py.
#
# Values are (ΔH° kcal/mol, ΔS° cal/(mol·K)).
R_DNA_BANERJEE2020: Dict[str, Tuple[float, float]] = {
    # Initiation is per terminal base pair (Biopython multiplies these by the
    # number of terminal A/T and G/C pairs). That reading is what reproduces the
    # paper's own claim that Sugimoto errs ≈ 4.9 °C at 100 mM; treating it as a
    # single per-duplex term does not.
    "init": (0.0, 0.0),
    "init_A/T": (0.0, -8.4),      # init. rA-dT / rU-dA
    "init_G/C": (0.0, -6.4),      # init. rG-dC / rC-dG
    "init_oneG/C": (0.0, 0.0),
    "init_allA/T": (0.0, 0.0),
    "init_5T/A": (0.0, 0.0),
    "sym": (0.0, 0.0),
    "AA/TT": (-7.8, -22.9),       # rAA/dTT
    "AC/TG": (-10.1, -27.7),      # rAC/dGT
    "AG/TC": (-9.4, -26.1),       # rAG/dCT
    "AT/TA": (-5.8, -17.4),       # rAU/dAT
    "CA/GT": (-9.8, -27.7),       # rCA/dTG
    "CC/GG": (-9.5, -25.1),       # rCC/dGG
    "CG/GC": (-9.0, -24.5),       # rCG/dCG
    "CT/GA": (-6.1, -18.4),       # rCU/dAG
    "GA/CT": (-8.6, -22.9),       # rGA/dTC
    "GC/CG": (-10.6, -27.7),      # rGC/dGC
    "GG/CC": (-13.3, -35.5),      # rGG/dCC
    "GT/CA": (-9.3, -25.5),       # rGU/dAC
    "TA/AT": (-6.6, -19.7),       # rUA/dTA
    "TC/AG": (-6.5, -16.4),       # rUC/dGA
    "TG/AC": (-8.9, -23.5),       # rUG/dCA
    "TT/AA": (-7.4, -24.5),       # rUU/dAA
}


@dataclass(frozen=True)
class NNModel:
    """A nearest-neighbor parameter set together with the salt it was measured at.

    Attributes:
        name: Stable identifier used in configs and cache keys.
        table: Biopython-format NN dict, ``{dimer: (ΔH, ΔS)}``.
        reference_monovalent_mM: Monovalent salt the parameters were derived at,
            or ``None`` for the 1 M convention Biopython's salt corrections
            assume. A number here means the correction must be applied
            *relative* to it — see :func:`melting_temperature`.
        citation: Where the numbers come from.
    """

    name: str
    table: Dict[str, Tuple[float, float]]
    reference_monovalent_mM: Optional[float]
    citation: str


DNA_DNA = NNModel(
    "santalucia1998", mt.DNA_NN4, None,
    "SantaLucia 1998 PNAS 95:1460 (unified); Allawi & SantaLucia 1997",
)
RNA_RNA = NNModel(
    "freier1986", mt.RNA_NN1, None,
    "Freier et al. 1986 PNAS 83:9373 (Xia 1998 is newer: mt.RNA_NN2)",
)
HYBRID_SUGIMOTO1995 = NNModel(
    "sugimoto1995", mt.R_DNA_NN1, None,
    "Sugimoto et al. 1995 Biochemistry 34:11211 — derived at 1 M NaCl",
)
HYBRID_BANERJEE2020 = NNModel(
    "banerjee2020", R_DNA_BANERJEE2020, 100.0,
    "Banerjee et al. 2020 NAR 48:12042, corrected NAR 2021 49:10796 — 100 mM NaCl",
)

#: Selectable DNA:RNA hybrid models. Sugimoto stays the default: it is the
#: software default everywhere (Biopython, IDT, MELTING), so it keeps our
#: numbers comparable with colleagues'. Banerjee is the more accurate set at
#: physiological salt (Tm error 1.1 °C vs Sugimoto's 4.9 °C at 100 mM) and is
#: opt-in via ReactionConditions.hybrid_nn_model.
HYBRID_MODELS: Dict[str, NNModel] = {
    HYBRID_SUGIMOTO1995.name: HYBRID_SUGIMOTO1995,
    HYBRID_BANERJEE2020.name: HYBRID_BANERJEE2020,
}

DEFAULT_HYBRID_MODEL = HYBRID_SUGIMOTO1995.name

#: Chemistry label -> (probe strand type, target strand type).
CHEMISTRY_STRANDS = {
    "drna": ("DNA", "RNA"),
    "ilock": ("DNA", "RNA"),
    "cdna": ("DNA", "DNA"),
    "rna": ("RNA", "RNA"),
}


def nn_model_for(
    sequence_type: str = "DNA",
    target_type: str = "RNA",
    hybrid_model: str = DEFAULT_HYBRID_MODEL,
) -> NNModel:
    """The NN model for a probe/target strand-type pair.

    Args:
        sequence_type: Probe strand, ``"DNA"`` or ``"RNA"``.
        target_type: Target strand, ``"DNA"`` or ``"RNA"``.
        hybrid_model: Which DNA:RNA set to use; ignored for the non-hybrid
            combinations, which have one canonical table each.

    Raises:
        ValueError: if ``hybrid_model`` is not registered — fail fast rather
            than silently falling back to a different set of physics.
    """
    if sequence_type == "DNA" and target_type == "RNA":
        try:
            return HYBRID_MODELS[hybrid_model]
        except KeyError:
            raise ValueError(
                f"unknown hybrid_nn_model {hybrid_model!r}; "
                f"known: {sorted(HYBRID_MODELS)}"
            ) from None
    if sequence_type == "RNA" and target_type == "RNA":
        return RNA_RNA
    return DNA_DNA


def nn_model_for_chemistry(
    chemistry: str, hybrid_model: str = DEFAULT_HYBRID_MODEL
) -> NNModel:
    """The NN model for a chemistry label (``dRNA`` / ``iLock`` / ``cDNA``)."""
    seq_t, tgt_t = CHEMISTRY_STRANDS.get(str(chemistry).lower(), ("DNA", "RNA"))
    return nn_model_for(seq_t, tgt_t, hybrid_model)


def melting_temperature(seq: str, model: NNModel, reaction) -> float:
    """Tm (°C) of ``seq`` under ``reaction``, honouring the model's salt reference.

    ``seq`` must already be the strand the table expects — for a DNA:RNA hybrid
    that is the RNA-sense strand (see ``chemistry.dna_revcomp_to_rna``).

    Co-solvent depression is **not** applied here; callers apply
    ``reaction.apply_solvents`` so the raw NN Tm stays inspectable.

    For a 1 M-referenced model this is plain ``Tm_NN``. For a model measured at
    a finite salt, the parameters already embody that condition, so the salt
    model is used only for the shift from the reference to the actual buffer:
    the (over-)correction from 1 M is common to both terms and cancels.
    """
    kwargs = reaction.tm_nn_kwargs()
    if model.reference_monovalent_mM is None:
        return mt.Tm_NN(seq, nn_table=model.table, **kwargs)

    at_reference_salt = mt.Tm_NN(
        seq, nn_table=model.table, **{**kwargs, "saltcorr": 0}
    )
    reference_kwargs = {
        **kwargs,
        "Na": 0.0,
        "K": model.reference_monovalent_mM,
        "Tris": 0.0,
        "Mg": 0.0,
        "dNTPs": 0.0,
    }
    shift = (
        mt.Tm_NN(seq, nn_table=model.table, **kwargs)
        - mt.Tm_NN(seq, nn_table=model.table, **reference_kwargs)
    )
    return at_reference_salt + shift


def duplex_melting_temperature(
    seq: str, c_seq: str, model: NNModel, reaction
) -> float:
    """Tm (°C) of a duplex held at a **fixed alignment**, mismatches included.

    :func:`melting_temperature` assumes a perfect duplex and lets the caller's
    single strand imply its partner. The cross-ligation screen needs the other
    thing: two sequences whose pairing register is already fixed by geometry,
    scored as they actually sit — mismatches and all. Letting an aligner choose
    the register instead is exactly the defect audited on 2026-09-02.

    Args:
        seq: one strand, 5'→3'.
        c_seq: the partner base-for-base, ``c_seq[i]`` pairing ``seq[i]`` (i.e.
            the partner read 3'→5'). This is Biopython's own ``c_seq``
            convention — *not* the reverse complement. Must match ``seq`` in
            length; clip both to the paired window before calling.

    Biopython's internal- and terminal-mismatch tables (``DNA_IMM1`` /
    ``DNA_TMM1``, its ``Tm_NN`` defaults) supply the mismatch terms.
    ``strict=False`` lets an unparameterised neighbour drop out of the sum
    rather than raise: a register carrying such a context is mismatch-rich and
    nowhere near the threshold either way, and raising would silently discard a
    candidate register — the very failure mode this function exists to prevent.
    """
    if len(seq) != len(c_seq):
        raise ValueError(
            f"fixed-alignment Tm needs equal lengths; got {len(seq)} vs "
            f"{len(c_seq)} — clip both to the paired window first"
        )
    kwargs = reaction.tm_nn_kwargs()
    if model.reference_monovalent_mM is None:
        return mt.Tm_NN(
            seq, c_seq=c_seq, nn_table=model.table, strict=False, **kwargs
        )

    at_reference_salt = mt.Tm_NN(
        seq, c_seq=c_seq, nn_table=model.table, strict=False,
        **{**kwargs, "saltcorr": 0},
    )
    reference_kwargs = {
        **kwargs,
        "Na": 0.0, "K": model.reference_monovalent_mM,
        "Tris": 0.0, "Mg": 0.0, "dNTPs": 0.0,
    }
    shift = (
        mt.Tm_NN(seq, c_seq=c_seq, nn_table=model.table, strict=False, **kwargs)
        - mt.Tm_NN(
            seq, c_seq=c_seq, nn_table=model.table, strict=False,
            **reference_kwargs,
        )
    )
    return at_reference_salt + shift


__all__ = [
    "CHEMISTRY_STRANDS",
    "DEFAULT_HYBRID_MODEL",
    "DNA_DNA",
    "HYBRID_BANERJEE2020",
    "HYBRID_MODELS",
    "HYBRID_SUGIMOTO1995",
    "NNModel",
    "R_DNA_BANERJEE2020",
    "RNA_RNA",
    "duplex_melting_temperature",
    "melting_temperature",
    "nn_model_for",
    "nn_model_for_chemistry",
]
