"""Configuration dataclasses for the mutation pipeline.

The CLI parses flags into a :class:`MutationConfig`, which the orchestrator
in ``pipeline.py`` consumes. A YAML file with the same shape is also
accepted via :meth:`MutationConfig.from_yaml`.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Public chemistry catalog. Per-experiment overrides should mutate the dict
# returned by :func:`default_chemistry_params` (a fresh copy each call) — never
# the module-level constant.
ALL_CHEMISTRIES: Tuple[str, ...] = ("iLock", "dRNA", "cDNA")

# Legacy chemistry tokens kept as input aliases. The CLI / YAML still accepts
# them but they are normalized to the new canonical names before any pipeline
# logic runs.
_LEGACY_CHEMISTRY_ALIASES: Dict[str, str] = {
    "mRNA_noiLock": "dRNA",
}


def _canonicalize_chemistry(name: str) -> str:
    """Map legacy chemistry tokens to the canonical name; pass others through."""
    if name in _LEGACY_CHEMISTRY_ALIASES:
        new = _LEGACY_CHEMISTRY_ALIASES[name]
        warnings.warn(
            f"Chemistry name {name!r} is deprecated; use {new!r} instead.",
            DeprecationWarning, stacklevel=3,
        )
        return new
    return name


@dataclass
class ChemistryParams:
    """Per-chemistry knobs. Defaults match the 2026-05 canonical chemistry set."""
    target_type: str            # "RNA" | "DNA"
    tm_model: str               # "R_DNA_NN1" | "DNA_NN4"
    tm_range: Tuple[float, float]
    check_rna_structure: bool
    mutation_position: str      # "5end_and_3tip" | "3end"
    designer: str               # "invader" | "3end" | "default"


def default_chemistry_params() -> Dict[str, ChemistryParams]:
    """Canonical 2026-05+ catalog.

    Default iLock Tm range is **50-70 C** (down from templates' 55-75).
    Per the 2026-05-10 audit decision, 50-70 rescues otherwise-unmappable
    invader-overlap mutations with no specificity loss.

    Chemistry keys use canonical names: ``iLock``, ``dRNA``, ``cDNA``.
    The legacy ``mRNA_noiLock`` token is normalized to ``dRNA`` on input.
    """
    return {
        "iLock": ChemistryParams(
            target_type="RNA", tm_model="R_DNA_NN1",
            tm_range=(50.0, 70.0),
            check_rna_structure=True,
            mutation_position="5end_and_3tip",
            designer="invader",
        ),
        "dRNA": ChemistryParams(
            target_type="RNA", tm_model="R_DNA_NN1",
            tm_range=(55.0, 75.0),
            check_rna_structure=True,
            mutation_position="3end",
            designer="3end",
        ),
        "cDNA": ChemistryParams(
            target_type="DNA", tm_model="DNA_NN4",
            tm_range=(45.0, 75.0),
            check_rna_structure=False,
            mutation_position="3end",
            designer="default",
        ),
    }


@dataclass
class MutationConfig:
    """Driver config for ``run_mutation_pipeline``.

    Required fields are passed positionally; everything else has a default.
    Use :meth:`from_yaml` to load from a YAML file.
    """
    # Required
    input_xlsx: Path
    output_dir: Path
    backbone_file: Path

    # Chemistries
    chemistries: List[str] = field(
        default_factory=lambda: list(ALL_CHEMISTRIES)
    )
    chemistry_params: Dict[str, ChemistryParams] = field(
        default_factory=default_chemistry_params
    )

    # Probe geometry knobs (shared across chemistries via FilterConfig)
    arm_len_range: Tuple[int, int] = (10, 30)
    ini_arm_len: int = 20
    # Schema-v2 (Phase 0, 2026-05-14): mutation pipeline gates on GC%; the
    # legacy g_content_range field is kept for back-compat but does not gate.
    # Mutation arms are length-variable (10–30 nt) and tightly constrained
    # by mutation context, so the default GC window mirrors the legacy
    # G-only permissiveness (essentially "anything goes"). Tightened only
    # when a specific panel needs it via YAML override.
    gc_content_range: Tuple[float, float] = (0.10, 0.90)
    g_content_range: Tuple[float, float] = (0.1, 0.9)  # DEPRECATED — informational only
    g_consecutive_max: int = 7
    mfe_min: float = -15.0
    max_tm_diff: float = 20.0
    pad: int = 50

    # Numbering / backbone
    start_no: Optional[int] = None
    last_no_from: Optional[Path] = None

    # iLock flap (default = HuangLab 2026+; legacy = TATATCCCTATAT)
    ilock_flap: str = "CGTTGCTGTGGCG"

    # RT primer (cDNA chemistry)
    rt_primer_gap: int = 15
    rt_primer_target_len: int = 20
    rt_primer_min_len: int = 15
    rt_primer_max_len: int = 30
    rt_primer_tm_range: Tuple[float, float] = (50.0, 65.0)

    # Phase toggles
    skip_blast: bool = False

    # Genome FASTA (default: <repo>/data/genome/GRCh38.fa)
    genome_path: Optional[Path] = None

    # Codebook name (e.g. "SP369"). Goes into every padlock probe_name and the
    # `codebook` column. If None, resolved from the backbone filename.
    codebook: Optional[str] = None

    def __post_init__(self) -> None:
        # Coerce path-likes
        self.input_xlsx = Path(self.input_xlsx)
        self.output_dir = Path(self.output_dir)
        self.backbone_file = Path(self.backbone_file)
        if self.last_no_from is not None:
            self.last_no_from = Path(self.last_no_from)
        if self.genome_path is not None:
            self.genome_path = Path(self.genome_path)

        # Normalize legacy chemistry tokens (e.g. mRNA_noiLock -> dRNA) so the
        # rest of the pipeline only sees canonical names. Also re-key the
        # per-chemistry params dict to match.
        self.chemistries = [_canonicalize_chemistry(c) for c in self.chemistries]
        self.chemistry_params = {
            _canonicalize_chemistry(k): v for k, v in self.chemistry_params.items()
        }

        # Validate chemistries
        unknown = [c for c in self.chemistries if c not in ALL_CHEMISTRIES]
        if unknown:
            raise ValueError(
                f"Unknown chemistry in chemistries: {unknown}. "
                f"Valid: {list(ALL_CHEMISTRIES)}"
            )
        if not self.chemistries:
            raise ValueError("chemistries must not be empty")

        # Validate Tm ranges (low < high)
        for chem, params in self.chemistry_params.items():
            lo, hi = params.tm_range
            if lo >= hi:
                raise ValueError(
                    f"Invalid tm_range for {chem}: {params.tm_range} "
                    f"(low must be strictly less than high)"
                )
        a_lo, a_hi = self.arm_len_range
        if a_lo >= a_hi:
            raise ValueError(f"Invalid arm_len_range: {self.arm_len_range}")
        g_lo, g_hi = self.g_content_range
        if g_lo >= g_hi:
            raise ValueError(f"Invalid g_content_range: {self.g_content_range}")
        rt_lo, rt_hi = self.rt_primer_tm_range
        if rt_lo >= rt_hi:
            raise ValueError(f"Invalid rt_primer_tm_range: {self.rt_primer_tm_range}")

        # Numbering: exactly one of start_no / last_no_from must be set
        if (self.start_no is None) == (self.last_no_from is None):
            raise ValueError(
                "Exactly one of start_no or last_no_from must be provided "
                f"(got start_no={self.start_no!r}, last_no_from={self.last_no_from!r})"
            )

    def resolve_start_no(self) -> int:
        """Return the No. to assign to the first probe in the panel.

        If ``start_no`` is set, returns it directly. Otherwise reads
        ``last_no_from`` and returns ``int(text) + 1``.
        """
        if self.start_no is not None:
            return int(self.start_no)
        assert self.last_no_from is not None
        return int(self.last_no_from.read_text(encoding="utf-8").strip()) + 1

    def resolve_genome_path(self, repo_root: Path) -> Path:
        """Return the genome FASTA path; default to ``<repo>/data/genome/GRCh38.fa``."""
        if self.genome_path is not None:
            return self.genome_path
        return repo_root / "data" / "genome" / "GRCh38.fa"

    def resolve_codebook(self) -> str:
        """Return the codebook name for this run.

        Uses :func:`probe_designer.io.probe_schema.resolve_codebook` to
        either honor the explicit ``codebook`` field or fall back to a
        regex match on the backbone filename.
        """
        from probe_designer.io.probe_schema import resolve_codebook
        return resolve_codebook(self.codebook, self.backbone_file)

    @classmethod
    def from_yaml(cls, path: Path) -> "MutationConfig":
        """Load a YAML file mapping field names to values.

        Per-chemistry overrides go under a ``chemistry_params:`` key, e.g.::

            input_xlsx: input/BZ30_mutations.xlsx
            output_dir: output/
            backbone_file: ../SPRINTseq_369_add_A_backbone.xlsx
            chemistries: [iLock]
            chemistry_params:
              iLock:
                tm_range: [45, 80]
            arm_len_range: [8, 35]
            mfe_min: -25.0
            start_no: 200
            ilock_flap: CGTTGCTGTGGCG
        """
        import yaml  # imported lazily so config.py stays cheap to import

        path = Path(path)
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}

        # Pull and rebuild chemistry_params if user overrode any subset.
        params_overrides = data.pop("chemistry_params", {}) or {}
        params = default_chemistry_params()
        for chem, overrides in params_overrides.items():
            chem = _canonicalize_chemistry(chem)
            if chem not in params:
                raise ValueError(f"chemistry_params: unknown chemistry {chem!r}")
            current = params[chem]
            for key, value in (overrides or {}).items():
                if not hasattr(current, key):
                    raise ValueError(
                        f"chemistry_params.{chem}: unknown field {key!r}"
                    )
                if key == "tm_range":
                    value = (float(value[0]), float(value[1]))
                setattr(current, key, value)
        data["chemistry_params"] = params

        # Coerce list-style ranges to tuples
        for key in ("arm_len_range", "g_content_range", "rt_primer_tm_range"):
            if key in data and not isinstance(data[key], tuple):
                v = data[key]
                data[key] = (type(v[0])(v[0]), type(v[1])(v[1]))

        return cls(**data)
