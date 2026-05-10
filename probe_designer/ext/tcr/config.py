"""Configuration dataclass for the TCR pipeline.

The CLI parses flags into a :class:`TcrConfig`, which the orchestrator in
``pipeline.py`` consumes. A YAML file with the same shape is also accepted
via :meth:`TcrConfig.from_yaml`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple


ALL_CHEMISTRIES: Tuple[str, ...] = ("mRNA", "cDNA")


@dataclass
class TcrConfig:
    """Driver config for ``run_tcr_pipeline``.

    Required fields are passed positionally; everything else has a default.
    Use :meth:`from_yaml` to load from a YAML file.
    """
    # Required
    input_xlsx: Path
    output_dir: Path
    backbone_file: Path

    # Chemistries — default = direct mRNA only (cDNA is opt-in).
    # Whenever 'cDNA' is in the list, RT primer design runs automatically.
    chemistries: List[str] = field(default_factory=lambda: ["mRNA"])

    # Scan geometry
    bds_len: int = 40
    arm_len: int = 20
    sites_per_clone: int = 3

    # mRNA arm Tm filter (cDNA Tm is informational, not gated)
    tm_range: Tuple[float, float] = (50.0, 60.0)

    # Numbering / backbone
    start_no: Optional[int] = None
    last_no_from: Optional[Path] = None

    # RT primer (only used when 'cDNA' is in chemistries)
    rt_primer_gap: int = 15
    rt_primer_target_len: int = 20
    rt_primer_min_len: int = 15
    rt_primer_max_len: int = 30
    rt_primer_tm_range: Tuple[float, float] = (50.0, 65.0)

    # Phase toggles
    skip_blast: bool = False
    plot_tm_landscape: bool = True

    def __post_init__(self) -> None:
        # Coerce path-likes
        self.input_xlsx = Path(self.input_xlsx)
        self.output_dir = Path(self.output_dir)
        self.backbone_file = Path(self.backbone_file)
        if self.last_no_from is not None:
            self.last_no_from = Path(self.last_no_from)

        # Validate chemistries
        unknown = [c for c in self.chemistries if c not in ALL_CHEMISTRIES]
        if unknown:
            raise ValueError(
                f"Unknown chemistry in chemistries: {unknown}. "
                f"Valid: {list(ALL_CHEMISTRIES)}"
            )
        if not self.chemistries:
            raise ValueError("chemistries must not be empty")

        # Validate scan geometry
        if 2 * self.arm_len != self.bds_len:
            raise ValueError(
                f"arm_len ({self.arm_len}) * 2 must equal bds_len ({self.bds_len})"
            )
        if self.sites_per_clone < 1:
            raise ValueError(f"sites_per_clone must be >= 1, got {self.sites_per_clone}")

        # Validate Tm ranges
        for label, lo_hi in (("tm_range", self.tm_range),
                             ("rt_primer_tm_range", self.rt_primer_tm_range)):
            lo, hi = lo_hi
            if lo >= hi:
                raise ValueError(
                    f"Invalid {label}: {lo_hi} (low must be strictly less than high)"
                )

        # Numbering: exactly one of start_no / last_no_from must be set
        if (self.start_no is None) == (self.last_no_from is None):
            raise ValueError(
                "Exactly one of start_no or last_no_from must be provided "
                f"(got start_no={self.start_no!r}, last_no_from={self.last_no_from!r})"
            )

    def has_cdna(self) -> bool:
        """True if cDNA chemistry is in the selection (RT primer phase will run)."""
        return "cDNA" in self.chemistries

    def resolve_start_no(self) -> int:
        """Return the No. to assign to the first probe in the panel."""
        if self.start_no is not None:
            return int(self.start_no)
        assert self.last_no_from is not None
        return int(self.last_no_from.read_text(encoding="utf-8").strip()) + 1

    @classmethod
    def from_yaml(cls, path: Path) -> "TcrConfig":
        """Load a YAML file mapping field names to values.

        Example::

            input_xlsx: input/BZ30_tcr_clones.xlsx
            output_dir: output/
            backbone_file: ../SPRINTseq_369_add_A_backbone.xlsx
            chemistries: [mRNA, cDNA]
            bds_len: 40
            arm_len: 20
            tm_range: [55, 65]
            start_no: 200
            rt_primer_gap: 15
        """
        import yaml  # lazy import

        path = Path(path)
        with path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}

        # Coerce list-style ranges to tuples
        for key in ("tm_range", "rt_primer_tm_range"):
            if key in data and not isinstance(data[key], tuple):
                v = data[key]
                data[key] = (float(v[0]), float(v[1]))

        return cls(**data)
