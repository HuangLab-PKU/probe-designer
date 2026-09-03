"""``probe-design design`` — full padlock probe design pipeline.

Thin CLI wrapper around :class:`probe_designer.pipeline.Pipeline`.
"""
from __future__ import annotations

import json
import logging
from enum import Enum
from pathlib import Path
from typing import List, Optional

import typer

from probe_designer.config import ConfigManager
from probe_designer.genome import DefaultIsoformProvider
from probe_designer.pipeline import Pipeline, ProgressHook
from probe_designer.utils import load_gene_list


logger = logging.getLogger(__name__)


class Strategy(str, Enum):
    """User-facing strategy names."""

    single_sequence = "single_sequence"
    isoform_consensus = "isoform_consensus"
    isoform_specific = "isoform_specific"


class _LoggingProgress:
    """ProgressHook that logs every stage transition at INFO level."""

    def on_gene_start(self, gene: str, index: int, total: int) -> None:
        logger.info("[%d/%d] Starting %s", index + 1, total, gene)

    def on_stage(self, gene: str, stage: str, info: dict) -> None:
        if info:
            logger.debug("[%s] stage=%s %s", gene, stage, info)
        else:
            logger.debug("[%s] stage=%s", gene, stage)

    def on_gene_complete(self, gene: str, site_count: int) -> None:
        logger.info("[%s] %d site(s) selected", gene, site_count)

    def on_pipeline_complete(self, gene_count: int) -> None:
        logger.info("Pipeline done (%d genes)", gene_count)


def _apply_buffer_overrides(
    filter_cfg, *, monovalent_mm, mg_mm, dntp_mm, strand_nm,
    formamide_pct, lab_temp_c, tm_margin_c, saltcorr,
) -> None:
    """Apply CLI reaction-buffer overrides onto a FilterConfig in place.

    Only non-None values override the config/defaults. If ``lab_temp_c`` or
    ``tm_margin_c`` change, the Tm window is re-anchored to the reaction
    temperature: min_tm = lab_temp_c + tm_margin_c, max_tm = lab_temp_c + 25.
    The buffer is validated eagerly (fail-fast) before the run starts.
    """
    overrides = {
        "monovalent_mM": monovalent_mm, "mg_mM": mg_mm, "dntp_mM": dntp_mm,
        "strand_nM": strand_nm, "formamide_pct": formamide_pct,
        "lab_temp_c": lab_temp_c, "tm_margin_c": tm_margin_c, "saltcorr": saltcorr,
    }
    for field, value in overrides.items():
        if value is not None:
            setattr(filter_cfg, field, value)
    if lab_temp_c is not None or tm_margin_c is not None:
        filter_cfg.min_tm = filter_cfg.lab_temp_c + filter_cfg.tm_margin_c
        filter_cfg.max_tm = filter_cfg.lab_temp_c + 25.0
    filter_cfg.reaction_conditions()  # validate the composition, fail-fast


def design(
    genes_file: Path = typer.Option(
        ..., "--genes-file", exists=True, dir_okay=False, readable=True,
        help="Text file with one gene symbol per line."
    ),
    strategy: Strategy = typer.Option(
        Strategy.isoform_consensus, "--strategy", case_sensitive=False,
        help="Search strategy. single_sequence = scan one mRNA/isoform FASTA; "
             "isoform_consensus = shared across isoforms; isoform_specific = unique to an isoform."
    ),
    config_path: Optional[Path] = typer.Option(
        None, "--config", exists=True, dir_okay=False, readable=True,
        help="YAML config. Defaults to configs/config_<strategy>.yaml bundled with the package."
    ),
    species: Optional[str] = typer.Option(
        None, "--species", help="Override species (e.g., human, mouse)."
    ),
    blast_species: Optional[List[str]] = typer.Option(
        None, "--blast-species", help="Target organism(s) for BLAST specificity filter."
    ),
    genome_fasta: Optional[Path] = typer.Option(
        None, "--genome-fasta", help="Local genome FASTA path (overrides config)."
    ),
    gtf_path: Optional[Path] = typer.Option(
        None, "--gtf", exists=True, dir_okay=False, readable=True,
        help="Local GTF path for offline isoform lookup (isoform strategies only)."
    ),
    output_dir: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output directory (overrides config.output.output_dir)."
    ),
    top_n: int = typer.Option(
        3, "--top-n", "-n", min=1,
        help="Final probes to select per gene after scoring and peak-ranking."
    ),
    min_gap: int = typer.Option(
        40, "--min-gap", min=0,
        help="Minimum bp between selected site midpoints."
    ),
    skip_blast: bool = typer.Option(
        False, "--skip-blast",
        help="Skip BLAST + specificity step (faster; only thermal filter applied)."
    ),
    isoform_biotype: Optional[List[str]] = typer.Option(
        None, "--isoform-biotype",
        help="Only consider isoforms with this Ensembl biotype during "
             "isoform_consensus / isoform_specific. Pass multiple times to "
             "keep multiple biotypes. Common: 'protein_coding'. Without this "
             "flag, all isoforms (including retained_intron, NMD, "
             "processed_transcript) participate, which can dilute consensus."
    ),
    # --- Isoform crediting (isoform_consensus) ---
    min_isoform_arm_nt: Optional[int] = typer.Option(
        None, "--min-isoform-arm-nt", min=0,
        help="An isoform counts as targeted only when a contiguous run of its "
             "own mRNA spans the ligation junction with at least this many nt "
             "on each side (the junction-centred core rule). Default 16. "
             "0 credits any isoform the junction is paired on."
    ),
    allow_noncoding_only: bool = typer.Option(
        False, "--allow-noncoding-only",
        help="Keep sites that bind no protein-coding transcript of their own "
             "gene. Off by default: such sites cannot report the gene's "
             "expression. Genes with no coding transcript are exempt anyway."
    ),
    # --- Reaction buffer overrides (chemistry.ReactionConditions) ---
    monovalent_mm: Optional[float] = typer.Option(
        None, "--monovalent-mm", help="Reaction monovalent K+ (mM). Default 75."
    ),
    mg_mm: Optional[float] = typer.Option(
        None, "--mg-mm", help="Reaction Mg2+ (mM). Default 10."
    ),
    dntp_mm: Optional[float] = typer.Option(
        None, "--dntp-mm", help="Reaction dNTP (mM). Default 0."
    ),
    strand_nm: Optional[float] = typer.Option(
        None, "--strand-nm", help="Probe/oligo concentration (nM). Default 100."
    ),
    formamide_pct: Optional[float] = typer.Option(
        None, "--formamide-pct", help="Formamide (percent v/v); depresses Tm. Default 20."
    ),
    lab_temp_c: Optional[float] = typer.Option(
        None, "--lab-temp-c",
        help="Hyb anneal temperature (C). Re-anchors Tm window "
             "(min_tm = lab_temp_c + tm_margin_c). Default 45."
    ),
    tm_margin_c: Optional[float] = typer.Option(
        None, "--tm-margin-c",
        help="Arm-Tm margin above lab_temp_c (C). Re-anchors min_tm. Default 5."
    ),
    saltcorr: Optional[int] = typer.Option(
        None, "--saltcorr", min=1, max=7,
        help="Biopython salt method (5=SantaLucia'98+vonAhsen Mg, 7=Owczarzy'08)."
    ),
    annotations_dir: Optional[Path] = typer.Option(
        None, "--annotations-dir",
        help="If set, write per-transcript Tm + accessibility bedGraph tracks "
             "(at the hyb conditions) here as a byproduct of the run, e.g. a NAS "
             "path beside the genome."
    ),
) -> None:
    """Run the end-to-end design pipeline and write scored/selected results."""
    if config_path is None:
        pkg_root = Path(__file__).resolve().parent.parent.parent  # designer/
        default_name = {
            Strategy.single_sequence: "config_single_sequence.yaml",
            Strategy.isoform_consensus: "config_consensus.yaml",
            Strategy.isoform_specific: "config_specific.yaml",
        }[strategy]
        config_path = pkg_root / "configs" / default_name
        if not config_path.exists():
            raise typer.BadParameter(
                f"No --config provided and default not found: {config_path}"
            )
    logger.info("Using config: %s", config_path)

    cfg = ConfigManager(str(config_path), species=species)
    cfg.search.search_strategy = strategy.value

    if min_isoform_arm_nt is not None:
        cfg.search.min_isoform_arm_nt = min_isoform_arm_nt
    if allow_noncoding_only:
        cfg.search.require_protein_coding = False
    if blast_species:
        cfg.blast.species = list(blast_species)
    if genome_fasta:
        cfg.genome.genome_fasta_path = str(genome_fasta)
    if output_dir:
        cfg.output.output_dir = str(output_dir)

    _apply_buffer_overrides(
        cfg.filter,
        monovalent_mm=monovalent_mm, mg_mm=mg_mm, dntp_mm=dntp_mm,
        strand_nm=strand_nm, formamide_pct=formamide_pct,
        lab_temp_c=lab_temp_c, tm_margin_c=tm_margin_c, saltcorr=saltcorr,
    )

    errors = cfg.validate_config()
    if errors:
        for err in errors:
            typer.echo(f"Config error: {err}", err=True)
        raise typer.Exit(code=2)

    run_dir = Path(cfg.output.output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    genes = load_gene_list(str(genes_file))
    logger.info("Processing %d genes with strategy=%s", len(genes), strategy.value)

    isoform_provider = None
    if strategy in (Strategy.isoform_consensus, Strategy.isoform_specific):
        isoform_provider = DefaultIsoformProvider(
            gtf_path=gtf_path,
            keep_biotypes=isoform_biotype or None,
        )

    pipeline = Pipeline(
        cfg,
        progress=_LoggingProgress(),
        isoform_provider=isoform_provider,
        output_dir=run_dir,
        annotations_dir=annotations_dir,
    )
    result = pipeline.run(
        genes,
        strategy=strategy.value,
        top_n=top_n,
        min_gap=min_gap,
        skip_blast=skip_blast,
    )

    # Persist results
    selected_json = run_dir / "selected_binding_sites.json"
    selected_json.write_text(
        json.dumps(result.sites_by_gene(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    config_snapshot_file = run_dir / "config_used.json"
    config_snapshot_file.write_text(
        json.dumps(result.config_snapshot, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    gene_report = run_dir / "gene_report.json"
    gene_report.write_text(
        json.dumps(
            {g: r.to_dict() for g, r in result.per_gene.items()},
            indent=2, ensure_ascii=False, default=str,
        ),
        encoding="utf-8",
    )

    typer.echo(
        f"Done. {result.total_sites} sites selected across {result.gene_count} genes "
        f"-> {selected_json}"
    )
