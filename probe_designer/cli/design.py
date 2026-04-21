"""``probe-design design`` — full padlock probe design pipeline.

Phase 1 port of scripts/find_consensus_probes.py into a single Typer command
that handles all three strategies AND auto-invokes scoring + peak_rank +
top-N selection at the end (parity with webapp; closes the "CLI returns
unranked probes" gap).

Strategy name mapping (Phase 1 internal):
    single_sequence  -> BruteForceStrategy  (renamed to SingleSequenceStrategy in Phase 5)
    isoform_consensus -> IsoformConsensusStrategy
    isoform_specific  -> IsoformSpecificStrategy
"""
from __future__ import annotations

import json
import logging
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

import typer

from probe_designer.config import ConfigManager
from probe_designer.database import DatabaseInterface
from probe_designer.filtering import SequenceFilter
from probe_designer.scoring import (
    compute_target_score,
    peak_rank,
    select_top_n_with_gap,
)
from probe_designer.search_strategies import BindingSiteSearcher
from probe_designer.utils import load_gene_list


logger = logging.getLogger(__name__)


class Strategy(str, Enum):
    """User-facing strategy names. Internal class dispatch in ``_INTERNAL_STRATEGY``."""
    single_sequence = "single_sequence"
    isoform_consensus = "isoform_consensus"
    isoform_specific = "isoform_specific"


# Phase 1: new user-facing name -> current BindingSiteSearcher dispatch key.
# Phase 5 will align class names with user-facing names.
_INTERNAL_STRATEGY: Dict[Strategy, str] = {
    Strategy.single_sequence: "brute_force",
    Strategy.isoform_consensus: "isoform_consensus",
    Strategy.isoform_specific: "isoform_specific",
}


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
    progress: bool = typer.Option(
        False, "--progress", help="Show per-gene progress bars during search."
    ),
) -> None:
    """Run the end-to-end design pipeline and write scored/selected results."""
    if config_path is None:
        # Bundled default config per strategy
        pkg_root = Path(__file__).resolve().parent.parent.parent  # designer/
        default_name = {
            Strategy.single_sequence: "config_bruteforce.yaml",  # file rename deferred to Phase 5
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
    cfg.search.search_strategy = _INTERNAL_STRATEGY[strategy]

    if blast_species:
        cfg.blast.species = list(blast_species)
    if genome_fasta:
        cfg.genome.genome_fasta_path = str(genome_fasta)
    if output_dir:
        cfg.output.output_dir = str(output_dir)

    # Validate final resolved config
    errors = cfg.validate_config()
    if errors:
        for err in errors:
            typer.echo(f"Config error: {err}", err=True)
        raise typer.Exit(code=2)

    run_dir = Path(cfg.output.output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    genes = load_gene_list(str(genes_file))
    logger.info("Processing %d genes with strategy=%s", len(genes), strategy.value)

    db = DatabaseInterface(cfg.database)
    accessor = db.initialize_genome_accessor(cfg.genome)

    # Strategy-specific input assembly
    if strategy is Strategy.single_sequence:
        sequences = _fetch_sequences(db, genes)
        isoforms_map = None
    else:
        isoforms_map = _fetch_isoforms(db, genes, run_dir)
        sequences = None

    # Search
    searcher = BindingSiteSearcher(
        cfg.search, cfg.filter, genome_accessor=accessor, show_progress=progress
    )
    binding_sites = searcher.search_all_genes(sequences=sequences, isoforms=isoforms_map)

    sites_json = run_dir / "binding_sites.json"
    sites_json.write_text(
        json.dumps(binding_sites, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    logger.info(
        "Pre-filter: %d sites across %d genes saved to %s",
        sum(len(s) for s in binding_sites.values()), len(binding_sites), sites_json,
    )

    # Filter + optional BLAST
    seq_filter = SequenceFilter(cfg.filter, cfg.blast)
    filtered = seq_filter.pre_blast_filter(binding_sites)

    if skip_blast:
        logger.warning("--skip-blast: specificity filter skipped; off-target risk not assessed.")
    else:
        blast_dir = run_dir / "blast_results"
        blast_dir.mkdir(exist_ok=True)
        try:
            blast_results_file = seq_filter.run_blast(filtered, str(blast_dir))
            filtered = seq_filter.specificity_filter(filtered, blast_results_file)
        except Exception as exc:  # Network / local BLAST failures are common
            logger.error("BLAST failed (%s); keeping pre-BLAST thermal results only.", exc)

    # Score + rank + select  (the Phase 1 parity win for CLI)
    scored_by_gene: Dict[str, List[Dict[str, Any]]] = {}
    for gene, sites in filtered.items():
        total_isoforms = len(isoforms_map.get(gene, [])) if isoforms_map else 1
        for site in sites:
            site["score"] = compute_target_score(
                site,
                min_arm_tm=cfg.filter.min_tm,
                max_tm_diff=cfg.filter.max_tm_diff,
                total_isoforms=total_isoforms,
            )
        ranked = peak_rank(sites, region_size=80, min_gap=min_gap)
        for idx, site in enumerate(ranked):
            site["peak_rank"] = idx
        scored_by_gene[gene] = ranked

    scored_json = run_dir / "scored_binding_sites.json"
    scored_json.write_text(
        json.dumps(scored_by_gene, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    logger.info("Scored and ranked results written to %s", scored_json)

    selected_by_gene = {
        gene: select_top_n_with_gap(sites, top_n=top_n, min_gap=min_gap)
        for gene, sites in scored_by_gene.items()
    }
    selected_json = run_dir / "selected_binding_sites.json"
    selected_json.write_text(
        json.dumps(selected_by_gene, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    typer.echo(
        f"Done. {sum(len(v) for v in selected_by_gene.values())} sites selected "
        f"across {len(selected_by_gene)} genes → {selected_json}"
    )


def _fetch_sequences(db: DatabaseInterface, genes: List[str]) -> Dict[str, Any]:
    """Fetch one reference sequence per gene for single_sequence strategy."""
    logger.info("Fetching per-gene reference sequences for single_sequence strategy")
    out = db.get_gene_sequences(genes)
    return out or {}


def _fetch_isoforms(
    db: DatabaseInterface, genes: List[str], run_dir: Path
) -> Dict[str, List[Dict]]:
    """Fetch per-gene isoform lists, caching each gene to its own JSON file."""
    cache_dir = run_dir / "isoform_info"
    cache_dir.mkdir(exist_ok=True)
    isoforms_map: Dict[str, List[Dict]] = {}
    cached = 0
    for g in genes:
        cache_file = cache_dir / f"{g}.json"
        if cache_file.exists():
            isoforms_map[g] = json.loads(cache_file.read_text(encoding="utf-8"))
            cached += 1
        else:
            isoforms = db.get_isoform_info(g)
            cache_file.write_text(
                json.dumps(isoforms, indent=2, ensure_ascii=False), encoding="utf-8"
            )
            isoforms_map[g] = isoforms
    logger.info(
        "Fetched isoforms for %d genes (%d cache hits)", len(isoforms_map), cached
    )
    return isoforms_map
