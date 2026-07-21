"""Central :class:`Pipeline` shared by CLI and webapp.

Owns the unified post-processing chain:
    fetch -> search -> pre_blast_filter -> BLAST -> specificity_filter ->
    score -> peak_rank -> select_top_n_with_gap

DB-awareness is injected via hooks; no direct DB import here.
"""
from __future__ import annotations

import logging
import shutil
import tempfile
import uuid
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from probe_designer.config import ConfigManager
from probe_designer.database import DatabaseInterface
from probe_designer.filtering import SequenceFilter
from probe_designer.pipeline.hooks import (
    NULL_EXISTING,
    NULL_PROGRESS,
    ExistingTargetsHook,
    IsoformProvider,
    ProgressHook,
)
from probe_designer.pipeline.result import GeneResult, PipelineResult
from probe_designer.policy import ThermoPolicy
from probe_designer.scoring import (
    compute_target_score,
    select_score_peaks,
)
from probe_designer.search_strategies import BindingSiteSearcher


logger = logging.getLogger(__name__)


# User-facing strategy names are now identical to the BindingSiteSearcher
# dispatch keys after the Phase 5 rename; the mapping is retained for
# forward-compatibility (e.g., if future user-facing names diverge again).
_STRATEGY_TO_SEARCHER_KEY: Dict[str, str] = {
    "single_sequence": "single_sequence",
    "isoform_consensus": "isoform_consensus",
    "isoform_specific": "isoform_specific",
}

_ISOFORM_STRATEGIES = {"isoform_consensus", "isoform_specific"}

_SPECIES_TO_ORGANISM = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "rat": "Rattus norvegicus",
    "zebrafish": "Danio rerio",
}


class Pipeline:
    """End-to-end probe design orchestrator.

    Args:
        config: loaded ConfigManager with all fields populated.
        progress: per-stage progress sink (webapp updates a DB row; CLI logs).
        existing_targets: lookup for previously-designed targets (DB-backed in webapp).
        isoform_provider: source of per-gene isoform lists. Defaults to the
            library's own DatabaseInterface-backed fetcher.
        output_dir: if set, result JSONs are written here by the CLI writer.
            Pipeline itself does NOT write files; callers decide.
    """

    def __init__(
        self,
        config: ConfigManager,
        *,
        progress: Optional[ProgressHook] = None,
        existing_targets: Optional[ExistingTargetsHook] = None,
        isoform_provider: Optional[IsoformProvider] = None,
        output_dir: Optional[Path] = None,
        annotations_dir: Optional[Path] = None,
    ) -> None:
        self.config = config
        self.progress: ProgressHook = progress or NULL_PROGRESS
        self.existing_targets: ExistingTargetsHook = existing_targets or NULL_EXISTING
        self.isoform_provider = isoform_provider  # may be None -> use DB default
        self.output_dir = Path(output_dir) if output_dir else None
        # Opt-in: emit per-transcript thermodynamic annotation tracks (arm Tm +
        # accessibility) at the current hyb conditions, e.g. to a NAS path.
        self._annotations_dir = Path(annotations_dir) if annotations_dir else None

        self._db = DatabaseInterface(config.database)
        self._genome_accessor = None
        try:
            self._genome_accessor = self._db.initialize_genome_accessor(config.genome)
        except Exception as exc:  # genome not essential for single_sequence path
            logger.warning("Genome accessor unavailable: %s", exc)

        self._seq_filter = SequenceFilter(config.filter, config.blast)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        genes: List[str],
        *,
        strategy: Optional[str] = None,
        top_n: int = 3,
        min_gap: int = 40,
        skip_blast: bool = False,
    ) -> PipelineResult:
        """Run the full pipeline across all genes and return in-memory results."""
        strategy = strategy or self.config.search.search_strategy
        if strategy not in _STRATEGY_TO_SEARCHER_KEY:
            raise ValueError(
                f"Unknown strategy '{strategy}'. "
                f"Expected one of: {sorted(set(_STRATEGY_TO_SEARCHER_KEY.keys()))}"
            )

        species = self.config.database.organism
        target_organism = _SPECIES_TO_ORGANISM.get(species, species)

        # Configure the searcher once; strategy name dispatches internally.
        self.config.search.search_strategy = _STRATEGY_TO_SEARCHER_KEY[strategy]
        searcher = BindingSiteSearcher(
            self.config.search,
            self.config.filter,
            genome_accessor=self._genome_accessor,
            show_progress=False,
        )

        per_gene: Dict[str, GeneResult] = {}
        total = len(genes)

        # One BLAST working dir for the whole run — avoids per-gene tempdir
        # leaks when output_dir is not provided (CLI dry-run, webapp default).
        owned_blast_tmp: Optional[Path] = None
        if self.output_dir is not None:
            blast_dir = Path(self.output_dir) / "blast_results"
            blast_dir.mkdir(parents=True, exist_ok=True)
        else:
            owned_blast_tmp = Path(tempfile.mkdtemp(prefix="probe_design_blast_"))
            blast_dir = owned_blast_tmp

        try:
            for idx, gene in enumerate(genes):
                self.progress.on_gene_start(gene, idx, total)
                result = self._run_single(
                    gene,
                    strategy=strategy,
                    species=species,
                    searcher=searcher,
                    target_organism=target_organism,
                    top_n=top_n,
                    min_gap=min_gap,
                    skip_blast=skip_blast,
                    blast_dir=blast_dir,
                )
                per_gene[gene] = result
                self.progress.on_gene_complete(gene, len(result.sites))

            self.progress.on_pipeline_complete(len(per_gene))

            return PipelineResult(
                per_gene=per_gene,
                config_snapshot=self._config_snapshot(),
                run_id=uuid.uuid4().hex[:12],
                output_dir=self.output_dir,
            )
        finally:
            if owned_blast_tmp is not None:
                shutil.rmtree(owned_blast_tmp, ignore_errors=True)

    # ------------------------------------------------------------------
    # Per-gene pipeline
    # ------------------------------------------------------------------

    def _run_single(
        self,
        gene: str,
        *,
        strategy: str,
        species: str,
        searcher: BindingSiteSearcher,
        target_organism: str,
        top_n: int,
        min_gap: int,
        skip_blast: bool,
        blast_dir: Path,
    ) -> GeneResult:
        """Run all post-processing stages for one gene."""
        result = GeneResult(gene=gene)
        result.existing_targets = self.existing_targets.lookup(gene, species)

        # 1. Sequence / isoform input for the chosen strategy
        sequences: Optional[Dict[str, Any]] = None
        isoforms: Optional[Dict[str, List[Dict[str, Any]]]] = None

        if strategy in _ISOFORM_STRATEGIES:
            self.progress.on_stage(gene, "isoforms", {})
            iso_list = self._get_isoforms(gene, species)
            if iso_list:
                isoforms = {gene: iso_list}
                result.isoform_count = len(iso_list)
            else:
                result.errors.append("No isoforms found; skipping gene")
                return result
        else:
            self.progress.on_stage(gene, "sequence", {})
            sequences = self._get_sequences_for_gene(gene, result)
            if not sequences:
                result.errors.append("No reference sequence; skipping gene")
                return result

        # 1b. Emit reference thermodynamic annotations (opt-in) for this gene's
        #     transcript(s), at the current hyb conditions. Never fatal to a run.
        if self._annotations_dir is not None:
            self._emit_annotations(gene, isoforms, sequences)

        # 2. Search binding sites
        self.progress.on_stage(gene, "search", {})
        try:
            sites_by_gene = searcher.search_all_genes(
                sequences=sequences, isoforms=isoforms
            )
        except Exception as exc:
            result.errors.append(f"search failed: {exc}")
            return result

        sites = sites_by_gene.get(gene, [])
        if not sites:
            result.errors.append("No candidates after search")
            return result

        # 3. Pre-BLAST thermal filter (strategy-agnostic)
        #    Isolated per gene like the search and BLAST stages: the filter now
        #    lets a ViennaRNA failure propagate instead of recording it as an
        #    ordinary probe rejection (audit R8), and without a boundary here one
        #    degenerate target sequence would abort an entire multi-gene run.
        #    The error stays visible in result.errors rather than being
        #    downgraded to a per-candidate failed_check.
        self.progress.on_stage(gene, "pre_blast", {"n_in": len(sites)})
        try:
            filtered_by_gene = self._seq_filter.pre_blast_filter({gene: sites})
        except Exception as exc:
            result.errors.append(f"pre-BLAST filter failed: {exc}")
            logger.warning("[%s] pre-BLAST filter failed: %s", gene, exc)
            return result
        sites = filtered_by_gene.get(gene, [])
        self.progress.on_stage(gene, "pre_blast", {"n_out": len(sites)})

        # 4. BLAST + specificity (optional)
        if not skip_blast and sites:
            self.progress.on_stage(gene, "blast", {"n_in": len(sites)})
            try:
                blast_file = self._seq_filter.run_blast(
                    {gene: sites}, str(blast_dir), target_organisms=[target_organism]
                )
                sites = self._seq_filter.specificity_filter(
                    {gene: sites}, blast_file
                ).get(gene, [])
                self.progress.on_stage(gene, "blast", {"n_out": len(sites)})
            except Exception as exc:
                result.errors.append(f"BLAST failed: {exc}")
                logger.warning(
                    "[%s] BLAST failed, continuing with pre-BLAST results: %s",
                    gene, exc,
                )

        # 5. Score, then pick the score PEAKS (greedy local-maxima / NMS):
        #    the highest-scoring sites whose binding sites don't overlap
        #    (min_distance = min_gap). Replaces the older peak_rank +
        #    select_top_n_with_gap two-stage spread.
        self.progress.on_stage(gene, "score", {"n_in": len(sites)})
        total_isoforms = result.isoform_count or 1
        for site in sites:
            site["score"] = compute_target_score(
                site,
                policy=ThermoPolicy.resolve(self.config.filter),
                total_isoforms=total_isoforms,
            )
        selected = select_score_peaks(sites, min_distance=min_gap, max_n=top_n)
        for idx, site in enumerate(selected):
            site["peak_rank"] = idx
        result.sites = selected
        return result

    def _emit_annotations(
        self,
        gene: str,
        isoforms: Optional[Dict[str, List[Dict[str, Any]]]],
        sequences: Optional[Dict[str, Any]],
    ) -> None:
        """Write Tm + accessibility bedGraph tracks for this gene's transcript(s)
        into the annotations dir, at the config's hyb conditions. Opt-in and
        best-effort: a failed isoform is skipped, never aborting the run."""
        from probe_designer.annotate import emit_annotations_for_sequences
        from probe_designer.search_strategies import IsoformAwareness

        reaction = self.config.filter.reaction_conditions()
        folding = self.config.filter.folding_conditions()
        arm_len = max(1, int(self.config.search.binding_site_length) // 2)

        seqs: Dict[str, str] = {}
        if isoforms and isoforms.get(gene):
            if not self._genome_accessor:
                logger.warning("[%s] annotations: no genome accessor; skipped", gene)
                return
            awareness = IsoformAwareness(isoforms[gene], self._genome_accessor)
            for iso in isoforms[gene]:
                name = iso.get("display_name") or iso.get("id") or gene
                try:
                    seqs[name] = awareness.get_isoform_mrna(iso)
                except Exception as exc:
                    logger.warning("[%s] annotations: skip isoform %s: %s",
                                   gene, name, exc)
        elif sequences:
            for sid, seq in sequences.items():
                if isinstance(seq, str):
                    seqs[str(sid)] = seq

        written = emit_annotations_for_sequences(
            seqs, reaction, self._annotations_dir, arm_length=arm_len,
            folding=folding,
        )

        # Canonical-transcript genome projection (for IGV overlay on the
        # assembly). Only the representative transcript is projected.
        if isoforms and isoforms.get(gene) and seqs:
            from probe_designer.annotate import build_canonical_genome_annotations
            from probe_designer.genome_projection import pick_canonical_isoform
            canonical = pick_canonical_isoform(isoforms[gene])
            cname = canonical.get("display_name") or gene
            if cname in seqs and canonical.get("Exon"):
                try:
                    written += build_canonical_genome_annotations(
                        canonical, seqs[cname], reaction,
                        out_dir=self._annotations_dir / "genome",
                        arm_length=arm_len, folding=folding, gene=gene,
                    )
                except Exception as exc:
                    logger.warning("[%s] genome-projection annotation failed: %s",
                                   gene, exc)

        logger.info("[%s] wrote %d annotation track(s) to %s",
                    gene, len(written), self._annotations_dir)

    # ------------------------------------------------------------------
    # Isoform / sequence data access
    # ------------------------------------------------------------------

    def _get_isoforms(self, gene: str, species: str) -> Optional[List[Dict[str, Any]]]:
        if self.isoform_provider is not None:
            return self.isoform_provider.get_isoforms(gene, species)
        # Default: designer's own DB interface (Ensembl REST / cached)
        try:
            return self._db.get_isoform_info(gene)
        except Exception as exc:
            logger.warning("[%s] isoform fetch failed: %s", gene, exc)
            return None

    def _get_sequences_for_gene(
        self, gene: str, result: GeneResult
    ) -> Optional[Dict[str, Any]]:
        try:
            fetched = self._db.get_gene_sequences([gene])
        except Exception as exc:
            result.errors.append(f"sequence fetch failed: {exc}")
            return None
        sequences = (fetched or {}).get("sequences", {})
        seq_data = sequences.get(gene)
        if not seq_data or not seq_data.get("sequence"):
            return None
        result.sequence_length = len(seq_data["sequence"])
        result.sequence_source = self.config.database.database_type
        return {gene: seq_data}

    # ------------------------------------------------------------------
    # Infrastructure
    # ------------------------------------------------------------------

    def _config_snapshot(self) -> Dict[str, Any]:
        """Serializable view of the current config for run provenance."""
        return {
            "database": asdict(self.config.database),
            "search": asdict(self.config.search),
            "filter": asdict(self.config.filter),
            "blast": asdict(self.config.blast),
            "genome": asdict(self.config.genome),
            "output": asdict(self.config.output),
        }
