"""
Configuration management module
Manages all configurable parameters for DNA probe design.
"""

import dataclasses
import os
import warnings
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import json
import yaml

from probe_designer.chemistry import FoldingConditions, ReactionConditions


@dataclass
class DatabaseConfig:
    """Database configuration."""
    organism: str = "mouse"  # species
    database_type: str = "ensembl"  # supported: ensembl, ncbi
    coord_system_version: str = "GRCh38"  # coordinate system version
    max_retries: int = 3  # max retries for API calls


@dataclass
class SearchConfig:
    """Binding site search configuration."""
    binding_site_length: int = 40  # binding site length
    max_binding_sites: int = 999999  # effectively unlimited (scoring handles ranking)
    search_strategy: str = "isoform_consensus"  # single_sequence, isoform_consensus, isoform_specific, exon_junction
    step_size: Optional[int] = None  # step size for single_sequence strategy
    genes_file: Optional[str] = None  # path to gene list file (can be overridden by command line)


@dataclass
class FilterConfig:
    """Sequence filtering configuration.

    Schema-v2 (Phase 0, 2026-05-14) distinguishes G-only content from GC content.
    Historically `min_g_content` / `max_g_content` counted only G nucleotides,
    which under-rejects cDNA-chemistry arms (G is depleted by mRNA→cDNA strand
    flip) and over-rejects when GC is normal. The new gate is GC fraction
    (G+C), with defaults anchored to a 2026-05-13 audit of BZ23/TNBC panels
    (`experiments/20260513_gc_content_audit/FINDINGS.md`).

    The legacy G-only fields are accepted for backward-compatibility and emit
    a DeprecationWarning when used; new pipelines should set `min/max_gc_content`.
    """
    # Pre-BLAST rules (thermal properties)
    # --- GC-content gate (Schema-v2, Phase 0) ---
    min_gc_content: float = 0.35  # min G+C fraction per arm (audit-anchored compromise)
    max_gc_content: float = 0.65  # max G+C fraction per arm
    # --- legacy G-only gate (kept for backward compat; emits deprecation warning) ---
    min_g_content: float = 0.3  # DEPRECATED — counts only G nucleotides; prefer min_gc_content
    max_g_content: float = 0.7  # DEPRECATED — counts only G nucleotides; prefer max_gc_content
    # --- homopolymer gates: all four bases now checked ---
    max_consecutive_g: int = 5  # max consecutive G; default raised to 5 to align with A/T/C
    max_consecutive_a: int = 5  # max consecutive A
    max_consecutive_t: int = 5  # max consecutive T
    max_consecutive_c: int = 5  # max consecutive C
    # --- Tm + structure ---
    # Re-anchored 2026-07-17 to the reaction temperature: min = lab_temp_c +
    # tm_margin_c (45 + 5 = 50), max = lab_temp_c + 25 (70). Computed at the
    # reaction buffer incl. formamide. Impact on shipped marker arms: 0.4%
    # (experiments/20260717_tm_buffer_config_impl/output/impact_summary.txt).
    min_tm: float = 50.0  # min melting temperature (°C); = lab_temp_c + tm_margin_c
    max_tm: float = 70.0  # max melting temperature (°C); = lab_temp_c + 25
    # Two-arm balance. Tightened 10 -> 5 on 2026-07-21 (audit R4) to the
    # Larsson 2010 / Krzywkowski 2017 field practice. It is a SCORING reference,
    # not a cutoff: as a hard gate 5 C would have rejected 17.7% of 232 shipped
    # probes (19.7% of cDNA), while as a reference it costs nothing and sharpens
    # ranking exactly where the discrimination matters.
    max_tm_diff: float = 5.0  # max Tm difference between 3' and 5' halves
    # 2026-07-19: the hard Tm-range gate is OFF by default. Absolute arm Tm is
    # now a SOFT scoring term (scoring.compute_target_score tm_proximity, peaked
    # at the reaction-anchored target min_tm) rather than a pass/fail cutoff, so
    # candidates rank by score instead of being dropped on Tm. min_tm/max_tm
    # remain the scoring target/reference. Set True to restore hard rejection.
    # 2026-07-21: the balance rule joined it, for the same reason (audit R4).
    # Both are read through policy.ThermoPolicy, shared with the ranker.
    enforce_tm_gate: bool = False
    enforce_tm_diff_gate: bool = False
    min_free_energy: float = -10.0  # min free energy (kcal/mol)
    check_rna_structure: bool = False  # whether to check RNA secondary structure
    # --- Phase 1A: target-accessibility gate (RNAplfold) ---
    # When min_accessibility > 0, the search strategy computes per-base
    # unpaired probability via RNAplfold over the FULL transcript and
    # gates each candidate window on its mean accessibility. When
    # min_accessibility == 0 (default), the legacy check_rna_structure
    # self-fold path runs instead. See filtering/accessibility.py.
    # Geometry re-anchored 2026-07-21 (audit R9) to Lange 2012: W = L + 50 with
    # L ~ 100. The old W=70/L=40 could not represent a base pair spanning >40 nt,
    # so sites behind a long-range stem read as open (mean p_unpaired over-stated
    # by ~0.10 on real transcripts). Build a FoldingConditions via
    # folding_conditions(); see chemistry.FoldingConditions.
    min_accessibility: float = 0.0  # 0 disables; recommended 0.30 when on
    plfold_window: int = 150        # RNAplfold W (sliding window size)
    plfold_span: int = 100          # RNAplfold L (max bp span within window)
    # None = track the reaction (ReactionConditions.effective_celsius), so the
    # accessibility that filters candidates matches the reference annotation.
    plfold_temperature: Optional[float] = None

    # --- Reaction buffer (2026-07-17): the real hybridization conditions used
    #     for every Tm/ΔG calculation. Flat fields so the generic YAML loader
    #     round-trips them; build a ReactionConditions via reaction_conditions().
    #     Defaults follow protocol rca.md v5.3 (see chemistry.ReactionConditions).
    #     Before this change all Tm used Biopython defaults (Na=50, Mg=0, 50 nM),
    #     understating arm Tm by ~10 °C — see experiments/20260715_tm_deltag_methods_audit/.
    monovalent_mM: float = 75.0        # K+ (Ampligase 25 + added KCl 50)
    mg_mM: float = 10.0                # Mg2+ (Ampligase)
    dntp_mM: float = 0.0               # dNTP (chelates Mg2+)
    strand_nM: float = 100.0           # 0.1 µM probe
    formamide_pct: float = 20.0        # % formamide (depresses Tm)
    formamide_deg_per_pct: float = 0.5  # °C Tm depression per % formamide
    solvents: Dict[str, float] = field(default_factory=dict)  # extra co-solvents {name: %}; see chemistry.SOLVENT_TM_COEFF
    lab_temp_c: float = 45.0           # hyb anneal temperature; arm-Tm anchor
    tm_margin_c: float = 5.0           # min arm Tm = lab_temp_c + tm_margin_c
    saltcorr: int = 5                  # Biopython salt method (5 = SantaLucia'98 + von Ahsen Mg)

    # Post-BLAST rules
    max_alignments: int = 5  # max allowed alignments (specificity)
    require_specificity: bool = True  # enforce specificity in BLAST filter
    target_organisms: List[str] = None  # target organisms to allow

    def reaction_conditions(self) -> ReactionConditions:
        """Materialize the flat buffer fields into a ReactionConditions.

        Validation (fail-fast) happens in ReactionConditions.__post_init__.
        """
        return ReactionConditions(
            monovalent_mM=self.monovalent_mM,
            mg_mM=self.mg_mM,
            dntp_mM=self.dntp_mM,
            strand_nM=self.strand_nM,
            formamide_pct=self.formamide_pct,
            formamide_deg_per_pct=self.formamide_deg_per_pct,
            solvents=dict(self.solvents),
            lab_temp_c=self.lab_temp_c,
            tm_margin_c=self.tm_margin_c,
            saltcorr=self.saltcorr,
        )

    def folding_conditions(self) -> FoldingConditions:
        """Materialize the flat RNAplfold fields into a FoldingConditions.

        Validation (fail-fast) happens in FoldingConditions.__post_init__.
        """
        return FoldingConditions(
            window=self.plfold_window,
            span=self.plfold_span,
            temperature_c=self.plfold_temperature,
        )


@dataclass
class ProbeConfig:
    """Probe assembly configuration."""
    backbone_file: Optional[str] = None  # Backbone Excel file path (must contain 'No.' and 'Sequence' columns)
    codebook: Optional[str] = None  # Codebook name (e.g. SP369); falls back to backbone-filename regex.


@dataclass
class BlastConfig:
    """BLAST configuration."""
    blast_type: str = "local"  # local, online
    database: str = "refseq_rna"  # BLAST database
    task: str = "megablast"  # BLAST task (megablast, blastn, etc.)
    evalue: float = 1e-5  # e-value threshold
    hitlist_size: int = 50  # number of hits to return
    alignments: int = 500  # number of alignments to show
    descriptions: int = 500  # number of descriptions to show
    megablast: bool = True  # use megablast algorithm (for blastn)
    short_query: bool = False  # adjust parameters for short queries
    filter: str = "none"  # filtering option
    format_type: str = "XML"  # output format
    service: str = "plain"  # BLAST service type
    # Batching and concurrency for online BLAST
    batch_size: int = 100  # number of sequences per online request
    concurrency: int = 3   # number of concurrent online requests
    species: Optional[List[str]] = None  # target species for BLAST search (if None, use config species)


@dataclass
class OutputConfig:
    """Output configuration."""
    output_dir: str = "results"  # output directory
    create_timestamp: bool = True  # create timestamped subdir
    save_intermediate: bool = True  # save intermediate files
    # file_formats: List[str] = None  # output formats
    save_fasta: bool = True
    save_json: bool = True
    save_excel: bool = True


@dataclass
class GenomeConfig:
    """Genome access configuration."""
    genome_fasta_path: Optional[str] = None  # path to local genome FASTA
    accessor_priority: Optional[List[str]] = None  # priority list: ['local', 'ensembl'] or ['ensembl', 'local']
    
    def __post_init__(self):
        """Set default priority if not specified."""
        if self.accessor_priority is None:
            self.accessor_priority = ['local', 'ensembl']


@dataclass
class SpeciesConfig:
    """Species-specific configuration."""
    organism: str
    coord_system_version: str
    genome_fasta_path: str
    display_name: str
    taxonomy_id: int


class ConfigManager:
    """Config manager for the application."""
    
    def __init__(self, config_file: Optional[str] = None, species: Optional[str] = None):
        self.database = DatabaseConfig()
        self.search = SearchConfig()
        self.filter = FilterConfig()
        self.probe = ProbeConfig()
        self.blast = BlastConfig()
        self.output = OutputConfig()
        self.genome = GenomeConfig()
        self.species_config: Optional[SpeciesConfig] = None  # Will be set by set_species()
        
        # if self.output.file_formats is None:
        #     self.output.file_formats = ["xlsx", "fasta", "json"]
        
        # load species configuration first
        self._load_species_config()
        
        # apply species-specific settings
        if species:
            self.set_species(species)

        # load main configuration
        if config_file and os.path.exists(config_file):
            # Load the main configuration first
            self.load_config(config_file)
            
            # Check if species is specified in config file and apply it after load_config
            # This ensures species-specific settings (like genome_fasta_path) take precedence
            # over any settings that might have been loaded from the config file
            config_data = None
            with open(config_file, 'r') as f:
                if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                    config_data = yaml.safe_load(f)
                else:
                    config_data = json.load(f)
            
            # If species is specified in config file but not as parameter, apply it
            if config_data and 'species' in config_data and not species:
                self.set_species(config_data['species'])
        elif config_file:
            print(f"Configuration file not found: {config_file}, using default config")

    
    def _load_species_config(self):
        """Load species configuration from species_config.json."""
        # __file__ is probe_designer/config/__init__.py; go up two levels to reach designer/
        package_root = os.path.dirname(os.path.dirname(__file__))  # .../probe_designer
        designer_root = os.path.dirname(package_root)               # .../designer
        possible_paths = [
            # Dev layout: designer/configs/
            os.path.join(designer_root, 'configs', 'species_config.json'),
            # Docker: /designer/configs/ (mounted read-only volume)
            os.path.join('/designer', 'configs', 'species_config.json'),
            # local/configs/ under designer root
            os.path.join(designer_root, 'local', 'configs', 'species_config.json'),
            # Relative to current working directory
            os.path.join('local', 'configs', 'species_config.json'),
            os.path.join('configs', 'species_config.json'),
        ]
        
        species_config_path = None
        for path in possible_paths:
            abs_path = os.path.abspath(path)
            if os.path.exists(abs_path):
                species_config_path = abs_path
                break
        
        if species_config_path:
            with open(species_config_path, 'r', encoding='utf-8') as f:
                species_data = json.load(f)
                self._species_data = species_data
        else:
            print(f"No species configuration file found. Tried: {possible_paths}")
            self._species_data = {"species": {}, "default_species": "mouse"}
    
    def set_species(self, species_name: str):
        """Set species-specific configuration."""
        if species_name in self._species_data["species"]:
            species_info = self._species_data["species"][species_name]
            self.species_config = SpeciesConfig(**species_info)
            
            # Update database and genome config
            self.database.organism = species_info["organism"]
            self.database.coord_system_version = species_info["coord_system_version"]
            self.genome.genome_fasta_path = species_info["genome_fasta_path"]
            
            print(f"Set species to: {species_info['display_name']} ({species_name})")
            print(f"  - Organism: {self.database.organism}")
            print(f"  - Genome FASTA: {self.genome.genome_fasta_path}")
        else:
            available_species = list(self._species_data["species"].keys())
            raise ValueError(f"Unknown species '{species_name}'. Available: {available_species}")
    
    def get_available_species(self) -> List[str]:
        """Get list of available species."""
        return list(self._species_data["species"].keys())
    
    def get_species_info(self, species_name: str) -> Optional[Dict]:
        """Get information for a specific species."""
        return self._species_data["species"].get(species_name)
    
    def load_config(self, config_file: str):
        """Load configuration from JSON or YAML file."""
        with open(config_file, 'r') as f:
            if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                config_data = yaml.safe_load(f)
            else:
                config_data = json.load(f)
        
        # Update configurations with loaded data
        # Note: database.organism may be overridden here, but species-specific settings
        # (like genome_fasta_path) should be set via set_species() after load_config
        if 'database' in config_data:
            for key, value in config_data['database'].items():
                if hasattr(self.database, key):
                    setattr(self.database, key, value)
        
        if 'search' in config_data:
            for key, value in config_data['search'].items():
                if hasattr(self.search, key):
                    setattr(self.search, key, value)
        
        if 'filter' in config_data:
            # Warn rather than silently discard: an unrecognised key is either a
            # typo (`formamide_pc`) or a knob that never existed, and both used
            # to vanish without a word. `final_probes_per_gene` lived in three
            # shipped YAMLs that way while the real control was --top-n (R7).
            unknown = [k for k in config_data['filter'] if not hasattr(self.filter, k)]
            if unknown:
                warnings.warn(
                    f"{config_file}: ignoring unknown filter option(s) "
                    f"{sorted(unknown)} — not a FilterConfig field. Check for a "
                    "typo; probes-per-gene is set with --top-n, not in the config.",
                    UserWarning,
                    stacklevel=2,
                )
            for key, value in config_data['filter'].items():
                if hasattr(self.filter, key):
                    setattr(self.filter, key, value)
        
        if 'blast' in config_data:
            for key, value in config_data['blast'].items():
                if hasattr(self.blast, key):
                    setattr(self.blast, key, value)
            
            # Handle blast species parameter - if not set, use species from config
            if 'species' not in config_data['blast'] or config_data['blast']['species'] is None:
                # Use species from species_config based on current species setting
                current_species = getattr(self, 'species_config', None)
                if current_species and hasattr(current_species, 'display_name'):
                    self.blast.species = [current_species.display_name]
        
        
        if 'output' in config_data:
            for key, value in config_data['output'].items():
                if hasattr(self.output, key):
                    setattr(self.output, key, value)
        
        if 'genome' in config_data:
            for key, value in config_data['genome'].items():
                if hasattr(self.genome, key):
                    setattr(self.genome, key, value)
        
        # Handle species from config file (but don't override if already set via parameter)
        # Note: species-specific settings like genome_fasta_path are handled by set_species()
        if 'species' in config_data:
            # Only set if not already set via constructor parameter
            # This is handled in __init__ after load_config
            pass
    
    def get_organism_gene_format(self, gene_name: str) -> str:
        """Format gene name based on organism conventions."""
        if self.database.organism == 'mouse':
            return gene_name.capitalize()
        elif self.database.organism == 'human':
            return gene_name.upper()
        else:
            return gene_name
    
    def save_config(self, output_file: str):
        """Save current configuration to YAML file."""
        config_data = {
            'database': {
                'organism': self.database.organism,
                'database_type': self.database.database_type,
                'coord_system_version': self.database.coord_system_version,
                'max_retries': self.database.max_retries
            },
            'search': {
                'search_strategy': self.search.search_strategy,
                'binding_site_length': self.search.binding_site_length,
                'max_binding_sites': self.search.max_binding_sites,
                'window_size': getattr(self.search, 'window_size', 50),
                'step_size': self.search.step_size
            },
            # Generated from the dataclass, never hand-listed: this block used
            # to mirror FilterConfig's fields by hand and silently dropped every
            # field added after it was written (enforce_tm_gate, the whole
            # RNAplfold block, max_alignments, target_organisms), so a
            # save -> load round trip lost settings (audit R7).
            'filter': dataclasses.asdict(self.filter),
            'blast': {
                'blast_type': self.blast.blast_type,
                'database': self.blast.database,
                'task': self.blast.task,
                'evalue': self.blast.evalue,
                'hitlist_size': self.blast.hitlist_size,
                'alignments': self.blast.alignments,
                'descriptions': self.blast.descriptions,
                'megablast': self.blast.megablast,
                'short_query': self.blast.short_query,
                'filter': self.blast.filter,
                'format_type': self.blast.format_type,
                'service': self.blast.service,
                'batch_size': self.blast.batch_size,
                'concurrency': self.blast.concurrency,
                'species': self.blast.species
            },
            'output': {
                'output_dir': self.output.output_dir,
                'save_fasta': getattr(self.output, 'save_fasta', True),
                'save_json': getattr(self.output, 'save_json', True),
                'save_excel': getattr(self.output, 'save_excel', True)
            },
            'genome': {
                'genome_fasta_path': self.genome.genome_fasta_path,
                'accessor_priority': self.genome.accessor_priority
            },
            'species': self.database.organism
        }
        
        with open(output_file, 'w', encoding='utf-8') as f:
            yaml.dump(config_data, f, default_flow_style=False, allow_unicode=True)
    
    def validate_config(self) -> List[str]:
        """Validate configuration values and return error list."""
        errors = []
        
        # search
        if self.search.binding_site_length <= 0:
            errors.append("binding_site_length must be > 0")
        if self.search.max_binding_sites <= 0:
            errors.append("max_binding_sites must be > 0")
        
        # filters: GC% gate (Schema-v2)
        if not (0 <= self.filter.min_gc_content <= 1):
            errors.append("min_gc_content must be between 0 and 1")
        if not (0 <= self.filter.max_gc_content <= 1):
            errors.append("max_gc_content must be between 0 and 1")
        if self.filter.min_gc_content >= self.filter.max_gc_content:
            errors.append("min_gc_content must be < max_gc_content")
        # filters: legacy G-only gate (kept for backward compat)
        if not (0 <= self.filter.min_g_content <= 1):
            errors.append("min_g_content must be between 0 and 1")
        if not (0 <= self.filter.max_g_content <= 1):
            errors.append("max_g_content must be between 0 and 1")
        if self.filter.min_g_content >= self.filter.max_g_content:
            errors.append("min_g_content must be < max_g_content")
        if self.filter.min_tm >= self.filter.max_tm:
            errors.append("min_tm must be < max_tm")
        
        return errors


# default config instance
default_config = ConfigManager()
