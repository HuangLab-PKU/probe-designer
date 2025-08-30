"""
Configuration management module
Manages all configurable parameters for DNA probe design.
"""

import os
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Any, Optional
import json
import yaml


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
    max_binding_sites: int = 30  # max sites per gene
    search_strategy: str = "exon_junction"  # exon_junction, brute_force, isoform_specific
    step_size: Optional[int] = None  # step size for brute_force


@dataclass
class FilterConfig:
    """Sequence filtering configuration."""
    # Pre-BLAST rules (thermal properties)
    min_g_content: float = 0.4  # min G fraction
    max_g_content: float = 0.7  # max G fraction
    max_consecutive_g: int = 4  # max consecutive Gs allowed
    min_tm: float = 45.0  # min melting temperature
    max_tm: float = 65.0  # max melting temperature
    max_tm_diff: float = 10.0  # max Tm difference between 3' and 5' halves
    min_free_energy: float = -10.0  # min free energy (kcal/mol)
    check_rna_structure: bool = False  # whether to check RNA secondary structure
    
    # Post-BLAST rules
    max_alignments: int = 5  # max allowed alignments (specificity)
    require_specificity: bool = True  # enforce specificity in BLAST filter
    target_organisms: List[str] = None  # target organisms to allow


@dataclass
class ProbeConfig:
    """Probe assembly configuration."""
    panel_type: str = "PRISM"  # PRISM, SPRINTseq, custom
    barcode_file: Optional[str] = None  # barcode Excel file path
    primer_left: Optional[str] = None  # left primer sequence
    primer_right: Optional[str] = None  # right primer sequence


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
    file_formats: List[str] = None  # output formats
    save_fasta: bool = True
    save_json: bool = True
    save_excel: bool = True


@dataclass
class GenomeConfig:
    """Genome access configuration."""
    genome_fasta_path: Optional[str] = None  # path to local genome FASTA
    use_local_first: bool = True  # prefer local over online access


@dataclass
class SpeciesConfig:
    """Species-specific configuration."""
    organism: str = "mouse"
    coord_system_version: str = "GRCm39"
    genome_fasta_path: str = ""
    display_name: str = "Mus musculus"
    taxonomy_id: int = 10090


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
        self.species_config = SpeciesConfig()
        
        # defaults
        if self.filter.target_organisms is None:
            self.filter.target_organisms = ["Mus musculus", "Homo sapiens"]
        
        if self.output.file_formats is None:
            self.output.file_formats = ["xlsx", "fasta", "json"]
        
        # load species configuration first
        self._load_species_config()
        
        # load main configuration
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
        
        # apply species-specific settings
        if species:
            self.set_species(species)
    
    def _load_species_config(self):
        """Load species configuration from species_config.json."""
        species_config_path = os.path.join(os.path.dirname(__file__), '..', 'configs', 'species_config.json')
        if os.path.exists(species_config_path):
            with open(species_config_path, 'r', encoding='utf-8') as f:
                species_data = json.load(f)
                self._species_data = species_data
        else:
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
        if 'database' in config_data:
            for key, value in config_data['database'].items():
                if hasattr(self.database, key):
                    setattr(self.database, key, value)
        
        if 'search' in config_data:
            for key, value in config_data['search'].items():
                if hasattr(self.search, key):
                    setattr(self.search, key, value)
        
        if 'filter' in config_data:
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
                else:
                    # Default fallback
                    self.blast.species = ["Mus musculus"]
        
        if 'output' in config_data:
            for key, value in config_data['output'].items():
                if hasattr(self.output, key):
                    setattr(self.output, key, value)
        
        if 'genome' in config_data:
            for key, value in config_data['genome'].items():
                if hasattr(self.genome, key):
                    setattr(self.genome, key, value)
    
    def save_config(self, config_file: str):
        """Save configuration to JSON or YAML file."""
        config_data = {
            'database': self.database.__dict__,
            'search': self.search.__dict__,
            'filter': self.filter.__dict__,
            'probe': self.probe.__dict__,
            'blast': self.blast.__dict__,
            'output': self.output.__dict__,
            'genome': self.genome.__dict__
        }
        
        with open(config_file, 'w', encoding='utf-8') as f:
            if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                yaml.dump(config_data, f, default_flow_style=False, allow_unicode=True)
            else:
                json.dump(config_data, f, indent=2, ensure_ascii=False)
    
    def get_organism_gene_format(self, gene_name: str) -> str:
        """Format gene name based on organism conventions."""
        if self.database.organism == 'mouse':
            return gene_name.capitalize()
        elif self.database.organism == 'human':
            return gene_name.upper()
        else:
            return gene_name
    
    def validate_config(self) -> List[str]:
        """Validate configuration values and return error list."""
        errors = []
        
        # search
        if self.search.binding_site_length <= 0:
            errors.append("binding_site_length must be > 0")
        if self.search.max_binding_sites <= 0:
            errors.append("max_binding_sites must be > 0")
        
        # filters
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
