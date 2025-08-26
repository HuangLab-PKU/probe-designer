"""
Configuration management module
Manages all configurable parameters for DNA probe design.
"""

import os
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Any, Optional
import json


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
    # Pre-BLAST rules
    min_g_content: float = 0.4  # min G fraction
    max_g_content: float = 0.7  # max G fraction
    max_consecutive_g: int = 4  # max consecutive Gs allowed
    min_tm: float = 45.0  # min melting temperature
    max_tm: float = 65.0  # max melting temperature
    min_free_energy: float = -10.0  # min free energy (kcal/mol)
    
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
    task: str = "megablast"  # BLAST task
    evalue: float = 1e-5  # e-value threshold


@dataclass
class OutputConfig:
    """Output configuration."""
    output_dir: str = "results"  # output directory
    create_timestamp: bool = True  # create timestamped subdir
    save_intermediate: bool = True  # save intermediate files
    file_formats: List[str] = None  # output formats


class ConfigManager:
    """Config manager for the application."""
    
    def __init__(self, config_file: Optional[str] = None):
        self.database = DatabaseConfig()
        self.search = SearchConfig()
        self.filter = FilterConfig()
        self.probe = ProbeConfig()
        self.blast = BlastConfig()
        self.output = OutputConfig()
        
        # defaults
        if self.filter.target_organisms is None:
            self.filter.target_organisms = ["Mus musculus", "Homo sapiens"]
        
        if self.output.file_formats is None:
            self.output.file_formats = ["xlsx", "fasta", "json"]
        
        # load from file
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
    
    def load_config(self, config_file: str):
        """Load configuration from JSON file."""
        with open(config_file, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
        
        # update sections
        for section, data in config_data.items():
            if hasattr(self, section):
                section_obj = getattr(self, section)
                for key, value in data.items():
                    if hasattr(section_obj, key):
                        setattr(section_obj, key, value)
    
    def save_config(self, config_file: str):
        """Save configuration to JSON file."""
        config_data = {
            'database': self.database.__dict__,
            'search': self.search.__dict__,
            'filter': self.filter.__dict__,
            'probe': self.probe.__dict__,
            'blast': self.blast.__dict__,
            'output': self.output.__dict__
        }
        
        with open(config_file, 'w', encoding='utf-8') as f:
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
