"""
Utility functions for the probe design system.
"""

import os
import json
import time
import random
import string
from typing import Dict, List, Any, Optional
from pathlib import Path
import pandas as pd
from tqdm import tqdm

# Global progress bar manager to avoid conflicts
class ProgressManager:
    """Manage progress bars to avoid conflicts and duplicate displays."""
    
    def __init__(self):
        self.current_position = 0
        self.active_bars = {}
    
    def create_bar(self, desc: str, total: int = None, leave: bool = True, position: int = None) -> tqdm:
        """Create a progress bar with proper positioning."""
        if position is None:
            position = self.current_position
            self.current_position += 1
        
        bar = tqdm(
            total=total,
            desc=desc,
            leave=leave,
            position=position,
            ncols=80,
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]'
        )
        
        if not leave:
            self.active_bars[position] = bar
        
        return bar
    
    def close_bar(self, position: int):
        """Close a specific progress bar."""
        if position in self.active_bars:
            self.active_bars[position].close()
            del self.active_bars[position]

# Global instance
progress_manager = ProgressManager()

def create_progress_bar(desc: str, total: int = None, leave: bool = True, position: int = None) -> tqdm:
    """Create a progress bar using the global manager."""
    return progress_manager.create_bar(desc, total, leave, position)

def close_progress_bar(position: int):
    """Close a progress bar using the global manager."""
    progress_manager.close_bar(position)


def create_output_directory(base_dir: str, create_timestamp: bool = True) -> str:
    """Create the output directory, optionally with a timestamp subdir."""
    if create_timestamp:
        timestamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        output_dir = os.path.join(base_dir, timestamp)
    else:
        output_dir = base_dir
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def load_gene_list(file_path: str, gene_column: str = "gene_name") -> List[str]:
    """Load gene list from .xlsx/.csv/.txt and return a clean list of strings."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Gene list file not found: {file_path}")
    file_ext = Path(file_path).suffix.lower()
    if file_ext in ('.xlsx', '.xls'):
        df = pd.read_excel(file_path)
    elif file_ext == '.csv':
        df = pd.read_csv(file_path)
    elif file_ext == '.txt':
        df = pd.read_csv(file_path, sep='\t', header=None, names=[gene_column])
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")
    if gene_column not in df.columns:
        available_columns = list(df.columns)
        raise ValueError(f"Column '{gene_column}' not found. Available: {available_columns}")
    gene_list = []
    for gene in df[gene_column].dropna():
        gene_str = str(gene).strip()
        if gene_str and gene_str != '0' and gene_str != 'nan':
            gene_list.append(gene_str)
    return gene_list


def save_gene_list(gene_list: List[str], output_file: str):
    """Save gene list to file by extension."""
    df = pd.DataFrame({'gene_name': gene_list})
    file_ext = Path(output_file).suffix.lower()
    if file_ext == '.xlsx':
        df.to_excel(output_file, index=False)
    elif file_ext == '.csv':
        df.to_csv(output_file, index=False)
    elif file_ext == '.txt':
        df.to_csv(output_file, sep='\t', index=False, header=False)
    else:
        df.to_csv(output_file, sep='\t', index=False, header=False)


def format_gene_name(gene_name: str, organism: str) -> str:
    """Format gene name according to organism conventions."""
    if organism == 'mouse':
        return gene_name.capitalize()
    elif organism == 'human':
        return gene_name.upper()
    else:
        return gene_name


def adjust_gene_name(gene_name: str, gene_list: List[str]) -> str:
    """Adjust gene name to match list (e.g., strip numeric suffix like -1, -2)."""
    gene_list_upper = [x.upper() for x in gene_list]
    match = re.search(r'(.+)-(\d+)$', gene_name)
    if match:
        base_gene_name = match.group(1)
        if base_gene_name.upper() in gene_list_upper:
            return base_gene_name
        elif gene_name.upper() in gene_list_upper:
            return gene_name
    return gene_name


def merge_results_from_directories(results_dir: str, output_file: str) -> pd.DataFrame:
    """Merge probes_wanted.xlsx across subdirectories under results_dir."""
    all_results = []
    for dir_name in os.listdir(results_dir):
        dir_path = os.path.join(results_dir, dir_name)
        if not os.path.isdir(dir_path):
            continue
        probes_file = os.path.join(dir_path, "probes_wanted.xlsx")
        if os.path.exists(probes_file):
            try:
                # Read without forcing an index so that 'gene_name' remains a column
                df = pd.read_excel(probes_file)
                # If previous files saved with 'gene_name' as index, bring it back as a column
                if 'gene_name' not in df.columns and getattr(df.index, 'name', None) == 'gene_name':
                    df = df.reset_index()
                all_results.append(df)
            except Exception as e:
                print(f"Failed to read {probes_file}: {e}")
    if not all_results:
        raise ValueError("No result files found")
    merged_df = pd.concat(all_results, ignore_index=True)
    merged_df.drop_duplicates(subset=["bds"], keep="first", inplace=True)
    merged_df.sort_values(["gene_name", "pos"], inplace=True)
    merged_df.to_excel(output_file)
    return merged_df


def find_missing_genes(gene_list: List[str], results_df: pd.DataFrame, gene_column: str = "gene_name") -> List[str]:
    """Return genes not present in results_df[gene_column] (case-insensitive)."""
    found_genes = set(results_df[gene_column].str.upper())
    missing_genes = []
    for gene in gene_list:
        if gene.upper() not in found_genes:
            missing_genes.append(gene)
    return missing_genes


def save_missing_genes(missing_genes: List[str], output_file: str):
    """Write missing genes to a text file (one per line)."""
    with open(output_file, 'w', encoding='utf-8') as f:
        for gene in missing_genes:
            f.write(f"{gene}\n")


def calculate_statistics(binding_sites: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Any]:
    """Compute aggregate statistics for binding sites dict."""
    total_genes = len(binding_sites)
    total_sites = sum(len(sites) for sites in binding_sites.values())
    genes_with_sites = len([g for g, sites in binding_sites.items() if sites])
    all_g_contents = []
    all_tms = []
    for sites in binding_sites.values():
        for site in sites:
            all_g_contents.append(site.get('g_content', 0))
            all_tms.append(site.get('tm', 0))
    stats = {
        'total_genes': total_genes,
        'total_binding_sites': total_sites,
        'genes_with_sites': genes_with_sites,
        'avg_sites_per_gene': total_sites / total_genes if total_genes > 0 else 0,
        'avg_g_content': sum(all_g_contents) / len(all_g_contents) if all_g_contents else 0,
        'avg_tm': sum(all_tms) / len(all_tms) if all_tms else 0,
        'min_g_content': min(all_g_contents) if all_g_contents else 0,
        'max_g_content': max(all_g_contents) if all_g_contents else 0,
        'min_tm': min(all_tms) if all_tms else 0,
        'max_tm': max(all_tms) if all_tms else 0
    }
    return stats


def save_statistics(stats: Dict[str, Any], output_file: str):
    """Save statistics dict to JSON file."""
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)


def validate_file_path(file_path: str, file_type: str = "file") -> bool:
    """Validate that file or directory exists based on file_type."""
    if file_type == "file":
        return os.path.isfile(file_path)
    elif file_type == "directory":
        return os.path.isdir(file_path)
    else:
        return os.path.exists(file_path)


def create_log_file(output_dir: str, log_name: str = "probe_design.log") -> str:
    """Return log file path inside output_dir."""
    log_file = os.path.join(output_dir, log_name)
    return log_file


def log_message(message: str, log_file: Optional[str] = None, print_to_console: bool = True):
    """Append a timestamped log message to file and/or console."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    log_entry = f"[{timestamp}] {message}"
    if print_to_console:
        print(log_entry)
    if log_file:
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(log_entry + '\n')


def check_dependencies() -> Dict[str, bool]:
    """Check presence of common Python dependencies used by this project."""
    dependencies = {
        'pandas': False,
        'numpy': False,
        'biopython': False,
        'requests': False,
        'tqdm': False
    }
    try:
        import pandas  # noqa: F401
        dependencies['pandas'] = True
    except ImportError:
        pass
    try:
        import numpy  # noqa: F401
        dependencies['numpy'] = True
    except ImportError:
        pass
    try:
        import Bio  # noqa: F401
        dependencies['biopython'] = True
    except ImportError:
        pass
    try:
        import requests  # noqa: F401
        dependencies['requests'] = True
    except ImportError:
        pass
    try:
        import tqdm  # noqa: F401
        dependencies['tqdm'] = True
    except ImportError:
        pass
    return dependencies


def print_dependency_status():
    """Print dependency status to console."""
    deps = check_dependencies()
    print("Dependency check:")
    for package, installed in deps.items():
        status = "✓" if installed else "✗"
        print(f"  {package}: {status}")
    missing = [pkg for pkg, installed in deps.items() if not installed]
    if missing:
        print(f"\nMissing: {', '.join(missing)}")
        print("Install via:")
        print(f"pip install {' '.join(missing)}")


def create_config_template(output_file: str):
    """Create a JSON config template on disk."""
    config_template = {
        "database": {
            "organism": "mouse",
            "database_type": "ensembl",
            "coord_system_version": "GRCh38",
            "max_retries": 3
        },
        "search": {
            "binding_site_length": 40,
            "max_binding_sites": 30,
            "search_strategy": "exon_junction",
            "step_size": None
        },
        "filter": {
            "min_g_content": 0.4,
            "max_g_content": 0.7,
            "max_consecutive_g": 4,
            "min_tm": 45.0,
            "max_tm": 65.0,
            "min_free_energy": -10.0,
            "max_alignments": 5,
            "require_specificity": True,
            "target_organisms": ["Mus musculus", "Homo sapiens"]
        },
        "probe": {
            "panel_type": "PRISM",
            "barcode_file": None,
            "primer_left": None,
            "primer_right": None
        },
        "blast": {
            "blast_type": "local",
            "database": "refseq_rna",
            "task": "megablast",
            "evalue": 1e-5
        },
        "output": {
            "output_dir": "results",
            "create_timestamp": True,
            "save_intermediate": True,
            "file_formats": ["xlsx", "fasta", "json"]
        }
    }
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(config_template, f, indent=2, ensure_ascii=False)
    print(f"Config template saved: {output_file}")


def get_file_size_mb(file_path: str) -> float:
    """Get file size in MB (0.0 if missing)."""
    if os.path.exists(file_path):
        return os.path.getsize(file_path) / (1024 * 1024)
    return 0.0


def format_file_size(size_bytes: int) -> str:
    """Human-readable file size string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"
