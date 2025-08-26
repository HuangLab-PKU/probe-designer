"""
Probe assembly module
Implements probe stitching for different panel types.
"""

import os
import json
import pandas as pd
from typing import List, Dict, Any, Optional
from tqdm import tqdm

from .config import ProbeConfig


class ProbeAssembler:
    """Assemble probes from binding sites for various panels."""
    
    def __init__(self, probe_config: ProbeConfig):
        self.probe_config = probe_config
        self.barcode_data = None
        self._load_barcode_data()
    
    def _load_barcode_data(self):
        """Load barcode data from file or fallback to defaults."""
        if self.probe_config.barcode_file and os.path.exists(self.probe_config.barcode_file):
            self.barcode_data = pd.read_excel(self.probe_config.barcode_file, index_col=0)
        else:
            self.barcode_data = self._get_default_barcodes()
    
    def _get_default_barcodes(self) -> pd.DataFrame:
        """Provide placeholder default barcodes when none specified."""
        if self.probe_config.panel_type == "PRISM":
            barcodes = {}
            for i in range(1, 32):
                barcode_seq = f"Barcode_{i:02d}"
                barcodes[f"PRISM_{i}"] = {
                    'Sequence': barcode_seq,
                    'Barcode (82bp)': barcode_seq
                }
            return pd.DataFrame.from_dict(barcodes, orient='index')
        elif self.probe_config.panel_type == "SPRINTseq":
            barcodes = {}
            for i in range(1, 370):
                barcode_seq = f"Barcode_{i:03d}"
                barcodes[f"SPRINTseq_{i}"] = {
                    'Barcode sequence': barcode_seq,
                    'Barcode(70bp)': self._create_sprintseq_probe(barcode_seq)
                }
            return pd.DataFrame.from_dict(barcodes, orient='index')
        else:
            return pd.DataFrame()
    
    def _create_sprintseq_probe(self, barcode_seq: str) -> str:
        """Build SPRINTseq composite barcode with primers."""
        primer_l = self.probe_config.primer_left or 'TCCCTACACGACGCTCTTCCGATCT'
        primer_r = self.probe_config.primer_right or 'CATTCCTGCTGAACCGCTCTTCCGA'
        return primer_l + barcode_seq + primer_r + barcode_seq
    
    def assemble_probes(self, binding_sites: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
        """Assemble probes from binding sites according to panel type."""
        probe_df = pd.DataFrame()
        if self.probe_config.panel_type == "PRISM":
            probe_df = self._assemble_prism_probes(binding_sites)
        elif self.probe_config.panel_type == "SPRINTseq":
            probe_df = self._assemble_sprintseq_probes(binding_sites)
        else:
            probe_df = self._assemble_custom_probes(binding_sites)
        return probe_df
    
    def _assemble_prism_probes(self, binding_sites: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
        """Assemble PRISM probes: right arm + BARCODE + left arm."""
        probe_df = pd.DataFrame()
        available_barcodes = list(self.barcode_data.index)
        for i, (gene_name, sites) in enumerate(binding_sites.items()):
            if not sites:
                continue
            site = sites[0]
            binding = site['sequence']
            binding_l = binding[:20].lower()
            binding_r = binding[20:].lower()
            if i < len(available_barcodes):
                barcode_name = available_barcodes[i]
                barcode = self.barcode_data.loc[barcode_name, "Barcode (82bp)"]
            else:
                barcode_name = available_barcodes[i % len(available_barcodes)]
                barcode = self.barcode_data.loc[barcode_name, "Barcode (82bp)"]
            probe = binding_r + barcode.upper() + binding_l
            probe_info = pd.DataFrame({
                "PRISM": [barcode_name],
                "gene": [gene_name],
                "probe": [probe],
                "barcode": [barcode],
                "binding": [binding],
                "binding_left": [binding_l],
                "binding_right": [binding_r],
                "position": [site['position']],
                "g_content": [site['g_content']],
                "tm": [site['tm']]
            })
            if len(probe_df) == 0:
                probe_df = probe_info
            else:
                probe_df = pd.concat([probe_df, probe_info], ignore_index=True)
        return probe_df.set_index('PRISM')
    
    def _assemble_sprintseq_probes(self, binding_sites: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
        """Assemble SPRINTseq probes: binding site + composite barcode."""
        probe_df = pd.DataFrame()
        available_barcodes = list(self.barcode_data.index)
        for i, (gene_name, sites) in enumerate(binding_sites.items()):
            if not sites:
                continue
            site = sites[0]
            binding = site['sequence']
            if i < len(available_barcodes):
                barcode_name = available_barcodes[i]
                barcode = self.barcode_data.loc[barcode_name, "Barcode(70bp)"]
            else:
                barcode_name = available_barcodes[i % len(available_barcodes)]
                barcode = self.barcode_data.loc[barcode_name, "Barcode(70bp)"]
            probe = binding + barcode
            probe_info = pd.DataFrame({
                "SPRINTseq": [barcode_name],
                "gene": [gene_name],
                "probe": [probe],
                "barcode": [barcode],
                "binding": [binding],
                "position": [site['position']],
                "g_content": [site['g_content']],
                "tm": [site['tm']]
            })
            if len(probe_df) == 0:
                probe_df = probe_info
            else:
                probe_df = pd.concat([probe_df, probe_info], ignore_index=True)
        return probe_df.set_index('SPRINTseq')
    
    def _assemble_custom_probes(self, binding_sites: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
        """Assemble custom probes using raw binding sequences."""
        probe_df = pd.DataFrame()
        for gene_name, sites in binding_sites.items():
            if not sites:
                continue
            site = sites[0]
            binding = site['sequence']
            probe = binding
            probe_info = pd.DataFrame({
                "gene": [gene_name],
                "probe": [probe],
                "binding": [binding],
                "position": [site['position']],
                "g_content": [site['g_content']],
                "tm": [site['tm']]
            })
            if len(probe_df) == 0:
                probe_df = probe_info
            else:
                probe_df = pd.concat([probe_df, probe_info], ignore_index=True)
        return probe_df
    
    def save_probes(self, probe_df: pd.DataFrame, output_dir: str):
        """Save probes to Excel, FASTA and JSON, plus stats."""
        os.makedirs(output_dir, exist_ok=True)
        excel_file = os.path.join(output_dir, f"{self.probe_config.panel_type}_probes.xlsx")
        probe_df.to_excel(excel_file)
        fasta_file = os.path.join(output_dir, f"{self.probe_config.panel_type}_probes.fasta")
        with open(fasta_file, 'w') as f:
            for idx, row in probe_df.iterrows():
                gene = row['gene']
                probe = row['probe']
                f.write(f">{gene}|{idx}\n")
                f.write(f"{probe}\n")
        json_file = os.path.join(output_dir, f"{self.probe_config.panel_type}_probes.json")
        probe_dict = probe_df.to_dict('index')
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(probe_dict, f, indent=2, ensure_ascii=False)
        stats = {
            'panel_type': self.probe_config.panel_type,
            'total_probes': len(probe_df),
            'unique_genes': probe_df['gene'].nunique(),
            'avg_probe_length': probe_df['probe'].str.len().mean(),
            'avg_g_content': probe_df['g_content'].mean(),
            'avg_tm': probe_df['tm'].mean()
        }
        stats_file = os.path.join(output_dir, f"{self.probe_config.panel_type}_stats.json")
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(stats, f, indent=2, ensure_ascii=False)
        return probe_df


class ProbeValidator:
    """Simple probe sequence and set validator."""
    
    def __init__(self):
        pass
    
    def validate_probe_sequence(self, probe_sequence: str) -> Dict[str, Any]:
        """Validate a single probe sequence with basic checks."""
        validation_result = {
            'is_valid': True,
            'errors': [],
            'warnings': []
        }
        if len(probe_sequence) < 50:
            validation_result['warnings'].append("Probe sequence is short")
        if len(probe_sequence) > 200:
            validation_result['warnings'].append("Probe sequence is long")
        valid_bases = set('ATCGN')
        sequence_bases = set(probe_sequence.upper())
        if not sequence_bases.issubset(valid_bases):
            invalid_bases = sequence_bases - valid_bases
            validation_result['errors'].append(f"Invalid bases: {invalid_bases}")
            validation_result['is_valid'] = False
        gc_content = (probe_sequence.upper().count('G') + probe_sequence.upper().count('C')) / len(probe_sequence)
        if gc_content < 0.3 or gc_content > 0.7:
            validation_result['warnings'].append(f"Abnormal GC content: {gc_content:.2f}")
        if self._has_long_repeats(probe_sequence):
            validation_result['warnings'].append("Contains long repeats")
        return validation_result
    
    def _has_long_repeats(self, sequence: str, min_repeat_length: int = 8) -> bool:
        """Detect if there are any long repeated substrings."""
        for i in range(len(sequence) - min_repeat_length + 1):
            repeat = sequence[i:i + min_repeat_length]
            if sequence.count(repeat) > 1:
                return True
        return False
    
    def validate_probe_set(self, probe_df: pd.DataFrame) -> Dict[str, Any]:
        """Validate a set of probes and return a report."""
        validation_result = {
            'is_valid': True,
            'errors': [],
            'warnings': [],
            'gene_coverage': {},
            'sequence_quality': {}
        }
        gene_counts = probe_df['gene'].value_counts()
        validation_result['gene_coverage'] = {
            'total_genes': len(gene_counts),
            'genes_with_probes': len(gene_counts[gene_counts > 0]),
            'avg_probes_per_gene': gene_counts.mean()
        }
        for idx, row in probe_df.iterrows():
            probe_sequence = row['probe']
            gene = row['gene']
            seq_validation = self.validate_probe_sequence(probe_sequence)
            validation_result['sequence_quality'][gene] = seq_validation
            if not seq_validation['is_valid']:
                validation_result['is_valid'] = False
                validation_result['errors'].extend([f"{gene}: {error}" for error in seq_validation['errors']])
            validation_result['warnings'].extend([f"{gene}: {warning}" for warning in seq_validation['warnings']])
        return validation_result
