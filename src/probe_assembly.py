"""
Probe assembly module
Simplified probe assembly from binding sites with backbones.
"""

import os
import json
import pandas as pd
from typing import List, Dict, Any, Optional, Union
from .config import ProbeConfig


class ProbeAssembler:
    """Assemble probes from binding sites with backbones."""
    
    def __init__(self, probe_config: ProbeConfig):
        self.probe_config = probe_config
        self.backbone_data = None
        if probe_config.backbone_file:
            self._load_backbone_data()
    
    def _load_backbone_data(self):
        """Load backbone data from Excel file."""
        if not os.path.exists(self.probe_config.backbone_file):
            raise FileNotFoundError(f"Backbone file not found: {self.probe_config.backbone_file}")
        
        # Read Excel file, expecting 'No.' and 'Sequence' columns
        df = pd.read_excel(self.probe_config.backbone_file)
        
        # Validate required columns
        if 'No.' not in df.columns or 'Sequence' not in df.columns:
            raise ValueError(f"Backbone file must contain 'No.' and 'Sequence' columns. Found: {df.columns.tolist()}")
        
        # Create mapping from No. to Sequence (convert No. to string for consistency)
        self.backbone_data = {str(no): seq for no, seq in zip(df['No.'], df['Sequence'])}
        
        # Validate no duplicate No.
        if len(self.backbone_data) != len(df):
            raise ValueError(f"Duplicate 'No.' values found in backbone file")
    
    def assemble_probes(
        self,
        binding_sites: Dict[str, List[Dict[str, Any]]],
        gene_info_file: str
    ) -> pd.DataFrame:
        """
        Assemble probes from binding sites using gene/bds info from Excel file.
        
        Args:
            binding_sites: Dict mapping gene_name to list of binding site dicts.
                          Each site must have: arm_3prime, arm_5prime, st, en, g_content, tm
            gene_info_file: Path to Excel file with 'gene_name' and 'No.' columns.
                          - Per-gene mode: one row per gene (same No. for all sites of that gene)
                          - Per-bds mode: one row per binding site (different No. for each site)
        
        Returns:
            DataFrame with columns: gene_name, probe_id, probe_seq, arm_5prime, arm_3prime,
                                   backbone, st, en, g_content, tm
        """
        if self.backbone_data is None:
            raise ValueError("Backbone data not loaded. Please provide backbone_file in ProbeConfig.")
        
        # Load gene_info from Excel file
        if not os.path.exists(gene_info_file):
            raise FileNotFoundError(f"Gene info file not found: {gene_info_file}")
        
        gene_info_df = pd.read_excel(gene_info_file)
        
        # Validate required columns
        if 'gene_name' not in gene_info_df.columns or 'No.' not in gene_info_df.columns:
            raise ValueError(f"Gene info file must contain 'gene_name' and 'No.' columns. Found: {gene_info_df.columns.tolist()}")
        
        # Normalize gene names: strip whitespace and convert to string
        gene_info_df['gene_name'] = gene_info_df['gene_name'].astype(str).str.strip()
        
        # Determine assembly mode: check if each gene has multiple rows
        gene_counts = gene_info_df['gene_name'].value_counts()
        is_per_bds_mode = (gene_counts > 1).any()
        
        # Build mapping from gene_name to No.(s) - use normalized names
        if is_per_bds_mode:
            # Per-bds mode: group by gene_name and create list of No.s
            gene_info = {}
            for gene_name, group in gene_info_df.groupby('gene_name'):
                gene_info[gene_name] = group['No.'].tolist()
        else:
            # Per-gene mode: one No. per gene
            gene_info = dict(zip(gene_info_df['gene_name'], gene_info_df['No.']))
        
        # Create a case-insensitive lookup map for better matching
        gene_info_lower = {k.strip().upper(): (k, v) for k, v in gene_info.items()}
        
        probe_records = []
        missing_genes = []
        
        for gene_name, sites in binding_sites.items():
            if not sites:
                continue
            
            # Normalize gene name for matching
            gene_name_normalized = gene_name.strip()
            gene_name_upper = gene_name_normalized.upper()
            
            # Try exact match first
            if gene_name_normalized in gene_info:
                probe_ids = gene_info[gene_name_normalized]
            # Try case-insensitive match
            elif gene_name_upper in gene_info_lower:
                original_name, probe_ids = gene_info_lower[gene_name_upper]
                print(f"警告: 基因名大小写不匹配，使用 '{original_name}' 匹配 '{gene_name}'")
            else:
                missing_genes.append(gene_name)
                continue
            
            probe_ids = gene_info[gene_name]
            
            # Check if it's per-gene or per-site mode
            if isinstance(probe_ids, (str, int, float)):
                # Per-gene mode: one No. for all sites
                probe_id = str(probe_ids)  # Convert to string for consistency
                if probe_id not in self.backbone_data:
                    raise ValueError(f"Probe ID '{probe_id}' not found in backbone file. Available IDs: {list(self.backbone_data.keys())}")
                
                backbone = self.backbone_data[probe_id]
                
                # Use same backbone for all sites of this gene
                for site in sites:
                    probe_seq = self._assemble_probe_sequence(
                        site['arm_5prime'].lower(),
                        site['arm_3prime'].lower(),
                        backbone
                    )
                    record = {
                        'gene_name': gene_name,
                        'probe_id': probe_id,
                        'probe_seq': probe_seq,
                        'arm_5prime': site['arm_5prime'],
                        'arm_3prime': site['arm_3prime'],
                        'backbone': backbone,
                        'st': site['st'],
                        'en': site['en'],
                        'g_content': site['g_content'],
                        'tm': site['tm'],
                        'tm_3prime': site.get('tm_3prime'),
                        'tm_5prime': site.get('tm_5prime')
                    }
                    if 'isoform_overlap_num' in site:
                        record['isoform_overlap_num'] = site['isoform_overlap_num']
                    probe_records.append(record)
            
            elif isinstance(probe_ids, list):
                # Per-site mode: one No. per binding site
                if len(probe_ids) != len(sites):
                    raise ValueError(f"Number of probe IDs ({len(probe_ids)}) does not match number of binding sites ({len(sites)}) for gene '{gene_name}'")
                
                for site, probe_id in zip(sites, probe_ids):
                    probe_id = str(probe_id)  # Convert to string for consistency
                    if probe_id not in self.backbone_data:
                        raise ValueError(f"Probe ID '{probe_id}' not found in backbone file. Available IDs: {list(self.backbone_data.keys())}")
                    
                    backbone = self.backbone_data[probe_id]
                    probe_seq = self._assemble_probe_sequence(
                        site['arm_5prime'],
                        site['arm_3prime'],
                        backbone
                    )
                    record = {
                        'gene_name': gene_name,
                        'probe_id': probe_id,
                        'probe_seq': probe_seq,
                        'arm_5prime': site['arm_5prime'],
                        'arm_3prime': site['arm_3prime'],
                        'backbone': backbone,
                        'st': site['st'],
                        'en': site['en'],
                        'g_content': site['g_content'],
                        'tm': site['tm'],
                        'tm_3prime': site.get('tm_3prime'),
                        'tm_5prime': site.get('tm_5prime')
                    }
                    if 'isoform_overlap_num' in site:
                        record['isoform_overlap_num'] = site['isoform_overlap_num']
                    probe_records.append(record)
            else:
                raise ValueError(f"Invalid gene_info format for gene '{gene_name}'. Expected str/int or list, got {type(probe_ids)}")
        
        # Report missing genes if any
        if missing_genes:
            available_genes = sorted(set(gene_info_df['gene_name'].tolist()))
            error_msg = (
                f"以下 {len(missing_genes)} 个基因在 gene_info 文件中未找到:\n"
                f"  缺失的基因: {', '.join(missing_genes)}\n"
                f"  gene_info 文件中的基因: {', '.join(available_genes[:20])}"
            )
            if len(available_genes) > 20:
                error_msg += f" ... (共 {len(available_genes)} 个)"
            raise ValueError(error_msg)
        
        return pd.DataFrame(probe_records)
    
    def _assemble_probe_sequence(self, arm_5prime: str, arm_3prime: str, backbone: str) -> str:
        """
        Assemble probe sequence: arm_5prime + backbone + arm_3prime
        
        Args:
            arm_5prime: 5' arm sequence
            arm_3prime: 3' arm sequence
            backbone: Backbone sequence
        
        Returns:
            Complete probe sequence
        """
        return arm_5prime + backbone + arm_3prime
    
    def add_ilock_modification(self, probe_df: pd.DataFrame, position: str = "3prime") -> pd.DataFrame:
        """
        Add iLock modification to probes.
        
        Args:
            probe_df: DataFrame with 'probe_seq' column
            position: Where to add iLock - "3prime", "5prime", or "both"
        
        Returns:
            DataFrame with modified 'probe_seq' and new 'modification' column
        """
        # iLock sequence (example - adjust as needed)
        ilock_seq = "CCTTCCTTCCTTCCTTCCTT"  # Typical iLock sequence
        
        probe_df = probe_df.copy()
        probe_df['modification'] = 'iLock'
        
        if position == "3prime":
            probe_df['probe_seq'] = probe_df['probe_seq'] + ilock_seq
        elif position == "5prime":
            probe_df['probe_seq'] = ilock_seq + probe_df['probe_seq']
        elif position == "both":
            probe_df['probe_seq'] = ilock_seq + probe_df['probe_seq'] + ilock_seq
        else:
            raise ValueError(f"Invalid position '{position}'. Must be '3prime', '5prime', or 'both'")
        
        return probe_df
    
    def save_probes(self, probe_df: pd.DataFrame, output_dir: str):
        """Save probes to Excel, FASTA and JSON, plus stats."""
        os.makedirs(output_dir, exist_ok=True)
        
        excel_file = os.path.join(output_dir, "probes.xlsx")
        probe_df.to_excel(excel_file, index=False)
        
        fasta_file = os.path.join(output_dir, "probes.fasta")
        with open(fasta_file, 'w') as f:
            for _, row in probe_df.iterrows():
                gene = row['gene_name']
                probe_id = row['probe_id']
                probe = row['probe_seq']
                f.write(f">{gene}|{probe_id}\n")
                f.write(f"{probe}\n")
        
        json_file = os.path.join(output_dir, "probes.json")
        probe_dict = probe_df.to_dict('records')
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(probe_dict, f, indent=2, ensure_ascii=False)
        
        stats = {
            'total_probes': len(probe_df),
            'unique_genes': probe_df['gene_name'].nunique(),
            'unique_probe_ids': probe_df['probe_id'].nunique(),
            'avg_probe_length': probe_df['probe_seq'].str.len().mean(),
            'avg_g_content': probe_df['g_content'].mean(),
            'avg_tm': probe_df['tm'].mean()
        }
        stats_file = os.path.join(output_dir, "probe_stats.json")
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(stats, f, indent=2, ensure_ascii=False)
        
        return probe_df
