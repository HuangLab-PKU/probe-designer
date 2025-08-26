"""
Sequence filtering module
Implements pre- and post-BLAST filtering rules.
"""

import os
import json
import pandas as pd
from typing import List, Dict, Any, Optional
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqUtils import MeltingTemp as mt
from tqdm import tqdm

from .config import FilterConfig, BlastConfig


class SequenceFilter:
    """Sequence filter for candidate binding sites."""
    
    def __init__(self, filter_config: FilterConfig, blast_config: BlastConfig):
        self.filter_config = filter_config
        self.blast_config = blast_config
    
    def pre_blast_filter(self, binding_sites: Dict[str, List[Dict[str, Any]]]) -> Dict[str, List[Dict[str, Any]]]:
        """Filter candidates before running BLAST."""
        filtered_sites = {}
        
        for gene_name, sites in tqdm(binding_sites.items(), desc="Pre-BLAST filtering"):
            filtered_sites[gene_name] = []
            for site in sites:
                if self._check_pre_blast_criteria(site):
                    filtered_sites[gene_name].append(site)
        
        return filtered_sites
    
    def _check_pre_blast_criteria(self, site: Dict[str, Any]) -> bool:
        """Check pre-BLAST criteria for a single site."""
        sequence = site['sequence']
        
        # G content
        g_content = sequence.count('G') / len(sequence)
        if not (self.filter_config.min_g_content <= g_content <= self.filter_config.max_g_content):
            return False
        
        # consecutive Gs
        if "G" * (self.filter_config.max_consecutive_g + 1) in sequence:
            return False
        
        # Tm
        try:
            tm = mt.Tm_NN(sequence, nn_table=mt.DNA_NN4)
            if not (self.filter_config.min_tm <= tm <= self.filter_config.max_tm):
                return False
        except:
            return False
        
        # RNA secondary structure (optional)
        if hasattr(self, '_check_rna_structure'):
            if not self._check_rna_structure(sequence):
                return False
        
        return True
    
    def run_blast(self, binding_sites: Dict[str, List[Dict[str, Any]]], output_dir: str) -> str:
        """Run BLAST analysis and return the path to XML results."""
        fasta_file = os.path.join(output_dir, "blast_input.fasta")
        self._prepare_blast_fasta(binding_sites, fasta_file)
        
        blast_results_file = os.path.join(output_dir, "blast_results.xml")
        
        if self.blast_config.blast_type == "local":
            self._run_local_blast(fasta_file, blast_results_file)
        else:
            self._run_online_blast(fasta_file, blast_results_file)
        
        return blast_results_file
    
    def _prepare_blast_fasta(self, binding_sites: Dict[str, List[Dict[str, Any]]], fasta_file: str):
        """Write a FASTA file for BLAST from candidate binding sites."""
        with open(fasta_file, 'w') as f:
            for gene_name, sites in binding_sites.items():
                for i, site in enumerate(sites):
                    f.write(f">{gene_name}_{i+1}|pos={site['position']}|g_content={site['g_content']:.2f}\n")
                    f.write(f"{site['sequence']}\n")
    
    def _run_local_blast(self, input_file: str, output_file: str):
        """Execute local BLAST."""
        blastn_cline = NcbiblastnCommandline(
            query=input_file,
            db=self.blast_config.database,
            task=self.blast_config.task,
            evalue=self.blast_config.evalue,
            outfmt=5,
            out=output_file
        )
        stdout, stderr = blastn_cline()
        if stderr:
            print(f"BLAST warning: {stderr}")
    
    def _run_online_blast(self, input_file: str, output_file: str):
        """Execute BLAST on NCBI servers."""
        from Bio.Blast import NCBIWWW
        with open(input_file, 'r') as f:
            fasta_string = f.read()
        print("Submitting BLAST queries to NCBI...")
        handle = NCBIWWW.qblast(
            program='blastn',
            database=self.blast_config.database,
            sequence=fasta_string,
            url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi',
            format_object='Alignment',
            format_type='Xml'
        )
        with open(output_file, 'w') as f:
            f.write(handle.read())
    
    def post_blast_filter(self, binding_sites: Dict[str, List[Dict[str, Any]]], blast_results_file: str) -> Dict[str, List[Dict[str, Any]]]:
        """Filter candidates after BLAST based on specificity and alignment metrics."""
        blast_data = self._parse_blast_results(blast_results_file)
        filtered_sites = {}
        for gene_name, sites in tqdm(binding_sites.items(), desc="Post-BLAST filtering"):
            filtered_sites[gene_name] = []
            for i, site in enumerate(sites):
                site_id = f"{gene_name}_{i+1}"
                if site_id in blast_data:
                    blast_info = blast_data[site_id]
                    if self._check_post_blast_criteria(blast_info):
                        site.update({
                            'blast_alignments': blast_info['alignments'],
                            'blast_evalue': blast_info['evalue'],
                            'blast_identity': blast_info['identity']
                        })
                        filtered_sites[gene_name].append(site)
        return filtered_sites
    
    def _parse_blast_results(self, blast_results_file: str) -> Dict[str, Dict[str, Any]]:
        """Parse BLAST XML results into a structured dictionary."""
        blast_data = {}
        with open(blast_results_file, 'r') as f:
            blast_records = NCBIXML.parse(f)
            for blast_record in blast_records:
                query_id = blast_record.query.split('|')[0]
                alignments = []
                evalues = []
                identities = []
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        alignments.append({
                            'hit_id': alignment.hit_id,
                            'hit_def': alignment.hit_def,
                            'evalue': hsp.expect,
                            'identity': hsp.identities / hsp.align_length,
                            'score': hsp.score,
                            'frame': hsp.frame
                        })
                        evalues.append(hsp.expect)
                        identities.append(hsp.identities / hsp.align_length)
                blast_data[query_id] = {
                    'alignments': alignments,
                    'evalue': min(evalues) if evalues else float('inf'),
                    'identity': max(identities) if identities else 0.0,
                    'alignment_count': len(alignments)
                }
        return blast_data
    
    def _check_post_blast_criteria(self, blast_info: Dict[str, Any]) -> bool:
        """Evaluate post-BLAST criteria for specificity and off-target risk."""
        if blast_info['alignment_count'] > self.filter_config.max_alignments:
            return False
        if self.filter_config.require_specificity:
            non_target_alignments = 0
            for alignment in blast_info['alignments']:
                hit_def = alignment['hit_def'].lower()
                is_target_organism = any(org.lower() in hit_def for org in self.filter_config.target_organisms)
                if not is_target_organism:
                    non_target_alignments += 1
            if non_target_alignments > 0:
                return False
        return True
    
    def save_filtered_results(self, filtered_sites: Dict[str, List[Dict[str, Any]]], output_dir: str):
        """Save filtered binding sites to Excel and JSON, plus stats."""
        os.makedirs(output_dir, exist_ok=True)
        all_sites = []
        for gene_name, sites in filtered_sites.items():
            for site in sites:
                site_data = {
                    'gene_name': gene_name,
                    'sequence': site['sequence'],
                    'position': site['position'],
                    'g_content': site['g_content'],
                    'tm': site['tm'],
                    'strategy': site['strategy']
                }
                if 'blast_alignments' in site:
                    site_data['blast_evalue'] = site['blast_evalue']
                    site_data['blast_identity'] = site['blast_identity']
                    site_data['alignment_count'] = len(site['blast_alignments'])
                all_sites.append(site_data)
        df = pd.DataFrame(all_sites)
        excel_file = os.path.join(output_dir, "filtered_binding_sites.xlsx")
        df.to_excel(excel_file, index=False)
        json_file = os.path.join(output_dir, "filtered_binding_sites.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(filtered_sites, f, indent=2, ensure_ascii=False)
        stats = {
            'total_genes': len(filtered_sites),
            'total_sites': sum(len(sites) for sites in filtered_sites.values()),
            'genes_with_sites': len([g for g, sites in filtered_sites.items() if sites]),
            'avg_sites_per_gene': sum(len(sites) for sites in filtered_sites.values()) / len(filtered_sites) if filtered_sites else 0
        }
        stats_file = os.path.join(output_dir, "filtering_stats.json")
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(stats, f, indent=2, ensure_ascii=False)
        return df


class RNAStructureFilter:
    """Optional RNA secondary structure filter using ViennaRNA if available."""
    
    def __init__(self, min_free_energy: float = -10.0):
        self.min_free_energy = min_free_energy
    
    def check_rna_structure(self, sequence: str) -> bool:
        """Check RNA secondary structure free energy threshold."""
        try:
            import RNA
            structure, energy = RNA.fold(sequence)
            return energy >= self.min_free_energy
        except ImportError:
            print("Warning: RNAfold not installed, skipping RNA structure check")
            return True
        except Exception as e:
            print(f"RNA structure check failed: {e}")
            return True
