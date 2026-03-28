"""
Unified database interface module
Provides access to Ensembl and NCBI for sequence retrieval.
"""

import os
import time
import requests
from typing import List, Dict, Any, Optional, Tuple, Callable
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from math import ceil

from .config import DatabaseConfig, GenomeConfig


class DatabaseInterface:
    """Unified database interface."""
    
    def __init__(self, config: DatabaseConfig):
        self.config = config
        self.server = "https://rest.ensembl.org"
        self.headers = {"Content-Type": "application/json"}
        
        # Set up NCBI Entrez credentials
        Entrez.email = "1418767067@qq.com"
        Entrez.api_key = '010eacb785458478918b0cb14bea9f9df609'
    
    def ensembl_genome_accessor(self, species: str = None, coord_system_version: str = None):
        """Create Ensembl REST API genome accessor function."""
        if species is None:
            species = self.config.organism
        if coord_system_version is None:
            coord_system_version = self.config.coord_system_version
        
        def fetch(seq_region_name: str, start: int, end: int) -> str:
            chrom = str(seq_region_name).replace("chr", "")
            ext = f"/sequence/region/{species}/{chrom}:{start}..{end}:1?coord_system_version={coord_system_version}"
            try:
                response = requests.get(self.server + ext, headers=self.headers)
                response.raise_for_status()
                return response.json().get('seq', '')
            except Exception as e:
                print(f"Failed to fetch genome sequence from Ensembl: {e}")
                return ''
        
        return fetch
    
    def local_genome_accessor(self, genome_fasta_path: str):
        """Create local FASTA genome accessor function using pyfaidx."""
        try:
            from pyfaidx import Fasta
            genome = Fasta(genome_fasta_path)
            
            def fetch(seq_region_name: str, start: int, end: int) -> str:
                chrom = str(seq_region_name)
                if chrom not in genome:
                    raise KeyError(f"Chromosome {chrom} not found in FASTA file")
                seq = genome[chrom][start-1:end]  # pyfaidx uses 0-based coordinates
                return str(seq)
            
            return fetch
        except ImportError:
            print("pyfaidx not available, falling back to Ensembl")
            return None
        except Exception as e:
            print(f"Failed to load local genome: {e}")
            return None
    
    def initialize_genome_accessor(self, genome_config: GenomeConfig) -> Callable:
        """
        Initialize genome accessor with fallback chain.
        
        Tries each accessor type in priority order (from genome_config.accessor_priority),
        returns the first successful one.
        
        Args:
            genome_config: GenomeConfig instance with accessor_priority and genome_fasta_path
            
        Returns:
            Genome accessor function (callable that takes seq_region_name, start, end)
            
        Raises:
            RuntimeError: If all accessors fail to initialize
        """
        accessor = None
        accessor_priority = getattr(genome_config, 'accessor_priority', ['local', 'ensembl'])
        accessor_used = None
        failed_accessors = []
        
        for accessor_type in accessor_priority:
            if accessor_type == 'local':
                fasta_path = getattr(genome_config, 'genome_fasta_path', None)
                if not fasta_path:
                    failed_accessors.append("local (genome_fasta_path not set)")
                    continue
                if not os.path.exists(fasta_path):
                    failed_accessors.append(f"local (file not found: {fasta_path})")
                    continue
                
                accessor = self.local_genome_accessor(fasta_path)
                if accessor:
                    accessor_used = f"local: {fasta_path}"
                    break
                else:
                    failed_accessors.append(f"local (initialization failed: {fasta_path})")
            
            elif accessor_type == 'ensembl':
                accessor = self.ensembl_genome_accessor()
                if accessor:
                    accessor_used = f"ensembl (species: {self.config.organism})"
                    break
                else:
                    failed_accessors.append("ensembl (initialization failed)")
        
        if accessor is None:
            error_msg = f"All genome accessors failed:\n  " + "\n  ".join(failed_accessors)
            print(f"ERROR: {error_msg}")
            raise RuntimeError(error_msg)
        
        if failed_accessors:
            print(f"Warning: Failed accessors: {'; '.join(failed_accessors)}")
        print(f"Using genome accessor: {accessor_used}")
        
        return accessor
    
    def get_isoform_info(self, gene_symbol: str) -> List[Dict]:
        """Retrieve detailed isoform information for a gene (Ensembl only)."""    
        if self.config.database_type == "ensembl":
            try:
                gene_info = self._ensembl_get_tx_info(gene_symbol)
                transcripts = gene_info['Transcript']
                return transcripts
            except Exception as e:
                print(f"Failed to get isoform info for gene {gene_symbol}: {e}")
                return []
        else:
            raise NotImplementedError("Isoform analysis supported for Ensembl only")

    def get_gene_sequences(self, gene_list: List[str], gene_ids: List[str] = None) -> Dict[str, Any]:
        """Retrieve sequences for a list of genes."""
        if self.config.database_type == "ensembl":
            return self._ensembl_get_gene_sequences(gene_list)
        elif self.config.database_type == "ncbi":
            return self._get_ncbi_sequences(gene_list, gene_ids=gene_ids)
        else:
            raise ValueError(f"Unsupported database type: {self.config.database_type}")
    
    def _ensembl_get_gene_sequences(self, gene_list: List[str]) -> Dict[str, Any]:
        """Fetch sequences from Ensembl."""
        sequences_of_all = {}
        error_messages = {gene: [] for gene in gene_list}
        
        for gene in tqdm(gene_list, desc="Fetch Ensembl sequences"):
            sequences_of_all[gene] = {}
            trial = 0
            while trial < self.config.max_retries:
                try:
                    # Lookup gene metadata
                    gene_data = self._ensembl_get_tx_info(gene)
                    if not gene_data:
                        error_messages[gene].append(f"Gene not found: {gene}")
                        break
                    
                    # List transcripts for gene
                    transcripts = gene_data.get('Transcript', [])
                    
                    # Choose transcript purely by actual sequence length
                    best = None
                    best_len = -1
                    for tx_info in transcripts:
                        seq = self._ensembl_get_cds(tx_info['id'])
                        tx_info['seq'] = seq
                        if tx_info and tx_info['seq']:
                            seq_len = len(tx_info['seq'])
                            if seq_len > best_len:
                                best = (tx_info['id'], tx_info)
                                best_len = seq_len
                    
                    if best is not None:
                        tx_id, tx_info = best
                        sequences_of_all[gene] = {
                            'gene_id': gene_data['id'],
                            'transcript_id': tx_id,
                            'sequence': tx_info['seq'],
                            'length': len(tx_info['seq']),
                            'source': 'ensembl',
                            'seq_region_name': tx_info.get('seq_region_name'),
                            'start': tx_info.get('start'),
                            'end': tx_info.get('end'),
                            'strand': tx_info.get('strand')
                        }
                        break
                    else:
                        # No transcript returned sequence this round, retry
                        error_messages[gene].append("No transcript returned sequence; retrying")
                        trial += 1
                        if trial < self.config.max_retries:
                            time.sleep(0.5)
                        continue
                except Exception as e:
                    error_messages[gene].append(f"Attempt {trial + 1} failed: {str(e)}")
                    trial += 1
                    if trial < self.config.max_retries:
                        time.sleep(1)
            # if no sequence after retries, sequences_of_all[gene] stays empty and errors captured
        
        return {
            'sequences': sequences_of_all,
            'errors': error_messages
        }
    
    def _ensembl_get_tx_info(self, gene_symbol: str) -> Optional[Dict]:
        """Lookup gene metadata by symbol."""
        lookup_ext = f"/lookup/symbol/{self.config.organism}/{gene_symbol}?expand=1"
        try:
            response = requests.get(self.server + lookup_ext, headers=self.headers, timeout=10)
            if response.ok:
                return response.json()
            else:
                return None
        except Exception:
            return None
    
    def _ensembl_get_tx_overlap(self, gene_id: str) -> List[Dict]:
        """List transcripts for the gene."""
        overlap_ext = f"/overlap/id/{gene_id}?feature=transcript"
        try:
            response = requests.get(self.server + overlap_ext, headers=self.headers, timeout=10)
            if response.ok:
                return response.json()
            else:
                return []
        except Exception:
            return []
    
    def _ensembl_get_cds(self, transcript_id: str) -> Optional[str]:
        """
        Fetch transcript sequence
        """
        # sequence type: cdna, cds, protein
        sequence_ext = f"/sequence/id/{transcript_id}?type=cds"
                
        try:
            response = requests.get(self.server + sequence_ext, headers=self.headers, timeout=10)
            response.raise_for_status() 
            data = response.json()            
            return data.get('seq', None)
        
        except (requests.exceptions.RequestException, requests.exceptions.JSONDecodeError) as e:
            return None
     

    def _get_ncbi_sequences(self, gene_list: List[str], gene_ids: List[str] = None) -> Dict[str, Any]:
        """Fetch sequences from NCBI using Entrez API."""
        sequences_of_all = {}
        error_messages = {gene: [] for gene in gene_list}

        # Step 1: Get gene IDs (skip search if preset IDs provided)
        if gene_ids is None:
            gene_ids = self._ncbi_get_gene_ids(gene_list)
        
        # Step 2: Fetch GenBank records in batches
        genbank_records = self._ncbi_fetch_genbank_records(gene_ids)
        
        # Step 3: Process GenBank records and extract sequences
        for gene, record in zip(gene_list, genbank_records):
            sequences_of_all[gene] = {}
            try:
                if record:
                    # Extract sequence and metadata from GenBank record
                    sequence, gene_info = SequenceProcessor.extract_coding_sequence(record)
                    
                    # Prefer genomic context from GenBank if record is genomic; fallback to Entrez
                    genomic_context = self._ncbi_get_genomic_context_from_genbank(record)
                    if genomic_context is None:
                        genomic_context = self._ncbi_get_genomic_context_from_nuccore(record.id)
                    
                    if sequence:
                        sequences_of_all[gene] = {
                            'gene_id': gene_info['gene_id'],
                            'gene_name': gene_info['gene_name'],
                            'sequence': sequence,
                            'length': len(sequence),
                            'organism': gene_info['organism'],
                            'mol_type': gene_info['mol_type'],
                            'source': 'ncbi',
                            'seq_region_name': genomic_context.get('seq_region_name') if genomic_context else None,
                            'start': genomic_context.get('start') if genomic_context else None,
                            'end': genomic_context.get('end') if genomic_context else None,
                            'strand': genomic_context.get('strand') if genomic_context else None,
                        }
                    else:
                        error_messages[gene].append("No coding sequence found in GenBank record")
                else:
                    error_messages[gene].append("No GenBank record found")
                    
            except Exception as e:
                error_messages[gene].append(f"Failed to process GenBank record: {str(e)}")
        
        return {
            'sequences': sequences_of_all,
            'errors': error_messages
        }
    
    def _ncbi_get_genomic_context_from_nuccore(self, nuccore_id: str) -> Optional[Dict[str, Any]]:
        """Map a nuccore (mRNA) ID to Gene and fetch GenomicInfo for coordinates and strand."""
        try:
            # Link nuccore to Gene
            link = Entrez.elink(dbfrom="nuccore", db="gene", id=nuccore_id)
            link_rec = Entrez.read(link)
            link.close()
            gene_ids = []
            for linkset in link_rec:
                for linksetdb in linkset.get('LinkSetDb', []):
                    if linksetdb.get('DbTo') == 'gene':
                        gene_ids.extend([l['Id'] for l in linksetdb.get('Link', [])])
            if not gene_ids:
                return None
            gene_id = gene_ids[0]
            # Gene esummary
            handle = Entrez.esummary(db="gene", id=gene_id)
            summary = Entrez.read(handle)
            handle.close()
            if not summary or not summary[0]:
                return None
            doc = summary[0]
            genomic_info = doc.get('GenomicInfo', []) or doc.get('GenomicInfoHC', [])
            if not genomic_info:
                return None
            # Pick the first locus (primary)
            gi = genomic_info[0]
            # Fields: ChrAccVer, ChrStart, ChrStop, Strand
            seq_region_name = gi.get('ChrAccVer') or gi.get('ChrLoc')
            start = int(gi.get('ChrStart')) if gi.get('ChrStart') is not None else None
            end = int(gi.get('ChrStop')) if gi.get('ChrStop') is not None else None
            strand = int(gi.get('Strand')) if gi.get('Strand') is not None else None
            return {
                'seq_region_name': seq_region_name,
                'start': start,
                'end': end,
                'strand': strand
            }
        except Exception:
            return None
    
    def _ncbi_get_genomic_context_from_genbank(self, genbank_record) -> Optional[Dict[str, Any]]:
        """Try to extract genomic context directly from a GenBank record when it is a genomic (chromosome/contig) record.
        Returns None for mRNA records (e.g., NM_/XM_)."""
        try:
            acc = getattr(genbank_record, 'id', '') or ''
            is_genomic = acc.startswith(('NC_', 'NG_', 'NW_', 'NT_'))
            if not is_genomic:
                return None
            # Prefer mRNA feature; fallback to CDS
            target_feat = None
            for feat in genbank_record.features:
                if feat.type == 'mRNA':
                    target_feat = feat
                    break
            if target_feat is None:
                for feat in genbank_record.features:
                    if feat.type == 'CDS':
                        target_feat = feat
                        break
            if target_feat is None:
                return None
            loc = target_feat.location
            start = int(loc.start)
            end = int(loc.end)
            strand = int(loc.strand) if loc.strand in (1, -1) else None
            seq_region_name = acc  # genomic accession
            # If chromosome qualifier exists on source feature, use it as human-readable
            for feat in genbank_record.features:
                if feat.type == 'source':
                    chrom = feat.qualifiers.get('chromosome', [None])[0]
                    if chrom:
                        seq_region_name = chrom
                        break
            return {
                'seq_region_name': seq_region_name,
                'start': start,
                'end': end,
                'strand': strand
            }
        except Exception:
            return None
    
    def _ncbi_get_gene_ids(self, gene_list: List[str]) -> List[str]:
        """Get NCBI gene IDs for a list of gene names."""
        id_list = []
        organism_mapping = {
            'mouse': 'Mus musculus',
            'human': 'Homo sapiens'
        }
        organism = organism_mapping.get(self.config.organism, 'Mus musculus')
        n_type = "mRNA"
        
        for gene in tqdm(gene_list, desc="Getting NCBI gene IDs"):
            try:
                # Create search term
                search_term = f"{gene}, {organism}, {n_type}"
                
                # Search for gene ID
                handle = Entrez.esearch(db="nuccore", term=search_term)
                record = Entrez.read(handle)
                handle.close()
                
                # Get the first (most relevant) result
                if record["IdList"]:
                    id_list.append(record["IdList"][0])
                else:
                    id_list.append(None)
                    
            except Exception as e:
                print(f"Failed to get ID for {gene}: {e}")
                id_list.append(None)
        
        return id_list
    
    def _ncbi_fetch_genbank_records(self, id_list: List[str]) -> List[Optional[SeqRecord]]:
        """Fetch GenBank records for a list of gene IDs."""
        records = []
        step = 3  # Batch size for fetching
        
        # Filter out None IDs
        valid_ids = [(i, id_val) for i, id_val in enumerate(id_list) if id_val is not None]
        
        for i in tqdm(range(ceil(len(valid_ids) / step)), desc="Fetching GenBank records"):
            batch_ids = valid_ids[i * step : (i + 1) * step]
            id_values = [id_val for _, id_val in batch_ids]
            
            try:
                # Fetch GenBank records for this batch
                handle = Entrez.efetch(
                    db="nucleotide",
                    strand=1,  # plus strand
                    id=id_values,
                    rettype="gbwithparts",
                    retmode="text",
                )
                
                # Parse GenBank records
                batch_records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                
                # Map records back to original positions
                for j, (original_idx, _) in enumerate(batch_ids):
                    if j < len(batch_records):
                        records.append((original_idx, batch_records[j]))
                    else:
                        records.append((original_idx, None))
                        
            except Exception as e:
                print(f"Failed to fetch batch {i}: {e}")
                for _, original_idx in batch_ids:
                    records.append((original_idx, None))
        
        # Sort records by original position
        records.sort(key=lambda x: x[0])
        return [record for _, record in records]
    
    def get_strand_by_symbol(self, gene_symbol: str, species: str = None) -> Optional[int]:
        """通过基因符号查询链信息"""
        if species is None:
            species = self.config.organism
        
        if self.config.database_type == "ensembl":
            return self._ensembl_get_strand_by_symbol(gene_symbol, species)
        else:
            # For NCBI, we would need to implement similar functionality
            # For now, return None to indicate not supported
            return None
   
    def _ensembl_get_strand_by_symbol(self, gene_symbol: str, species: str) -> Optional[int]:
        """通过Ensembl API查询基因的链信息"""
        ext = f"/lookup/symbol/{species}/{gene_symbol}"
        try:
            response = requests.get(self.server + ext, headers=self.headers, timeout=10)
            if response.ok: 
                return response.json().get('strand')
            else:
                print(f"Warning: {gene_symbol} not found (HTTP {response.status_code})")
                return None
        except Exception as e:
            print(f"Error for {gene_symbol}: {e}")
            return None

    def save_sequences_to_fasta(self, sequences: Dict[str, Any], output_file: str):
        """Save sequences to FASTA file."""
        with open(output_file, 'w') as f:
            for gene, data in sequences.items():
                if 'sequence' in data:
                    f.write(f">{gene}|{data.get('gene_id', 'unknown')}|{data.get('transcript_id', 'unknown')}\n")
                    f.write(f"{data['sequence']}\n")
    
    def save_sequences_to_json(self, sequences: Dict[str, Any], output_file: str):
        """Save sequences to JSON file."""
        import json
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(sequences, f, indent=2, ensure_ascii=False)
    

class SequenceProcessor:
    """Helpers for sequence extraction and transformation."""
    
    @staticmethod
    def extract_coding_sequence(genbank_record) -> Tuple[str, Dict]:
        """Extract CDS from a GenBank record and return sequence and metadata."""
        coding_sequence = ""
        gene_info = {
            'gene_id': genbank_record.id,
            'gene_name': 'Unknown',
            'mol_type': genbank_record.annotations.get("molecule_type", ""),
            'organism': genbank_record.annotations.get("organism", ""),
            'seq': str(genbank_record.seq)
        }
        if genbank_record.features:
            for feature in genbank_record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                    coding_sequence = feature.location.extract(genbank_record).seq
                    gene_info['gene_name'] = gene_name
                    break
        if coding_sequence:
            gene_info['seq'] = str(coding_sequence)
        return gene_info['seq'], gene_info
    
    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Return reverse complement of a DNA sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(complement[base] for base in reversed(sequence))
