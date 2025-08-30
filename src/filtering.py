"""
Sequence filtering module
Implements pre- and post-BLAST filtering rules.
"""

import os
import json
import pandas as pd
from typing import List, Dict, Any, Optional, Literal
from Bio.Blast import NCBIXML
from Bio.SeqUtils import MeltingTemp as mt
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from .config import FilterConfig, BlastConfig


class SequenceFilter:
    """Sequence filter for candidate binding sites."""
    
    def __init__(self, filter_config: FilterConfig, blast_config: BlastConfig):
        self.filter_config = filter_config
        self.blast_config = blast_config
    
    def thermal_filter(self, arm_3prime: str, arm_5prime: str,
                      sequence_type: Literal["DNA", "RNA"] = "DNA",
                      target_type: Literal["DNA", "RNA"] = "RNA",
                      target_sequence: Optional[str] = None,
                      **kwargs) -> Dict[str, Any]:
        """
        Thermal filter function for screening sequences based on thermodynamic properties
        
        Args:
            arm_3prime: 3' arm sequence
            arm_5prime: 5' arm sequence
            sequence_type: Type of sequence ("DNA" or "RNA")
            target_type: Type of target ("DNA" or "RNA")
            target_sequence: Target sequence for RNA secondary structure analysis
            **kwargs: Optional filter parameter overrides
            
        Returns:
            Dictionary containing filter results and thermodynamic parameters
        """
        # Get filter parameters, prioritize parameters in kwargs
        min_g_content = kwargs.get('min_g_content', self.filter_config.min_g_content)
        max_g_content = kwargs.get('max_g_content', self.filter_config.max_g_content)
        max_consecutive_g = kwargs.get('max_consecutive_g', self.filter_config.max_consecutive_g)
        min_tm = kwargs.get('min_tm', self.filter_config.min_tm)
        max_tm = kwargs.get('max_tm', self.filter_config.max_tm)
        max_tm_diff = kwargs.get('max_tm_diff', getattr(self.filter_config, 'max_tm_diff', 10.0))
        min_free_energy = kwargs.get('min_free_energy', self.filter_config.min_free_energy)
        check_rna_structure = kwargs.get('check_rna_structure', getattr(self.filter_config, 'check_rna_structure', False))
        
        # Combine arms to get full sequence
        full_sequence = arm_5prime + arm_3prime
        
        # Initialize result dictionary
        result = {
            'arm_3prime': arm_3prime,
            'arm_5prime': arm_5prime,
            'full_sequence': full_sequence,
            'target_sequence': target_sequence,
            'sequence_type': sequence_type,
            'target_type': target_type,
            'passed': False,
            'g_content': 0.0,
            'tm': 0.0,
            'tm_3prime': 0.0,
            'tm_5prime': 0.0,
            'tm_diff': 0.0,
            'free_energy': 0.0,
            'has_consecutive_g': False,
            'failed_checks': [],
            'arm_3prime_length': len(arm_3prime),
            'arm_5prime_length': len(arm_5prime)
        }
        
        # 1. Check G content for full sequence
        g_content = full_sequence.count('G') / len(full_sequence)
        result['g_content'] = g_content
        if not (min_g_content <= g_content <= max_g_content):
            result['failed_checks'].append('g_content')
        
        # 2. Check consecutive Gs in full sequence
        consecutive_g_count = max_consecutive_g + 1
        has_consecutive_g = "G" * consecutive_g_count in full_sequence
        result['has_consecutive_g'] = has_consecutive_g
        if has_consecutive_g:
            result['failed_checks'].append('consecutive_g')
        
        # 3. Calculate melting temperature for full sequence
        try:
            if sequence_type == "DNA" and target_type == "RNA":
                # DNA-RNA hybridization
                tm = mt.Tm_NN(full_sequence, nn_table=mt.R_DNA_NN1)
            elif sequence_type == "RNA" and target_type == "RNA":
                # RNA-RNA hybridization
                tm = mt.Tm_NN(full_sequence, nn_table=mt.RNA_NN1)
            else:
                # DNA-DNA hybridization
                tm = mt.Tm_NN(full_sequence, nn_table=mt.DNA_NN4)
            result['tm'] = tm
        except Exception as e:
            result['failed_checks'].append('tm_calculation_error')
            tm = 0.0
            result['tm'] = tm
        
        # 4. Check melting temperature range for full sequence
        if not (min_tm <= tm <= max_tm):
            result['failed_checks'].append('tm_range')
        
        # 5. Calculate melting temperature for 3' and 5' arms
        try:
            if sequence_type == "DNA" and target_type == "RNA":
                tm_3 = mt.Tm_NN(arm_3prime, nn_table=mt.R_DNA_NN1)
                tm_5 = mt.Tm_NN(arm_5prime, nn_table=mt.R_DNA_NN1)
            elif sequence_type == "RNA" and target_type == "RNA":
                tm_3 = mt.Tm_NN(arm_3prime, nn_table=mt.RNA_NN1)
                tm_5 = mt.Tm_NN(arm_5prime, nn_table=mt.RNA_NN1)
            else:
                tm_3 = mt.Tm_NN(arm_3prime, nn_table=mt.DNA_NN4)
                tm_5 = mt.Tm_NN(arm_5prime, nn_table=mt.DNA_NN4)
            
            result['tm_3prime'] = tm_3
            result['tm_5prime'] = tm_5
            result['tm_diff'] = abs(tm_5 - tm_3)
            
            # Check melting temperature for 3' and 5' arms
            if not (min_tm <= tm_3 <= max_tm):
                result['failed_checks'].append('tm_3prime_range')
            if not (min_tm <= tm_5 <= max_tm):
                result['failed_checks'].append('tm_5prime_range')
            
            # Check melting temperature difference
            if result['tm_diff'] > max_tm_diff:
                result['failed_checks'].append('tm_diff')
                
        except Exception as e:
            result['failed_checks'].append('tm_half_calculation_error')
        
        # 6. Check RNA secondary structure (if enabled and target sequence provided)
        if check_rna_structure and target_sequence and target_type == "RNA":
            try:
                import RNA
                _, mfe = RNA.fold(target_sequence)
                result['free_energy'] = mfe
                if mfe < min_free_energy:
                    result['failed_checks'].append('rna_structure')
            except ImportError:
                # RNAfold not available, skip this check
                pass
            except Exception as e:
                result['failed_checks'].append('rna_structure_error')
        
        # 7. Determine if all checks passed
        result['passed'] = len(result['failed_checks']) == 0
        
        return result
    
    def _optimize_subsequence(self, positions: List[int], length: int, min_gap: int = 1, 
                            better_gap: int = 1, gene: str = "") -> List[int]:
        """Optimize subsequence selection to maximize distance between sites."""
        if len(positions) == 0:
            print(f"Gene {gene}: No valid positions found, please loosen threshold conditions.")
            return []
        
        if len(positions) <= length:
            return positions
        
        positions.sort()
        
        def is_valid(min_difference: int, target_length: int, min_gap: int) -> bool:
            count = 1
            current_min = positions[0]
            
            for i in range(1, len(positions)):
                if positions[i] - current_min >= min_difference:
                    count += 1
                    current_min = positions[i]
            
            return count >= target_length and min_difference > min_gap
        
        # Binary search for optimal spacing
        left, right = 0, positions[-1] - positions[0]
        result = []
        
        while left <= right:
            mid = (left + right) // 2
            if is_valid(mid, length, min_gap):
                result = [positions[0]]
                current_min = positions[0]
                
                for i in range(1, len(positions)):
                    if positions[i] - current_min >= mid:
                        result.append(positions[i])
                        current_min = positions[i]
                        if len(result) >= length:
                            break
                
                left = mid + 1
            else:
                right = mid - 1
        
        if not result:
            print(f"Gene {gene}: Not enough positions for {length} binding sites.")
            result = positions[:length]
        
        if mid < better_gap:
            print(f"Gene {gene}: Conditions too harsh, consider loosening parameters for better results")
        
        return result
    
    def pre_blast_filter(self, binding_sites: Dict[str, List[Dict[str, Any]]]) -> Dict[str, List[Dict[str, Any]]]:
        """Filter candidates before running BLAST using thermal filter."""
        filtered_sites = {}
        
        print(f"Applying pre-BLAST filters to {len(binding_sites)} genes...")
        for gene_name, sites in binding_sites.items():
            filtered_sites[gene_name] = []
            for site in sites:
                # Use thermal_filter for screening
                thermal_result = self.thermal_filter(
                    site['arm_3prime'],
                    site['arm_5prime'],
                    sequence_type="DNA",
                    target_type="RNA",
                    target_sequence=site['target_sequence']
                )
                
                if thermal_result['passed']:
                    # Update site dictionary with thermodynamic parameters
                    site.update({
                        'g_content': thermal_result['g_content'],
                        'tm': thermal_result['tm'],
                        'tm_3prime': thermal_result['tm_3prime'],
                        'tm_5prime': thermal_result['tm_5prime'],
                        'tm_diff': thermal_result['tm_diff'],
                        'free_energy': thermal_result['free_energy'],
                        'thermal_filter_passed': True
                    })
                    filtered_sites[gene_name].append(site)
        
        total_original = sum(len(sites) for sites in binding_sites.values())
        total_filtered = sum(len(sites) for sites in filtered_sites.values())
        print(f"Pre-BLAST filtering: {total_original} -> {total_filtered} sites")
        
        return filtered_sites
    

    
    def run_blast(self, binding_sites: Dict[str, List[Dict[str, Any]]], output_dir: str, target_organisms: Optional[List[str]] = None) -> str:
        """Run BLAST analysis and return the path to results (file or directory)."""
        os.makedirs(output_dir, exist_ok=True)
        fasta_file = os.path.join(output_dir, "blast_input.fasta")
        self._prepare_blast_fasta(binding_sites, fasta_file)
        
        # Use blast species from config if not provided, fallback to filter target_organisms
        if target_organisms is None:
            if self.blast_config.species is not None:
                target_organisms = self.blast_config.species
            else:
                target_organisms = self.filter_config.target_organisms
        
        if self.blast_config.blast_type == "local":
            blast_results_file = os.path.join(output_dir, "blast_results.xml")
            self._run_local_blast(fasta_file, blast_results_file, target_organisms)
            return blast_results_file
        else:
            # Online qblast with batching and concurrency
            results_dir = os.path.join(output_dir, "blast_batches")
            os.makedirs(results_dir, exist_ok=True)
            self._run_online_blast_batched(binding_sites, results_dir, target_organisms)
            return results_dir
    
    def _prepare_blast_fasta(self, binding_sites: Dict[str, List[Dict[str, Any]]], fasta_file: str):
        """Write a FASTA file for BLAST from candidate binding sites."""
        with open(fasta_file, 'w') as f:
            for gene_name, sites in binding_sites.items():
                for i, site in enumerate(sites):
                    f.write(f">{gene_name}_{i+1}|pos={site['position']}|g_content={site['g_content']:.2f}\n")
                    f.write(f"{site['sequence']}\n")
    
    def _run_local_blast(self, input_file: str, output_file: str, target_organisms: List[str]):
        """Execute local BLAST with fallback to online BLAST if local fails."""
        try:
            # Use subprocess instead of Bio.Application to avoid deprecation warning
            import subprocess
            
            cmd = [
                'blastn',
                '-query', input_file,
                '-db', self.blast_config.database,
                '-task', self.blast_config.task,
                '-evalue', str(self.blast_config.evalue),
                '-outfmt', '5',
                '-out', output_file
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            if result.stderr:
                print(f"BLAST warning: {result.stderr}")
            
            print(f"Local BLAST completed successfully with target organisms: {target_organisms}")
            
        except subprocess.CalledProcessError as e:
            print(f"Local BLAST failed: {e}")
            print("Falling back to online BLAST...")
            
            # Switch to online BLAST
            self.blast_config.blast_type = "online"
            self._run_online_blast(input_file, output_file)
            print(f"Online BLAST completed with target organisms: {target_organisms}")
        except FileNotFoundError:
            print("BLAST+ not found in PATH. Falling back to online BLAST...")
            
            # Switch to online BLAST
            self.blast_config.blast_type = "online"
            self._run_online_blast(input_file, output_file)
            print(f"Online BLAST completed with target organisms: {target_organisms}")
    
    def _run_online_blast(self, input_file: str, output_file: str):
        """Execute BLAST on NCBI servers."""
        from Bio.Blast import NCBIWWW
        import time
        
        with open(input_file, 'r') as f:
            fasta_string = f.read()
        
        print("Submitting BLAST queries to NCBI...")
        print("Note: Online BLAST may take several minutes and may have limitations.")
        print("Consider using local BLAST for better performance.")
        
        try:
            # Use correct parameter names for NCBIWWW.qblast
            handle = NCBIWWW.qblast(
                program='blastn',
                database=self.blast_config.database,
                sequence=fasta_string,
                url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                expect=self.blast_config.evalue,
                hitlist_size=self.blast_config.hitlist_size,
                alignments=self.blast_config.alignments,
                descriptions=self.blast_config.descriptions,
                megablast=self.blast_config.megablast,
                short_query=self.blast_config.short_query,
                filter=self.blast_config.filter,
                format_type=self.blast_config.format_type,
                service=self.blast_config.service
            )
            
            # Read the result
            result = handle.read()
            
            # Check if result is valid XML
            if result.strip().startswith('<?xml'):
                with open(output_file, 'w') as f:
                    f.write(result)
                print("BLAST results saved successfully.")
            else:
                print("Warning: BLAST returned non-XML result. This might be due to:")
                print("1. NCBI service limitations")
                print("2. Query too large")
                print("3. Network issues")
                print("Creating empty BLAST results file for now.")
                # Create an empty XML file
                with open(output_file, 'w') as f:
                    f.write('<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>')
                
        except Exception as e:
            print(f"BLAST online query failed: {e}")
            print("Creating empty BLAST results file.")
            with open(output_file, 'w') as f:
                f.write('<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>')
    
    def _run_online_blast_batched(self, binding_sites: Dict[str, List[Dict[str, Any]]], results_dir: str, target_organisms: List[str]):
        """Run online BLAST in batches with concurrency. Writes multiple XML files in results_dir."""
        from Bio.Blast import NCBIWWW
        
        # Prepare all FASTA entries in-memory to preserve headers
        fasta_entries: List[str] = []
        for gene_name, sites in binding_sites.items():
            for i, site in enumerate(sites):
                header = f">{gene_name}_{i+1}|pos={site['position']}|g_content={site['g_content']:.2f}"
                seq = site['sequence']
                fasta_entries.append(f"{header}\n{seq}\n")
        
        if not fasta_entries:
            return
        
        # Chunk into batches
        batch_size = max(1, getattr(self.blast_config, 'batch_size', 100))
        batches = ["".join(fasta_entries[i:i+batch_size]) for i in range(0, len(fasta_entries), batch_size)]
        concurrency = max(1, getattr(self.blast_config, 'concurrency', 3))
        
        print(f"Submitting {len(batches)} BLAST batch(es) with concurrency={concurrency}, batch_size={batch_size}")
        print(f"Target organisms: {target_organisms}")
        
        def submit_batch(batch_index: int, fasta_string: str) -> str:
            try:
                handle = NCBIWWW.qblast(
                    program='blastn',
                    database=self.blast_config.database,
                    sequence=fasta_string,
                    url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                    expect=self.blast_config.evalue,
                    hitlist_size=self.blast_config.hitlist_size,
                    alignments=self.blast_config.alignments,
                    descriptions=self.blast_config.descriptions,
                    megablast=self.blast_config.megablast,
                    short_query=self.blast_config.short_query,
                    filter=self.blast_config.filter,
                    format_type=self.blast_config.format_type,
                    service=self.blast_config.service
                )
                result = handle.read()
                out_path = os.path.join(results_dir, f"blast_results_{batch_index:04d}.xml")
                if result.strip().startswith('<?xml'):
                    with open(out_path, 'w') as f:
                        f.write(result)
                else:
                    # Write minimal xml to keep parser happy
                    with open(out_path, 'w') as f:
                        f.write('<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>')
                return out_path
            except Exception as e:
                print(f"BLAST batch {batch_index} failed: {e}")
                out_path = os.path.join(results_dir, f"blast_results_{batch_index:04d}.xml")
                with open(out_path, 'w') as f:
                    f.write('<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>')
                return out_path
        
        # Run with thread pool
        with ThreadPoolExecutor(max_workers=concurrency) as executor:
            futures = {executor.submit(submit_batch, idx, batch): idx for idx, batch in enumerate(batches)}
            for future in tqdm(as_completed(futures), total=len(futures), desc="BLAST batches"):
                _ = future.result()

    def specificity_filter(self, binding_sites: Dict[str, List[Dict[str, Any]]], blast_results_file: str) -> Dict[str, List[Dict[str, Any]]]:
        """Filter candidates based on specificity and complementarity with reference sequences."""
        blast_data = self._parse_blast_results(blast_results_file)
        filtered_sites = {}
        
        # If no BLAST data available, return empty results
        if not blast_data:
            print("No BLAST data available - specificity filtering failed.")
            print("No sequences will be returned due to BLAST failure.")
            return {}
        
        print(f"Applying specificity filters to {len(binding_sites)} genes...")
        # Specificity filtering
        for gene_name, sites in binding_sites.items():
            filtered_sites[gene_name] = []
            for i, site in enumerate(sites):
                site_id = f"{gene_name}_{i+1}"
                if site_id in blast_data:
                    blast_info = blast_data[site_id]
                    if self._check_specificity_criteria(blast_info, site):
                        site.update({
                            'blast_alignments': blast_info['alignments'],
                            'blast_evalue': blast_info['evalue'],
                            'blast_identity': blast_info['identity'],
                            'specificity_status': 'specificity_filtered'
                        })
                        filtered_sites[gene_name].append(site)
                else:
                    # If no BLAST data for this sequence, skip it
                    print(f"Warning: No BLAST data found for {site_id}")
        
        total_original = sum(len(sites) for sites in binding_sites.values())
        total_filtered = sum(len(sites) for sites in filtered_sites.values())
        print(f"Specificity filtering: {total_original} -> {total_filtered} sites")
        
        return filtered_sites
    
    def _parse_blast_results(self, blast_results_path: str) -> Dict[str, Dict[str, Any]]:
        """Parse BLAST XML results (single file or directory of files) into a dictionary."""
        blast_data: Dict[str, Dict[str, Any]] = {}
        xml_files: List[str] = []
        
        if os.path.isdir(blast_results_path):
            # Collect all XML files from directory (batches)
            xml_files = [
                os.path.join(blast_results_path, f)
                for f in sorted(os.listdir(blast_results_path)) if f.lower().endswith('.xml')
            ]
            if not xml_files:
                print("Warning: No BLAST XML files found in directory.")
                return {}
        else:
            xml_files = [blast_results_path]
        
        from io import StringIO
        total_records = 0
        for xml_path in xml_files:
            try:
                with open(xml_path, 'r') as f:
                    content = f.read().strip()
                if not content or content == '<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>':
                    continue
                records = NCBIXML.parse(StringIO(content))
                for blast_record in records:
                    total_records += 1
                    query_def = blast_record.query
                    query_id = query_def.split('|')[0] if '|' in query_def else query_def
                    alignments = []
                    evalues = []
                    identities = []
                    if hasattr(blast_record, 'alignments') and blast_record.alignments:
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
            except Exception as e:
                print(f"Error parsing BLAST XML {xml_path}: {e}")
                continue
        print(f"Parsed {total_records} BLAST record(s) from {len(xml_files)} file(s)")
        return blast_data
    
    def _check_specificity_criteria(self, blast_info: Dict[str, Any], site: Dict[str, Any]) -> bool:
        """Evaluate specificity criteria including organism specificity and complementarity."""
        if self.filter_config.require_specificity:
            # Check organism specificity and gene specificity
            target_organism_alignments = 0
            non_target_gene_alignments = 0
            complementarity_found = False
            
            for alignment in blast_info['alignments']:
                hit_def = alignment['hit_def'].lower()
                is_target_organism = any(org.lower() in hit_def for org in self.filter_config.target_organisms)
                
                if is_target_organism:
                    target_organism_alignments += 1
                    # Check if this is the same gene (case-insensitive)
                    gene_name = site.get('gene_name', '').lower()
                    is_same_gene = gene_name in hit_def or any(gene_part in hit_def for gene_part in gene_name.split('_'))
                    
                    if not is_same_gene:
                        non_target_gene_alignments += 1
                    else:
                        # Check complementarity with target sequence
                        if self._check_complementarity(site['sequence'], alignment):
                            complementarity_found = True
            
            # Must have target organism alignment AND complementarity AND no non-target gene alignments
            if target_organism_alignments == 0 or not complementarity_found or non_target_gene_alignments > 0:
                return False
        
        return True
    
    def _check_complementarity(self, probe_sequence: str, alignment: Dict[str, Any]) -> bool:
        """Check if probe sequence is complementary to the aligned reference sequence."""
        try:
            # Check complementarity based on BLAST frame information
            # frame[1] = -1 indicates complementary strand
            if 'frame' in alignment and len(alignment['frame']) >= 2:
                frame_1 = alignment['frame'][1]
                if frame_1 == -1:
                    return True  # Complementary strand
                elif frame_1 == 1:
                    return False  # Same strand
                else:
                    # Unknown frame, use identity as fallback
                    if alignment['identity'] > 0.95:
                        return True
            
            # Fallback: if identity is very high (>95%), assume it's complementary
            if alignment['identity'] > 0.95:
                return True
            
            return False
        except Exception as e:
            print(f"Error checking complementarity: {e}")
            return False
    
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
