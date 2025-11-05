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
# from tqdm import tqdm  # Imported when needed
# from concurrent.futures import ThreadPoolExecutor, as_completed  # Imported when needed

try:
    from .config import FilterConfig, BlastConfig
except ImportError:
    # Fallback for when running as standalone
    from config import FilterConfig, BlastConfig


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
        
        
        # Initialize result dictionary
        result = {
            'arm_3prime': arm_3prime,
            'arm_5prime': arm_5prime,
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
        
        # 1. Check G content per arm (padlock arms evaluated分别)
        def frac_g(seq: str) -> float:
            return (seq.upper().count('G') / len(seq)) if seq else 0.0
        g_content_arm3 = frac_g(arm_3prime)
        g_content_arm5 = frac_g(arm_5prime)
        # Backward-compat summary as均值
        result['g_content'] = (g_content_arm3 + g_content_arm5) / 2.0
        result['g_content_arm3'] = g_content_arm3
        result['g_content_arm5'] = g_content_arm5
        if not (min_g_content <= g_content_arm3 <= max_g_content):
            result['failed_checks'].append('g_content_arm3')
        if not (min_g_content <= g_content_arm5 <= max_g_content):
            result['failed_checks'].append('g_content_arm5')
        
        # 2. Check consecutive Gs (G-quadruplex risk) on each component separately
        def has_g_run(seq: str, threshold: int) -> bool:
            if not seq:
                return False
            run = 0
            for ch in seq.upper():
                if ch == 'G':
                    run += 1
                    if run >= threshold:
                        return True
                else:
                    run = 0
            return False

        has_g_arm3 = has_g_run(arm_3prime, max_consecutive_g)
        has_g_arm5 = has_g_run(arm_5prime, max_consecutive_g)
        has_g_target = has_g_run(target_sequence, max_consecutive_g) if target_sequence else False

        result['has_consecutive_g'] = any([has_g_arm3, has_g_arm5, has_g_target])
        result['has_consecutive_g_arm3'] = has_g_arm3
        result['has_consecutive_g_arm5'] = has_g_arm5
        result['has_consecutive_g_target'] = has_g_target
        if result['has_consecutive_g']:
            result['failed_checks'].append('consecutive_g')
        
        # 3-5. Compute padlock thermodynamic parameters (full and arm Tm) via helper
        thermo = self._compute_padlock_thermo(
            arm_3prime=arm_3prime,
            arm_5prime=arm_5prime,
            sequence_type=sequence_type,
            target_type=target_type,
            target_sequence=target_sequence
        )
        result['tm'] = thermo.get('tm', 0.0)
        tm_3 = thermo.get('tm_3prime', 0.0)
        tm_5 = thermo.get('tm_5prime', 0.0)
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
            except Exception:
                result['failed_checks'].append('rna_structure_error')
        
        # 7. Determine if all checks passed
        result['passed'] = len(result['failed_checks']) == 0
        
        return result

    def _compute_padlock_thermo(
        self,
        arm_3prime: str,
        arm_5prime: str,
        sequence_type: Literal["DNA", "RNA"],
        target_type: Literal["DNA", "RNA"],
        target_sequence: Optional[str]
    ) -> Dict[str, float]:
        """Compute Tm for full padlock and its two arms under given hybridization model.
        For DNA–RNA, we compute the RNA reverse-complement strand of DNA and use R_DNA_NN tables.
        """
        def dna_revcomp_to_rna(seq: str) -> str:
            comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
            return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

        # Assemble full probe sequence from arms
        full_sequence = arm_3prime + arm_5prime

        tm_full = 0.0
        tm_3 = 0.0
        tm_5 = 0.0
        
        try:
            if sequence_type == "DNA" and target_type == "RNA":
                rna_full = dna_revcomp_to_rna(full_sequence)
                tm_full = mt.Tm_NN(rna_full, nn_table=mt.R_DNA_NN1)
            elif sequence_type == "RNA" and target_type == "RNA":
                tm_full = mt.Tm_NN(full_sequence, nn_table=mt.RNA_NN1)
            else:
                tm_full = mt.Tm_NN(full_sequence, nn_table=mt.DNA_NN4)
        except Exception:
            tm_full = 0.0
        
        try:
            if sequence_type == "DNA" and target_type == "RNA":
                rna_3 = dna_revcomp_to_rna(arm_3prime)
                rna_5 = dna_revcomp_to_rna(arm_5prime)
                tm_3 = mt.Tm_NN(rna_3, nn_table=mt.R_DNA_NN1)
                tm_5 = mt.Tm_NN(rna_5, nn_table=mt.R_DNA_NN1)
            elif sequence_type == "RNA" and target_type == "RNA":
                tm_3 = mt.Tm_NN(arm_3prime, nn_table=mt.RNA_NN1)
                tm_5 = mt.Tm_NN(arm_5prime, nn_table=mt.RNA_NN1)
            else:
                tm_3 = mt.Tm_NN(arm_3prime, nn_table=mt.DNA_NN4)
                tm_5 = mt.Tm_NN(arm_5prime, nn_table=mt.DNA_NN4)
        except Exception:
            tm_3 = 0.0
            tm_5 = 0.0
        
        return {
            'tm': tm_full,
            'tm_3prime': tm_3,
            'tm_5prime': tm_5
        }
    
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
            # Check if we switched to online mode during execution
            if self.blast_config.blast_type == "online":
                # Return the batch directory that was created
                return os.path.join(output_dir, "blast_batches")
            else:
                return blast_results_file
        else:
            # Online qblast with batching and concurrency
            results_dir = os.path.join(output_dir, "blast_batches")
            os.makedirs(results_dir, exist_ok=True)
            self._run_online_blast_batched(binding_sites, results_dir, target_organisms)
            return results_dir
    
    def _generate_query_id(self, gene_name: str, site_index: int, position: int, g_content: float) -> str:
        """Generate a consistent query ID for BLAST operations.
        
        Args:
            gene_name: Name of the gene
            site_index: 0-based index of the site within the gene
            position: Position of the binding site
            g_content: GC content of the sequence
            
        Returns:
            Query ID in format: gene_name_index|pos=position|g_content=content
        """
        return f"{gene_name}_{site_index+1}|pos={position}|g_content={g_content:.2f}"
    
    def _prepare_blast_fasta(self, binding_sites: Dict[str, List[Dict[str, Any]]], fasta_file: str):
        """Write a FASTA file for BLAST from candidate binding sites."""
        with open(fasta_file, 'w') as f:
            for gene_name, sites in binding_sites.items():
                for i, site in enumerate(sites):
                    query_id = self._generate_query_id(gene_name, i, site['position'], site['g_content'])
                    f.write(f">{query_id}\n")
                    f.write(f"{site['sequence']}\n")
    
    def _load_binding_sites_from_fasta(self, fasta_file: str) -> Dict[str, List[Dict[str, Any]]]:
        """Load binding sites from FASTA file for use with batched online BLAST."""
        binding_sites = {}
        
        with open(fasta_file, 'r') as f:
            current_header = None
            current_sequence = None
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous entry if exists
                    if current_header and current_sequence:
                        self._add_sequence_to_binding_sites(binding_sites, current_header, current_sequence)
                    
                    # Start new entry
                    current_header = line[1:]  # Remove '>'
                    current_sequence = ""
                else:
                    if current_sequence is not None:
                        current_sequence += line
            
            # Process last entry
            if current_header and current_sequence:
                self._add_sequence_to_binding_sites(binding_sites, current_header, current_sequence)
        
        return binding_sites
    
    def _add_sequence_to_binding_sites(self, binding_sites: Dict[str, List[Dict[str, Any]]], 
                                      header: str, sequence: str):
        """Add a sequence to binding_sites dictionary based on FASTA header."""
        # Parse header: gene_name_index|pos=position|g_content=content
        parts = header.split('|')
        if len(parts) >= 1:
            gene_part = parts[0]
            # Extract gene name and index
            if '_' in gene_part:
                gene_name = '_'.join(gene_part.split('_')[:-1])
                # index = int(gene_part.split('_')[-1]) - 1  # Convert to 0-based (not used)
            else:
                gene_name = gene_part
            
            # Extract position and g_content from header
            position = 0
            g_content = 0.5
            for part in parts[1:]:
                if part.startswith('pos='):
                    position = int(part.split('=')[1])
                elif part.startswith('g_content='):
                    g_content = float(part.split('=')[1])
            
            # Create site dictionary
            site = {
                'sequence': sequence,
                'position': position,
                'g_content': g_content,
                'arm_3prime': sequence[:len(sequence)//2],  # Approximate split
                'arm_5prime': sequence[len(sequence)//2:],
                'target_sequence': sequence
            }
            
            # Add to binding_sites
            if gene_name not in binding_sites:
                binding_sites[gene_name] = []
            binding_sites[gene_name].append(site)
    
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
            print("Falling back to online BLAST with batching...")
            
            # Switch to online BLAST with batching
            self.blast_config.blast_type = "online"
            # Create batch directory for online results
            batch_dir = os.path.join(os.path.dirname(output_file), "blast_batches")
            os.makedirs(batch_dir, exist_ok=True)
            # Load binding sites from the input file to use with batched online BLAST
            binding_sites = self._load_binding_sites_from_fasta(input_file)
            self._run_online_blast_batched(binding_sites, batch_dir, target_organisms)
            print(f"Online BLAST with batching completed with target organisms: {target_organisms}")
        
        except FileNotFoundError:
            print("BLAST+ not found in PATH. Falling back to online BLAST with batching...")
            
            # Switch to online BLAST with batching
            self.blast_config.blast_type = "online"
            # Create batch directory for online results
            batch_dir = os.path.join(os.path.dirname(output_file), "blast_batches")
            os.makedirs(batch_dir, exist_ok=True)
            # Load binding sites from the input file to use with batched online BLAST
            binding_sites = self._load_binding_sites_from_fasta(input_file)
            self._run_online_blast_batched(binding_sites, batch_dir, target_organisms)
            print(f"Online BLAST with batching completed with target organisms: {target_organisms}")
    
    def _run_online_blast(self, input_file: str, output_file: str):
        """Execute BLAST on NCBI servers."""
        from Bio.Blast import NCBIWWW
        
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
                
        except Exception as err:
            print(f"BLAST online query failed: {err}")
            print("Creating empty BLAST results file.")
            with open(output_file, 'w') as f:
                f.write('<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>')
    
    def _run_online_blast_batched(self, binding_sites: Dict[str, List[Dict[str, Any]]], results_dir: str, target_organisms: List[str]):
        """Run online BLAST in batches with concurrency, caching, and error handling. Writes multiple XML files in results_dir."""
        import requests
        import concurrent.futures
        import time
        
        # Prepare all FASTA entries in-memory to preserve headers
        fasta_entries: List[str] = []
        for gene_name, sites in binding_sites.items():
            for i, site in enumerate(sites):
                header = f">{gene_name}_{i+1}|pos={site['position']}|g_content={site['g_content']:.2f}"
                seq = site['sequence']
                fasta_entries.append(f"{header}\n{seq}\n")
        
        if not fasta_entries:
            return
        
        # Chunk into smaller batches to avoid "Process size limit exceeded" error
        batch_size = max(1, getattr(self.blast_config, 'batch_size', 100))
        batches = ["".join(fasta_entries[i:i+batch_size]) for i in range(0, len(fasta_entries), batch_size)]
        concurrency = max(1, getattr(self.blast_config, 'concurrency', 3))
        
        print(f"Submitting {len(batches)} BLAST batch(es) with concurrency={concurrency}, batch_size={batch_size}")
        print(f"Target organisms: {target_organisms}")
        
        # Check for existing batch results (cache mechanism)
        existing_batches = {}
        for i in range(len(batches)):
            batch_file = os.path.join(results_dir, f"blast_results_{i+1:04d}.xml")
            if os.path.exists(batch_file) and os.path.getsize(batch_file) > 100:  # Check if file exists and has content
                existing_batches[i] = batch_file
                # print(f"Found existing results for batch {i + 1}")  # Reduced output
        
        # Get all problematic batches (missing, empty, or corrupted)
        problematic_batches = self._get_problematic_batches(results_dir, len(batches))
        
        if problematic_batches:
            print(f"Processing {len(problematic_batches)} missing/corrupted batches: {[i+1 for i in problematic_batches]}")
            
            def submit_batch(batch_index: int, fasta_string: str) -> str:
                """Submit a single batch using NCBI BLAST web service with proper error handling."""
                try:
                    # NCBI BLAST web service parameters
                    blast_params = {
                        'CMD': 'Put',
                        'PROGRAM': 'blastn',
                        'DATABASE': self.blast_config.database,
                        'QUERY': fasta_string,
                        'HITLIST_SIZE': str(self.blast_config.hitlist_size),
                        'EXPECT': str(self.blast_config.evalue),
                        'WORD_SIZE': '7',
                        'GAPCOSTS': '5 2',
                        'PENALTY': '-3',
                        'REWARD': '2',
                        'FORMAT_TYPE': 'XML'
                    }
                    
                    print(f"Submitting batch {batch_index+1}/{len(batches)} with {len(fasta_string.split('>'))-1} sequences...")
                    
                    # Submit BLAST job
                    response = requests.post('https://blast.ncbi.nlm.nih.gov/Blast.cgi', data=blast_params, timeout=30)
                    response.raise_for_status()
                    
                    # Extract RID (Request ID) from response
                    rid = None
                    for line in response.text.split('\n'):
                        if 'RID =' in line:
                            rid = line.split('RID =')[1].strip()
                            break
                    
                    if not rid:
                        print(f"Failed to get RID from BLAST submission for batch {batch_index+1}")
                        return None
                    
                    print(f"Batch {batch_index+1} submitted with RID: {rid}")
                    
                    # Poll for results with longer timeout
                    max_attempts = 20  # 30 minutes max per batch
                    for attempt in range(max_attempts):
                        time.sleep(90)  # Wait 180 seconds between checks
                        
                        status_params = {
                            'CMD': 'Get',
                            'RID': rid,
                            'FORMAT_TYPE': 'XML'
                        }
                        
                        status_response = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', params=status_params, timeout=30)
                        status_response.raise_for_status()
                        
                        # Check if we got XML results (completed)
                        if '<?xml version=' in status_response.text and '<BlastOutput>' in status_response.text:
                            print(f"Batch {batch_index+1} completed successfully")
                            
                            # Save batch results - use 1-based numbering for user-friendly file names
                            batch_output = os.path.join(results_dir, f"blast_results_{batch_index+1:04d}.xml")
                            with open(batch_output, 'w', encoding='utf-8') as f:
                                f.write(status_response.text)
                            
                            print(f"Saved batch results to: {batch_output}")
                            return batch_output
                        elif 'Status=WAITING' in status_response.text:
                            print(f"Batch {batch_index+1} still running... (attempt {attempt + 1}/{max_attempts})")
                            continue
                        elif 'Status=READY' in status_response.text:
                            print(f"Batch {batch_index+1} ready, fetching results...")
                            continue
                        else:
                            print(f"BLAST status for batch {batch_index+1}: {status_response.text[:100]}")
                            continue
                    
                    print(f"Batch {batch_index+1} timed out")
                    return None
                    
                except requests.RequestException as e:
                    print(f"Error during BLAST search for batch {batch_index+1}: {e}")
                    return None
                except Exception as e:
                    print(f"Unexpected error during BLAST search for batch {batch_index+1}: {e}")
                    return None
            
            # Process missing batches concurrently (limit to avoid overwhelming NCBI)
            max_workers = min(concurrency, len(problematic_batches))
            completed_batches = 0
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                # Submit only missing batches
                future_to_batch = {executor.submit(submit_batch, i, batches[i]): i for i in problematic_batches}
                
                # Collect results in order
                batch_results = {}
                for future in concurrent.futures.as_completed(future_to_batch):
                    batch_idx = future_to_batch[future]
                    try:
                        result = future.result()
                        if result:
                            batch_results[batch_idx] = result
                            completed_batches += 1
                            print(f"Completed batch {batch_idx+1}/{len(batches)}")
                        else:
                            print(f"Batch {batch_idx+1} failed")
                    except Exception as e:
                        print(f"Batch {batch_idx+1} generated an exception: {e}")
                
                # Add existing batch results
                for batch_idx, batch_file in existing_batches.items():
                    batch_results[batch_idx] = batch_file
        else:
            print("All batches already exist, using cached results")
            batch_results = existing_batches
            completed_batches = len(batch_results)
        
        print(f"BLAST search completed: {completed_batches}/{len(batches)} batches successful")
        
        # Check for any remaining missing batches and retry up to 3 times
        max_retries = 3
        for retry in range(max_retries):
            still_problematic = self._get_problematic_batches(results_dir, len(batches))
            
            if not still_problematic:
                print(f"All batches completed successfully after {retry + 1} attempts")
                break
            
            if retry < max_retries - 1:
                print(f"Retry {retry + 1}: {len(still_problematic)} batches still missing: {[i+1 for i in still_problematic]}")
                # Re-run missing batches
                with concurrent.futures.ThreadPoolExecutor(max_workers=min(concurrency, len(still_problematic))) as retry_executor:
                    retry_futures = {retry_executor.submit(submit_batch, i, batches[i]): i for i in still_problematic}
                    for future in concurrent.futures.as_completed(retry_futures):
                        batch_idx = retry_futures[future]
                        try:
                            result = future.result()
                            if result:
                                print(f"Retry completed batch {batch_idx+1}")
                        except Exception as e:
                            print(f"Retry failed for batch {batch_idx+1}: {e}")
            else:
                print(f"Warning: {len(still_problematic)} batches still missing after {max_retries} attempts: {[i+1 for i in still_problematic]}")
        
        if completed_batches == 0:
            print("No BLAST batches completed successfully")
            return
        
        print("All BLAST batch files saved successfully")
    
    def specificity_filter(self, binding_sites: Dict[str, List[Dict[str, Any]]], blast_results_file: str) -> Dict[str, List[Dict[str, Any]]]:
        """Filter candidates based on specificity and complementarity with reference sequences."""
        blast_data = self._parse_blast_results(blast_results_file, binding_sites)
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
                # Create the query ID using the same format as FASTA generation
                query_id = self._generate_query_id(gene_name, i, site['position'], site['g_content'])
                
                if query_id in blast_data:
                    blast_info = blast_data[query_id]
                    
                    # Check if this is missing BLAST data
                    if blast_info.get('missing', False):
                        print(f"Warning: Missing BLAST data for {query_id}")
                        continue
                    
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
                    print(f"Warning: No BLAST data found for {query_id}")
        
        total_original = sum(len(sites) for sites in binding_sites.values())
        total_filtered = sum(len(sites) for sites in filtered_sites.values())
        print(f"Specificity filtering: {total_original} -> {total_filtered} sites")
        
        return filtered_sites

    def _parse_blast_results(self, blast_results_path: str, binding_sites: Dict[str, List[Dict[str, Any]]] = None) -> Dict[str, Dict[str, Any]]:
        """Parse BLAST XML results (single file or directory of files) into a dictionary.
        
        For batch processing, this method ensures alignment between BLAST results and binding sites
        by inferring batch size and handling missing files properly.
        """
        blast_data: Dict[str, Dict[str, Any]] = {}
        
        if os.path.isdir(blast_results_path):
            # Batch processing mode - need to ensure alignment with binding sites
            if binding_sites is None:
                print("Warning: binding_sites required for batch processing alignment")
                return {}
            
            # Get total number of binding sites
            total_sites = sum(len(sites) for sites in binding_sites.values())
            print(f"Total binding sites to process: {total_sites}")
            
            # Find all batch files and sort by batch number (1-based)
            all_files = [f for f in os.listdir(blast_results_path) if f.lower().endswith('.xml')]
            batch_files = []
            for f in all_files:
                if f.startswith('blast_results_') and f.endswith('.xml'):
                    try:
                        batch_num = int(f.replace('blast_results_', '').replace('.xml', ''))
                        batch_files.append((batch_num, f))
                    except ValueError:
                        print(f"Warning: Cannot parse batch number from {f}")
                        continue
            
            if not batch_files:
                print("Warning: No BLAST XML files found in directory.")
                return {}
            
            # Sort by batch number (1-based)
            batch_files.sort(key=lambda x: x[0])
            print(f"Found batch files: {[f[0] for f in batch_files]}")
            
            # Infer batch size from first few valid files
            batch_size = self._infer_batch_size(blast_results_path, batch_files, total_sites)
            if batch_size is None:
                print("Warning: Could not infer batch size")
                return {}
            
            print(f"Inferred batch size: {batch_size}")
            
            # Calculate expected number of batches
            expected_batches = (total_sites + batch_size - 1) // batch_size  # Ceiling division
            print(f"Expected number of batches: {expected_batches}")
            
            # Identify problematic batches using the unified function
            problematic_batches = self._get_problematic_batches(blast_results_path, expected_batches)
            
            if problematic_batches:
                print(f"Problematic batches: {[i+1 for i in problematic_batches]}")
                print("Binding sites in problematic batches will be marked as missing BLAST data")
            
            # First, build a mapping from global site index to (gene_name, gene_internal_index, site)
            site_mapping = []
            for gene_name, sites in binding_sites.items():
                for gene_internal_index, site in enumerate(sites):
                    site_mapping.append({
                        'gene_name': gene_name,
                        'gene_internal_index': gene_internal_index,
                        'site': site
                    })
            
            # Process each batch file and map to binding sites
            for batch_num, filename in batch_files:
                xml_path = os.path.join(blast_results_path, filename)
                batch_blast_data = self._parse_single_batch_file(xml_path, batch_num)
                
                # Map this batch's results to the corresponding binding sites
                batch_start = (batch_num - 1) * batch_size
                batch_end = min(batch_start + batch_size, total_sites)
                
                # print(f"Processing batch {batch_num}: sites {batch_start+1}-{batch_end}")  # Reduced output
                
                # Process sites in this batch range
                for global_site_index in range(batch_start, batch_end):
                    if global_site_index < len(site_mapping):
                        site_info = site_mapping[global_site_index]
                        gene_name = site_info['gene_name']
                        gene_internal_index = site_info['gene_internal_index']
                        site = site_info['site']
                        
                        # Create query ID using the standard function with gene-internal index
                        query_id = self._generate_query_id(gene_name, gene_internal_index, site['position'], site['g_content'])
                        
                        if query_id in batch_blast_data:
                            blast_data[query_id] = batch_blast_data[query_id]
                        else:
                            # Mark as missing BLAST data
                            blast_data[query_id] = {
                                'alignments': [],
                                'evalue': float('inf'),
                                'identity': 0.0,
                                'alignment_count': 0,
                                'missing': True
                            }
            
            # Handle missing batches - mark all sites in missing batches as missing
            for problematic_batch_num in problematic_batches:
                batch_start = (problematic_batch_num - 1) * batch_size
                batch_end = min(batch_start + batch_size, total_sites)
                print(f"Marking sites {batch_start+1}-{batch_end} as problematic (batch {problematic_batch_num})")
                
                # Process sites in this batch range using site_mapping
                for global_site_index in range(batch_start, batch_end):
                    if global_site_index < len(site_mapping):
                        site_info = site_mapping[global_site_index]
                        gene_name = site_info['gene_name']
                        gene_internal_index = site_info['gene_internal_index']
                        site = site_info['site']
                        
                        # Create query ID using the standard function
                        query_id = self._generate_query_id(gene_name, gene_internal_index, site['position'], site['g_content'])
                        
                        blast_data[query_id] = {
                            'alignments': [],
                            'evalue': float('inf'),
                            'identity': 0.0,
                            'alignment_count': 0,
                            'missing': True
                        }
            
        else:
            # Single file mode
            blast_data = self._parse_single_batch_file(blast_results_path, 1)
        
        print(f"Parsed BLAST data for {len(blast_data)} queries")
        return blast_data
    
    def _parse_single_batch_file(self, xml_path: str, batch_num: int) -> Dict[str, Dict[str, Any]]:
        """Parse a single BLAST XML file and return the results."""
        blast_data = {}
        
        try:
            with open(xml_path, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            
            if not content or content == '<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>':
                print(f"Warning: Empty or minimal BLAST results in batch {batch_num}")
                return {}
            
            # Check for error messages in the XML
            if '<Iteration_message>' in content and 'Error:' in content:
                print(f"Warning: BLAST error found in batch {batch_num}")
                return {}
            
            from io import StringIO
            from Bio.Blast import NCBIXML
            records = NCBIXML.parse(StringIO(content))
            batch_records = 0
            
            for blast_record in records:
                batch_records += 1
                # Get query definition from the blast record
                # Use the raw query string to preserve exact formatting
                query_def = blast_record.query
                # The query_def should already be in the format: gene_name_index|pos=position|g_content=content
                # This matches what we expect in our binding sites
                query_id = query_def
                
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
            
            # print(f"Parsed {batch_records} records from batch {batch_num}")  # Reduced output
            
        except Exception as e:
            print(f"Error parsing BLAST XML batch {batch_num}: {e}")
            return {}
        
        return blast_data
    
    def _get_problematic_batches(self, results_dir: str, total_batches: int) -> List[int]:
        """Get list of batch indices that are missing, empty, or corrupted."""
        problematic_batches = []
        for i in range(total_batches):
            batch_file = os.path.join(results_dir, f"blast_results_{i+1:04d}.xml")
            if not os.path.exists(batch_file) or os.path.getsize(batch_file) <= 100:
                problematic_batches.append(i)
            else:
                # Check if file is corrupted
                try:
                    with open(batch_file, 'r', encoding='utf-8') as f:
                        content = f.read().strip()
                        if not content or content == '<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>':
                            problematic_batches.append(i)
                        elif '<Iteration_message>' in content and 'Error:' in content:
                            problematic_batches.append(i)
                except Exception:
                    problematic_batches.append(i)
        return problematic_batches

    def _infer_batch_size(self, results_dir: str, batch_files: List[tuple], total_sites: int) -> int:
        """Infer batch size from the first few valid batch files."""
        batch_sizes = []
        
        # Check first 3 valid files to infer batch size
        for batch_num, filename in batch_files[:3]:
            xml_path = os.path.join(results_dir, filename)
            try:
                with open(xml_path, 'r', encoding='utf-8') as f:
                    content = f.read().strip()
                
                if not content or content == '<?xml version="1.0"?>\n<BlastOutput>\n</BlastOutput>':
                    continue
                
                if '<Iteration_message>' in content and 'Error:' in content:
                    continue
                
                # Count iterations in this batch
                from io import StringIO
                from Bio.Blast import NCBIXML
                records = NCBIXML.parse(StringIO(content))
                batch_size = sum(1 for _ in records)
                batch_sizes.append(batch_size)
                # print(f"Batch {batch_num}: {batch_size} records")  # Reduced output
                
            except Exception as e:
                print(f"Error reading batch {batch_num}: {e}")
                continue
        
        if not batch_sizes:
            return None
        
        # Use the most common batch size (should be consistent)
        from collections import Counter
        size_counts = Counter(batch_sizes)
        inferred_size = size_counts.most_common(1)[0][0]
        
        print(f"Batch sizes found: {batch_sizes}")
        print(f"Inferred batch size: {inferred_size}")
        
        return inferred_size
    
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
                    'position': site['position'],  # local position
                    'g_content': site['g_content'],
                    'tm': site['tm'],
                    'strategy': site['strategy']
                }
                
                # Add genomic position information if available
                if 'genomic_start' in site:
                    site_data['genomic_start'] = site['genomic_start']
                if 'genomic_end' in site:
                    site_data['genomic_end'] = site['genomic_end']
                
                # Add strand information if available
                if 'strand' in site:
                    site_data['strand'] = site['strand']
                elif 'isoform_id' in site:
                    # Try to get strand from isoform info
                    strand = self._get_strand_from_isoform(site['isoform_id'], gene_name, output_dir)
                    if strand is not None:
                        site_data['strand'] = strand
                
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
    
    def _get_strand_from_isoform(self, isoform_id: str, gene_name: str, output_dir: str) -> Optional[int]:
        """Get strand information from isoform info file."""
        try:
            import json
            isoform_info_file = os.path.join(output_dir, "isoform_info", f"{gene_name}.json")
            if os.path.exists(isoform_info_file):
                with open(isoform_info_file, 'r', encoding='utf-8') as f:
                    isoform_data = json.load(f)
                    for isoform in isoform_data:
                        if isoform.get('id') == isoform_id or isoform.get('display_name') == isoform_id:
                            return isoform.get('strand')
        except Exception as e:
            print(f"Warning: Could not get strand for isoform {isoform_id}: {e}")
        return None


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
