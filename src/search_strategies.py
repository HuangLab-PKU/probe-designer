"""
Binding site search strategies
Implement different methods to search binding sites on sequences.
"""

import os
import json
import random
from typing import List, Dict, Any, Tuple, Optional
from Bio.SeqUtils import MeltingTemp as mt
from tqdm import tqdm

from .config import SearchConfig, FilterConfig


class SearchStrategy:
    """Base class for search strategies."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig):
        self.search_config = search_config
        self.filter_config = filter_config
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        """Search binding sites (must be implemented by subclasses)."""
        raise NotImplementedError
    
    def _pre_filter_sequence(self, sequence: str) -> bool:
        """Run quick pre-filters before more expensive checks."""
        g_content = sequence.count('G') / len(sequence)
        if not (self.filter_config.min_g_content <= g_content <= self.filter_config.max_g_content):
            return False
        if "G" * (self.filter_config.max_consecutive_g + 1) in sequence:
            return False
        return True
    
    def _calculate_tm(self, sequence: str) -> float:
        """Calculate melting temperature (Tm)."""
        try:
            return mt.Tm_NN(sequence, nn_table=mt.DNA_NN4)
        except:
            return 0.0
    
    def _check_tm_constraints(self, sequence: str) -> bool:
        """Check Tm constraints."""
        tm = self._calculate_tm(sequence)
        return self.filter_config.min_tm <= tm <= self.filter_config.max_tm


class BruteForceStrategy(SearchStrategy):
    """Search across the whole sequence with comprehensive filtering and optimization."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, genomic_context: Optional[Dict[str, Any]] = None):
        super().__init__(search_config, filter_config)
        # genomic_context: {'seq_region_name': str, 'start': int, 'end': int, 'strand': int}
        self.genomic_context = genomic_context or {}
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        """Search binding sites using brute force approach with comprehensive filtering."""
        if len(sequence) < self.search_config.binding_site_length:
            print(f"Gene {gene_name}: sequence too short ({len(sequence)} bp), skipping.")
            return []
        
        # Calculate search parameters
        bds_len = self.search_config.binding_site_length
        pin_gap = 0.1  # Skip 10% from start and end
        seq_gap = int(len(sequence) * pin_gap)
        
        # Generate candidate positions
        positions = list(range(seq_gap, len(sequence) - seq_gap - bds_len))
        
        # Find all valid candidates
        candidates = []
        for pos in tqdm(positions, desc=f"Searching {gene_name}", leave=False):
            target_seq = sequence[pos:pos + bds_len]
            
            # Apply comprehensive filtering
            if self._is_valid_binding_site_comprehensive(target_seq):
                tm_3 = self._calculate_tm_3_prime(target_seq)
                tm_5 = self._calculate_tm_5_prime(target_seq)
                tm_full = self._calculate_tm(target_seq)
                
                # Check Tm constraints for 3' and 5' halves
                if (self.filter_config.min_tm <= tm_3 <= self.filter_config.max_tm and 
                    self.filter_config.min_tm <= tm_5 <= self.filter_config.max_tm):
                    
                    # Check Tm difference between 3' and 5' halves
                    tm_diff = abs(tm_5 - tm_3)
                    if tm_diff <= 5:  # Tm difference threshold
                        
                        # Check RNA secondary structure (if available)
                        if self._check_rna_structure(target_seq):
                            candidates.append({
                                'pos': pos,
                                'sequence': target_seq,
                                'tm': tm_full,
                                'tm_3': tm_3,
                                'tm_5': tm_5,
                                'tm_diff': tm_diff,
                                'g_content': target_seq.count('G') / len(target_seq)
                            })
        
        # Optimize selection to get well-distributed sites
        selected_positions = self._optimize_subsequence(
            [c['pos'] for c in candidates], 
            self.search_config.max_binding_sites,
            min_gap=40,
            better_gap=80,
            gene=gene_name
        )
        
        # Build final results with genomic coordinates when context available
        binding_sites = []
        strand = self.genomic_context.get('strand')
        g_start = self.genomic_context.get('start')
        g_end = self.genomic_context.get('end')
        for candidate in candidates:
            if candidate['pos'] in selected_positions:
                rel_pos = candidate['pos']
                genomic_start = None
                genomic_end = None
                if strand in (1, -1) and isinstance(g_start, int) and isinstance(g_end, int):
                    if strand == 1:
                        genomic_start = g_start + rel_pos
                        genomic_end = genomic_start + bds_len
                    else:
                        genomic_end = g_end - rel_pos
                        genomic_start = genomic_end - bds_len
                binding_sites.append({
                    'gene_name': gene_name,
                    'sequence': candidate['sequence'],
                    'position': rel_pos,
                    'genomic_start': genomic_start,
                    'genomic_end': genomic_end,
                    'length': bds_len,
                    'g_content': candidate['g_content'],
                    'tm': candidate['tm'],
                    'tm_3': candidate['tm_3'],
                    'tm_5': candidate['tm_5'],
                    'tm_diff': candidate['tm_diff'],
                    'strategy': 'brute_force'
                })
        
        return binding_sites
    
    def _is_valid_binding_site_comprehensive(self, sequence: str) -> bool:
        """Comprehensive validation of binding site candidate."""
        # Check consecutive Gs
        if "G" * (self.filter_config.max_consecutive_g + 1) in sequence:
            return False
        
        # Check G content
        g_content = sequence.count('G') / len(sequence)
        if not (self.filter_config.min_g_content <= g_content <= self.filter_config.max_g_content):
            return False
        
        return True
    
    def _calculate_tm_3_prime(self, sequence: str) -> float:
        """Calculate Tm for 3' half of sequence."""
        try:
            half_len = len(sequence) // 2
            return mt.Tm_NN(sequence[half_len:], nn_table=mt.R_DNA_NN1)
        except:
            return 0.0
    
    def _calculate_tm_5_prime(self, sequence: str) -> float:
        """Calculate Tm for 5' half of sequence."""
        try:
            half_len = len(sequence) // 2
            return mt.Tm_NN(sequence[:half_len], nn_table=mt.R_DNA_NN1)
        except:
            return 0.0
    
    def _check_rna_structure(self, sequence: str) -> bool:
        """Check RNA secondary structure free energy."""
        try:
            import RNA
            _, mfe = RNA.fold(sequence)
            return mfe >= self.filter_config.min_free_energy
        except ImportError:
            # RNAfold not available, skip this check
            return True
        except Exception:
            return True
    
    def _optimize_subsequence(self, positions: List[int], length: int, min_gap: int, 
                            better_gap: int, gene: str) -> List[int]:
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


class ExonJunctionStrategy(SearchStrategy):
    """Search around exon-exon junctions (simplified placeholder)."""
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        binding_sites = []
        # In real use, exon coordinates are needed. Here we use a simple heuristic.
        seq_len = len(sequence)
        junction_positions = [seq_len // 3, 2 * seq_len // 3]
        for pos in junction_positions:
            sites = self._search_around_position(sequence, pos, gene_name)
            binding_sites.extend(sites)
        return binding_sites[:self.search_config.max_binding_sites]
    
    def _search_around_position(self, sequence: str, center_pos: int, gene_name: str) -> List[Dict[str, Any]]:
        sites = []
        search_range = 100
        bds_len = self.search_config.binding_site_length
        start_pos = max(0, center_pos - search_range)
        end_pos = min(len(sequence) - bds_len, center_pos + search_range)
        for pos in range(start_pos, end_pos, 5):
            if pos + bds_len > len(sequence):
                break
            candidate = sequence[pos:pos + bds_len]
            if self._is_valid_binding_site(candidate):
                sites.append({
                    'gene_name': gene_name,
                    'sequence': candidate,
                    'position': pos,
                    'length': bds_len,
                    'g_content': candidate.count('G') / len(candidate),
                    'tm': self._calculate_tm(candidate),
                    'strategy': 'exon_junction'
                })
        return sites
        

class IsoformSpecificStrategy(SearchStrategy):
    """Design probes unique to single isoforms using exon/intron boundary analysis and genome sequence access."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, isoforms: List[Dict], genome_accessor=None):
        super().__init__(search_config, filter_config)
        # Expect isoforms as list of dicts with exon coordinates and genome info
        # Each isoform should have: {'id', 'external_name', 'exons': [{'start', 'end'}], 'strand', 'seq_region_name', 'start', 'end'}
        self.isoforms = isoforms
        self.genome_accessor = genome_accessor  # Function to fetch genome sequence
        
        # Pre-process exon splice regions and counts
        if self.isoforms and self.genome_accessor:
            self._preprocess_exon_regions()
    
    def _preprocess_exon_regions(self):
        """Pre-process exon regions and calculate overlap counts."""
        # Collect all exon boundaries from all isoforms
        all_exon_boundaries = set()
        for iso in self.isoforms:
            for exon in iso.get('exons', []):
                all_exon_boundaries.add(exon['start'])
                all_exon_boundaries.add(exon['end'])
        
        split_points = sorted(all_exon_boundaries)
        
        # Calculate exon splice counts for each isoform
        exon_splice_counts = {}
        for iso in self.isoforms:
            iso_name = iso.get('external_name', iso.get('id'))
            exons = sorted(iso.get('exons', []), key=lambda x: x['start'])
            
            for exon in exons:
                exon_start, exon_end = exon['start'], exon['end']
                # Split exon at all boundaries
                cur_splits = [exon_start] + [p for p in split_points if exon_start < p < exon_end] + [exon_end]
                
                for j in range(len(cur_splits) - 1):
                    exon_splice = (cur_splits[j], cur_splits[j + 1])
                    if exon_splice not in exon_splice_counts:
                        exon_splice_counts[exon_splice] = []
                    exon_splice_counts[exon_splice].append(iso_name)
        
        self.exon_splice_counts = exon_splice_counts
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        """Find binding sites unique to single isoforms using exon/intron boundary analysis.
        
        Process:
        1. Use sliding window approach across each isoform
        2. For each window, determine which exon regions it overlaps
        3. Check if the overlapping regions are unique to a single isoform
        4. Apply comprehensive filtering and return probes with genomic positions
        """
        if not self.isoforms or not self.genome_accessor:
            # Fallback to brute force
            return BruteForceStrategy(self.search_config, self.filter_config).search_binding_sites(sequence, gene_name)
        
        # Get global coordinates for the gene locus
        global_start = min(iso['start'] for iso in self.isoforms)
        global_end = max(iso['end'] for iso in self.isoforms)
        
        # Fetch genome sequence for the entire gene region
        sample_isoform = self.isoforms[0]
        seq_region = sample_isoform.get('seq_region_name', '1')  # Default chromosome
        genome_seq = self.genome_accessor(seq_region, global_start, global_end)
        
        all_probes = []
        bds_len = self.search_config.binding_site_length
        tolerance_len = 10  # Tolerance length for isoform specificity check
        
        # Process each isoform
        for isoform in self.isoforms:
            isoform_probes = []
            isoform_start = isoform['start']
            isoform_end = isoform['end']
            
            # Get exon splices for this isoform
            exon_splices = list(self.exon_splice_counts.keys())
            
            # Slide window across the isoform
            current_pos = isoform_start
            while current_pos + bds_len <= isoform_end:
                # Find the end position considering exon boundaries
                end_pos = self._find_end_position(current_pos, bds_len, exon_splices)
                if end_pos >= isoform_end:
                    break
                
                # Get the exon regions that this window overlaps
                window_regions = self._get_overlapping_regions(current_pos, end_pos, exon_splices)
                
                # Check isoform specificity based on overlapping regions
                target_isoforms = self._get_target_isoforms(window_regions)
                
                # Only proceed if this window is specific to a single isoform
                if len(target_isoforms) == 1 and target_isoforms[0] == isoform.get('external_name', isoform.get('id')):
                    # Extract sequence from genome
                    local_start = current_pos - global_start
                    local_end = end_pos - global_start
                    window_seq = genome_seq[local_start:local_end]
                    
                    if isoform.get('strand', 1) == -1:
                        window_seq = self._reverse_complement(window_seq)
                    
                    # Apply comprehensive filtering
                    if self._is_valid_binding_site_comprehensive(window_seq):
                        tm_3 = self._calculate_tm_3_prime(window_seq)
                        tm_5 = self._calculate_tm_5_prime(window_seq)
                        tm_full = self._calculate_tm(window_seq)
                        
                        if (self.filter_config.min_tm <= tm_3 <= self.filter_config.max_tm and 
                            self.filter_config.min_tm <= tm_5 <= self.filter_config.max_tm):
                            
                            tm_diff = abs(tm_5 - tm_3)
                            if tm_diff <= 5 and self._check_rna_structure(window_seq):
                                isoform_probes.append({
                                    'gene_name': gene_name,
                                    'sequence': window_seq,
                                    'position': current_pos - isoform_start,  # Local position
                                    'genomic_start': current_pos,
                                    'genomic_end': end_pos,
                                    'length': bds_len,
                                    'g_content': window_seq.count('G') / len(window_seq),
                                    'tm': tm_full,
                                    'tm_3': tm_3,
                                    'tm_5': tm_5,
                                    'tm_diff': tm_diff,
                                    'strategy': 'isoform_specific',
                                    'target_isoforms': target_isoforms,
                                    'isoform_overlap_num': len(target_isoforms),
                                    'isoform_id': isoform.get('id'),
                                    'isoform_name': isoform.get('external_name', isoform.get('id')),
                                    'overlapping_regions': window_regions
                                })
                
                # Move to next position
                current_pos = self._find_next_position(current_pos, exon_splices)
            
            all_probes.extend(isoform_probes)
        
        # Optimize selection across all isoforms
        if len(all_probes) > self.search_config.max_binding_sites:
            # Sort by genomic position and optimize spacing
            all_probes.sort(key=lambda x: x['genomic_start'])
            positions = [p['genomic_start'] for p in all_probes]
            selected_positions = self._optimize_subsequence(
                positions,
                self.search_config.max_binding_sites,
                min_gap=40,
                better_gap=80,
                gene=gene_name
            )
            all_probes = [p for p in all_probes if p['genomic_start'] in selected_positions]
        
        return all_probes
    
    def _find_end_position(self, start_pos: int, length: int, exon_splices: List[Tuple[int, int]]) -> int:
        """Find the end position for a window of given length starting from start_pos."""
        current_pos = start_pos
        remaining_length = length
        
        for splice_start, splice_end in exon_splices:
            if splice_start <= current_pos < splice_end:
                # Current position is within this splice
                available_in_splice = splice_end - current_pos
                if remaining_length <= available_in_splice:
                    return current_pos + remaining_length
                else:
                    remaining_length -= available_in_splice
                    current_pos = splice_end
        
        return start_pos + length
    
    def _get_overlapping_regions(self, start_pos: int, end_pos: int, exon_splices: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """Get exon regions that overlap with the given window."""
        overlapping = []
        for splice_start, splice_end in exon_splices:
            if (splice_start < end_pos and splice_end > start_pos):
                overlapping.append((splice_start, splice_end))
        return overlapping
    
    def _get_target_isoforms(self, regions: List[Tuple[int, int]]) -> List[str]:
        """Get the list of isoforms that share all the given regions."""
        if not regions:
            return []
        
        # Start with isoforms that share the first region
        target_isoforms = set(self.exon_splice_counts.get(regions[0], []))
        
        # Intersect with isoforms that share other regions
        for region in regions[1:]:
            region_isoforms = set(self.exon_splice_counts.get(region, []))
            target_isoforms = target_isoforms.intersection(region_isoforms)
        
        return list(target_isoforms)
    
    def _find_next_position(self, current_pos: int, exon_splices: List[Tuple[int, int]]) -> int:
        """Find the next position to start the next window."""
        # Simple approach: move by 1 position
        return current_pos + 1
    
    def _reverse_complement(self, sequence: str) -> str:
        """Return reverse complement of a DNA sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(complement[base] for base in reversed(sequence))


class BindingSiteSearcher:
    """High-level facade to run a strategy across many genes."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, genome_accessor=None):
        self.search_config = search_config
        self.filter_config = filter_config
        self.genome_accessor = genome_accessor
        self.strategies = {
            'exon_junction': ExonJunctionStrategy,
            'brute_force': BruteForceStrategy,
            'isoform_specific': IsoformSpecificStrategy
        }
    
    def create_strategy(self, strategy_name: str, **kwargs) -> SearchStrategy:
        if strategy_name not in self.strategies:
            raise ValueError(f"Unsupported strategy: {strategy_name}")
        strategy_class = self.strategies[strategy_name]
        
        # Add genome_accessor for isoform_specific strategy
        if strategy_name == 'isoform_specific' and self.genome_accessor:
            kwargs['genome_accessor'] = self.genome_accessor
        
        # Forward genomic_context for brute_force if provided via kwargs
        return strategy_class(self.search_config, self.filter_config, **kwargs)
    
    def search_all_genes(self, sequences: Dict[str, Any], isoforms: Optional[Dict] = None) -> Dict[str, List[Dict[str, Any]]]:
        results = {}
        
        for gene_name, gene_data in tqdm(sequences.items(), desc="Search binding sites"):
            if 'sequence' not in gene_data:
                continue
            
            # Get sequence directly
            sequence = gene_data['sequence']
            
            # Create appropriate strategy
            if self.search_config.search_strategy == 'isoform_specific':
                if isoforms and gene_name in isoforms:
                    strategy = self.create_strategy('isoform_specific', isoforms=isoforms[gene_name])
                else:
                    strategy = self.create_strategy('brute_force')
            else:
                if self.search_config.search_strategy == 'brute_force':
                    # Assemble genomic context if available in gene_data
                    genomic_context = {
                        'seq_region_name': gene_data.get('seq_region_name'),
                        'start': gene_data.get('start'),
                        'end': gene_data.get('end'),
                        'strand': gene_data.get('strand')
                    }
                    strategy = self.create_strategy('brute_force', genomic_context=genomic_context)
                else:
                    strategy = self.create_strategy(self.search_config.search_strategy)
            
            # Search binding sites
            binding_sites = strategy.search_binding_sites(sequence, gene_name)
            results[gene_name] = binding_sites
        
        return results
    
    def save_binding_sites(self, binding_sites: Dict[str, List[Dict[str, Any]]], output_dir: str):
        os.makedirs(output_dir, exist_ok=True)
        fasta_file = os.path.join(output_dir, "total_binding_sites.fasta")
        with open(fasta_file, 'w') as f:
            for gene_name, sites in binding_sites.items():
                for i, site in enumerate(sites):
                    pos_field = f"pos={site['position']}"
                    gpos_field = f"gpos={site.get('genomic_start')}-{site.get('genomic_end')}"
                    f.write(f">{gene_name}_{i+1}|{pos_field}|{gpos_field}|g_content={site['g_content']:.2f}|tm={site['tm']:.1f}\n")
                    f.write(f"{site['sequence']}\n")
        json_file = os.path.join(output_dir, "binding_sites.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(binding_sites, f, indent=2, ensure_ascii=False)
        stats_file = os.path.join(output_dir, "binding_sites_stats.json")
        stats = {}
        for gene_name, sites in binding_sites.items():
            stats[gene_name] = {
                'total_sites': len(sites),
                'avg_g_content': sum(s['g_content'] for s in sites) / len(sites) if sites else 0,
                'avg_tm': sum(s['tm'] for s in sites) / len(sites) if sites else 0
            }
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(stats, f, indent=2, ensure_ascii=False)


def _is_valid_binding_site(self, sequence: str) -> bool:
    """Predicate: is the candidate a valid binding site?"""
    if not self._pre_filter_sequence(sequence):
        return False
    if not self._check_tm_constraints(sequence):
        return False
    return True


SearchStrategy._is_valid_binding_site = _is_valid_binding_site
