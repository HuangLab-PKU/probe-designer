"""
Binding site search strategies
Implement different methods to search binding sites on sequences.
"""

import os
import sys
import json
import random  # noqa: F401 (reserved for future stochastic strategies)
from typing import List, Dict, Any, Tuple, Optional
from tqdm import tqdm  # noqa: F401 (progress bars may be enabled later)

from .config import SearchConfig, FilterConfig
from .filtering import SequenceFilter


class SearchStrategy:
    """Base class for search strategies."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, show_progress: bool = False):
        self.search_config = search_config
        self.filter_config = filter_config
        self.show_progress = show_progress
        # Create filter instance for thermal screening
        from .config import BlastConfig
        self.filter = SequenceFilter(filter_config, BlastConfig())
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        """Search binding sites (must be implemented by subclasses)."""
        raise NotImplementedError
    
    def _reverse_complement(self, sequence: str) -> str:
        """Return reverse complement of a DNA sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(complement[base] for base in reversed(sequence))
    
    def _optimize_subsequence(self, positions: List[int], num_required: int, min_gap: int, 
                            better_gap: int, gene: str) -> List[int]:
        """Optimize subsequence selection to maximize distance between sites."""
        if len(positions) == 0:
            print(f"Gene {gene}: No valid positions found, please loosen threshold conditions.")
            return []
        
        if len(positions) <= num_required:
            return positions
        
        positions.sort()
        
        def is_valid(min_difference: int, num_required: int, min_gap: int) -> bool:
            count = 1
            current_min = positions[0]
            
            for i in range(1, len(positions)):
                if positions[i] - current_min >= min_difference:
                    count += 1
                    current_min = positions[i]
            
            return count >= num_required and min_difference > min_gap
        
        # Binary search for optimal spacing
        left, right = 0, positions[-1] - positions[0]
        result = []
        
        while left <= right:
            mid = (left + right) // 2
            if is_valid(mid, num_required, min_gap):
                result = [positions[0]]
                current_min = positions[0]
                
                for i in range(1, len(positions)):
                    if positions[i] - current_min >= mid:
                        result.append(positions[i])
                        current_min = positions[i]
                        if len(result) >= num_required:
                            break
                
                left = mid + 1
            else:
                right = mid - 1
        
        if not result:
            print(f"Gene {gene}: Not enough positions for {num_required} binding sites.")
            result = positions[:num_required]
        
        if mid < better_gap:
            print(f"Gene {gene}: Conditions too harsh, consider loosening parameters for better results")
        
        return result


class BruteForceStrategy(SearchStrategy):
    """Search across the whole sequence with comprehensive filtering and optimization."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, genomic_context: Optional[Dict[str, Any]] = None, show_progress: bool = False):
        super().__init__(search_config, filter_config, show_progress)
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
        
        # Setup progress bar if requested
        use_progress = self.show_progress
        # Create iterator with or without progress bar
        if use_progress:
            from tqdm import tqdm
            position_iterator = tqdm(positions, desc=f"Searching {gene_name}", unit="pos")
        else:
            position_iterator = positions
        
        for pos in position_iterator:
            target_seq = sequence[pos:pos + bds_len]
            # Get reverse complement for probe design (DNA probe binding to RNA target)
            probe_seq = self._reverse_complement(target_seq)
            
            # Split probe sequence into arms (3' arm is front, 5' arm is back for padlock probe)
            arm_3prime = probe_seq[:bds_len//2]  # Front half
            arm_5prime = probe_seq[bds_len//2:]  # Back half
            
            # Use thermal_filter for screening
            thermal_result = self.filter.thermal_filter(
                arm_3prime,
                arm_5prime,
                sequence_type="DNA",
                target_type="RNA",
                target_sequence=target_seq
            )
            
            if thermal_result['passed']:
                candidates.append({
                    'pos': pos,
                    'sequence': probe_seq,  # Store the probe sequence (reverse complement)
                    'target_sequence': target_seq,  # Store original target sequence
                    'arm_3prime': thermal_result['arm_3prime'],
                    'arm_5prime': thermal_result['arm_5prime'],
                    'tm': thermal_result['tm'],
                    'tm_3': thermal_result['tm_3prime'],
                    'tm_5': thermal_result['tm_5prime'],
                    'tm_diff': thermal_result['tm_diff'],
                    'g_content': thermal_result['g_content'],
                    'free_energy': thermal_result['free_energy']
                })
        
        # Optimize selection using subsequence optimization
        if candidates:
            positions = [c['pos'] for c in candidates]
            optimized_positions = self._optimize_subsequence(
                positions, 
                self.search_config.max_binding_sites,
                min_gap=1,  # Updated gap parameters
                better_gap=1,
                gene=gene_name
            )
            
            # Build final results with genomic coordinates when context available
            binding_sites = []
            strand = self.genomic_context.get('strand')
            g_start = self.genomic_context.get('start')
            g_end = self.genomic_context.get('end')
            
            for candidate in candidates:
                if candidate['pos'] in optimized_positions:
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
                        'st': genomic_start,
                        'en': genomic_end,
                        'length': bds_len,
                        'g_content': candidate['g_content'],
                        'tm': candidate['tm'],
                        'tm_3': candidate['tm_3'],
                        'tm_5': candidate['tm_5'],
                        'tm_diff': candidate['tm_diff'],
                        'free_energy': candidate['free_energy'],
                        'strategy': 'brute_force'
                    })
            
            return binding_sites
        
        return []


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
            target_seq = sequence[pos:pos + bds_len]
            # Get reverse complement for probe design (DNA probe binding to RNA target)
            probe_seq = self._reverse_complement(target_seq)
            
            # Split probe sequence into arms (3' arm is front, 5' arm is back for padlock probe)
            arm_3prime = probe_seq[:bds_len//2]  # Front half
            arm_5prime = probe_seq[bds_len//2:]  # Back half
            
            # Use thermal_filter for screening
            thermal_result = self.filter.thermal_filter(
                arm_3prime,
                arm_5prime,
                sequence_type="DNA",
                target_type="RNA",
                target_sequence=target_seq
            )
            
            if thermal_result['passed']:
                sites.append({
                    'gene_name': gene_name,
                    'sequence': probe_seq,  # Store the probe sequence (reverse complement)
                    'target_sequence': target_seq,  # Store original target sequence
                    'arm_3prime': thermal_result['arm_3prime'],
                    'arm_5prime': thermal_result['arm_5prime'],
                    'length': bds_len,
                    'g_content': thermal_result['g_content'],
                    'tm': thermal_result['tm'],
                    'tm_3': thermal_result['tm_3prime'],
                    'tm_5': thermal_result['tm_5prime'],
                    'tm_diff': thermal_result['tm_diff'],
                    'free_energy': thermal_result['free_energy'],
                    'strategy': 'exon_junction'
                })
        return sites
        

class IsoformAwareStrategy(SearchStrategy):
    """Unified strategy for isoform-aware probe design.
    
    Supports two modes:
    - 'specific': Design probes unique to single isoforms (uses isoforms_intersection)
    - 'consensus': Design probes that bind as many isoforms as possible (uses isoforms_union)
    """
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, isoforms: List[Dict], 
                 genome_accessor=None, mode: str = 'consensus', show_progress: bool = False):
        super().__init__(search_config, filter_config, show_progress)
        if mode not in ('specific', 'consensus'):
            raise ValueError(f"Mode must be 'specific' or 'consensus', got '{mode}'")
        self.mode = mode
        self.isoforms = isoforms
        self.genome_accessor = genome_accessor
        self.awareness = IsoformAwareness(self.isoforms, self.genome_accessor) if (self.isoforms and self.genome_accessor) else None
        if not self.awareness:
            raise ValueError(f"IsoformAwareStrategy requires isoforms and genome_accessor")
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        """Search binding sites using isoform-aware approach."""
        # Initialize genome sequence (may take time on first call)
        _ = self.awareness.get_genome_seq()
        
        bds_len = self.search_config.binding_site_length
        candidates: List[Dict[str, Any]] = []
        seen_positions = set() if self.mode == 'consensus' else None  # Only track duplicates for consensus mode
        
        # Setup progress bars if requested
        if self.show_progress:
            isoform_pbar = tqdm(self.isoforms, desc=f"{gene_name} (isoforms)", unit="isoform", leave=True, dynamic_ncols=True, file=sys.stdout)
        else:
            isoform_pbar = None
        
        for isoform in (isoform_pbar if isoform_pbar else self.isoforms):
            iso_start = isoform['start']
            iso_end = isoform['end']
            strand = isoform['strand']
            isoform_display_name = isoform['display_name']
            
            # Calculate positions for this isoform (only in exon regions, excluding introns)
            if self.show_progress:
                # Sum all exon lengths (exons use inclusive coordinates, so length = end - start + 1)
                total_exon_length = sum(exon['end'] - exon['start'] + 1 for exon in isoform['Exon'])
                isoform_positions = max(1, total_exon_length - bds_len + 1)  # Ensure at least 1 for progress bar
                isoform_pbar.set_postfix_str(f"Processing {isoform_display_name}")
                position_pbar = tqdm(total=isoform_positions, desc=f"  {isoform_display_name}", unit="pos", leave=False, dynamic_ncols=True, file=sys.stdout)
            else:
                isoform_positions = None
                position_pbar = None
            
            # Get sorted exons for this isoform
            exons = sorted(isoform['Exon'], key=lambda x: x['start'])
            current_exon_idx = 0
            current_pos = None
            
            # Find the first exon that contains or is after iso_start
            for idx, exon in enumerate(exons):
                if exon['end'] >= iso_start:
                    current_exon_idx = idx
                    current_pos = max(iso_start, exon['start'])
                    break
            
            # If no valid exon found, skip this isoform
            if current_pos is None:
                if position_pbar:
                    position_pbar.close()
                continue
            
            while current_exon_idx < len(exons):
                current_exon = exons[current_exon_idx]
                exon_start = current_exon['start']
                exon_end = current_exon['end']
                
                # Process positions within current exon
                while current_pos <= exon_end - bds_len + 1:
                    # For specific mode, constrain to current isoform's exons; for consensus, use union
                    # Pass current_exon_idx to start searching from the current exon
                    target_regions = self.awareness.find_target_region(current_pos, bds_len, isoform, exon_idx=current_exon_idx)
                    if not target_regions:
                        # If find_target_region returns empty, we've reached the end of valid positions
                        break
                    
                    # Get the last position of the target region for boundary checking
                    # target_regions[-1][1] is inclusive
                    end_pos_inclusive = target_regions[-1][1]
                    if end_pos_inclusive > iso_end:
                        break
                    
                    # Get overlapping regions for all segments in target_regions
                    regions = []
                    for reg_start, reg_end in target_regions:
                        # reg_start and reg_end are both inclusive
                        regions.extend(self.awareness.get_overlapping_regions(reg_start, reg_end))
                    
                    # Determine target isoforms based on mode
                    if self.mode == 'specific':
                        target_isoforms = self.awareness.isoforms_intersection(regions)
                        # Only proceed if this position is unique to the current isoform
                        if len(target_isoforms) != 1 or target_isoforms[0] != isoform_display_name:
                            current_pos += 1
                            if position_pbar:
                                position_pbar.update(1)
                            continue
                        covered = target_isoforms
                    elif self.mode == 'consensus':  # consensus mode
                        covered = set(self.awareness.isoforms_union(regions))
                        if not covered:
                            current_pos += 1
                            if position_pbar: 
                                position_pbar.update(1)
                            continue
                    
                    # Get target sequence from regions
                    target_seq = self.awareness.get_target_seq(target_regions, strand)
                    
                    # Get reverse complement for probe design (DNA probe binding to RNA target)
                    probe_seq = self._reverse_complement(target_seq)
                    
                    # Split probe sequence into arms (3' arm is front, 5' arm is back for padlock probe)
                    arm_3prime = probe_seq[:bds_len//2]  # Front half
                    arm_5prime = probe_seq[bds_len//2:]  # Back half
                    
                    # Use thermal_filter for screening
                    thermal_result = self.filter.thermal_filter(
                        arm_3prime,
                        arm_5prime,
                        sequence_type="DNA",
                        target_type="RNA",
                        target_sequence=target_seq
                    )
                    
                    if not thermal_result['passed']:
                        current_pos += 1
                        if position_pbar:
                            position_pbar.update(1)
                        continue
                    
                    # Calculate genomic position
                    # end_pos_inclusive is the last inclusive position of the target region
                    genomic_start = current_pos if strand == 1 else end_pos_inclusive - bds_len + 1
                    genomic_end = current_pos + bds_len - 1 if strand == 1 else end_pos_inclusive
                    
                    # Skip if we've already seen this genomic position (consensus mode only)
                    if seen_positions is not None:
                        if (genomic_start, genomic_end) in seen_positions:
                            current_pos += 1
                            if position_pbar:
                                position_pbar.update(1)
                            continue
                        seen_positions.add((genomic_start, genomic_end))
                    
                    # Build probe dictionary
                    probe_dict = {
                        'gene_name': gene_name,
                        'sequence': probe_seq,  # Store the probe sequence (reverse complement)
                        'target_sequence': target_seq,  # Store original target sequence
                        'arm_3prime': thermal_result['arm_3prime'],
                        'arm_5prime': thermal_result['arm_5prime'],
                        'st': genomic_start,
                        'en': genomic_end,
                        'target_regions': target_regions,
                        'length': bds_len,
                        'g_content': thermal_result['g_content'],
                        'tm': thermal_result['tm'],
                        'tm_3': thermal_result['tm_3prime'],
                        'tm_5': thermal_result['tm_5prime'],
                        'tm_diff': thermal_result['tm_diff'],
                        'free_energy': thermal_result['free_energy'],
                        'strategy': f'isoform_{self.mode}',
                        'isoform_id': isoform['id'],
                        'isoform_name': isoform['display_name'],
                        'overlapping_regions': regions
                    }
                    
                    # Add mode-specific fields
                    if self.mode == 'specific':
                        probe_dict['target_isoforms'] = covered
                        probe_dict['isoform_overlap_num'] = len(covered)
                        probe_dict['segments'] = target_regions  # Alias for target_regions
                    else:  # consensus
                        probe_dict['target_isoforms'] = sorted(list(covered))
                        probe_dict['isoform_overlap_num'] = len(covered)
                    
                    candidates.append(probe_dict)
                    current_pos += 1
                    if position_pbar:
                        position_pbar.update(1)
                        position_pbar.set_postfix_str(f"Found {len(candidates)} candidates")
                
                # Move to next exon
                current_exon_idx += 1
                if current_exon_idx < len(exons):
                    current_pos = exons[current_exon_idx]['start']
                else:
                    break
            
            if position_pbar:
                position_pbar.close()
            
            if isoform_pbar:
                isoform_pbar.set_postfix_str(f"Found {len(candidates)} candidates total")
        
        if isoform_pbar:
            isoform_pbar.close()
        
        # Optimize selection
        if self.mode == 'consensus':
            # Prefer higher coverage then spacing optimization
            candidates.sort(key=lambda x: (-x['isoform_overlap_num'], x['st']))
            
            # If we have many candidates, pre-filter to keep only the highest coverage ones
            # This ensures that spacing optimization (which ignores coverage) operates only on high-quality candidates
            if len(candidates) > self.search_config.max_binding_sites:
                # Keep top 3x candidates to give optimizer room for spacing
                # This prevents low-coverage well-spaced sites from displacing high-coverage sites
                pool_size = min(len(candidates), self.search_config.max_binding_sites * 3)
                candidates = candidates[:pool_size]
                
        else:
            # For specific mode, just sort by genomic position
            candidates.sort(key=lambda x: x['st'])
        
        if len(candidates) > self.search_config.max_binding_sites:
            # Optimization needs positions sorted by coordinate
            candidates.sort(key=lambda x: x['st'])
            positions = [c['st'] for c in candidates]
            selected = self._optimize_subsequence(positions, self.search_config.max_binding_sites, 0, 0, gene_name)
            
            # Filter candidates to keep only selected positions
            # Note: We need to handle potential duplicate positions if they exist (unlikely in same gene same strand)
            selected_set = set(selected)
            candidates = [c for c in candidates if c['st'] in selected_set]
            
            # If consensus mode, re-sort by coverage for final output
            if self.mode == 'consensus':
                candidates.sort(key=lambda x: (-x['isoform_overlap_num'], x['st']))
        
        return candidates


# Backward compatibility aliases
class IsoformSpecificStrategy(IsoformAwareStrategy):
    """Design probes unique to single isoforms (backward compatibility alias)."""
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, isoforms: List[Dict], 
                 genome_accessor=None, show_progress: bool = False):
        super().__init__(search_config, filter_config, isoforms, genome_accessor, mode='specific', show_progress=show_progress)


class IsoformConsensusStrategy(IsoformAwareStrategy):
    """Design probes that bind as many isoforms as possible (backward compatibility alias)."""
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, isoforms: List[Dict], 
                 genome_accessor=None, show_progress: bool = False):
        super().__init__(search_config, filter_config, isoforms, genome_accessor, mode='consensus', show_progress=show_progress)


class BindingSiteSearcher:
    """High-level facade to run a strategy across many genes."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, genome_accessor=None, show_progress: bool = False):
        self.search_config = search_config
        self.filter_config = filter_config
        self.genome_accessor = genome_accessor
        self.show_progress = show_progress
        self.strategies = {
            'exon_junction': ExonJunctionStrategy,
            'brute_force': BruteForceStrategy,
            'isoform_specific': IsoformSpecificStrategy,
            'isoform_consensus': IsoformConsensusStrategy
        }
    
    def create_strategy(self, strategy_name: str, **kwargs) -> SearchStrategy:
        if strategy_name not in self.strategies:
            raise ValueError(f"Unsupported strategy: {strategy_name}")
        strategy_class = self.strategies[strategy_name]
        
        # Add genome_accessor for isoform-based strategies
        if strategy_name in ('isoform_specific', 'isoform_consensus') and self.genome_accessor:
            kwargs['genome_accessor'] = self.genome_accessor
        
        # Forward show_progress to strategy
        if 'show_progress' not in kwargs:
            kwargs['show_progress'] = self.show_progress
        
        # Forward genomic_context for brute_force if provided via kwargs
        return strategy_class(self.search_config, self.filter_config, **kwargs)
    
    def search_all_genes(self, sequences: Optional[Dict[str, Any]] = None, isoforms: Optional[Dict] = None) -> Dict[str, List[Dict[str, Any]]]:
        """
        Search binding sites for genes.
        
        Args:
            sequences: Optional dict of {gene_name: {'sequence': str, ...}} for strategies that need direct sequence access (e.g., brute_force)
            isoforms: Optional dict of {gene_name: [isoform_dicts]} for strategies that use isoform information (e.g., isoform_consensus)
        
        Returns:
            Dict of {gene_name: [binding_site_dicts]}
        """
        results = {}
        
        # Determine which genes to process based on strategy
        if self.search_config.search_strategy in ('isoform_consensus', 'isoform_specific'):
            # For isoform-based strategies, use isoforms dict
            if isoforms is None:
                raise ValueError(f"Strategy '{self.search_config.search_strategy}' requires isoforms parameter")
            genes_to_process = list(isoforms.keys())
        else:
            # For sequence-based strategies (brute_force, exon_junction), use sequences dict
            if sequences is None:
                raise ValueError(f"Strategy '{self.search_config.search_strategy}' requires sequences parameter")
            genes_to_process = list(sequences.keys())
        
        total_genes = len(genes_to_process)
        
        # Don't show progress bar at this level - let individual strategies handle it
        # This avoids conflicts with nested progress bars in strategies
        
        # Iterate through genes
        for gene_name in genes_to_process:
            # Create appropriate strategy
            if self.search_config.search_strategy == 'isoform_specific':
                if isoforms and gene_name in isoforms:
                    strategy = self.create_strategy('isoform_specific', isoforms=isoforms[gene_name])
                else:
                    # Fallback to brute_force if no isoforms
                    if sequences and gene_name in sequences and 'sequence' in sequences[gene_name]:
                        sequence = sequences[gene_name]['sequence']
                        strategy = self.create_strategy('brute_force')
                    else:
                        if not self.show_progress:
                            print(f"Skipping {gene_name}: no isoforms and no sequence")
                        continue
            elif self.search_config.search_strategy == 'isoform_consensus':
                if isoforms and gene_name in isoforms:
                    strategy = self.create_strategy('isoform_consensus', isoforms=isoforms[gene_name], genome_accessor=self.genome_accessor)
                else:
                    # Fallback to brute_force if no isoforms
                    if sequences and gene_name in sequences and 'sequence' in sequences[gene_name]:
                        sequence = sequences[gene_name]['sequence']
                        strategy = self.create_strategy('brute_force')
                    else:
                        if not self.show_progress:
                            print(f"Skipping {gene_name}: no isoforms and no sequence")
                        continue
            else:
                # Sequence-based strategies (brute_force, exon_junction)
                if sequences is None or gene_name not in sequences:
                    if not self.show_progress:
                        print(f"Skipping {gene_name}: no sequence data")
                    continue
                
                gene_data = sequences[gene_name]
                if 'sequence' not in gene_data:
                    if not self.show_progress:
                        print(f"Skipping {gene_name}: no sequence")
                    continue
                
                sequence = gene_data['sequence']
                
                if self.search_config.search_strategy == 'brute_force':
                    # Assemble genomic context if available in gene_data
                    genomic_context = {
                        'seq_region_name': gene_data['seq_region_name'],
                        'start': gene_data['start'],
                        'end': gene_data['end'],
                        'strand': gene_data['strand']
                    }
                    strategy = self.create_strategy('brute_force', genomic_context=genomic_context)
                else:
                    strategy = self.create_strategy(self.search_config.search_strategy)
            
            # Search binding sites (progress bars are handled inside strategies)
            # For isoform-based strategies, sequence is not needed (will be fetched via genome_accessor)
            # For sequence-based strategies, pass the sequence
            if self.search_config.search_strategy in ('isoform_consensus', 'isoform_specific'):
                # For isoform strategies, pass placeholder - sequence will be fetched internally via genome_accessor
                binding_sites = strategy.search_binding_sites('', gene_name)
            else:
                # For sequence-based strategies, pass the actual sequence
                binding_sites = strategy.search_binding_sites(sequence, gene_name)
            
            results[gene_name] = binding_sites
        
        return results
    
    def save_binding_sites(self, binding_sites: Dict[str, List[Dict[str, Any]]], output_dir: str):
        os.makedirs(output_dir, exist_ok=True)
        fasta_file = os.path.join(output_dir, "total_binding_sites.fasta")
        with open(fasta_file, 'w') as f:
            for gene_name, sites in binding_sites.items():
                for i, site in enumerate(sites):
                    gpos_field = f"gpos={site['st']}-{site['en']}"
                    f.write(f">{gene_name}_{i+1}|{gpos_field}|g_content={site['g_content']:.2f}|tm={site['tm']:.1f}\n")
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


# --- Shared isoform awareness helpers ---

class IsoformAwareness:
    """Reusable isoform-aware preprocessing over exon splices and genome locus."""
    
    def __init__(self, isoforms: List[Dict], genome_accessor):
        self.isoforms = isoforms
        self.genome_accessor = genome_accessor
        # locus
        self.global_start = min(iso['start'] for iso in self.isoforms)
        self.global_end = max(iso['end'] for iso in self.isoforms)
        self.seq_region_name = self.isoforms[0]['seq_region_name']
        # splices
        self.exon_splice_to_isoforms = self._build_exon_splice_to_isoforms()
        self._genome_seq_cache: Optional[str] = None
    
    def _build_exon_splice_to_isoforms(self) -> Dict[Tuple[int, int], List[str]]:
        """Build mapping from exon splice (start,end) to list of isoform names that contain it."""
        all_exon_boundaries = set()
        for iso in self.isoforms:
            for exon in iso['Exon']:
                all_exon_boundaries.add(exon['start'])
                all_exon_boundaries.add(exon['end'])
        split_points = sorted(all_exon_boundaries)

        exon_splice_to_isoforms: Dict[Tuple[int, int], List[str]] = {}
        for iso in self.isoforms:
            isoform_name = iso['display_name']
            exons = sorted(iso['Exon'], key=lambda x: x['start'])
            for exon in exons:
                exon_start, exon_end = exon['start'], exon['end']
                cur_splits = [exon_start] + [p for p in split_points if exon_start < p < exon_end] + [exon_end]
                for j in range(len(cur_splits) - 1):
                    exon_splice = (cur_splits[j], cur_splits[j + 1])
                    exon_splice_to_isoforms.setdefault(exon_splice, []).append(isoform_name)
        return exon_splice_to_isoforms
    
    def get_genome_seq(self) -> str:
        if self._genome_seq_cache is None:
            self._genome_seq_cache = self.genome_accessor(self.seq_region_name, self.global_start, self.global_end)
        return self._genome_seq_cache

    def find_target_region(self, start_pos: int, length: int, isoform: Optional[Dict[str, Any]] = None, exon_idx: Optional[int] = None) -> List[Tuple[int, int]]:
        """Find probe target region as a list of continuous segments after traversing 'length' bases along exonic regions starting at 'start_pos'.
        Returns a list of tuples (start, end) representing the target region, which may be discontinuous due to alternative splicing.
        Both start and end are inclusive (Ensembl-style coordinates).
        If 'isoform' is provided, constrain traversal to that isoform's exons only; otherwise, use union of exons.
        
        Args:
            start_pos: Starting genomic position
            length: Length of target region to find
            isoform: Optional isoform dict to constrain search to its exons
            exon_idx: Optional index to start from in the sorted exon list (0-based). If provided, traversal starts from this exon.
        """
        exons = sorted(isoform['Exon'], key=lambda x: x['start'])
    
        current_pos = start_pos
        remain = length
        target_regions = []  # List of (start, end) tuples, both inclusive
        
        # If exon_idx is provided, start from that exon; otherwise start from the beginning
        exon_iter = exons[exon_idx:] if exon_idx is not None and 0 <= exon_idx < len(exons) else exons
        
        for exon in exon_iter:
            seg_start, seg_end = exon['start'], exon['end']
            # seg_start and seg_end are inclusive (from Ensembl)
            if seg_end < current_pos:
                continue
            # If current_pos is in an intron before this exon, jump to exon start
            if seg_start > current_pos:
                current_pos = seg_start
            
            # Now consume within exon segment
            # seg_end is inclusive, so available length = seg_end - current_pos + 1
            avail = max(0, seg_end - current_pos + 1)
            if avail <= 0:
                continue
            
            if remain <= avail:
                # This segment completes the target region
                # Store as (start, end) where both are inclusive
                target_regions.append((current_pos, current_pos + remain - 1))
                return target_regions
            else:
                # Consume the entire segment and continue
                # Store seg_end as inclusive
                target_regions.append((current_pos, seg_end))
                remain -= avail
                current_pos = seg_end + 1  # Move to next position after this segment

        # Fallback: if exons exhausted, return naive region (start_pos, start_pos + length - 1)
        if not target_regions:
            target_regions.append((start_pos, start_pos + length - 1))
        return target_regions
    
    def get_target_seq(self, regions: List[Tuple[int, int]], strand: int) -> str:
        """Get target sequence from given genomic regions.
        
        Args:
            regions: List of (start, end) tuples representing genomic segments (both inclusive)
            strand: Strand direction (1 for positive, -1 for negative)
            
        Returns:
            Target sequence from the regions (reverse complemented if negative strand)
        """
        if not regions:
            return ""
        
        genome_seq = self.get_genome_seq()
        seq_parts: List[str] = []
        
        # Validate genome sequence length
        # Ensembl/genome accessors use inclusive coordinates, so length = end - start + 1
        expected_length = self.global_end - self.global_start + 1
        if len(genome_seq) != expected_length:
            raise ValueError(f"Genome sequence length mismatch: expected {expected_length}, got {len(genome_seq)}")
        
        # Extract sequence from each region
        # regions are stored as (start, end) where both are inclusive
        # Convert to Python slicing (exclusive end) for genome_seq access
        for reg_start, reg_end in regions:
            ls = reg_start - self.global_start
            # reg_end is inclusive, so for Python slicing we need reg_end + 1
            le = reg_end - self.global_start + 1
            
            # Boundary check
            if ls < 0 or le > len(genome_seq):
                raise IndexError(f"Index out of range: ls={ls}, le={le}, genome_seq_len={len(genome_seq)}, "
                               f"reg_start={reg_start}, reg_end={reg_end}, global_start={self.global_start}")
            
            seq_parts.append(genome_seq[ls:le])
        
        seq = ''.join(seq_parts)
        
        if strand == -1:
            # Reverse complement for negative strand
            comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
            seq = ''.join(comp.get(b, 'N') for b in reversed(seq))
        
        return seq
    
    def get_overlapping_regions(self, start_pos: int, end_pos: int) -> List[Tuple[int, int]]:
        """
        Find overlapping exon splice regions.
        
        Args:
            start_pos: Start position (inclusive)
            end_pos: End position (inclusive)
        
        Returns:
            List of (start, end) tuples from exon_splice_to_isoforms (inclusive coordinates)
        """
        # exon_splice_to_isoforms keys are (start, end) in inclusive coordinates
        # target_regions are (start, end) in inclusive coordinates
        # Overlap condition: [s, e] (inclusive) overlaps [start_pos, end_pos] (inclusive)
        # if s <= end_pos and e >= start_pos
        return [(s, e) for (s, e) in self.exon_splice_to_isoforms.keys() if s <= end_pos and e >= start_pos]
    
    def isoforms_intersection(self, regions: List[Tuple[int, int]]) -> List[str]:
        if not regions:
            return []
        inter = set(self.exon_splice_to_isoforms.get(regions[0], []))
        for r in regions[1:]:
            inter &= set(self.exon_splice_to_isoforms.get(r, []))
        return list(inter)
    
    def isoforms_union(self, regions: List[Tuple[int, int]]) -> List[str]:
        covered = set()
        for r in regions:
            covered.update(self.exon_splice_to_isoforms.get(r, []))
        return list(covered)
