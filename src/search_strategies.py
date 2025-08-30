"""
Binding site search strategies
Implement different methods to search binding sites on sequences.
"""

import os
import json
import random
from typing import List, Dict, Any, Tuple, Optional
from tqdm import tqdm

from .config import SearchConfig, FilterConfig
from .filtering import SequenceFilter


class SearchStrategy:
    """Base class for search strategies."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig):
        self.search_config = search_config
        self.filter_config = filter_config
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
        print(f"Searching {len(positions)} positions for {gene_name}...")
        for pos in positions:
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
                        'position': rel_pos,
                        'genomic_start': genomic_start,
                        'genomic_end': genomic_end,
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
                    'position': pos,
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
        

class IsoformSpecificStrategy(SearchStrategy):
    """Design probes unique to single isoforms using exon/intron boundary analysis and genome sequence access."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, isoforms: List[Dict], genome_accessor=None):
        super().__init__(search_config, filter_config)
        self.isoforms = isoforms
        self.genome_accessor = genome_accessor
        self.awareness = IsoformAwareness(self.isoforms, self.genome_accessor) if (self.isoforms and self.genome_accessor) else None
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        if not self.awareness:
            return BruteForceStrategy(self.search_config, self.filter_config).search_binding_sites(sequence, gene_name)
        
        genome_seq = self.awareness.get_genome_seq()
        bds_len = self.search_config.binding_site_length
        all_probes: List[Dict[str, Any]] = []
        
        for isoform in self.isoforms:
            isoform_probes = []
            iso_start = isoform['start']
            iso_end = isoform['end']
            strand = isoform.get('strand', 1)
            current_pos = iso_start
            while current_pos + bds_len <= iso_end:
                end_pos = self.awareness.find_end_position(current_pos, bds_len)
                if end_pos > iso_end:
                    break
                regions = self.awareness.get_overlapping_regions(current_pos, end_pos)
                target_isoforms = self.awareness.isoforms_intersection(regions)
                if len(target_isoforms) == 1 and target_isoforms[0] == isoform.get('external_name', isoform.get('id')):
                    ls = current_pos - self.awareness.global_start
                    le = end_pos - self.awareness.global_start
                    target_seq = genome_seq[ls:le]
                    if strand == -1:
                        target_seq = self._reverse_complement(target_seq)
                    
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
                        all_probes.append({
                            'gene_name': gene_name,
                            'sequence': probe_seq,  # Store the probe sequence (reverse complement)
                            'target_sequence': target_seq,  # Store original target sequence
                            'arm_3prime': thermal_result['arm_3prime'],
                            'arm_5prime': thermal_result['arm_5prime'],
                            'position': current_pos - iso_start,
                            'genomic_start': current_pos if strand == 1 else end_pos - bds_len,
                            'genomic_end': current_pos + bds_len if strand == 1 else end_pos,
                            'length': bds_len,
                            'g_content': thermal_result['g_content'],
                            'tm': thermal_result['tm'],
                            'tm_3': thermal_result['tm_3prime'],
                            'tm_5': thermal_result['tm_5prime'],
                            'tm_diff': thermal_result['tm_diff'],
                            'free_energy': thermal_result['free_energy'],
                            'strategy': 'isoform_specific',
                            'target_isoforms': target_isoforms,
                            'isoform_overlap_num': len(target_isoforms),
                            'isoform_id': isoform.get('id'),
                            'isoform_name': isoform.get('external_name', isoform.get('id')),
                            'overlapping_regions': regions
                        })
                current_pos += 1
        
        # global downselect if necessary
        if len(all_probes) > self.search_config.max_binding_sites:
            all_probes.sort(key=lambda x: x['genomic_start'])
            positions = [p['genomic_start'] for p in all_probes]
            selected = self._optimize_subsequence(positions, self.search_config.max_binding_sites, 1, 1, gene_name)
            all_probes = [p for p in all_probes if p['genomic_start'] in selected]
        return all_probes


class IsoformConsensusStrategy(SearchStrategy):
    """Design probes that bind as many isoforms as possible (consensus sites)."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, isoforms: List[Dict], genome_accessor=None):
        super().__init__(search_config, filter_config)
        self.isoforms = isoforms
        self.genome_accessor = genome_accessor
        self.awareness = IsoformAwareness(self.isoforms, self.genome_accessor) if (self.isoforms and self.genome_accessor) else None
    
    def search_binding_sites(self, sequence: str, gene_name: str) -> List[Dict[str, Any]]:
        if not self.awareness:
            return BruteForceStrategy(self.search_config, self.filter_config).search_binding_sites(sequence, gene_name)
        
        genome_seq = self.awareness.get_genome_seq()
        bds_len = self.search_config.binding_site_length
        candidates: List[Dict[str, Any]] = []
        seen_positions = set()  # Track seen genomic positions to avoid duplicates
        
        for isoform in self.isoforms:
            iso_start = isoform['start']
            iso_end = isoform['end']
            strand = isoform.get('strand', 1)
            current_pos = iso_start
            while current_pos + bds_len <= iso_end:
                end_pos = self.awareness.find_end_position(current_pos, bds_len)
                if end_pos > iso_end:
                    break
                regions = self.awareness.get_overlapping_regions(current_pos, end_pos)
                covered = set(self.awareness.isoforms_union(regions))
                if not covered:
                    current_pos += 1; continue
                ls = current_pos - self.awareness.global_start
                le = end_pos - self.awareness.global_start
                target_seq = genome_seq[ls:le]
                if strand == -1:
                    target_seq = self._reverse_complement(target_seq)
                
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
                    current_pos += 1; continue
                
                # Calculate genomic position
                genomic_start = current_pos if strand == 1 else end_pos - bds_len
                genomic_end = current_pos + bds_len if strand == 1 else end_pos
                
                # Skip if we've already seen this genomic position
                if (genomic_start, genomic_end) in seen_positions:
                    current_pos += 1; continue
                
                seen_positions.add((genomic_start, genomic_end))
                
                candidates.append({
                    'gene_name': gene_name,
                    'sequence': probe_seq,  # Store the probe sequence (reverse complement)
                    'target_sequence': target_seq,  # Store original target sequence
                    'arm_3prime': thermal_result['arm_3prime'],
                    'arm_5prime': thermal_result['arm_5prime'],
                    'position': current_pos - iso_start,
                    'genomic_start': genomic_start,
                    'genomic_end': genomic_end,
                    'length': bds_len,
                    'g_content': thermal_result['g_content'],
                    'tm': thermal_result['tm'],
                    'tm_3': thermal_result['tm_3prime'],
                    'tm_5': thermal_result['tm_5prime'],
                    'tm_diff': thermal_result['tm_diff'],
                    'free_energy': thermal_result['free_energy'],
                    'strategy': 'isoform_consensus',
                    'target_isoforms': sorted(list(covered)),
                    'isoform_overlap_num': len(covered),
                    'isoform_id': isoform.get('id'),
                    'isoform_name': isoform.get('external_name', isoform.get('id')),
                    'overlapping_regions': regions
                })
                current_pos += 1
        
        # Prefer higher coverage then spacing optimization
        candidates.sort(key=lambda x: (-x['isoform_overlap_num'], x['genomic_start']))
        if len(candidates) > self.search_config.max_binding_sites:
            positions = [c['genomic_start'] for c in candidates]
            selected = self._optimize_subsequence(positions, self.search_config.max_binding_sites, 1, 1, gene_name)
            candidates = [c for c in candidates if c['genomic_start'] in selected]
        return candidates


class BindingSiteSearcher:
    """High-level facade to run a strategy across many genes."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig, genome_accessor=None):
        self.search_config = search_config
        self.filter_config = filter_config
        self.genome_accessor = genome_accessor
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
        
        # Add genome_accessor for isoform_specific strategy
        if strategy_name == 'isoform_specific' and self.genome_accessor:
            kwargs['genome_accessor'] = self.genome_accessor
        
        # Forward genomic_context for brute_force if provided via kwargs
        return strategy_class(self.search_config, self.filter_config, **kwargs)
    
    def search_all_genes(self, sequences: Dict[str, Any], isoforms: Optional[Dict] = None) -> Dict[str, List[Dict[str, Any]]]:
        results = {}
        
        print(f"Searching binding sites for {len(sequences)} genes...")
        for gene_name, gene_data in sequences.items():
            if 'sequence' not in gene_data:
                continue
            
            print(f"Processing gene: {gene_name}")
            
            # Get sequence directly
            sequence = gene_data['sequence']
            
            # Create appropriate strategy
            if self.search_config.search_strategy == 'isoform_specific':
                if isoforms and gene_name in isoforms:
                    strategy = self.create_strategy('isoform_specific', isoforms=isoforms[gene_name])
                else:
                    strategy = self.create_strategy('brute_force')
            elif self.search_config.search_strategy == 'isoform_consensus':
                if isoforms and gene_name in isoforms:
                    strategy = self.create_strategy('isoform_consensus', isoforms=isoforms[gene_name], genome_accessor=self.genome_accessor)
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
            print(f"Found {len(binding_sites)} binding sites for {gene_name}")
        
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


# --- Shared isoform awareness helpers ---

def build_exon_splice_counts(isoforms: List[Dict]) -> Dict[Tuple[int, int], List[str]]:
    """Build mapping from exon splice (start,end) to list of isoform names that contain it."""
    all_exon_boundaries = set()
    for iso in isoforms:
        for exon in iso.get('exons', []):
            all_exon_boundaries.add(exon['start'])
            all_exon_boundaries.add(exon['end'])
    split_points = sorted(all_exon_boundaries)

    exon_splice_counts: Dict[Tuple[int, int], List[str]] = {}
    for iso in isoforms:
        iso_name = iso.get('external_name', iso.get('id'))
        exons = sorted(iso.get('exons', []), key=lambda x: x['start'])
        for exon in exons:
            exon_start, exon_end = exon['start'], exon['end']
            cur_splits = [exon_start] + [p for p in split_points if exon_start < p < exon_end] + [exon_end]
            for j in range(len(cur_splits) - 1):
                exon_splice = (cur_splits[j], cur_splits[j + 1])
                exon_splice_counts.setdefault(exon_splice, []).append(iso_name)
    return exon_splice_counts


class IsoformAwareness:
    """Reusable isoform-aware preprocessing over exon splices and genome locus."""
    
    def __init__(self, isoforms: List[Dict], genome_accessor):
        self.isoforms = isoforms
        self.genome_accessor = genome_accessor
        # locus
        self.global_start = min(iso['start'] for iso in self.isoforms)
        self.global_end = max(iso['end'] for iso in self.isoforms)
        sample_isoform = self.isoforms[0]
        self.seq_region_name = sample_isoform.get('seq_region_name', '1')
        # splices
        self.exon_splice_counts = build_exon_splice_counts(self.isoforms)
        self._genome_seq_cache: Optional[str] = None
    
    def get_genome_seq(self) -> str:
        if self._genome_seq_cache is None:
            self._genome_seq_cache = self.genome_accessor(self.seq_region_name, self.global_start, self.global_end)
        return self._genome_seq_cache
    
    def find_end_position(self, start_pos: int, length: int) -> int:
        exon_splices = list(self.exon_splice_counts.keys())
        current_pos = start_pos
        remaining_length = length
        for splice_start, splice_end in exon_splices:
            if splice_start <= current_pos < splice_end:
                avail = splice_end - current_pos
                if remaining_length <= avail:
                    return current_pos + remaining_length
                remaining_length -= avail
                current_pos = splice_end
        return start_pos + length
    
    def get_overlapping_regions(self, start_pos: int, end_pos: int) -> List[Tuple[int, int]]:
        return [(s, e) for (s, e) in self.exon_splice_counts.keys() if s < end_pos and e > start_pos]
    
    def isoforms_intersection(self, regions: List[Tuple[int, int]]) -> List[str]:
        if not regions:
            return []
        inter = set(self.exon_splice_counts.get(regions[0], []))
        for r in regions[1:]:
            inter &= set(self.exon_splice_counts.get(r, []))
        return list(inter)
    
    def isoforms_union(self, regions: List[Tuple[int, int]]) -> List[str]:
        covered = set()
        for r in regions:
            covered.update(self.exon_splice_counts.get(r, []))
        return list(covered)
