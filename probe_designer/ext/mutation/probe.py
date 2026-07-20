"""
Mutation probe design module
Designs padlock probes (iLock probes) for mutation detection based on WES data.

Key design principles:
1. Ref/Alt are relative to reference genome sequence
2. For Ref padlock: 5' and 3' positions are genomic absolute positions
3. For mutation padlock: 5' and 3' positions reflect actual positions in mutated sequence
4. Strand handling:
   - strand=1: mRNA = reference genome sequence = target seq
   - strand=-1: mRNA = reverse complement of reference genome = target seq
5. Padlock arms are reverse complementary to target sequence
6. Binding sequence = plp arm 3' + plp arm 5' (reverse complementary to target seq)
"""

import re
from typing import List, Dict, Any, Optional, Tuple
import logging
from tqdm import tqdm
from copy import deepcopy

from probe_designer.chemistry import reverse_complement
from probe_designer.config import SearchConfig, FilterConfig, BlastConfig
from probe_designer.filtering import SequenceFilter


class MutationProbeDesigner:
    """Designs padlock probes for mutation detection."""
    
    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig,
                 target_type: str = "RNA"):
        if target_type not in ("RNA", "DNA"):
            raise ValueError(f"target_type must be 'RNA' or 'DNA', got {repr(target_type)}")
        self.search_config = search_config
        self.filter_config = filter_config
        self.target_type = target_type  # "RNA" for mRNA, "DNA" for cDNA
        self.filter = SequenceFilter(filter_config, BlastConfig())
    
    def _global_to_local_pos(self, global_pos: int, global_start: int) -> int:
        """Convert global genomic position to local sequence position (strand=1)."""
        return global_pos - global_start
    
    def _local_to_global_pos(self, local_pos: int, global_start: int, strand: int) -> int:
        """Convert local sequence position to global genomic position."""
        if strand == 1:
            return local_pos + global_start
        else:  # strand == -1
            # For negative strand, we need to reverse the position
            # This is only for recording purposes, not for sequence extraction
            return global_start + local_pos
    
    def _get_padlock_local_positions(self, mut_start: int, mut_end: int, mut_type: str, 
                                   strand: int, is_mutation_probe: bool = False) -> Tuple[int, int]:
        """
        Calculate padlock local positions based on mutation type and strand.
        
        Key principle: pos5 should always be at the mutation position for proper detection.
        
        Args:
            mut_start: Local mutation start position (strand=1 coordinates)
            mut_end: Local mutation end position (strand=1 coordinates)
            mut_type: Type of mutation (SNP, DEL, INS, MNP)
            strand: Strand direction (1 or -1)
            is_mutation_probe: Whether this is for mutation probe (True) or reference probe (False)
            
        Returns:
            pos5_local in local coordinates (strand=1)
        """
        # For SNP, pos5 should be at the SNP position
        if mut_type == "SNP": pos5_local = mut_start
                
        # For MNP, pos5 should be at the mutation site
        elif mut_type == "MNP":
            if strand == 1:
                pos5_local = int(-(-(mut_start+mut_end)//2))
            elif strand == -1:
                pos5_local = int((mut_start+mut_end)//2)
 
        elif mut_type == "DEL":
            if is_mutation_probe:
                # For deletion mutation probes, position at deletion site
                # here the mut_start and mut_end are in the mutated sequence
                if strand == 1:
                    pos5_local = mut_start
                elif strand == -1:
                    pos5_local = mut_end + 1
            else:
                # For reference probes, span the entire deletion region
                if strand == 1:
                    pos5_local = mut_start-1
                elif strand == -1:
                    pos5_local = mut_end+1
               
        elif mut_type == "INS":
            # For insertion, pos5 should be at the insertion site
            # here the mut_start and mut_end are in the mutated sequence
            if is_mutation_probe:
                if strand == 1:
                    pos5_local = mut_start
                elif strand == -1:
                    pos5_local = mut_end
            else: 
                if strand == 1:
                    pos5_local = mut_start
                elif strand == -1:
                    pos5_local = mut_end + 1
        else:
            raise ValueError(f"Unknown mutation type: {mut_type}")
            
        return pos5_local
    
    
    def _mutate_genomic_sequence(self, genomic_seq: str, mut_start: int, mut_end: int, 
                                ref: str, alt: str, strand: int) -> Tuple[str, int, int]:
        """
        Apply mutation to genomic sequence (always strand=1 direction).
        
        Args:
            genomic_seq: Genomic sequence in strand=1 direction
            mut_start: Local mutation start position
            mut_end: Local mutation end position
            ref: Reference allele
            alt: Alternative allele
            strand: Strand direction (1 or -1)
            
        Returns:
            Tuple of (mutated_sequence, new_mut_start, new_mut_end) where positions
            are in the mutated sequence
        """
        # Ref/Alt are always relative to the genomic sequence (strand=1 direction)
        # We don't need to convert them here, the conversion happens in _extract_padlock_arms
        ref_conv = ref
        alt_conv = alt
        
        if ref_conv == "-":  # Insertion
            mutated = genomic_seq[:mut_end+1] + alt_conv + genomic_seq[mut_end+1:]
            # In mutated sequence, insertion region is from mut_end to mut_end+len(alt)
            new_mut_start = mut_end + 1
            new_mut_end = mut_end + len(alt_conv)
        elif alt_conv == "-":  # Deletion
            mutated = genomic_seq[:mut_start] + genomic_seq[mut_end+1:]
            # In mutated sequence, deletion region is collapsed to mut_start position
            new_mut_start = mut_start - 1
            new_mut_end = mut_start - 1
        else:  # Substitution (SNP or MNP)
            mutated = genomic_seq[:mut_start] + alt_conv + genomic_seq[mut_end+1:]
            # In mutated sequence, substitution region is from mut_start to mut_start+len(alt)
            new_mut_start = mut_start
            new_mut_end = mut_start + len(alt_conv) - 1
        
        return mutated, new_mut_start, new_mut_end
    
    def _get_mutation_type(self, ref: str, alt: str, start: int, end: int) -> str:
        """Determine mutation type based on reference and alternative alleles."""
        if alt == "-":
            return "DEL"
        elif ref == "-":
            return "INS"
        elif len(ref) == 1 and len(alt) == 1:
            return "SNP"
        else:
            return "MNP"
    
    def _extract_padlock_arms(self, plus_seq: str, pos5_local: int, arm5_len: int, 
                             arm3_len: int, strand: int) -> Dict[str, Any]:
        """
        Extract padlock arms from genomic sequence based on strand direction.
        
        Key principles:
        1. Input plus_seq is always in strand=1 direction
        2. pos5_local is in strand=1 coordinates
        3. Function performs sequence conversion and coordinate transformation based on strand
        4. pos3 is automatically generated based on converted pos5
        5. For strand=1: mRNA = plus_seq, padlock arms = reverse complement of mRNA
        6. For strand=-1: mRNA = reverse complement of plus_seq, padlock arms = reverse complement of mRNA
        7. target_sequence = arm5_seq + arm3_seq (continuous in mRNA)
        
        Args:
            plus_seq: Genomic sequence in strand=1 direction
            pos5_local: Local position for 5' arm extraction (strand=1 coordinates)
            arm5_len: Length of 5' arm
            arm3_len: Length of 3' arm
            strand: Strand direction (1 or -1)
            
        Returns:
            Dictionary containing arm sequences and padlock sequences
        """
        # Step 1: Convert sequence based on strand
        if strand == 1:
            # For positive strand, use plus_seq directly
            target_seq = plus_seq
            # pos5 and pos3 remain the same in local coordinates
            pos5_target = pos5_local
            pos3_target = pos5_local + 1 # pos3 equals pos5 for Python coordinates (left-closed, right-open)
        else:  # strand == -1
            # For negative strand, convert plus_seq to mRNA (reverse complement)
            target_seq = reverse_complement(plus_seq)
            # Convert coordinates: for negative strand, we need to reverse the positions
            # pos5_local in plus_seq corresponds to len(plus_seq) - 1 - pos5_local in target_seq
            pos5_target = len(plus_seq) - 1 - pos5_local
            pos3_target = pos5_target + 1  # pos3 equals pos5 for Python coordinates
        
        # Step 2: Extract arm sequences from target sequence
        arm5_seq = target_seq[pos5_target+1-arm5_len:pos5_target+1]  # Left arm from target
        arm3_seq = target_seq[pos3_target:pos3_target+arm3_len]  # Right arm from target
        
        # Step 3: Generate padlock arms based on target type
        if self.target_type == "RNA":
            # mRNA targeting: padlock arms = RC(mRNA region)
            probe_arm5 = reverse_complement(arm5_seq)
            probe_arm3 = reverse_complement(arm3_seq)
        else:
            # cDNA targeting: cDNA = RC(mRNA), so padlock arms = mRNA sequence
            # Arms are swapped because cDNA runs antiparallel to mRNA
            probe_arm5 = arm3_seq
            probe_arm3 = arm5_seq
        
        # Step 4: Create binding sequence (plp arm 3' + plp arm 5')
        binding_seq = probe_arm3 + probe_arm5
        
        # Step 5: Create target_sequence (the actual target strand the probe binds)
        if self.target_type == "RNA":
            target_sequence = arm5_seq + arm3_seq  # mRNA sequence
        else:
            target_sequence = reverse_complement(arm5_seq + arm3_seq)  # cDNA sequence

        return {
            'probe_arm5': probe_arm5,  # Padlock 5' arm
            'probe_arm3': probe_arm3,  # Padlock 3' arm
            'target_sequence': target_sequence,  # Actual target sequence (mRNA or cDNA)
            'pos5_target': pos5_target,  # pos5 in target sequence coordinates
            'pos3_target': pos3_target   # pos3 in target sequence coordinates
        }
    
    def _evaluate_probe_arms(self, arm5: str, arm3: str, target_seq: str) -> Dict[str, Any]:
        """Evaluate probe arms using thermal filter."""
        thermal_result = self.filter.thermal_filter(
            arm3,  # 3' arm
            arm5,  # 5' arm
            sequence_type="DNA",
            target_type=self.target_type,
            target_sequence=target_seq
        )
        
        return thermal_result
    
    def _generate_probe_candidates(self, plus_seq: str, pos5_local: int, strand: int, 
                                 plp_params: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate probe candidates for a given position."""
        candidates = []
        
        # Get arm length parameters
        ini_arm_len = plp_params.get('ini_arm_len', 20)
        min_arm_len = plp_params.get('min_arm_len', 15)
        max_arm_len = plp_params.get('max_arm_len', 25)
        
        # Generate all possible arm length combinations
        for arm5_len in range(min_arm_len, max_arm_len + 1):
            for arm3_len in range(min_arm_len, max_arm_len + 1):
                # Check if we have enough sequence
                # For Python coordinates (left-closed, right-open), pos3 equals pos5
                seq_len = len(plus_seq)
                if strand == 1 and (pos5_local - arm5_len < 0 or pos5_local + arm3_len > seq_len): 
                    continue
                if strand == -1 and (pos5_local - arm3_len < 0 or pos5_local + arm5_len > seq_len): 
                    continue
                
                # Extract padlock arms
                arm_info = self._extract_padlock_arms(plus_seq, pos5_local, arm5_len, arm3_len, strand)
                
                # Evaluate using thermal filter
                thermal_result = self._evaluate_probe_arms(
                    arm_info['probe_arm3'],  # 3' arm (comes first in thermal filter)
                    arm_info['probe_arm5'],  # 5' arm (comes second in thermal filter)
                    arm_info['target_sequence']
                )
                
                # Keep all candidates for diagnostics, filter later
                # Debug: record failure reasons for failed probes
                if not thermal_result['passed']:
                    if len(candidates) < 5:
                        logging.debug("Failed thermal filter: %s", thermal_result.get('failed_checks', []))
                
                # Calculate score (lower is better)
                tm_diff = abs(thermal_result['tm_5prime'] - thermal_result['tm_3prime'])
                arm_len_diff = abs(arm5_len - ini_arm_len) + abs(arm3_len - ini_arm_len)
                mfe_penalty = -thermal_result['free_energy']  # Convert to positive penalty
                score = tm_diff + arm_len_diff + mfe_penalty
                
                candidate = {
                    'arm5_len': arm5_len,
                    'arm3_len': arm3_len,
                    # 'arm5_seq': arm_info['arm5_seq'],
                    # 'arm3_seq': arm_info['arm3_seq'],
                    'probe_arm5': arm_info['probe_arm5'],
                    'probe_arm3': arm_info['probe_arm3'],
                    'target_sequence': arm_info['target_sequence'],
                    'score': score,
                    'tm': thermal_result['tm'],
                    'tm_5prime': thermal_result['tm_5prime'],
                    'tm_3prime': thermal_result['tm_3prime'],
                    'tm_diff': thermal_result['tm_diff'],
                    'g_content': thermal_result['g_content'],
                    'gc_content': thermal_result.get('gc_content'),  # Schema-v2
                    'free_energy': thermal_result['free_energy'],
                    'passed_thermal_filter': thermal_result['passed'],
                    'failed_thermal_checks': thermal_result.get('failed_checks', [])
                }
                
                candidates.append(candidate)
        
        # Sort by score (lower is better)
        candidates.sort(key=lambda x: x['score'])
        
        return candidates
    
    def design_mutation_probes(self, mut_info_list: List[Dict[str, Any]], 
                             genome_accessor, plp_params: Dict[str, Any], 
                             pad: int = 100) -> List[Dict[str, Any]]:
        """
        Design padlock probes for mutation detection.
        
        Args:
            mut_info_list: List of mutation information dictionaries
            genome_accessor: Function to access genome sequence
            plp_params: Probe design parameters
            pad: Padding around mutation site
            
        Returns:
            List of mutation information with designed probes
        """
        results = deepcopy(mut_info_list)
        
        for mut_info in tqdm(results, desc="Designing mutation probes"):
            # Get mutation coordinates
            chrom = mut_info['Chr'].replace('chr', '')
            start = mut_info['Start']
            end = mut_info['End']
            ref = mut_info['Ref']
            alt = mut_info['Alt']
            strand = mut_info.get('strand', 1)
            
            # Determine mutation type
            mut_type = self._get_mutation_type(ref, alt, start, end)
            mut_info['mut_type'] = mut_type
            
            # Get genomic sequence around mutation
            # End coordinate is inclusive in this dataset
            global_start = max(0, start - pad)
            global_end = end + pad
            
            try:
                local_seq = genome_accessor(chrom, global_start, global_end)
                if not local_seq:
                    logging.debug("Could not fetch sequence for %s", mut_info.get('Gene.refGene.x', 'unknown'))
                    mut_info['plp_ref'] = []
                    mut_info['plp_mut'] = []
                    continue
            except Exception as e:
                logging.debug("Error fetching sequence for %s: %s", mut_info.get('Gene.refGene.x', 'unknown'), e)
                mut_info['plp_ref'] = []
                mut_info['plp_mut'] = []
                continue
            
            # Convert to local coordinates
            # End coordinate is inclusive, so we need to add 1 for proper slicing
            local_mut_start = start - global_start
            local_mut_end = end - global_start
            
            # Design probes for reference sequence (local coordinates)
            plp_ref_local = self._design_reference_probes(
                local_seq, local_mut_start, local_mut_end, mut_type, strand, plp_params
            )
            
            # Convert reference probe local positions to genomic coordinates
            for probe in plp_ref_local:
                if 'plp_pos5_local' in probe and 'plp_pos3_local' in probe:
                    # Convert local positions to global positions (considering strand for recording)
                    probe['plp_pos5'] = self._local_to_global_pos(probe['plp_pos5_local'], global_start, strand)
                    probe['plp_pos3'] = self._local_to_global_pos(probe['plp_pos3_local'], global_start, strand)
                probe['coord_system'] = 'genomic'
            mut_info['plp_ref'] = plp_ref_local
            
            # Design probes for mutated sequence (local coordinates)
            plp_mut_local = self._design_mutated_probes(
                local_seq, local_mut_start, local_mut_end, ref, alt, mut_type, strand, plp_params
            )
            
            # Convert mutation probe local positions to genomic coordinates
            for probe in plp_mut_local:
                if 'plp_pos5_local' in probe and 'plp_pos3_local' in probe:
                    # For mutation probes, we need to account for the mutation effect on coordinates
                    # The local positions are in the mutated sequence, but we need to map them
                    # back to the original genomic coordinates
                    
                    # Apply mutation to get the mutated sequence and new positions
                    mut_seq, mut_new_start, mut_new_end = self._mutate_genomic_sequence(
                        local_seq, local_mut_start, local_mut_end, ref, alt, strand
                    )
                    
                    # Map positions in mutated sequence back to genomic coordinates
                    # The local positions are already correctly set in _design_mutated_probes
                    # We just need to convert them to global coordinates (considering strand for recording)
                    probe['plp_pos5'] = self._local_to_global_pos(probe['plp_pos5_local'], global_start, strand)
                    probe['plp_pos3'] = self._local_to_global_pos(probe['plp_pos3_local'], global_start, strand)
                
                probe['coord_system'] = 'genomic'
            mut_info['plp_mut'] = plp_mut_local
        
        return results
    
    def _design_reference_probes(self, genomic_seq: str, local_mut_start: int, 
                               local_mut_end: int, mut_type: str, strand: int, 
                               plp_params: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Design probes for reference sequence.
        
        Args:
            genomic_seq: Genomic sequence in strand=1 direction
            local_mut_start: Local mutation start position
            local_mut_end: Local mutation end position
            mut_type: Type of mutation (SNP, DEL, INS, MNP)
            strand: Strand direction (1 or -1)
            plp_params: Padlock parameters
            
        Returns:
            List of reference probe candidates
        """
        probes = []
        
        # Get padlock local positions for reference probes
        pos5_local = self._get_padlock_local_positions(
            local_mut_start, local_mut_end, mut_type, strand, is_mutation_probe=False
        )
        
        # Generate probe candidates
        candidates = self._generate_probe_candidates(genomic_seq, pos5_local, strand, plp_params)
        
        # Add position information to candidates
        for candidate in candidates:
            candidate['plp_pos5_local'] = pos5_local
            # candidate['plp_pos3_local'] = pos3_local
            candidate['probe_type'] = 'reference'
        
        probes.extend(candidates)
        
        # Sort by score and return
        probes.sort(key=lambda x: x['score'])
        return probes
    
    def _design_mutated_probes(self, genomic_seq: str, local_mut_start: int, 
                             local_mut_end: int, ref: str, alt: str, mut_type: str,
                             strand: int, plp_params: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Design probes for mutated sequence.
        
        Args:
            genomic_seq: Genomic sequence in strand=1 direction
            local_mut_start: Local mutation start position
            local_mut_end: Local mutation end position
            ref: Reference allele
            alt: Alternative allele
            mut_type: Type of mutation (SNP, DEL, INS, MNP)
            strand: Strand direction (1 or -1)
            plp_params: Padlock parameters
            
        Returns:
            List of mutation probe candidates
        """
        probes = []
        
        # Apply mutation to genomic sequence
        mut_genomic_seq, mut_new_start, mut_new_end = self._mutate_genomic_sequence(
            genomic_seq, local_mut_start, local_mut_end, ref, alt, strand
        )
        
        # Get padlock local positions for mutation probes
        pos5_local= self._get_padlock_local_positions(
            mut_new_start, mut_new_end, mut_type, strand, is_mutation_probe=True
        )
        
        # Generate probe candidates from mutated sequence
        candidates = self._generate_probe_candidates(mut_genomic_seq, pos5_local, strand, plp_params)
        
        # Add position information to candidates
        for candidate in candidates:
            candidate['plp_pos5_local'] = pos5_local
            # candidate['plp_pos3_local'] = pos3_local
            candidate['mutation_type'] = mut_type.lower()
            candidate['probe_type'] = 'mutation'
        
        probes.extend(candidates)
        
        # Sort by score and return
        probes.sort(key=lambda x: x['score'])
        return probes


class MutationProbeDesignerInvader(MutationProbeDesigner):
    """iLock invader-overlap designer (canonical Krzywkowski 2017 chemistry).

    arm5_binding and arm3_binding overlap by ONE nucleotide on the target — the
    SNP/discrimination position. As a result probe_arm5[0] == probe_arm3[-1] ==
    RC(SNP nt), i.e. the same SNP-complement carried at both the probe's 5' body
    start and its 3' tip. With the standard iLock stitch
    (linker = probe_arm3[-1]), the duplicated invader linker carries the
    SNP-complement and competes with arm3's 3' tip for the same target nt;
    FEN1 cleavage is then strictly conditional on SNP base pairing.

    The parent class places arm3 starting at SNP+1 (no overlap, gap=0), so the
    linker does not in general equal the SNP-complement and invader chemistry
    only "works" by coincidence (when SNP nt happens to equal the next nt).
    Use this subclass for the iLock chemistry only — the cDNA/DNA pass keeps
    the parent's arm-swap behavior and must NOT use this subclass.
    """

    def __init__(self, search_config: SearchConfig, filter_config: FilterConfig,
                 target_type: str = "RNA"):
        if target_type != "RNA":
            raise ValueError(
                f"MutationProbeDesignerInvader requires target_type='RNA', got "
                f"{target_type!r}. cDNA chemistry uses the parent's arm-swap logic."
            )
        super().__init__(search_config, filter_config, target_type="RNA")

    def _extract_padlock_arms(self, plus_seq: str, pos5_local: int, arm5_len: int,
                              arm3_len: int, strand: int) -> Dict[str, Any]:
        """Same as parent except pos3_target == pos5_target (overlap by 1 at SNP).

        target_sequence is the continuous mRNA region with the duplicate SNP nt
        deduplicated, so its length is `arm5_len + arm3_len - 1`.
        """
        if strand == 1:
            target_seq = plus_seq
            pos5_target = pos5_local
        else:
            target_seq = reverse_complement(plus_seq)
            pos5_target = len(plus_seq) - 1 - pos5_local
        # KEY DIFFERENCE vs parent: arm3 starts AT pos5 (overlap by 1 nt at SNP)
        # rather than at pos5+1. Both arms include the SNP at their adjacent end.
        pos3_target = pos5_target

        arm5_seq = target_seq[pos5_target + 1 - arm5_len : pos5_target + 1]
        arm3_seq = target_seq[pos3_target : pos3_target + arm3_len]

        probe_arm5 = reverse_complement(arm5_seq)
        probe_arm3 = reverse_complement(arm3_seq)

        # Continuous mRNA region, deduplicated (the SNP nt at arm5_seq[-1] ==
        # arm3_seq[0] is shared on the target).
        target_sequence = arm5_seq + arm3_seq[1:]

        return {
            'probe_arm5': probe_arm5,
            'probe_arm3': probe_arm3,
            'target_sequence': target_sequence,
            'pos5_target': pos5_target,
            'pos3_target': pos3_target,
        }


class MutationProbeDesigner3End(MutationProbeDesigner):
    """Places the mutation at the probe's 3' tip for non-iLock direct-RNA padlocks.

    The default `MutationProbeDesigner` sets `pos5_local = mut_start` (SNP), so
    `arm5_seq` ends AT the mutation (mutation ends up at the probe's 5' body
    start after RC). For non-iLock padlocks we want the mutation at the probe's
    3' tip — i.e., at `arm3_seq`'s FIRST base. This requires `arm3_seq` to start
    at the mutation in `target_seq` (= mRNA), which means shifting `pos5_local`
    by one in the strand-appropriate direction:

      strand +1: target_seq = plus_seq; shift `pos5_local` by -1.
      strand -1: target_seq = RC(plus_seq); shift `pos5_local` by +1.

    Only applied when `target_type == "RNA"`. The DNA (cDNA) pass keeps the
    default pos5, because its arm-swap convention already places the mutation
    at the probe's 3' tip.
    """

    def _get_padlock_local_positions(self, mut_start, mut_end, mut_type, strand,
                                      is_mutation_probe=False):
        pos5 = super()._get_padlock_local_positions(
            mut_start, mut_end, mut_type, strand, is_mutation_probe
        )
        if self.target_type != "RNA":
            return pos5
        return pos5 - 1 if strand == 1 else pos5 + 1


def verify_iLock_probe(probe_seq: str, plp_arm5: str, plp_arm3: str,
                       flap: str, *, check_invader: bool = True) -> None:
    """Hard runtime guard for iLock probe stitching.

    Raises ValueError if the assembled probe does not match the canonical
    iLock invader layout `flap + plp_arm5 + backbone + plp_arm3`.

    The two failure modes the lab has actually shipped (and that this guard
    must catch before any future panel goes to a synthesis vendor):

    1. **Swapped stitch geometry** (BZ09 / BZ09_sup / BZ10 / BZ12 era):
       the probe was assembled as `flap + plp_arm3 + bb + plp_arm5`,
       i.e. plp_arm3 at the 5' body and plp_arm5 at the 3' end. Because
       mutation_probe.py defines plp_arm5 = "arm binding mRNA 5' side" (whose
       5' tip lands at the ligation junction on target), this swap places
       both probe tips at the FAR ENDS of the binding region — the probe
       binds but cannot ligate. ALL such probes are non-functional.

    2. **Invader invariant violation**: plp_arm5[0] != plp_arm3[-1]. With
       invader-overlap chemistry the same SNP-complement base sits at both
       arm tips so they compete for the target SNP nucleotide; if they
       differ, FEN1 cleavage is no longer SNP-discriminating. To fix at
       the design level, step1 must use `MutationProbeDesignerInvader`.

    Layout note (2026-06-08 update): the canonical iLock no longer has a
    separate 1-nt linker between flap and arm5. The earlier convention
    duplicated arm3[-1] right after the flap, which left an extra
    SNP-complement base at the junction (probe ended up 1 nt longer with
    XX at flap–arm5 boundary). The invader designer's overlap already
    puts the SNP-complement at arm5[0], so the linker insertion was
    redundant and pushed the FEN1 cleavage site out of register.

    Parameters
    ----------
    probe_seq : str
        The fully-assembled iLock probe sequence (5'->3', case-insensitive).
    plp_arm5 : str
        The arm binding mRNA's 5' side (= mutation_probe.py's plp_arm5 column).
        MUST end up at the probe's 5' body, immediately after the flap.
    plp_arm3 : str
        The arm binding mRNA's 3' side (= plp_arm3 column). MUST end up at
        the probe's 3' end, after the backbone.
    flap : str
        The iLock invader flap. Must match the probe's 5' prefix.
        Default-ish values: 'CGTTGCTGTGGCG' (HuangLab 2026+) or
        'TATATCCCTATAT' (BZ09 legacy).
    check_invader : bool, optional
        If True (default), additionally enforce the invader-overlap invariant
        plp_arm5[0] == plp_arm3[-1]. Set False only if you know the panel
        is intentionally a flap-padlock without invader chemistry.
    """
    p = probe_seq.upper()
    f = flap.upper()
    a5 = plp_arm5.upper(); a3 = plp_arm3.upper()

    if not p.startswith(f):
        raise ValueError(
            f"iLock probe does not start with the expected flap {flap!r}. "
            f"Got prefix {probe_seq[:len(flap)]!r}. Did you pass the right flap?"
        )
    rest = p[len(f):]
    if len(rest) < len(a5):
        raise ValueError("iLock probe too short — no arm5 body after flap.")

    # Geometry guards — the key check.
    body_arm = rest[:len(a5)]
    if body_arm != a5:
        # Detect specifically the BZ09-style swap so the error message points
        # the user at the actual bug.
        if rest.startswith(a3) and rest.endswith(a5):
            raise ValueError(
                "iLock STITCH GEOMETRY SWAPPED (BZ09-style bug): probe places "
                "plp_arm3 at the 5' body and plp_arm5 at the 3' end. "
                "mutation_probe.py defines plp_arm5 = 'arm binding mRNA 5'-side' "
                "and that arm MUST sit at the probe's 5' body so its 5' tip "
                "lands at the ligation junction on the target. Swap them: use "
                "`flap + plp_arm5 + backbone + plp_arm3`. With this "
                "stitch the resulting probe binds the target but cannot "
                "circularize/ligate."
            )
        raise ValueError(
            f"iLock probe 5'-body arm does not match plp_arm5. Expected "
            f"{plp_arm5!r}, got {body_arm!r}. The probe layout must be "
            f"`flap + plp_arm5 + backbone + plp_arm3`."
        )
    if not p.endswith(a3):
        raise ValueError(
            f"iLock probe 3'-end arm does not match plp_arm3. The probe must "
            f"end with plp_arm3={plp_arm3!r}; got tail {probe_seq[-len(a3):]!r}."
        )

    if check_invader:
        if a5[:1] != a3[-1:]:
            raise ValueError(
                f"iLock invader invariant violated: plp_arm5[0]={plp_arm5[0]!r} "
                f"!= plp_arm3[-1]={plp_arm3[-1]!r}. Both must be RC(SNP_alt) so "
                f"the SNP-complement base sits at both arm tips and competes "
                f"for the target SNP nucleotide. Step1 must use "
                f"`MutationProbeDesignerInvader` to make arm5 and arm3 overlap "
                f"by 1 nt at the SNP. (Pass check_invader=False only for "
                f"intentional flap-padlock panels with no invader chemistry.)"
            )


def assemble_ilock(arm5: str, arm3: str, backbone: str, flap: str) -> str:
    """Assemble a canonical iLock probe.

    Layout: ``flap + arm5.lower() + backbone.upper() + arm3.lower()``.

    No separate 1-nt linker between flap and arm5: the invader-overlap
    designer (:class:`MutationProbeDesignerInvader`) already shifts arm5 so
    that ``arm5[0] == arm3[-1] == RC(SNP)``. ``arm5[0]`` itself sits right
    after the flap and is the SNP-complement base that competes with
    ``arm3[-1]`` for the target SNP nucleotide — FEN1 cleaves the flap +
    arm5[0] block off when the pairing is correct. Inserting an additional
    linker would double the SNP-complement at the flap–arm5 boundary and
    push the FEN1 cleavage site out of register.

    Calls :func:`verify_iLock_probe` on the result; raises ``ValueError`` if
    the assembled probe violates canonical layout or the invader invariant.
    """
    if not arm3:
        raise ValueError("assemble_ilock: arm3 must not be empty")
    arm5_l = arm5.lower()
    arm3_l = arm3.lower()
    probe = flap + arm5_l + backbone.upper() + arm3_l
    verify_iLock_probe(probe, plp_arm5=arm5_l, plp_arm3=arm3_l, flap=flap)
    return probe


def design_mutation_probes(mut_info_list: List[Dict[str, Any]],
                          genome_accessor, plp_params: Dict[str, Any],
                          pad: int = 50, filter_config: FilterConfig = None,
                          target_type: str = "RNA") -> List[Dict[str, Any]]:
    """
    Convenience function to design mutation probes.

    Args:
        mut_info_list: List of mutation information dictionaries
        genome_accessor: Function to access genome sequence
        plp_params: Probe design parameters
        pad: Padding around mutation site
        filter_config: Optional filter configuration
        target_type: "RNA" for mRNA targeting, "DNA" for cDNA targeting

    Returns:
        List of mutation information with designed probes
    """
    # Create default configurations
    search_config = SearchConfig()
    if filter_config is None:
        filter_config = FilterConfig()

    # Create designer
    designer = MutationProbeDesigner(search_config, filter_config, target_type=target_type)

    # Design probes
    return designer.design_mutation_probes(mut_info_list, genome_accessor, plp_params, pad)
