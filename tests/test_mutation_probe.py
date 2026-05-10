#!/usr/bin/env python3
"""Tests for the mutation probe design module (probe_designer.ext.mutation_probe)."""

import pytest

from probe_designer.ext.mutation_probe import MutationProbeDesigner, design_mutation_probes
from probe_designer.config import SearchConfig, FilterConfig, DatabaseConfig
from probe_designer.database import DatabaseInterface


# ---------------------------------------------------------------------------
# Test: imports
# ---------------------------------------------------------------------------

class TestImports:
    def test_core_imports(self):
        """Core config and database classes are importable."""
        assert SearchConfig is not None
        assert FilterConfig is not None
        assert DatabaseConfig is not None
        assert DatabaseInterface is not None

    def test_mutation_module_imports(self):
        """MutationProbeDesigner and helper function are importable."""
        assert MutationProbeDesigner is not None
        assert design_mutation_probes is not None


# ---------------------------------------------------------------------------
# Test: MutationProbeDesigner basics
# ---------------------------------------------------------------------------

class TestMutationProbeDesigner:

    @pytest.fixture
    def designer(self):
        return MutationProbeDesigner(SearchConfig(), FilterConfig())

    def test_mutation_type_snp(self, designer):
        assert designer._get_mutation_type("A", "T", 100, 101) == "SNP"

    def test_mutation_type_del(self, designer):
        assert designer._get_mutation_type("AT", "-", 100, 102) == "DEL"

    def test_mutation_type_ins(self, designer):
        assert designer._get_mutation_type("-", "AT", 100, 100) == "INS"


# ---------------------------------------------------------------------------
# Test: genomic sequence mutation
# ---------------------------------------------------------------------------

class TestGenomicSequenceMutation:

    @pytest.fixture
    def designer(self):
        return MutationProbeDesigner(SearchConfig(), FilterConfig())

    def test_snp_mutation(self, designer):
        seq = "ATCGATCGATCG"
        # SNP: replace position 2..2 (C -> G)
        mutated, new_start, new_end = designer._mutate_genomic_sequence(
            seq, 2, 2, "C", "G", strand=1
        )
        assert mutated == "ATGGATCGATCG" or mutated[2] == "G"

    def test_deletion(self, designer):
        seq = "ATCGATCGATCG"
        # DEL: remove positions 2..3 (CG -> -)
        mutated, new_start, new_end = designer._mutate_genomic_sequence(
            seq, 2, 3, "CG", "-", strand=1
        )
        assert len(mutated) == len(seq) - 2

    def test_insertion(self, designer):
        seq = "ATCGATCGATCG"
        # INS: insert GG after position 2
        mutated, new_start, new_end = designer._mutate_genomic_sequence(
            seq, 2, 2, "-", "GG", strand=1
        )
        assert len(mutated) == len(seq) + 2


# ---------------------------------------------------------------------------
# Test: padlock arm extraction
# ---------------------------------------------------------------------------

class TestPadlockArmExtraction:

    @pytest.fixture
    def designer(self):
        return MutationProbeDesigner(SearchConfig(), FilterConfig())

    def test_arm_keys_present(self, designer):
        target_seq = "ATCGATCGATCGATCGATCG"
        arm_info = designer._extract_padlock_arms(target_seq, 10, 5, 5, 1)
        for key in ("probe_arm5", "probe_arm3", "target_sequence"):
            assert key in arm_info


# ---------------------------------------------------------------------------
# Test: DatabaseInterface creation
# ---------------------------------------------------------------------------

class TestDatabaseInterface:

    def test_create_ensembl_interface(self):
        config = DatabaseConfig(
            database_type="ensembl",
            organism="human",
            coord_system_version="GRCh38",
            max_retries=3,
        )
        db = DatabaseInterface(config)
        assert db is not None


# ---------------------------------------------------------------------------
# Test: MutationProbeDesignerInvader (iLock invader-overlap subclass)
# ---------------------------------------------------------------------------
#
# Canonical iLock invader chemistry requires the discriminating SNP nucleotide
# to be carried at BOTH (a) the probe's 5' body start (= probe_arm5[0]) and
# (b) the probe's 3' tip (= probe_arm3[-1]). The 3' tip + duplicated invader
# linker compete for the same target nt, gating FEN1 cleavage on SNP pairing.
#
# This requires arm5_binding and arm3_binding on the target to OVERLAP by one
# nucleotide at the SNP — i.e. pos3_target == pos5_target (vs the parent's
# pos5_target + 1).

class TestMutationProbeDesignerInvader:

    @pytest.fixture
    def designer(self):
        from probe_designer.ext.mutation_probe import MutationProbeDesignerInvader
        return MutationProbeDesignerInvader(SearchConfig(), FilterConfig())

    def test_subclass_imports(self):
        from probe_designer.ext.mutation_probe import MutationProbeDesignerInvader
        assert issubclass(MutationProbeDesignerInvader, MutationProbeDesigner)

    def test_rna_only_enforcement(self):
        """Subclass refuses target_type='DNA' (cDNA pass uses parent's swap logic)."""
        from probe_designer.ext.mutation_probe import MutationProbeDesignerInvader
        with pytest.raises(ValueError, match=r"RNA"):
            MutationProbeDesignerInvader(SearchConfig(), FilterConfig(), target_type="DNA")

    def test_pos3_overlaps_pos5(self, designer):
        """Invader overrides parent's pos3 = pos5 + 1 → pos3 = pos5 (1-nt overlap)."""
        info = designer._extract_padlock_arms(
            "AAAAAAAAAAGCCCCCCCCCC", pos5_local=10, arm5_len=5, arm3_len=5, strand=1
        )
        assert info["pos3_target"] == info["pos5_target"]

    def test_arm_geometry_strand_plus(self, designer):
        """strand=+1: arm5 ends at SNP, arm3 starts at SNP. Both share that nt.
        plus_seq[10] = 'G' (the SNP). RC of G = C, so probe arm5[0] = arm3[-1] = C."""
        plus_seq = "AAAAAAAAAAGCCCCCCCCCC"  # length 21, position 10 = G
        info = designer._extract_padlock_arms(plus_seq, 10, 5, 5, strand=1)
        # arm5 region on target = plus_seq[6:11] = "AAAAG"; probe_arm5 = RC = "CTTTT"
        assert info["probe_arm5"] == "CTTTT"
        # arm3 region on target = plus_seq[10:15] = "GCCCC"; probe_arm3 = RC = "GGGGC"
        assert info["probe_arm3"] == "GGGGC"
        # Invader invariant
        assert info["probe_arm5"][0] == info["probe_arm3"][-1] == "C"

    def test_arm_geometry_strand_minus(self, designer):
        """strand=-1: target_seq = RC(plus_seq); same overlap logic in target coords."""
        plus_seq = "AAAAAAAAAACGGGGGGGGGG"  # length 21, position 10 = C on plus_seq
        # target_seq = RC(plus_seq) = "CCCCCCCCCCGTTTTTTTTTT"; mRNA SNP at target pos 10 = G
        info = designer._extract_padlock_arms(plus_seq, 10, 5, 5, strand=-1)
        # arm5 on target = target_seq[6:11] = "CCCCG"; probe_arm5 = RC = "CGGGG"
        assert info["probe_arm5"] == "CGGGG"
        # arm3 on target = target_seq[10:15] = "GTTTT"; probe_arm3 = RC = "AAAAC"
        assert info["probe_arm3"] == "AAAAC"
        assert info["probe_arm5"][0] == info["probe_arm3"][-1] == "C"

    def test_target_sequence_dedup(self, designer):
        """Continuous mRNA region length = arm5_len + arm3_len - 1 (overlap by 1)."""
        plus_seq = "AAAAAAAAAAGCCCCCCCCCC"
        info = designer._extract_padlock_arms(plus_seq, 10, 5, 5, strand=1)
        assert len(info["target_sequence"]) == 9  # 5 + 5 - 1
        assert info["target_sequence"] == "AAAAGCCCC"

    def test_invader_invariant_holds_for_arbitrary_inputs(self, designer):
        """For any plus_seq + pos5 + arm lengths, arm5[0] == arm3[-1] = RC(SNP nt)."""
        cases = [
            ("AAAAAAAAAAGCCCCCCCCCC", 10, 6, 7, 1),
            ("CTAGCATCGAGCTAGCATCGA", 10, 5, 6, 1),
            ("TTTTTAAAAAACCCCCCGGGG", 12, 7, 5, 1),
            ("GCATGCATGCATGCATGCATG", 11, 6, 6, -1),
        ]
        for plus_seq, pos5, a5, a3, strand in cases:
            info = designer._extract_padlock_arms(plus_seq, pos5, a5, a3, strand)
            assert info["probe_arm5"][0] == info["probe_arm3"][-1], (
                f"Invader invariant violated for plus_seq={plus_seq}, pos5={pos5}, "
                f"strand={strand}: arm5={info['probe_arm5']}, arm3={info['probe_arm3']}"
            )

    def test_parent_pos3_unchanged(self):
        """Sanity: parent class still uses pos3 = pos5 + 1 (no leak from subclass)."""
        parent = MutationProbeDesigner(SearchConfig(), FilterConfig())
        info = parent._extract_padlock_arms(
            "AAAAAAAAAAGCCCCCCCCCC", pos5_local=10, arm5_len=5, arm3_len=5, strand=1
        )
        assert info["pos3_target"] == info["pos5_target"] + 1


# ---------------------------------------------------------------------------
# Test: verify_iLock_probe — runtime guard against stitch geometry swaps
# ---------------------------------------------------------------------------
#
# `verify_iLock_probe` is the hard-fail runtime check that step3 stitch
# templates call after assembling each iLock probe. It guards against the two
# bugs found in 2026-05:
#
#   (a) STITCH GEOMETRY: BZ09-era stitch placed `plp_arm3` at the 5' body and
#       `plp_arm5` at the 3' end. This is WRONG: `plp_arm5` (= RC of mRNA 5'
#       side, per mutation_probe.py convention) must sit at the probe 5' body
#       so its 5' tip lands at the ligation junction on the target. If swapped,
#       the probe binds but cannot ligate (5' tip and 3' tip land at the FAR
#       ENDS of the binding region, ~40 nt apart).
#   (b) INVADER OVERLAP: the duplicated 1-nt linker after the flap must equal
#       both `plp_arm5[0]` and `plp_arm3[-1]` (= RC of the SNP). If not, FEN1
#       cleavage isn't strictly SNP-conditional.

class TestVerifyILockProbe:

    @pytest.fixture
    def verify(self):
        from probe_designer.ext.mutation_probe import verify_iLock_probe
        return verify_iLock_probe

    @pytest.fixture
    def good_probe(self):
        """A geometrically-correct, invader-overlap-correct iLock probe.
        arm5 = RC(mRNA 5' side), arm3 = RC(mRNA 3' side), arm5[0] == arm3[-1]."""
        flap = 'CGTTGCTGTGGCG'
        arm5 = 'CTTTT'         # arm5[0] = C
        arm3 = 'GGGGC'         # arm3[-1] = C; arm3[-1] == arm5[0] ✓
        backbone = 'TCCCTACACGAC'
        linker = arm3[-1].upper()    # 'C'
        probe = flap + linker + arm5 + backbone + arm3
        return {'probe': probe, 'flap': flap, 'arm5': arm5, 'arm3': arm3}

    def test_passes_for_correct_probe(self, verify, good_probe):
        # No exception — silent pass.
        verify(good_probe['probe'], good_probe['arm5'], good_probe['arm3'],
               flap=good_probe['flap'])

    def test_rejects_swapped_geometry_BZ09_style(self, verify, good_probe):
        """BZ09 stitch swap: arm3 placed at 5' body, arm5 placed at 3' end.
        Must raise ValueError with a 'geometry' / 'swap' clue in the message."""
        flap = good_probe['flap']
        arm5 = good_probe['arm5']; arm3 = good_probe['arm3']
        bad_linker = arm5[-1].upper()     # BZ09 uses arm5[-1] (= last nt of arm5)
        bad_probe = flap + bad_linker + arm3 + 'TCCCTACACGAC' + arm5
        with pytest.raises(ValueError, match=r"(?i)swap|geom|wrong"):
            verify(bad_probe, arm5, arm3, flap=flap)

    def test_rejects_missing_flap(self, verify, good_probe):
        """Probe must start with the flap; otherwise reject."""
        no_flap = good_probe['probe'][len(good_probe['flap']):]  # strip flap
        with pytest.raises(ValueError, match=r"(?i)flap"):
            verify(no_flap, good_probe['arm5'], good_probe['arm3'],
                   flap=good_probe['flap'])

    def test_rejects_invader_invariant_violation(self, verify):
        """arm5[0] != arm3[-1] -> linker isn't SNP_compl -> reject."""
        flap = 'CGTTGCTGTGGCG'
        arm5 = 'AAAAA'         # arm5[0] = A
        arm3 = 'GGGGT'         # arm3[-1] = T (≠ A)
        backbone = 'TCCCTACACGAC'
        linker = arm3[-1].upper()
        probe = flap + linker + arm5 + backbone + arm3
        with pytest.raises(ValueError, match=r"(?i)invader|snp"):
            verify(probe, arm5, arm3, flap=flap)

    def test_legacy_flap_accepted(self, verify):
        """Legacy BZ09 flap TATATCCCTATAT must be accepted when explicitly passed."""
        flap = 'TATATCCCTATAT'
        arm5 = 'CTTTT'; arm3 = 'GGGGC'
        backbone = 'TCCCTACACGAC'
        probe = flap + arm3[-1].upper() + arm5 + backbone + arm3
        verify(probe, arm5, arm3, flap=flap)  # no exception

    def test_does_not_silently_match_arm_swap(self, verify):
        """If arm5 and arm3 happen to be each other's reverse, the swap might
        accidentally parse — guard MUST still reject based on column identity."""
        flap = 'CGTTGCTGTGGCG'
        # Construct arm5 / arm3 such that swap places different sequences at body/end
        arm5 = 'CAAAA'   # arm5[0] = C
        arm3 = 'TTTTC'   # arm3[-1] = C; invader OK
        backbone = 'TCCCTACACGAC'
        # BZ09 swap probe — body is arm3, end is arm5
        bad_probe = flap + arm5[-1].upper() + arm3 + backbone + arm5
        with pytest.raises(ValueError):
            verify(bad_probe, arm5, arm3, flap=flap)
