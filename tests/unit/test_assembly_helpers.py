"""Tests for the chemistry-agnostic + mutation-specific assembly helpers."""
from __future__ import annotations

import pytest

from probe_designer.ext.mutation.probe import assemble_ilock, verify_iLock_probe
from probe_designer.probe_assembly import assemble_plain_padlock


# ---------------------------------------------------------------------------
# assemble_plain_padlock — chemistry-agnostic primitive
# ---------------------------------------------------------------------------


def test_assemble_plain_padlock_concatenates_in_order():
    # Layout: arm5 + backbone + arm3
    out = assemble_plain_padlock(arm5="aaaa", arm3="tttt", backbone="BBBBB")
    assert out == "aaaaBBBBBtttt"


def test_assemble_plain_padlock_preserves_input_casing():
    # Caller is responsible for casing — function does not normalize.
    out = assemble_plain_padlock(arm5="Arm5", arm3="Arm3", backbone="Backbone")
    assert out == "Arm5BackboneArm3"


# ---------------------------------------------------------------------------
# assemble_ilock — geometry + invariant must hold; verify_iLock_probe runs
# ---------------------------------------------------------------------------


@pytest.fixture
def valid_invader_arms():
    """Arms satisfying the invader invariant arm5[0] == arm3[-1].

    Used as a happy-path fixture for assemble_ilock.
    """
    # SNP complement = "G". arm5 starts with G; arm3 ends with G.
    return {
        "arm5": "GAAATTTGCC",   # length 10, starts with G
        "arm3": "TTACCAAAAG",   # length 10, ends with G
        "backbone": "AAAAAACCCCCCTTTTTTGGGGGG",
        "flap": "CGTTGCTGTGGCG",
    }


def test_assemble_ilock_canonical_layout(valid_invader_arms):
    probe = assemble_ilock(**valid_invader_arms)
    # Layout: flap + linker(arm3[-1].upper()) + arm5.lower() + backbone.upper() + arm3.lower()
    expected = (
        valid_invader_arms["flap"]
        + valid_invader_arms["arm3"][-1].upper()
        + valid_invader_arms["arm5"].lower()
        + valid_invader_arms["backbone"].upper()
        + valid_invader_arms["arm3"].lower()
    )
    assert probe == expected


def test_assemble_ilock_runs_verify(valid_invader_arms):
    """verify_iLock_probe must accept the assembled probe."""
    probe = assemble_ilock(**valid_invader_arms)
    # Should not raise.
    verify_iLock_probe(
        probe,
        plp_arm5=valid_invader_arms["arm5"].lower(),
        plp_arm3=valid_invader_arms["arm3"].lower(),
        flap=valid_invader_arms["flap"],
    )


def test_assemble_ilock_rejects_invader_invariant_violation(valid_invader_arms):
    """If arm5[0] != arm3[-1], verify_iLock_probe inside assemble_ilock fires."""
    bad = dict(valid_invader_arms)
    # arm5 starts with A; arm3 still ends with G — invariant broken.
    bad["arm5"] = "AAAATTTGCC"
    with pytest.raises(ValueError, match="invader"):
        assemble_ilock(**bad)


def test_assemble_ilock_rejects_empty_arm3():
    with pytest.raises(ValueError, match="arm3 must not be empty"):
        assemble_ilock(arm5="GAAA", arm3="", backbone="TTTT", flap="CGTTGCTGTGGCG")
