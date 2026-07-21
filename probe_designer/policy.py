"""Design rules for arm thermodynamics — the contract between gate and ranker.

``filtering`` decides whether a candidate is *acceptable*; ``scoring`` decides
whether it is *good*. Both questions are asked of the same two numbers (the arm
Tm target and the two-arm balance tolerance), but each subsystem used to keep
its own copy, and the copies had drifted: ``scoring.compute_target_score``
defaulted to ``min_arm_tm=50 / max_tm_diff=10`` regardless of what FilterConfig
said, and FilterConfig's own ``min_tm`` default (45) contradicted every shipped
YAML (50). A rule that lives in two places eventually disagrees with itself
(audit R7).

:class:`ThermoPolicy` holds those thresholds *and* the normalized 0-1 scoring
shapes derived from them. The composite weights stay in ``scoring`` — how much
a rule is worth relative to isoform coverage is a ranking decision, not a
thermodynamic one — but the rule itself has exactly one definition.

Gate vs score. Every rule here can be either a hard gate (reject) or a soft
preference (rank lower, survive). Soft is the default, following the Phase 1
decision to rank candidates rather than drop them on absolute Tm: a threshold
chosen from the literature is a statement about what is *preferable*, and
turning it into a cliff discards otherwise-good probes for being on the wrong
side of a number that carries a few degrees of model error anyway. Flip
``enforce_*`` to restore hard rejection.
"""
from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any, Optional

from probe_designer.chemistry import ReactionConditions

# Width of the two-sided Tm-proximity falloff (C). An arm this far from the
# target scores 0. Wide on purpose: hybrid-Tm model error is a few degrees, so
# a narrow falloff would rank on noise.
DEFAULT_TM_FALLOFF_C = 20.0

# Two-arm Tm balance tolerance (C). Larsson 2010 / Krzywkowski 2017 require the
# two arms to be comparable; field practice is <= 3-5 C, against the 10 (20 for
# mutation) this package used. Tightening it as a *gate* would have rejected
# 17.7% of 232 shipped probes, so it tightened as a scoring reference instead
# (audit R4; experiments/20260717_tm_buffer_config_impl/).
DEFAULT_MAX_TM_DIFF_C = 5.0

# Upper arm-Tm bound as an offset above the target, for specificity.
_MAX_TM_OFFSET_C = 20.0


@dataclass(frozen=True)
class ThermoPolicy:
    """Arm-Tm thresholds plus their normalized scoring shapes.

    Attributes:
        min_arm_tm: Target arm Tm (C) — the reaction-anchored
            ``lab_temp_c + tm_margin_c``. Doubles as the lower gate bound and
            the peak of the proximity score.
        max_arm_tm: Upper arm-Tm gate bound (C).
        max_tm_diff: Two-arm Tm balance tolerance (C); also the width over
            which the balance score decays to zero.
        enforce_tm_gate: Reject arms outside ``[min_arm_tm, max_arm_tm]``.
        enforce_tm_diff_gate: Reject probes whose arms differ by more than
            ``max_tm_diff``.
        tm_falloff_c: Width of the two-sided proximity falloff (C).
    """

    min_arm_tm: float = 50.0
    max_arm_tm: float = 70.0
    max_tm_diff: float = DEFAULT_MAX_TM_DIFF_C
    enforce_tm_gate: bool = False
    enforce_tm_diff_gate: bool = False
    tm_falloff_c: float = DEFAULT_TM_FALLOFF_C

    def __post_init__(self) -> None:
        if self.max_arm_tm < self.min_arm_tm:
            raise ValueError(
                f"max_arm_tm ({self.max_arm_tm}) must be >= min_arm_tm "
                f"({self.min_arm_tm})"
            )
        if self.tm_falloff_c <= 0:
            raise ValueError(f"tm_falloff_c must be > 0, got {self.tm_falloff_c}")
        if self.max_tm_diff < 0:
            raise ValueError(f"max_tm_diff must be >= 0, got {self.max_tm_diff}")

    # ---------------------------------------------------------- construction
    @classmethod
    def from_reaction(cls, reaction: ReactionConditions, **overrides: Any) -> "ThermoPolicy":
        """Anchor the Tm target to a reaction's ``lab_temp_c + tm_margin_c``.

        Raising the hybridization temperature raises what counts as a good arm,
        without anyone editing a threshold.
        """
        target = reaction.min_arm_tm
        base = {"min_arm_tm": target, "max_arm_tm": target + _MAX_TM_OFFSET_C}
        base.update({k: v for k, v in overrides.items() if v is not None})
        return cls(**base)

    # FilterConfig field name -> policy field name. The only place the two
    # vocabularies meet; FilterConfig speaks of a probe's Tm window, the policy
    # speaks of an arm's target.
    _CONFIG_ALIASES = {"min_tm": "min_arm_tm", "max_tm": "max_arm_tm"}

    @classmethod
    def resolve(cls, config: Any, **overrides: Any) -> "ThermoPolicy":
        """Build from a ``FilterConfig``, applying per-call keyword overrides.

        ``thermal_filter`` accepts one-off threshold overrides in ``**kwargs``;
        this folds them in so a caller-supplied value and a config value are
        never read through two different code paths. ``None`` overrides are
        ignored, so ``kwargs.get(name)`` on an absent key is harmless.
        """
        own = {f.name for f in fields(cls)}
        values = {}
        for cfg_name, policy_name in cls._CONFIG_ALIASES.items():
            if hasattr(config, cfg_name):
                values[policy_name] = getattr(config, cfg_name)
        for name in own:
            if name not in values and hasattr(config, name):
                values[name] = getattr(config, name)
        for name, value in overrides.items():
            if value is None:
                continue
            values[cls._CONFIG_ALIASES.get(name, name)] = value
        return cls(**{k: v for k, v in values.items() if k in own})

    # ------------------------------------------------------------- decisions
    def tm_range_ok(self, arm_tm: float) -> bool:
        """Does this arm Tm pass? Always True unless the hard gate is on."""
        if not self.enforce_tm_gate:
            return True
        return self.min_arm_tm <= arm_tm <= self.max_arm_tm

    def tm_balance_ok(self, tm_diff: float) -> bool:
        """Are the arms balanced enough? Always True unless the gate is on."""
        if not self.enforce_tm_diff_gate:
            return True
        return tm_diff <= self.max_tm_diff

    # ---------------------------------------------------------------- scores
    def tm_range_score(self, arm_tm: float) -> float:
        """Proximity of an arm Tm to the target, 1.0 at the target down to 0.0.

        Two-sided: an arm well below the target will not stay bound at the
        hybridization temperature, and one well above costs specificity and
        cross-probe Tm uniformity.
        """
        deviation = abs(arm_tm - self.min_arm_tm)
        return max(0.0, 1.0 - deviation / self.tm_falloff_c)

    def tm_balance_score(self, tm_diff: float) -> float:
        """Two-arm balance, 1.0 for identical arms down to 0.0 at the tolerance."""
        if self.max_tm_diff <= 0:
            return 0.0
        return max(0.0, 1.0 - tm_diff / self.max_tm_diff)


DEFAULT_POLICY = ThermoPolicy()

__all__ = [
    "DEFAULT_MAX_TM_DIFF_C",
    "DEFAULT_POLICY",
    "DEFAULT_TM_FALLOFF_C",
    "ThermoPolicy",
]
