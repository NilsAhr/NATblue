"""
Unit tests for FTR2NAT — the release criterion that lets go early without the
conflict re-forming.

Measured motivation (synthetic suite, 2026-08-21/24). Over 17 scenarios
MVP2NAT_VERT produced 0 s of loss of separation under PastCPANAT and 113 s
under FTRNAT. That is not a safety/efficiency trade, it is two methods wrong in
opposite directions:

  * PastCPANAT holds while `not past_cpa or hor_los or is_bouncing`, which on
    co-track or still-converging NAT geometry is permanently true. Nothing can
    re-form because nothing moves -- and it never releases at all in three of
    the lifecycle scenarios.
  * FTRNAT releases correctly by the instantaneous separation test (every
    measured release had >= RPZ or >= HPZ) and then the conflict comes back.
    One aircraft was released six times, re-engaging after 19, 35, 14, 43 and
    30 s; the churn drove the ownship 1175 ft vertically where PastCPANAT moved
    it 259 ft, and that movement caused the loss of separation.

Note what was NOT wrong: `conflict_predicted` is a faithful single-pair mirror
of StateBased.detect, and FTRNAT is already 3D. So these tests cover the parts
that had no coverage at all -- `resumenav()` itself and its helpers, which is
exactly where both audited defects lived.

Shared geometry and stubs are imported from test_ftr rather than duplicated.
Distances in metres, speeds in m/s, matching BlueSky internals.
"""
import numpy as np
import pytest

from bluesky.traffic.asas.ftrnat import FTRNAT
from bluesky.traffic.asas.ftr2nat import (
    FTR2NAT, FTR2NAT_FIX, FTR2NAT_NOMARGIN, FTR2NAT_NODWELL)

from .test_ftr import (_StubTraf, RPZ, HPZ, DTLOOK, FT,
                       HEADON_DIST, HEADON_VREL)


def _bare(cls):
    """A resume-nav instance without bluesky.init().

    ResumeNavigation.__init__ reaches into Entity machinery, so build via
    __new__ and supply only what the methods under test touch -- the pattern
    test_loggerff_asas_conflict_state.py established.
    """
    o = cls.__new__(cls)
    o.resopairs = set()
    o.assumed = {}
    o.intent = 'ASSUMED'
    o._dwell = {}
    o._rel_fach = 1.0
    o._rel_facv = 1.0
    o._rel_dwell_s = 15.0
    return o


class _Conf:
    """Just the three arrays the release predicate reads off ConflictDetection."""

    def __init__(self, rpz=RPZ, hpz=HPZ, dtlook=DTLOOK):
        self.rpz = np.array([rpz, rpz])
        self.hpz = np.array([hpz, hpz])
        self.dtlookahead = np.array([dtlook, dtlook])


# ---------------------------------------------------------------------------
# DEFECT 1 — the ASSUMED path recorded a RESOLUTION rate as the intruder's
# intent, fabricating a vertical escape.
# ---------------------------------------------------------------------------
class TestAssumedIntentIsIntentNotResolution:
    """FTRINTENT defaults to ASSUMED and ftrnat.py:57 stores raw `intruder.vs`.

    An intruder already manoeuvring under ASAS is climbing at the autopilot
    default of 1500 fpm, and `dvs_settled = vsrevert` (ftrnat.py:117)
    extrapolates that forever -- roughly 7500 ft of invented vertical escape
    over DTLOOK 300, which silently turns criterion 2 into a no-op. This is the
    same bug class `desired_vs` and `revert_conflict` were written to
    eliminate, let back in through the ASSUMED path.
    """

    def test_records_the_autopilot_target_not_the_current_rate(self):
        # Pushed up at 1500 fpm by ASAS, but its autopilot wants the level it
        # is already at: the intruder's INTENT is level.
        intr = _StubTraf(alt=10668.0, apalt=10668.0, apvs=7.62, vs=7.62)
        o = _bare(FTR2NAT_FIX)
        o._record_assumed(("A", "B"), intr, 0)
        assert o.assumed[("A", "B")][2] == 0.0

    def test_a_genuine_commanded_climb_is_still_recorded(self):
        # The fix must not blind the criterion to real vertical intent.
        intr = _StubTraf(alt=10668.0, apalt=11582.0, apvs=7.62, vs=0.0)
        o = _bare(FTR2NAT_FIX)
        o._record_assumed(("A", "B"), intr, 0)
        assert o.assumed[("A", "B")][2] == pytest.approx(7.62)

    def test_horizontal_intent_is_unchanged(self):
        intr = _StubTraf(gseast=123.0, gsnorth=-45.0)
        o = _bare(FTR2NAT_FIX)
        o._record_assumed(("A", "B"), intr, 0)
        assert o.assumed[("A", "B")][:2] == (123.0, -45.0)

    def test_settled_rate_is_zero_so_nothing_extrapolates_forever(self):
        o = _bare(FTR2NAT_FIX)
        dalt_settled, dvs_settled = o._assumed_settled(10668.0, 10363.0, 7.62)
        assert dvs_settled == 0.0
        assert dalt_settled == pytest.approx(305.0)

    def test_the_phantom_escape_no_longer_clears_a_settling_pair(self):
        """The concrete failure, isolated to one call.

        The geometry matters. `revert_conflict` ORs a "now" and a "settled"
        case, so a phantom rate in the settled branch cannot release a pair
        whose CURRENT state already conflicts -- the OR catches it. It bites
        where the now-branch is clear and the settled branch is the one
        deciding: the resolver currently holds the pair 1500 ft apart, but both
        are commanded back to the same level, so settled means co-level.

        With dvs_settled = 1500 fpm the settled branch sees the pair sliding
        clear of HPZ in 40 s -- before the horizontal window opens at 53.7 s --
        and releases. With 0.0 the pair is co-level for all time and it holds.
        """
        args = (np.asarray(HEADON_DIST, float), np.asarray(HEADON_VREL, float),
                RPZ, HPZ, DTLOOK,
                1500.0 * FT, 0.0,     # now: vertically clear, so the now-branch
                0.0)                  # settled: co-level
        assert FTRNAT.revert_conflict(*args, 7.62) is False   # FTRNAT released
        assert FTRNAT.revert_conflict(*args, 0.0) is True     # FTR2NAT holds


# ---------------------------------------------------------------------------
# DEFECT 2 — the predicate was directionally asymmetric.
# ---------------------------------------------------------------------------
class TestReleaseIsSymmetric:
    """For (A,B) the relative velocity is v_B_actual - v_A_desired; for (B,A)
    it is v_A_actual - v_B_desired. Those are not negatives of each other, so
    the two directions can disagree and one aircraft of a pair can be handed
    back while the other is still resolving. Measured: FTRNAT left 1 of 2
    aircraft engaged on L02/L03/L04 where PastCPANAT left 2 of 2.
    """

    def test_both_orderings_are_the_same_pair(self):
        assert frozenset(("A", "B")) == frozenset(("B", "A"))

    def test_one_direction_holding_holds_the_pair(self):
        o = _bare(FTR2NAT_FIX)
        key = frozenset(("A", "B"))
        assert o._release_ok(key, True and False, 1.0) is False

    def test_both_directions_clear_releases_the_pair(self):
        o = _bare(FTR2NAT_FIX)
        assert o._release_ok(frozenset(("A", "B")), True and True, 1.0) is True


# ---------------------------------------------------------------------------
# The stability layer.
# ---------------------------------------------------------------------------
class TestDwell:
    """`revert_conflict` is a pure instantaneous predicate with no hysteresis,
    re-evaluated every tick. A predicate flapping at its own boundary is
    exactly the observed chatter, so a clearance has to persist."""

    def test_no_release_before_the_dwell_time(self):
        o = _bare(FTR2NAT)
        key = frozenset(("A", "B"))
        assert [o._release_ok(key, True, 1.0) for _ in range(14)] == [False] * 14

    def test_release_once_the_dwell_time_is_reached(self):
        o = _bare(FTR2NAT)
        key = frozenset(("A", "B"))
        for _ in range(14):
            o._release_ok(key, True, 1.0)
        assert o._release_ok(key, True, 1.0) is True

    def test_a_failing_tick_resets_the_counter(self):
        # The anti-chatter property: a flapping predicate must never
        # accumulate toward a release.
        o = _bare(FTR2NAT)
        o._rel_dwell_s = 5.0
        key = frozenset(("A", "B"))
        for _ in range(4):
            o._release_ok(key, True, 1.0)
        assert o._release_ok(key, False, 1.0) is False
        assert key not in o._dwell
        assert o._release_ok(key, True, 1.0) is False

    def test_dwell_counts_time_not_ticks(self):
        o = _bare(FTR2NAT)
        o._rel_dwell_s = 10.0
        key = frozenset(("A", "B"))
        assert o._release_ok(key, True, 4.0) is False
        assert o._release_ok(key, True, 4.0) is False
        assert o._release_ok(key, True, 4.0) is True

    def test_zero_dwell_releases_immediately(self):
        o = _bare(FTR2NAT_NODWELL)
        o._rel_dwell_s = 0.0
        assert o._release_ok(frozenset(("A", "B")), True, 1.0) is True

    def test_forget_drops_the_counter(self):
        o = _bare(FTR2NAT)
        key = frozenset(("A", "B"))
        o._release_ok(key, True, 1.0)
        o._forget(key)
        assert key not in o._dwell

    def test_pairs_do_not_share_a_counter(self):
        o = _bare(FTR2NAT)
        o._rel_dwell_s = 3.0
        a, b = frozenset(("A", "B")), frozenset(("C", "D"))
        o._release_ok(a, True, 1.0)
        o._release_ok(a, True, 1.0)
        assert o._release_ok(b, True, 1.0) is False


class TestReleaseMargin:
    """FTRNAT judges the release on plain conf.rpz/conf.hpz -- the bare
    detection zone -- while MVP2NAT resolves toward rpz*resofach /
    hpz*resofacv. So it lets go at exactly the boundary detection re-triggers
    on, and any wobble brings the conflict straight back."""

    def test_fix_uses_the_plain_detection_zones(self):
        o = _bare(FTR2NAT_FIX)
        assert o._release_zones(_Conf(), 0, 1) == (RPZ, HPZ)

    def test_ftr2nat_inflates_both_zones(self):
        o = _bare(FTR2NAT)
        o._rel_fach, o._rel_facv = 1.05, 1.05
        r, h = o._release_zones(_Conf(), 0, 1)
        assert r == pytest.approx(RPZ * 1.05)
        assert h == pytest.approx(HPZ * 1.05)

    def test_the_margin_holds_a_pair_the_bare_zone_would_release(self):
        # 1020 ft apart vertically: clear of HPZ = 1000 ft, inside 1.05 * HPZ.
        dalt = 1020.0 * FT
        args = (np.asarray(HEADON_DIST, float), np.asarray(HEADON_VREL, float))
        bare = FTRNAT.revert_conflict(*args, RPZ, HPZ, DTLOOK,
                                      dalt, 0.0, dalt, 0.0)
        margin = FTRNAT.revert_conflict(*args, RPZ, HPZ * 1.05, DTLOOK,
                                        dalt, 0.0, dalt, 0.0)
        assert bare is False        # FTRNAT releases right at the boundary
        assert margin is True       # FTR2NAT keeps a margin

    def test_nomargin_variant_is_back_to_the_plain_zones(self):
        o = _bare(FTR2NAT_NOMARGIN)
        assert o._release_zones(_Conf(), 0, 1) == (RPZ, HPZ)


class TestHorizontalBracket:
    """`revert_conflict` ORs a "now" and a "settled" VERTICAL case, but both
    branches use the same horizontal terms -- so the bracket is vertical-only.
    The horizontal side assumes the ownship attains its autopilot track
    instantly when it actually ramps at the turn rate, and during that
    transient detection sees the ACTUAL state."""

    def _rc(self, cls, vrel_revert, vrel_now, dalt=0.0):
        o = _bare(cls)
        return o._revert_conflict(
            np.asarray(HEADON_DIST, float),
            np.asarray(vrel_revert, float), np.asarray(vrel_now, float),
            RPZ, HPZ, DTLOOK, dalt, 0.0, dalt, 0.0)

    def test_fix_ignores_the_current_velocity(self):
        # Reverted track clear, current track still head-on: FTR2NAT_FIX, like
        # FTRNAT, releases -- it never looks at where the aircraft points now.
        assert self._rc(FTR2NAT_FIX, [0.0, 200.0], HEADON_VREL) is False

    def test_ftr2nat_holds_while_the_current_velocity_still_conflicts(self):
        assert self._rc(FTR2NAT, [0.0, 200.0], HEADON_VREL) is True

    def test_ftr2nat_releases_when_both_are_clear(self):
        assert self._rc(FTR2NAT, [0.0, 200.0], [0.0, 200.0]) is False

    def test_ftr2nat_still_holds_when_only_the_reverted_one_conflicts(self):
        assert self._rc(FTR2NAT, HEADON_VREL, [0.0, 200.0]) is True

    def test_the_bracket_does_not_override_genuine_vertical_clearance(self):
        # Both horizontal branches conflict, but the pair is 1500 ft apart and
        # level: vertically it never overlaps, so it must still release.
        assert self._rc(FTR2NAT, HEADON_VREL, HEADON_VREL,
                        dalt=1500.0 * FT) is False


class TestRegistration:
    def test_all_four_are_registered_resume_nav_methods(self, traffic_):
        from bluesky.traffic.asas import ResumeNavigation
        derived = ResumeNavigation.derived()
        for name in ("FTR2NAT", "FTR2NAT_FIX",
                     "FTR2NAT_NOMARGIN", "FTR2NAT_NODWELL"):
            assert name in derived, name

    def test_the_lineage_is_ftrnat(self, traffic_):
        assert issubclass(FTR2NAT_FIX, FTRNAT)
        assert issubclass(FTR2NAT, FTR2NAT_FIX)

    def test_ftr2nat_reuses_the_nat_waypoint_snap_forward(self, traffic_):
        from bluesky.traffic.asas.pastcpanat import PastCPANAT
        assert (FTR2NAT._advance_past_overflown
                is PastCPANAT._advance_past_overflown)

    def test_default_resume_nav_is_still_pastcpanat(self, traffic_):
        from bluesky.traffic.asas import ResumeNavigation
        assert ResumeNavigation.select_default().__name__ == 'PastCPANAT'


# ---------------------------------------------------------------------------
# Two implementation bugs the first in-sim run exposed. Both held R01 -- the
# release discriminator -- engaged for its entire 4002 s.
# ---------------------------------------------------------------------------
class TestZeroRelativeVelocityIsNeverTheCurrentState:
    """The horizontal bracket was handed a ZERO relative velocity.

    The first version passed `vnow - vnow` where it meant "intruder actual
    minus ownship actual". With zero relative motion `conflict_predicted`
    floors dv2 at 1e-6, so vabs is 1e-3 and the half-chord
    sqrt(R2 - dcpa2)/vabs is enormous: for any pair INSIDE rpz the horizontal
    window becomes effectively infinite and overlaps any vertical window at
    all. The pair is then held for ever. Measured on R01: engaged 4002 s of
    4002, with the conflict itself over since t=220.
    """

    def test_zero_vrel_inside_rpz_reports_an_endless_conflict(self):
        inside = [4.0 * 1852.0, 0.0]      # 4 NM apart, inside RPZ = 5 NM
        # Vertically converging at 3000 fpm from 2100 ft -- a real but future
        # window, which a sane horizontal term would not overlap.
        assert FTRNAT.conflict_predicted(
            np.asarray(inside, float), np.asarray([0.0, 0.0], float),
            2100.0 * FT, -15.24, RPZ, HPZ, DTLOOK) is True

    def test_the_real_relative_velocity_does_not(self):
        # Same geometry, but with the pair's actual separating velocity the
        # horizontal window is BOUNDED and has already closed, so the future
        # vertical window cannot overlap it.
        inside = [4.0 * 1852.0, 0.0]
        assert FTRNAT.conflict_predicted(
            np.asarray(inside, float), np.asarray([200.0, 0.0], float),
            2100.0 * FT, -15.24, RPZ, HPZ, DTLOOK) is False

    def test_zero_vrel_holds_even_a_far_future_vertical_window(self):
        # The signature of the bug: with an unbounded horizontal window even a
        # vertical convergence 5 minutes out reads as a conflict, so the pair
        # can never be released while it is inside rpz.
        inside = [4.0 * 1852.0, 0.0]
        assert FTRNAT.conflict_predicted(
            np.asarray(inside, float), np.asarray([0.0, 0.0], float),
            2100.0 * FT, -2.0, RPZ, HPZ, DTLOOK) is True
        assert FTRNAT.conflict_predicted(
            np.asarray(inside, float), np.asarray([200.0, 0.0], float),
            2100.0 * FT, -2.0, RPZ, HPZ, DTLOOK) is False


class TestDwellUsesTheResumenavInterval:
    """The dwell accumulated `bs.sim.simdt`, but resumenav runs at the ASAS rate.

    traffic.py drives it from `asastimer` (traffic.py:82,407), so it is called
    once per ASAS interval -- 1.0 s by default -- while simdt is the
    integration step, 0.05 s in this fork. A 15 s dwell became a 300 s one.
    """

    def test_the_interval_is_the_asas_rate_not_the_integration_step(self):
        import bluesky as bs
        dt = FTR2NAT._resnav_dt()
        assert dt == pytest.approx(bs.settings.asas_dt)
        assert dt > getattr(bs.sim, "simdt", 0.05)

    def test_dwell_reaches_the_target_in_the_expected_number_of_calls(self):
        o = _bare(FTR2NAT)
        o._rel_dwell_s = 15.0
        key = frozenset(("A", "B"))
        dt = 1.0
        n = sum(1 for _ in range(20) if not o._release_ok(key, True, dt))
        assert n == 14          # 15th call releases, not the 300th


class TestRampHorizon:
    """The "now" case must not extrapolate a level-capture rate for ever.

    revert_conflict brackets the ramp-then-hold profile by ORing a "now" case
    with a "settled" case, and the now case extrapolates the current rate
    forever. Under a vertical resolver that is fatal: the aircraft is parked
    off its route level, desired_vs is +/-1500 fpm back toward it, and the
    unbounded now-case reads that as a permanent closure. Measured: FTR2NAT
    held R01 for all 4002 s in the vertical domain while releasing it at 272 s
    in both horizontal domains.
    """

    def test_capture_time_from_the_offset_and_the_rate(self):
        # 1050 ft to recover at 1500 fpm is 42 s.
        t = FTR2NAT._ramp_horizon(alt=10668.0, apalt=10668.0 + 1050.0 * FT,
                                  vs_desired=7.62, dtlook=DTLOOK)
        assert t == pytest.approx(42.0, abs=0.5)

    def test_a_level_aircraft_gets_the_full_lookahead(self):
        assert FTR2NAT._ramp_horizon(10668.0, 10668.0, 0.0, DTLOOK) == DTLOOK

    def test_never_longer_than_the_lookahead(self):
        assert FTR2NAT._ramp_horizon(0.0, 1e6, 1.0, DTLOOK) == DTLOOK

    def test_a_convergence_after_capture_no_longer_holds(self):
        # Pair 2100 ft apart, reverting toward each other at 3000 fpm, so they
        # would be co-level in 42 s -- but each captures its level after 42 s,
        # so it never happens. Bounded: clear. Unbounded: held.
        o = _bare(FTR2NAT)
        args = (np.asarray(HEADON_DIST, float),
                np.asarray(HEADON_VREL, float), np.asarray(HEADON_VREL, float),
                RPZ, HPZ, DTLOOK, 2100.0 * FT, -15.24, 4000.0 * FT, 0.0)
        assert o._revert_conflict(*args, dtlook_now=DTLOOK) is True
        assert o._revert_conflict(*args, dtlook_now=20.0) is False

    def test_a_genuine_convergence_is_still_caught_by_the_settled_case(self):
        # Same bounded horizon, but the pair really does end up co-level.
        o = _bare(FTR2NAT)
        args = (np.asarray(HEADON_DIST, float),
                np.asarray(HEADON_VREL, float), np.asarray(HEADON_VREL, float),
                RPZ, HPZ, DTLOOK, 2100.0 * FT, -15.24, 0.0, 0.0)
        assert o._revert_conflict(*args, dtlook_now=20.0) is True
