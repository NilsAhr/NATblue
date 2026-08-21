"""
Unit tests for the MVP2NAT vertical runaway guard and per-domain channel
isolation (paper-3 ticket 07, sections 2 and 2b).

``TestClampVertical`` follows ``test_ftr.py`` and avoids the session-scoped
``traffic_`` fixture entirely — ``_clamp_vertical`` is a pure static method over
numpy input, so those tests pin the *maths* and run in milliseconds.

``TestVsAnchors`` and ``TestChannelOwnership`` DO take ``traffic_``: a
``ConflictResolution`` is a ``TrafficArrays`` Entity whose ``__init__`` reaches
for ``TrafficArrays.root.ntraf``, so instantiating one without ``bluesky.init()``
leaves a half-built object. They still assert policy, not simulation behaviour.

Units follow BlueSky internals: altitudes in metres, vertical speeds in m/s.
"""
import numpy as np

from bluesky.traffic.asas_nat.mvp2nat import (
    MVP2NAT, MVP2NAT_VERT, MVP2NAT_SPD, MVP2NAT_HDG, MVP2NAT_BOTH,
    MVP2NAT_VERT_LEGACY, MVP2NAT_HYBRID3,
    _MAX_RESO_VS, _MAX_VNAV_DEV_M,
)

FT = 0.3048
FPM = FT / 60.0


# ---------------------------------------------------------------------------
# Dive guard — the maths
# ---------------------------------------------------------------------------
class TestClampVertical:
    """MVP2NAT._clamp_vertical(vs, alt, anchor) -> (vs, alt), both bounded."""

    def test_vs_within_cap_is_untouched(self, traffic_):
        vs = np.array([1000.0 * FPM, -1200.0 * FPM])
        alt = np.array([10000.0, 10000.0])
        anchor = np.array([10000.0, 10000.0])
        out_vs, _ = MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.allclose(out_vs, vs)

    def test_runaway_vs_is_capped_both_signs(self, traffic_):
        # the observed failure mode was a ~44 000 fpm commanded dive
        vs = np.array([44000.0 * FPM, -44000.0 * FPM])
        alt = np.array([10000.0, 10000.0])
        anchor = np.array([10000.0, 10000.0])
        out_vs, _ = MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.allclose(out_vs, [_MAX_RESO_VS, -_MAX_RESO_VS])

    def test_altitude_bounded_below_anchor(self, traffic_):
        # a 20 000 ft dive command against a FL380 anchor
        anchor = np.array([38000.0 * FT])
        alt = np.array([18000.0 * FT])
        _, out_alt = MVP2NAT._clamp_vertical(np.array([0.0]), alt, anchor)
        assert np.isclose(out_alt[0], anchor[0] - _MAX_VNAV_DEV_M)

    def test_altitude_bounded_above_anchor(self, traffic_):
        anchor = np.array([38000.0 * FT])
        alt = np.array([58000.0 * FT])
        _, out_alt = MVP2NAT._clamp_vertical(np.array([0.0]), alt, anchor)
        assert np.isclose(out_alt[0], anchor[0] + _MAX_VNAV_DEV_M)

    def test_altitude_inside_the_band_is_untouched(self, traffic_):
        anchor = np.array([38000.0 * FT])
        alt = np.array([36000.0 * FT])
        _, out_alt = MVP2NAT._clamp_vertical(np.array([0.0]), alt, anchor)
        assert np.isclose(out_alt[0], alt[0])

    def test_nan_anchor_leaves_altitude_unclamped(self, traffic_):
        # no anchor yet (aircraft not under ASAS vertical control) -> the
        # altitude passes through, but the VS cap still applies
        anchor = np.array([np.nan])
        alt = np.array([18000.0 * FT])
        vs = np.array([44000.0 * FPM])
        out_vs, out_alt = MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.isclose(out_alt[0], alt[0])
        assert np.isclose(out_vs[0], _MAX_RESO_VS)

    def test_inputs_are_not_modified_in_place(self, traffic_):
        vs = np.array([44000.0 * FPM])
        alt = np.array([18000.0 * FT])
        anchor = np.array([38000.0 * FT])
        MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.isclose(vs[0], 44000.0 * FPM)
        assert np.isclose(alt[0], 18000.0 * FT)

    def test_guard_does_not_bind_the_hybrid3_hold(self, traffic_):
        # HYBRID3's co-track hold is bounded at 1500 fpm and targets roughly
        # 1.1 x HPZ (1100 ft) below the leader: both well inside the guard, so
        # the accepted iter3b behaviour must pass through unchanged.
        assert _MAX_RESO_VS > 1500.0 * FPM
        assert _MAX_VNAV_DEV_M > 1100.0 * FT


# ---------------------------------------------------------------------------
# Dive guard — the anchor bookkeeping that makes the bound CUMULATIVE
# ---------------------------------------------------------------------------
class TestVsAnchors:
    """MVP2NAT._update_vs_anchors keeps one anchor per aircraft under ASAS
    vertical control, and forgets it on release."""

    def test_anchor_set_on_first_engagement_and_held(self, traffic_):
        r = MVP2NAT()
        r._update_vs_anchors(["AC1"], {"AC1": 38000.0 * FT})
        assert r._vs_anchor["AC1"] == 38000.0 * FT
        # the aircraft has since descended; the anchor must NOT follow it --
        # that is exactly the per-tick re-anchoring that let the dive ratchet
        r._update_vs_anchors(["AC1"], {"AC1": 30000.0 * FT})
        assert r._vs_anchor["AC1"] == 38000.0 * FT

    def test_anchor_dropped_on_release_and_reset_on_re_engagement(self, traffic_):
        r = MVP2NAT()
        r._update_vs_anchors(["AC1"], {"AC1": 38000.0 * FT})
        r._update_vs_anchors([], {"AC1": 30000.0 * FT})
        assert "AC1" not in r._vs_anchor
        r._update_vs_anchors(["AC1"], {"AC1": 30000.0 * FT})
        assert r._vs_anchor["AC1"] == 30000.0 * FT

    def test_multiple_aircraft_are_tracked_independently(self, traffic_):
        r = MVP2NAT()
        r._update_vs_anchors(["AC1", "AC2"],
                             {"AC1": 38000.0 * FT, "AC2": 34000.0 * FT})
        r._update_vs_anchors(["AC2"], {"AC1": 20000.0 * FT, "AC2": 30000.0 * FT})
        assert "AC1" not in r._vs_anchor
        assert r._vs_anchor["AC2"] == 34000.0 * FT

    def test_unknown_callsign_is_ignored(self, traffic_):
        r = MVP2NAT()
        r._update_vs_anchors(["GHOST"], {})
        assert r._vs_anchor == {}

    def test_anchor_survives_the_conflict_ending_while_asas_still_owns_the_ac(
            self, traffic_):
        """REGRESSION (20250714 smoke cell, 271/568 aircraft ran away).

        `timesolveV` reverts to 1e9 the moment the pair leaves confpairs, but
        resume-nav keeps `active` True until it releases the aircraft, and
        update() keeps re-running resolve() while ANY pair is in conflict. If
        the anchor is dropped at that point the altitude clamp is skipped while
        the command is still applied, asasalttemp becomes a moving target below
        the CURRENT altitude, and the aircraft chases it into the ground
        (-23 308 ft observed).

        So: an aircraft that is still ASAS-active must KEEP its original
        anchor even after its vertical solve has gone away.
        """
        r = MVP2NAT()
        # tick 1: conflict present, anchor captured at FL317
        r._update_vs_anchors(["AC1"], {"AC1": 31703.0 * FT})
        assert r._vs_anchor["AC1"] == 31703.0 * FT
        # tick 2..N: conflict gone but resume-nav has not released the aircraft,
        # so it is still in the engaged set and the anchor must NOT move down
        for alt_ft in (27482.0, 14784.0, 2087.0, -2146.0):
            r._update_vs_anchors(["AC1"], {"AC1": alt_ft * FT})
            assert r._vs_anchor["AC1"] == 31703.0 * FT, alt_ft

    def test_the_clamp_bounds_the_descent_once_the_anchor_is_held(self,
                                                                 traffic_):
        """With the anchor held, the commanded altitude cannot fall further
        than _MAX_VNAV_DEV_M below the engagement altitude, so the aircraft
        levels off instead of descending indefinitely."""
        anchor = np.array([31703.0 * FT])
        runaway_target = np.array([-23308.0 * FT])      # the observed command
        _, out_alt = MVP2NAT._clamp_vertical(np.array([0.0]), runaway_target,
                                             anchor)
        assert np.isclose(out_alt[0], anchor[0] - _MAX_VNAV_DEV_M)
        assert out_alt[0] > 0.0                          # above ground


# ---------------------------------------------------------------------------
# Isolated-domain channel ownership
# ---------------------------------------------------------------------------
class TestChannelOwnership:
    """MVP2NAT._channel_owners() -> which channels the resolver commands.

    Isolated-domain semantics (ticket section 2b, decision 1): a domain commands
    ONLY its own channel; every other channel stays with the autopilot.
    """

    def test_vert_owns_only_the_vertical_channels(self, traffic_):
        assert MVP2NAT_VERT()._channel_owners() == {
            "hdg": False, "tas": False, "vs": True, "alt": True}

    def test_spd_owns_only_the_speed_channel(self, traffic_):
        assert MVP2NAT_SPD()._channel_owners() == {
            "hdg": False, "tas": True, "vs": False, "alt": False}

    def test_hdg_owns_only_the_heading_channel(self, traffic_):
        assert MVP2NAT_HDG()._channel_owners() == {
            "hdg": True, "tas": False, "vs": False, "alt": False}

    def test_both_owns_heading_and_speed_only(self, traffic_):
        assert MVP2NAT_BOTH()._channel_owners() == {
            "hdg": True, "tas": True, "vs": False, "alt": False}

    def test_the_four_domain_classes_isolate(self, traffic_):
        for cls in (MVP2NAT_VERT, MVP2NAT_SPD, MVP2NAT_HDG, MVP2NAT_BOTH):
            assert cls()._isolate_domain is True, cls.__name__

    def test_base_and_hybrid_do_not_isolate(self, traffic_):
        # The base and every hybrid subclass keep legacy channel behaviour, so
        # the accepted HYBRID3 / iter3b result stays byte-identical.
        assert MVP2NAT()._isolate_domain is False
        assert MVP2NAT_HYBRID3()._isolate_domain is False

    def test_legacy_variants_do_not_isolate(self, traffic_):
        assert MVP2NAT_VERT_LEGACY()._isolate_domain is False
        # ...but they keep their domain switches, so they are still "vertical"
        assert MVP2NAT_VERT_LEGACY().swresovert is True

    def test_isolated_returns_none_when_isolation_is_off(self, traffic_):
        # None is the sentinel that tells the channel properties to fall through
        # to legacy behaviour, so it must not be confused with an empty array.
        assert MVP2NAT()._isolated("tas") is None


# ---------------------------------------------------------------------------
# Ablation ladder + stock reference (synthetic suite, 2026-08-21)
#
# The ladder puts every mvp2nat mechanism behind a switch so the synthetic
# scenarios can attribute a behaviour change to ONE mechanism at a single git
# SHA. That is only sound if two things hold, and both are pinned here:
#   * the fully-enabled layer reproduces the original class exactly, and
#   * the default class still has every mechanism ON, so nothing that already
#     works changes.
# The at-scale version of the first check is compare_runs.py, which asserts
# MVP2NAT_VERT_L7 and MVP2NAT_VERT produce the same trajectories; these tests
# are the cheap version that runs in milliseconds.
# ---------------------------------------------------------------------------
from bluesky.traffic.asas_nat.mvp2nat import (            # noqa: E402
    MVP2NAT_STOCK, MVP2NAT_STOCK_VERT, MVP2NAT_VERT_L0, MVP2NAT_VERT_L4,
    MVP2NAT_VERT_L7, MVP2NAT_VERT_L8, MVP2NAT_VERT_L9, MVP2NAT_SPD_L9,
    MVP2NAT_HDG_L9, MVP2NAT_BOTH_L9,
    _ALL_SWITCHES, _LADDER, _DOMAINS)
import bluesky.traffic.asas_nat.mvp2nat as _m             # noqa: E402


def _switches(cls):
    obj = cls()
    return {s: getattr(obj, s) for s in _ALL_SWITCHES}


class TestAblationLadder:
    """The switch grid itself -- pure __init__ bookkeeping, no fixture needed."""

    def test_default_mvp2nat_has_every_mechanism_on(self, traffic_):
        # The whole refactor is safe only because the default is unchanged.
        sw = _switches(MVP2NAT)
        assert all(sw[s] for s in _ALL_SWITCHES if s != "_isolate_domain")
        assert sw["_isolate_domain"] is False       # base never isolates

    def test_l0_is_stock_every_switch_off(self, traffic_):
        assert all(v is False for v in _switches(MVP2NAT_VERT_L0).values())

    def test_ladder_is_cumulative(self, traffic_):
        # Each layer must be a strict superset of the previous one, otherwise
        # a layer could silently turn a mechanism back OFF and the attribution
        # of any observed change would be wrong.
        seen = set()
        for _, adds in _LADDER:
            assert not (set(adds) & seen), "layer re-enables %s" % adds
            seen |= set(adds)
        assert seen == set(_ALL_SWITCHES), "ladder does not cover every switch"

    def test_every_ladder_class_exists_for_every_domain(self, traffic_):
        for dom in _DOMAINS:
            for layer, _ in _LADDER:
                name = "MVP2NAT_%s_%s" % (dom, layer)
                assert hasattr(_m, name), "missing " + name

    def test_top_layer_matches_the_isolated_domain_class(self, traffic_):
        # The TOP layer is the production class. L8 became the top when the
        # BUG-01 divergence guard was added; L7 is deliberately left at the
        # pre-fix state so the before/after is measurable at one SHA.
        for lad, dom in ((MVP2NAT_VERT_L9, MVP2NAT_VERT),
                         (MVP2NAT_SPD_L9, MVP2NAT_SPD),
                         (MVP2NAT_HDG_L9, MVP2NAT_HDG),
                         (MVP2NAT_BOTH_L9, MVP2NAT_BOTH)):
            assert _switches(lad) == _switches(dom), lad.__name__
            a, b = lad(), dom()
            for f in ("swresohoriz", "swresospd", "swresohdg", "swresovert"):
                assert getattr(a, f) == getattr(b, f), (lad.__name__, f)

    def test_l7_is_the_pre_bug01_state(self, traffic_):
        # L7 must keep _sw_vert_diverge OFF and everything else ON, so that
        # "L7 vs L8" isolates exactly the BUG-01 fix.
        sw = _switches(MVP2NAT_VERT_L7)
        assert sw["_sw_vert_diverge"] is False and sw["_sw_vert_hold"] is False
        assert all(v for k, v in sw.items()
                   if k not in ("_sw_vert_diverge", "_sw_vert_hold"))

    def test_l8_is_bug01_only(self, traffic_):
        # L8 vs L9 isolates the BUG-04 cost fix from the BUG-01 safety fix.
        sw = _switches(MVP2NAT_VERT_L8)
        assert sw["_sw_vert_diverge"] is True and sw["_sw_vert_hold"] is False

    def test_symbreak_layer_is_where_the_breaker_turns_on(self, traffic_):
        # L4 is the layer the "head-on, both aircraft descend" report is
        # attributed to: below it the co-altitude sign term is stock's.
        assert _switches(MVP2NAT_VERT_L0)["_sw_symbreak"] is False
        assert _switches(MVP2NAT_VERT_L4)["_sw_symbreak"] is True

    def test_ladder_never_enables_a_hybrid_channel(self, traffic_):
        # The ladder isolates the BASE resolver; the hybrid controller is out
        # of scope, and a stray FL-hold would contaminate every layer above it.
        for dom in _DOMAINS:
            for layer, _ in _LADDER:
                o = getattr(_m, "MVP2NAT_%s_%s" % (dom, layer))()
                assert o.swfallbackhdg is False
                assert o.swspeedcotrack is False
                assert o.swaltholdcotrack is False


class TestStockReference:
    """MVP2NAT_STOCK reproduces stock MVP's algorithm, not its unit bug."""

    def test_stock_has_every_mechanism_off(self, traffic_):
        assert all(v is False for v in _switches(MVP2NAT_STOCK).values())

    def test_stock_uses_stock_domain_defaults(self, traffic_):
        # Stock MVP defaults to HORIZONTAL (mvp.py:12-18); mvp2nat inverted
        # that to vertical, so the reference has to invert it back.
        o = MVP2NAT_STOCK()
        assert (o.swresohoriz, o.swresospd, o.swresohdg, o.swresovert) == \
               (True, True, True, False)

    def test_stock_vert_is_the_rmethv_configuration(self, traffic_):
        o = MVP2NAT_STOCK_VERT()
        assert (o.swresovert, o.swresohoriz) == (True, False)
        assert all(v is False for v in _switches(o.__class__).values())


# ---------------------------------------------------------------------------
# Vertical divergence guard (BUG-01, found 2026-08-21 on synthetic L05)
#
# Observed: a co-altitude pair deadlocked in a PERIOD-2 limit cycle. Every tick
# both aircraft's commanded altitude reversed sign exactly (-499 ft, +499 ft,
# -499 ft, ...) at a commanded +/-3000 fpm, and after 2000 s of this they had
# separated by 51 ft. On the head-on scenarios the same cycle capped the split
# at ~668 ft against a 1000 ft HPZ, i.e. a residual intrusion.
#
# Mechanism, confirmed from the log and reproduced below:
#   1. Both level -> |vrel_z| < _VREL_VERT_FLOOR -> the symmetry breaker splits
#      the pair, one up one down.
#   2. That split IS a relative vertical rate (270 fpm = 1.37 m/s), which is
#      13x the 0.1 m/s floor, so the NEXT tick takes the normal MVP branch.
#   3. The normal branch is a convergence-NULLING rule: dv3 is set to oppose
#      the current relative vertical rate. Applied to a pair that is already
#      separating, it commands them back together.
#   4. tsolV = |drel_z / vrel_z| is tiny just after the split (11.6 s), so
#      dv3 = iV / tsolV saturates the 3000 fpm cap.
#   -> steps 1-4 repeat forever.
#
# MVP's vertical rule assumes the pair is vertically CONVERGING; a diverging
# pair is outside the model's domain. The guard detects that with the sign
# test drel_z * vrel_z > 0 (the gap is growing) and leaves the rate alone.
# ---------------------------------------------------------------------------
class TestVerticalDivergenceGuard:
    """MVP() must not reverse a pair that is already separating vertically."""

    class _Perf:
        def __init__(self, n):
            self.hmax = np.full(n, 42000.0 * FT)
            self.hmo = np.full(n, 43000.0 * FT)

    class _Ac:
        """Minimal ownship/intruder stand-in: MVP() only reads these fields."""
        def __init__(self, alt_ft, vs_fpm, gse=230.0, gsn=0.0):
            n = len(alt_ft)
            self.alt = np.array(alt_ft, dtype=float) * FT
            self.vs = np.array(vs_fpm, dtype=float) * FPM
            self.gseast = np.full(n, gse)
            self.gsnorth = np.full(n, gsn)
            self.perf = TestVerticalDivergenceGuard._Perf(n)

    class _Conf:
        def __init__(self, n, rpz_nm=5.0, hpz_ft=1000.0, dtlook=500.0):
            self.rpz = np.full(n, rpz_nm * 1852.0)
            self.hpz = np.full(n, hpz_ft * FT)
            self.dtlookahead = np.full(n, dtlook)

    def _dv3(self, cr, alt_ft, vs_fpm, dist_nm=4.0, tcpa=600.0, tLOS=600.0):
        own = self._Ac(alt_ft, vs_fpm)
        conf = self._Conf(len(alt_ft))
        dv, _, _ = cr.MVP(own, own, conf, qdr=90.0, dist=dist_nm * 1852.0,
                          tcpa=tcpa, tLOS=tLOS, idx1=0, idx2=1)
        return dv[2]

    def test_the_split_still_happens_when_both_are_level(self, traffic_):
        # Regime 1: co-altitude and level -> the symmetry breaker must fire,
        # otherwise the guard has broken the case it depends on.
        cr = MVP2NAT()
        assert self._dv3(cr, [35000.0, 35000.0], [0.0, 0.0]) != 0.0

    def test_a_diverging_pair_is_left_alone(self, traffic_):
        # Regime 2, the defect: ownship 51 ft BELOW and descending, intruder
        # climbing -- the exact state logged at L05 t=2001. The gap is growing,
        # so the vertical resolution must not command anything.
        cr = MVP2NAT()
        assert self._dv3(cr, [35000.0, 35051.0], [-132.0, 138.0]) == 0.0

    def test_the_mirror_ordering_is_also_left_alone(self, traffic_):
        # Both orderings of the pair are evaluated every tick; if only one is
        # guarded the cycle survives at half amplitude.
        cr = MVP2NAT()
        assert self._dv3(cr, [35051.0, 35000.0], [138.0, -132.0]) == 0.0

    def test_a_converging_pair_is_still_resolved(self, traffic_):
        # The guard must not disarm the resolver on the case it exists for:
        # intruder 2000 ft above and descending toward the ownship.
        cr = MVP2NAT()
        assert self._dv3(cr, [35000.0, 37000.0], [0.0, -500.0]) != 0.0

    def test_a_climb_through_from_below_is_still_resolved(self, traffic_):
        # S02's geometry: intruder 2000 ft below, climbing through.
        cr = MVP2NAT()
        assert self._dv3(cr, [35000.0, 33000.0], [0.0, 500.0]) != 0.0

    def test_the_guard_can_be_switched_off_for_the_ablation(self, traffic_):
        # L7 must keep the pre-fix behaviour so the before/after is measurable
        # at one SHA, exactly like the *_LEGACY classes.
        cr = MVP2NAT()
        cr._sw_vert_diverge = False
        assert self._dv3(cr, [35000.0, 35051.0], [-132.0, 138.0]) != 0.0


# ---------------------------------------------------------------------------
# Vertical hold (BUG-04, found 2026-08-21 on synthetic S01 after the BUG-01 fix)
#
# Observed: the pair needed 1050 ft of vertical separation (HPZ 1000 x RFACV
# 1.05) and built 4305 ft -- each aircraft flew ~2150 ft, four times the
# requirement, and stayed there until resume-nav let go.
#
# Mechanism: once the divergence guard sets dv3 = 0, the commanded altitude is
# asasalttemp = clip(vs * tsolV, +/- _MAX_DH_M) + alt, i.e. the aircraft is
# told to continue to a point _MAX_DH_M (2000 ft) away. Nothing in that
# expression knows how much separation the conflict actually required, so the
# manoeuvre size is set by a constant instead of by the geometry.
#
# The fix stops the vertical manoeuvre once every conflict of that aircraft is
# vertically satisfied -- separated by at least the required minimum AND not
# converging -- and holds the level reached until resume-nav releases it.
# Holding rather than returning to profile is deliberate: returning early would
# fly the aircraft back into the conflict it just solved.
# ---------------------------------------------------------------------------
class TestVerticalHold:
    """MVP2NAT._vert_hold_target: stop climbing once the requirement is met."""

    def test_no_hold_while_still_converging(self, traffic_):
        cr = MVP2NAT()
        # 3000 ft apart but closing: the resolver must stay engaged.
        assert cr._vert_satisfied(dalt=3000.0 * FT, dvs=-500.0 * FPM,
                                  hpz_req=1050.0 * FT) is False

    def test_no_hold_while_separation_is_insufficient(self, traffic_):
        cr = MVP2NAT()
        # Diverging, but only 400 ft apart against a 1050 ft requirement.
        assert cr._vert_satisfied(dalt=400.0 * FT, dvs=500.0 * FPM,
                                  hpz_req=1050.0 * FT) is False

    def test_hold_once_diverging_and_separated(self, traffic_):
        cr = MVP2NAT()
        assert cr._vert_satisfied(dalt=1100.0 * FT, dvs=500.0 * FPM,
                                  hpz_req=1050.0 * FT) is True

    def test_sign_symmetric(self, traffic_):
        # Both orderings of the pair see the same answer; if only one holds,
        # one aircraft levels off while the other keeps going.
        cr = MVP2NAT()
        assert cr._vert_satisfied(dalt=-1100.0 * FT, dvs=-500.0 * FPM,
                                  hpz_req=1050.0 * FT) is True

    def test_level_pair_at_sufficient_separation_holds(self, traffic_):
        # Diverging is "not converging", so a pair that is far enough apart and
        # steady needs no vertical action either.
        cr = MVP2NAT()
        assert cr._vert_satisfied(dalt=1100.0 * FT, dvs=0.0,
                                  hpz_req=1050.0 * FT) is True

    def test_hold_target_is_captured_once_and_held(self, traffic_):
        cr = MVP2NAT()
        cr._vert_hold.clear()
        cr._update_vert_hold(["AC1"], {"AC1": 35000.0 * FT})
        first = cr._vert_hold["AC1"]
        # A later tick at a different altitude must NOT move the target,
        # otherwise the hold ratchets exactly like the pre-guard dive did.
        cr._update_vert_hold(["AC1"], {"AC1": 34000.0 * FT})
        assert cr._vert_hold["AC1"] == first

    def test_hold_is_dropped_when_the_aircraft_is_no_longer_satisfied(self, traffic_):
        cr = MVP2NAT()
        cr._vert_hold.clear()
        cr._update_vert_hold(["AC1"], {"AC1": 35000.0 * FT})
        cr._update_vert_hold([], {"AC1": 35000.0 * FT})
        assert "AC1" not in cr._vert_hold

    def test_switch_off_disables_the_hold(self, traffic_):
        cr = MVP2NAT()
        cr._sw_vert_hold = False
        cr._vert_hold.clear()
        cr._update_vert_hold(["AC1"], {"AC1": 35000.0 * FT})
        assert cr._vert_hold == {}
