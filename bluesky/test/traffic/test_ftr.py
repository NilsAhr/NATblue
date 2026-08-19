"""
Unit tests for the Free-To-Revert (FTR) resume-navigation criteria.

These tests deliberately avoid the session-scoped ``traffic_`` fixture in
conftest.py (which runs a full ``bluesky.init()``). Every criterion under test
is a pure static method taking plain numpy input, so the tests run in
milliseconds and pin the *maths* rather than the simulation wiring.

Sign convention throughout (matches bluesky/traffic/asas/statebased.py, where
every relative quantity for a conflict pair is intruder-minus-ownship):

    dist  = p_intruder - p_ownship   [m]   horizontal, (east, north)
    vrel  = v_intruder - v_ownship   [m/s] horizontal, (east, north)
    dalt  = alt_intruder - alt_ownship     [m]
    dvs   = vs_intruder  - vs_ownship      [m/s]

``bool(...)`` wraps every predicate result because numpy comparisons yield
``np.bool_``, which is not identical to Python's ``True``/``False``.
"""
import numpy as np
import pytest

from bluesky.traffic.asas.ftr import FTR
from bluesky.traffic.asas.ftrnat import FTRNAT


NM = 1852.0
FT = 0.3048

RPZ = 5.0 * NM      # 9260 m  — ZONER 5
HPZ = 1000.0 * FT   # 304.8 m — ZONEDH 1000
DTLOOK = 300.0      # s       — DTLOOK 300


def _clears(dist, vrel, rpz=RPZ):
    return bool(FTR.clears(np.asarray(dist, dtype=float),
                           np.asarray(vrel, dtype=float), rpz))


# ---------------------------------------------------------------------------
# FTR.clears — the horizontal forward-CPA test (criterion kernel)
# ---------------------------------------------------------------------------
class TestClears:
    """FTR.clears(dist, vrel, rpz) -> True if reverting keeps CPA outside rpz."""

    def test_head_on_collision_course_does_not_clear(self):
        # Intruder 10 km east, closing head-on at 200 m/s: CPA distance is zero.
        assert _clears([10000.0, 0.0], [-200.0, 0.0]) is False

    def test_offset_pass_outside_rpz_clears(self):
        # Same closure but offset 10 km north: miss distance 10 km > 5 NM.
        assert _clears([10000.0, 10000.0], [-200.0, 0.0]) is True

    def test_diverging_pair_judged_on_present_distance(self):
        # tcpa < 0. The pair is opening and already outside rpz -> clear. This
        # is the branch that stops FTR being fooled by a CPA lying in the past.
        assert _clears([10000.0, 0.0], [200.0, 0.0]) is True

    def test_diverging_pair_still_inside_rpz_does_not_clear(self):
        # Opening, but currently in horizontal LOS -> must not release.
        assert _clears([5000.0, 0.0], [200.0, 0.0]) is False

    def test_zero_relative_velocity_uses_distance(self):
        assert _clears([10000.0, 0.0], [0.0, 0.0]) is True
        assert _clears([5000.0, 0.0], [0.0, 0.0]) is False


# ---------------------------------------------------------------------------
# FTR.desired_spd_trk — the wind triangle
# ---------------------------------------------------------------------------
class _StubTraf:
    """Minimal stand-in for bs.traf exposing only what desired_spd_trk reads."""

    class _Ap:
        def __init__(self, trk, tas):
            self.trk = np.array([trk])
            self.tas = np.array([tas])

    def __init__(self, trk=90.0, tas=200.0, windeast=0.0, windnorth=0.0,
                 gseast=0.0, gsnorth=0.0, alt=9144.0, apalt=9144.0,
                 apvs=7.62, vs=0.0):
        self.ap = self._Ap(trk, tas)
        self.ap.alt = np.array([apalt])
        self.ap.vs = np.array([apvs])
        self.windeast = np.array([windeast])
        self.windnorth = np.array([windnorth])
        self.gseast = np.array([gseast])
        self.gsnorth = np.array([gsnorth])
        self.alt = np.array([alt])
        self.vs = np.array([vs])


class TestDesiredSpdTrk:
    """The reason FTR is worth having in NATBlue: it is wind-aware."""

    def test_zero_wind_reduces_to_tas_along_track(self):
        traf = _StubTraf(trk=90.0, tas=200.0)
        assert FTR.desired_spd_trk(traf, 0) == pytest.approx([200.0, 0.0],
                                                             abs=1e-9)

    def test_pure_tailwind_adds_to_ground_speed(self):
        traf = _StubTraf(trk=90.0, tas=200.0, windeast=50.0)
        assert FTR.desired_spd_trk(traf, 0) == pytest.approx([250.0, 0.0],
                                                             abs=1e-9)

    def test_pure_headwind_subtracts_from_ground_speed(self):
        traf = _StubTraf(trk=90.0, tas=200.0, windeast=-50.0)
        assert FTR.desired_spd_trk(traf, 0) == pytest.approx([150.0, 0.0],
                                                             abs=1e-9)

    def test_crosswind_crabs_and_keeps_track_good(self):
        # 50 m/s from the south while tracking east: ground velocity stays due
        # east (track made good) but shrinks to sqrt(tas^2 - wperp^2).
        traf = _StubTraf(trk=90.0, tas=200.0, windnorth=50.0)
        expected = np.sqrt(200.0 ** 2 - 50.0 ** 2)
        assert FTR.desired_spd_trk(traf, 0) == pytest.approx([expected, 0.0],
                                                             abs=1e-9)

    def test_crosswind_exceeding_tas_falls_back_to_current_gs(self):
        traf = _StubTraf(trk=90.0, tas=200.0, windnorth=250.0,
                         gseast=11.0, gsnorth=22.0)
        assert FTR.desired_spd_trk(traf, 0) == pytest.approx([11.0, 22.0])

    def test_northerly_track_uses_east_north_ordering(self):
        # Guards the (sin, cos) ordering: tracking north must give +north.
        traf = _StubTraf(trk=0.0, tas=200.0)
        assert FTR.desired_spd_trk(traf, 0) == pytest.approx([0.0, 200.0],
                                                             abs=1e-9)


# ---------------------------------------------------------------------------
# Two criteria vs one: the discriminating geometry
# ---------------------------------------------------------------------------
def test_criterion2_holds_where_criterion1_would_release():
    """An overtake in which the intruder is *currently* speed-matched but will
    revert to a slower cruise speed.

    Criterion 1 (intruder holds its observed velocity) sees zero closure and
    releases. Criterion 2 (intruder also reverts to its autopilot velocity)
    sees the ownship close from 12 km to zero miss distance and holds.

    This is the whole point of the two-criterion method — a single-criterion
    release cannot see it.
    """
    dist = [12000.0, 0.0]                 # intruder 12 km ahead (east)
    vown = np.array([250.0, 0.0])         # ownship reverts to 250 m/s east
    vintr_observed = np.array([250.0, 0.0])   # currently speed-matched
    vintr_reverted = np.array([200.0, 0.0])   # will slow back to 200 m/s

    assert _clears(dist, vintr_observed - vown) is True
    assert _clears(dist, vintr_reverted - vown) is False


# ---------------------------------------------------------------------------
# FTRNAT.conflict_predicted — the 3D criterion (mirrors statebased detection)
# ---------------------------------------------------------------------------
def _predicted(dist, vrel, dalt, dvs, dtlook=DTLOOK):
    return bool(FTRNAT.conflict_predicted(np.asarray(dist, dtype=float),
                                          np.asarray(vrel, dtype=float),
                                          float(dalt), float(dvs),
                                          RPZ, HPZ, dtlook))


# A head-on closing geometry shared by several cases below. Intruder 20 km east
# closing at 200 m/s -> tcpa 100 s, zero miss distance, horizontal conflict
# window [53.7 s, 146.3 s].
HEADON_DIST = [20000.0, 0.0]
HEADON_VREL = [-200.0, 0.0]


class TestConflictPredicted:
    """Free to revert iff the reverted state is NOT predicted in conflict."""

    def test_vertical_escape_clears_a_converging_pair(self):
        # THE case that motivates FTRNAT: co-level, converging horizontally,
        # but the vertical resolution opens past HPZ (30.5 s) well before the
        # horizontal conflict window opens (53.7 s). Free to revert.
        assert _predicted(HEADON_DIST, HEADON_VREL, dalt=0.0, dvs=10.0) is False

    def test_slow_vertical_escape_still_in_conflict(self):
        # Same geometry, 1 m/s vertical opening: still inside HPZ (until
        # 304.8 s) while horizontally in conflict -> hold the resolution.
        assert _predicted(HEADON_DIST, HEADON_VREL, dalt=0.0, dvs=1.0) is True

    def test_co_level_zero_vs_reduces_to_horizontal(self):
        # dvs == 0 and |dalt| < HPZ: vertically in conflict for all time, so
        # the horizontal criterion alone decides.
        assert _predicted(HEADON_DIST, HEADON_VREL, dalt=0.0, dvs=0.0) is True

    def test_level_but_vertically_separated_is_never_in_conflict(self):
        # dvs == 0 and |dalt| > HPZ: no vertical overlap ever, regardless of
        # how badly the pair converges horizontally.
        assert _predicted(HEADON_DIST, HEADON_VREL, dalt=500.0, dvs=0.0) is False

    def test_horizontally_clear_is_never_in_conflict(self):
        # Miss distance 20 km > RPZ, co-level and co-altitude.
        assert _predicted([20000.0, 20000.0], HEADON_VREL,
                          dalt=0.0, dvs=0.0) is False

    def test_conflict_beyond_dtlook_is_not_predicted(self):
        # Horizontal window opens at 53.7 s; with DTLOOK 30 s it is out of range.
        assert _predicted(HEADON_DIST, HEADON_VREL, dalt=0.0, dvs=1.0,
                          dtlook=30.0) is False

    def test_already_past_the_conflict_is_not_predicted(self):
        # Diverging: the whole conflict window lies in the past (toutconf < 0).
        assert _predicted(HEADON_DIST, [200.0, 0.0], dalt=0.0, dvs=0.0) is False

    def test_vertical_closure_into_the_horizontal_window_is_a_conflict(self):
        # Starts vertically clear (600 m) but closes at 5 m/s, entering HPZ at
        # ~59 s — inside the horizontal window [53.7, 146.3].
        assert _predicted(HEADON_DIST, HEADON_VREL,
                          dalt=600.0, dvs=-5.0) is True

    def test_vertical_closure_after_the_horizontal_window_clears(self):
        # Same closure rate but starting 1500 m apart: enters HPZ at ~239 s,
        # after the pair is horizontally clear at 146.3 s.
        assert _predicted(HEADON_DIST, HEADON_VREL,
                          dalt=1500.0, dvs=-5.0) is False

    def test_zero_relative_velocity_co_altitude_is_a_conflict(self):
        # Formation-flying degenerate case: stationary relative to each other,
        # inside both zones.
        assert _predicted([5000.0, 0.0], [0.0, 0.0], dalt=0.0, dvs=0.0) is True

    def test_zero_relative_velocity_outside_rpz_is_not_a_conflict(self):
        assert _predicted([20000.0, 0.0], [0.0, 0.0], dalt=0.0, dvs=0.0) is False


def test_ftrnat_criterion2_holds_where_criterion1_would_release_in_3d():
    """The 3D analogue of the two-criterion discriminator.

    The intruder is currently climbing away, so criterion 1 clears on the
    vertical escape. But its autopilot will level off the moment ASAS releases
    it — criterion 2 sees no vertical escape at all and holds the resolution.
    """
    dist = np.asarray(HEADON_DIST)
    vrel = np.asarray(HEADON_VREL)

    # Criterion 1: intruder holds its observed 10 m/s climb -> clears.
    assert bool(FTRNAT.conflict_predicted(dist, vrel, 0.0, 10.0,
                                          RPZ, HPZ, DTLOOK)) is False
    # Criterion 2: intruder reverts to level flight -> conflict remains.
    assert bool(FTRNAT.conflict_predicted(dist, vrel, 0.0, 0.0,
                                          RPZ, HPZ, DTLOOK)) is True


# ---------------------------------------------------------------------------
# Wiring
# ---------------------------------------------------------------------------
def test_ftr_and_ftrnat_are_registered_resume_nav_methods():
    from bluesky.traffic.asas import ResumeNavigation
    derived = ResumeNavigation.derived()
    assert 'FTR' in derived
    assert 'FTRNAT' in derived
    assert issubclass(FTRNAT, FTR)


def test_ftrnat_reuses_the_nat_waypoint_snap_forward():
    # FTRNAT must reuse PastCPANAT's overflown-waypoint recovery rather than
    # the stock findact+direct, or a sustained vertical manoeuvre can leave the
    # FMS pointing backwards.
    from bluesky.traffic.asas.pastcpanat import PastCPANAT
    assert FTRNAT._advance_past_overflown is PastCPANAT._advance_past_overflown


def test_default_resume_nav_is_still_pastcpanat():
    # Adding FTR must not change what an existing scenario gets from RESO.
    from bluesky.traffic.asas import ResumeNavigation
    assert ResumeNavigation.select_default().__name__ == 'PastCPANAT'


# ---------------------------------------------------------------------------
# FTRNAT.desired_vs — ap.vs is a DEFAULT RATE, not an intent
# ---------------------------------------------------------------------------
class TestDesiredVs:
    """Autopilot.vs falls back to vsdef (1500 fpm, autopilot.py:101) whenever no
    vertical speed is selected, so a level aircraft still reports ap.vs = 7.62.

    Reading that as a real climb makes every co-level pair look like it escapes
    vertically in ~40 s, releasing the resolution almost immediately. The signed
    rate must come from the altitude error instead, mirroring
    Traffic.update_airspeed: target_vs = swaltsel * sign(delta_alt) * |ap.vs|.
    """

    def test_level_at_commanded_altitude_is_zero(self):
        traf = _StubTraf(alt=9144.0, apalt=9144.0, apvs=7.62, vs=0.0)
        assert FTRNAT.desired_vs(traf, 0) == 0.0

    def test_commanded_climb_is_positive_vsdef(self):
        traf = _StubTraf(alt=9144.0, apalt=10668.0, apvs=7.62, vs=0.0)
        assert FTRNAT.desired_vs(traf, 0) == pytest.approx(7.62)

    def test_commanded_descent_is_negative_vsdef(self):
        traf = _StubTraf(alt=10668.0, apalt=9144.0, apvs=7.62, vs=0.0)
        assert FTRNAT.desired_vs(traf, 0) == pytest.approx(-7.62)

    def test_sign_comes_from_altitude_error_not_from_ap_vs(self):
        # aporasas forces its vs positive (aporasas.py:43); a negative value
        # must still yield a climb when the target is above.
        traf = _StubTraf(alt=9144.0, apalt=10668.0, apvs=-7.62, vs=0.0)
        assert FTRNAT.desired_vs(traf, 0) == pytest.approx(7.62)

    def test_within_capture_band_counts_as_level(self):
        # 2 ft of residual altitude error is inside the 10 ft dead band.
        traf = _StubTraf(alt=9144.0, apalt=9144.0 + 2.0 * FT, apvs=7.62)
        assert FTRNAT.desired_vs(traf, 0) == 0.0

    def test_held_below_cruise_reverts_to_a_climb(self):
        # The FL-hold case: MVP2NAT parks the aircraft 1000 ft low via cr.alt,
        # but ap.alt still points at cruise, so reverting means climbing back.
        traf = _StubTraf(alt=9144.0 - 1000.0 * FT, apalt=9144.0, apvs=7.62,
                         vs=0.0)
        assert FTRNAT.desired_vs(traf, 0) == pytest.approx(7.62)


def test_co_level_pair_at_cruise_is_not_falsely_free_to_revert():
    """Regression for the ap.vs-as-intent bug.

    Two aircraft level at FL300 on a converging track have NO vertical escape.
    Feeding ap.vs (7.62 m/s) in as a relative rate fabricates one and releases
    the resolution; the correct relative rate is zero, which holds it.
    """
    own = _StubTraf(alt=9144.0, apalt=9144.0, apvs=7.62, vs=0.0)
    intr = _StubTraf(alt=9144.0, apalt=9144.0, apvs=7.62, vs=0.0)

    dvs = FTRNAT.desired_vs(intr, 0) - FTRNAT.desired_vs(own, 0)
    assert dvs == 0.0
    assert _predicted(HEADON_DIST, HEADON_VREL, dalt=0.0, dvs=dvs) is True
    # The bug: using raw ap.vs would have produced a phantom escape.
    assert _predicted(HEADON_DIST, HEADON_VREL, dalt=0.0, dvs=-7.62) is False


# ---------------------------------------------------------------------------
# FTRNAT.revert_conflict — ramp-then-hold bracketing
# ---------------------------------------------------------------------------
class TestRevertConflict:
    """conflict_predicted extrapolates a constant vertical rate forever, which
    is only exact until an aircraft captures its commanded level. revert_conflict
    brackets the real ramp-then-hold profile with the two linear cases."""

    def _rc(self, dalt_now, dvs_now, dalt_settled, dvs_settled,
            dist=None, vrel=None):
        return bool(FTRNAT.revert_conflict(
            np.asarray(dist if dist is not None else HEADON_DIST, dtype=float),
            np.asarray(vrel if vrel is not None else HEADON_VREL, dtype=float),
            RPZ, HPZ, DTLOOK,
            dalt_now, dvs_now, dalt_settled, dvs_settled))

    def test_phantom_climb_is_caught_by_the_settled_case(self):
        # THE regression. 11 ft below cruise with the autopilot default rate:
        # the "now" case sees a 1500 fpm escape and clears, but the aircraft
        # settles co-level, so the pair is still in conflict.
        assert self._rc(dalt_now=-3.4, dvs_now=7.62,
                        dalt_settled=0.0, dvs_settled=0.0) is True

    def test_genuine_vertical_separation_clears_both_cases(self):
        # A real resolution: 1500 ft apart now, and both commanded levels stay
        # 1500 ft apart, so reverting keeps them separated.
        assert self._rc(dalt_now=1500.0 * FT, dvs_now=0.0,
                        dalt_settled=1500.0 * FT, dvs_settled=0.0) is False

    def test_settled_co_level_holds_even_when_currently_separated(self):
        # Currently 1500 ft apart but both commanded to the same level:
        # reverting converges them, so this must NOT release.
        assert self._rc(dalt_now=1500.0 * FT, dvs_now=0.0,
                        dalt_settled=0.0, dvs_settled=0.0) is True

    def test_transient_conflict_holds_even_when_settled_state_is_clear(self):
        # Co-level now, climbing to a level 1500 ft apart, but so slowly that
        # the pair is still within HPZ through the horizontal window.
        assert self._rc(dalt_now=0.0, dvs_now=1.0,
                        dalt_settled=1500.0 * FT, dvs_settled=0.0) is True

    def test_clear_when_neither_case_conflicts(self):
        assert self._rc(dalt_now=0.0, dvs_now=0.0,
                        dalt_settled=0.0, dvs_settled=0.0,
                        dist=[20000.0, 20000.0]) is False
