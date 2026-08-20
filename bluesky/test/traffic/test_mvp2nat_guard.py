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

    def test_vs_within_cap_is_untouched(self):
        vs = np.array([1000.0 * FPM, -1200.0 * FPM])
        alt = np.array([10000.0, 10000.0])
        anchor = np.array([10000.0, 10000.0])
        out_vs, _ = MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.allclose(out_vs, vs)

    def test_runaway_vs_is_capped_both_signs(self):
        # the observed failure mode was a ~44 000 fpm commanded dive
        vs = np.array([44000.0 * FPM, -44000.0 * FPM])
        alt = np.array([10000.0, 10000.0])
        anchor = np.array([10000.0, 10000.0])
        out_vs, _ = MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.allclose(out_vs, [_MAX_RESO_VS, -_MAX_RESO_VS])

    def test_altitude_bounded_below_anchor(self):
        # a 20 000 ft dive command against a FL380 anchor
        anchor = np.array([38000.0 * FT])
        alt = np.array([18000.0 * FT])
        _, out_alt = MVP2NAT._clamp_vertical(np.array([0.0]), alt, anchor)
        assert np.isclose(out_alt[0], anchor[0] - _MAX_VNAV_DEV_M)

    def test_altitude_bounded_above_anchor(self):
        anchor = np.array([38000.0 * FT])
        alt = np.array([58000.0 * FT])
        _, out_alt = MVP2NAT._clamp_vertical(np.array([0.0]), alt, anchor)
        assert np.isclose(out_alt[0], anchor[0] + _MAX_VNAV_DEV_M)

    def test_altitude_inside_the_band_is_untouched(self):
        anchor = np.array([38000.0 * FT])
        alt = np.array([36000.0 * FT])
        _, out_alt = MVP2NAT._clamp_vertical(np.array([0.0]), alt, anchor)
        assert np.isclose(out_alt[0], alt[0])

    def test_nan_anchor_leaves_altitude_unclamped(self):
        # no anchor yet (aircraft not under ASAS vertical control) -> the
        # altitude passes through, but the VS cap still applies
        anchor = np.array([np.nan])
        alt = np.array([18000.0 * FT])
        vs = np.array([44000.0 * FPM])
        out_vs, out_alt = MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.isclose(out_alt[0], alt[0])
        assert np.isclose(out_vs[0], _MAX_RESO_VS)

    def test_inputs_are_not_modified_in_place(self):
        vs = np.array([44000.0 * FPM])
        alt = np.array([18000.0 * FT])
        anchor = np.array([38000.0 * FT])
        MVP2NAT._clamp_vertical(vs, alt, anchor)
        assert np.isclose(vs[0], 44000.0 * FPM)
        assert np.isclose(alt[0], 18000.0 * FT)

    def test_guard_does_not_bind_the_hybrid3_hold(self):
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
