"""
Behavioural tests for `PastCPANAT._advance_past_overflown` -- the NAT recovery
that hands a released aircraft to a waypoint AHEAD of it instead of letting the
FMS turn it around.

Why this file exists. Every NAT release method calls it on every release
transition -- PastCPANAT, FTRNAT (aliased) and FTR2NAT -- and until now the
only coverage was an identity assertion that FTR2NAT reuses the same function
object. Nothing tested what it does.

The case that matters for the paper is the LATERALLY DISPLACED aircraft. In the
heading and combined domains the resolver pushes aircraft up to 80 NM off route
(measured on L05 of the synthetic suite), and the walk-forward decides whether
a waypoint is "behind" from the bearing to it. Bearings degrade in exactly that
regime: 80 NM abeam a waypoint a few NM ahead, the bearing is dominated by the
lateral component and sits within a couple of degrees of perpendicular -- right
on the > 90 deg threshold. A false "behind" verdict does not fail loudly. It
hands the aircraft further down the route, cutting the corner and silently
understating flown distance and fuel in precisely the arms that deviate most.

So the invariant under test is not "which waypoint" but "never past one that is
still ahead": `_ahead_of` computes the first waypoint the aircraft has not yet
passed, independently of the code under test, and no case may exceed it.

Geometry: five waypoints along 50 N, 0.25 deg of longitude apart (~9.6 NM),
the aircraft tracking east, so along-track order is simply increasing
longitude.

Units: `Traffic.cre` is the internal API and takes SI -- altitude in metres and
speed as CAS in m/s. Passing 35000/250 as if they were feet and knots gives an
aircraft at Mach 3.5, and `Route.findact`'s turn-time heuristic then advances
the active waypoint on its own, which looks exactly like a snap-forward bug and
is not one.
"""
import numpy as np
import pytest

import bluesky as bs
from bluesky.traffic.asas.pastcpanat import PastCPANAT
from bluesky.traffic.route import Route


LAT = 50.0
LONS = [-30.0, -29.75, -29.5, -29.25, -29.0]      # W1..W5, index 0..4
NAMES = ["W1", "W2", "W3", "W4", "W5"]
FT = 0.3048
CRUISE_ALT_M = 35000.0 * FT
CRUISE_CAS_MS = 133.0                              # ~Mach 0.80 at FL350
EIGHTY_NM = 80.0 / 60.0                            # ~1.333 deg of latitude


def _ahead_of(lon, lons=None):
    """First waypoint the aircraft has NOT yet passed, computed independently.

    The route runs due east along a parallel, so along-track order is just
    increasing longitude and this needs no projection.
    """
    lons = LONS if lons is None else lons
    for i, lon_wp in enumerate(lons):
        if lon_wp > lon:
            return i
    return len(lons) - 1


def _place(acid, lat, lon, trk=90.0, lons=None, names=None):
    """One aircraft with the five-waypoint route, at a chosen position.

    The route is built first and the position imposed afterwards, so the test
    controls exactly where the aircraft sits relative to its own plan -- which
    is what a release from a resolved conflict looks like.
    """
    lons = LONS if lons is None else lons
    names = NAMES if names is None else names
    bs.traf.cre(acid, "B77W", LAT, lons[0], trk, CRUISE_ALT_M, CRUISE_CAS_MS)
    idx = bs.traf.id2idx(acid)
    acrte = bs.traf.ap.route[idx]
    for name, lon_wp in zip(names, lons):
        acrte.addwpt_simple(idx, name, Route.wplatlon, LAT, lon_wp)
    bs.traf.lat[idx] = lat
    bs.traf.lon[idx] = lon
    bs.traf.trk[idx] = trk
    bs.traf.coslat[idx] = np.cos(np.radians(lat))
    return idx, acrte


@pytest.fixture(autouse=True)
def _clean(traffic_):
    bs.traf.reset()
    yield
    bs.traf.reset()


class TestAircraftOnItsRoute:
    """The no-op case: nothing has been overflown, so nothing may be skipped."""

    def test_a_waypoint_still_ahead_is_kept(self):
        lon = LONS[2] - 0.1                      # 0.1 deg short of W3
        idx, acrte = _place("ON1", LAT, lon)
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.wpname[acrte.iactwp] == "W3"
        assert acrte.iactwp <= _ahead_of(lon)

    def test_it_matches_plain_findact_plus_direct(self):
        # The first call mirrors the base sequence, so an aircraft that has not
        # overflown anything must land where the base method would put it.
        idx, acrte = _place("ON2", LAT, LONS[1] - 0.05)
        expected = acrte.wpname[acrte.findact(idx)]
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.wpname[acrte.iactwp] == expected


class TestOverflownWaypoint:
    """The case the function was written for."""

    def test_one_overflown_waypoint_is_walked_past(self):
        # Just past W3: the nearest waypoint is still W3, but it is behind.
        lon = LONS[2] + 0.05
        idx, acrte = _place("OF1", LAT, lon)
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.wpname[acrte.iactwp] == "W4"
        assert acrte.iactwp <= _ahead_of(lon)

    def test_it_never_walks_past_the_last_waypoint(self):
        # Past W5 entirely. There is nothing to advance to, and the guard has
        # to return rather than index off the end of the route.
        idx, acrte = _place("OF2", LAT, LONS[4] + 0.05)
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.iactwp == len(NAMES) - 1

    def test_an_empty_route_is_survivable(self):
        bs.traf.cre("OF3", "B77W", LAT, LONS[0], 90.0, CRUISE_ALT_M,
                    CRUISE_CAS_MS)
        idx = bs.traf.id2idx("OF3")
        PastCPANAT._advance_past_overflown(idx)      # findact returns -1
        assert bs.traf.ap.route[idx].nwp == 0


class TestLaterallyDisplacedAircraft:
    """80 NM off route -- the heading and combined domains.

    A displaced aircraft is the regime where the bearing test is weakest, and a
    false "behind" verdict here shortens the route rather than raising an
    error.
    """

    def test_a_waypoint_ahead_is_not_skipped_from_80_nm_abeam(self):
        lon = LONS[2] - 0.1                      # still short of W3
        idx, acrte = _place("LD1", LAT + EIGHTY_NM, lon)
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.wpname[acrte.iactwp] == "W3"
        assert acrte.iactwp <= _ahead_of(lon)

    def test_the_route_is_not_shortcut_from_80_nm_abeam(self):
        lon = LONS[2] + 0.05                     # displaced and past W3
        idx, acrte = _place("LD2", LAT + EIGHTY_NM, lon)
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.iactwp <= _ahead_of(lon)

    def test_displacement_alone_never_advances_the_active_waypoint(self):
        # Same along-track position, on route and 80 NM off it. Lateral
        # displacement is not evidence that a waypoint has been overflown, so
        # the two must agree.
        lon = LONS[1] + 0.05
        idx_on, rte_on = _place("LD3", LAT, lon)
        PastCPANAT._advance_past_overflown(idx_on)
        on_route = rte_on.wpname[rte_on.iactwp]

        idx_off, rte_off = _place("LD4", LAT + EIGHTY_NM, lon)
        PastCPANAT._advance_past_overflown(idx_off)
        assert rte_off.wpname[rte_off.iactwp] == on_route
        assert rte_off.iactwp <= _ahead_of(lon)


class TestRealNatWaypointSpacing:
    """The synthetic suite spaces waypoints 10 NM apart; real NAT plans space
    them by whole degrees of longitude, 75-200 NM per leg.

    That matters here because the "is this waypoint behind me" test compares
    the bearing to the waypoint against the bearing of the LEG THAT ENDS AT IT,
    measured from its start. Over a 200 NM great-circle leg at 50 N the leg
    bearing swings by about 3 deg end to end, so the reference the test uses is
    not the local track. The threshold has 90 deg of room, but the amount of
    room shrinks with leg length and it is worth pinning.
    """

    LONS_2DEG = [-30.0, -28.0, -26.0, -24.0, -22.0]      # ~77 NM legs
    LONS_5DEG = [-40.0, -35.0, -30.0, -25.0, -20.0]      # ~193 NM legs
    NAT_NAMES = ["N1", "N2", "N3", "N4", "N5"]

    @pytest.mark.parametrize("lons", [LONS_2DEG, LONS_5DEG])
    @pytest.mark.parametrize("offset_nm", [0.0, 80.0])
    def test_a_waypoint_ahead_is_never_skipped(self, lons, offset_nm):
        # A third of the way into the leg that ends at N3, i.e. genuinely
        # short of it, on route and 80 NM off it.
        lon = lons[1] + (lons[2] - lons[1]) / 3.0
        idx, acrte = _place("NS1", LAT + offset_nm / 60.0, lon,
                            lons=lons, names=self.NAT_NAMES)
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.iactwp <= _ahead_of(lon, lons)

    @pytest.mark.parametrize("lons", [LONS_2DEG, LONS_5DEG])
    @pytest.mark.parametrize("offset_nm", [0.0, 80.0])
    def test_an_overflown_waypoint_is_still_walked_past(self, lons, offset_nm):
        # Just past N3 on a long leg: the walk must still fire, otherwise the
        # aircraft is handed back to a waypoint behind it and turns around --
        # the whole reason the function exists.
        lon = lons[2] + (lons[3] - lons[2]) / 20.0
        idx, acrte = _place("NS2", LAT + offset_nm / 60.0, lon,
                            lons=lons, names=self.NAT_NAMES)
        PastCPANAT._advance_past_overflown(idx)
        assert acrte.wpname[acrte.iactwp] == "N4"
