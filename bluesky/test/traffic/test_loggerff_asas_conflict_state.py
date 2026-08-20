"""
Unit tests for the per-EVENT conflict bookkeeping in ``loggerff_asas``.

Background (paper-3, 2026-08-20). ``update()`` ended a conflict by deleting
only ``self.duration[pair]``. Everything else the event accumulates is written
under an ``if cpf not in self.<dict>`` guard, so a pair that left and re-entered
``confpairs`` logged its SECOND row carrying the FIRST event's ``tinconf``, an
inherited ``intrusion`` flag, and a ``dist`` minimised across both events.

That corrupts the study's primary metric three ways:
  * ``n_conflicts``  -- the postprocessing dedup key is (pair, tinconf), so
    re-detections sharing a stale tinconf collapse into one row.
  * ``n_intrusions`` -- the surviving row carries the inherited flag.
  * ``gap_to_prev_resolution_s = tinconf_k - toutconf_{k-1}`` -- goes negative
    for exactly the re-detections it exists to measure.

Like ``test_ftr.py`` and ``test_mvp2nat_guard.py``, these tests avoid the
session-scoped ``traffic_`` fixture: the two methods under test are pure
dictionary bookkeeping over ``self``, so a bare instance pins the semantics and
runs in milliseconds.

Distances are in metres and times in seconds, matching BlueSky internals.
"""
import pytest

from bluesky.plugins.loggerff_asas import LoggerffAsas


PAIR = frozenset(("AC001", "AC002"))
OTHER = frozenset(("AC003", "AC004"))


@pytest.fixture
def log():
    """A bare logger carrying only the conflict-state dictionaries.

    ``LoggerffAsas`` is a ``TrafficArrays`` Entity whose ``__init__`` reaches
    into ``TrafficArrays.root``, so it is constructed via ``__new__`` and the
    dictionaries the methods under test touch are supplied directly.
    """
    o = LoggerffAsas.__new__(LoggerffAsas)
    o.tinconf = {}
    o.toutconf = {}
    o.dist = {}
    o.tLOS = {}
    o.intrusion_occurred = {}
    o.dcpa = {}
    o.dalt = {}
    o.tcpa = {}
    o.qdr = {}
    return o


# ---------------------------------------------------------------------------
# Within a single event: the write-once / running-minimum semantics must hold
# ---------------------------------------------------------------------------
class TestWithinOneEvent:
    def test_tinconf_is_write_once(self, log):
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=100.0, in_los=False)
        log._note_conflict_state(PAIR, dist_now=8000.0, simt=110.0, in_los=False)
        assert log.tinconf[PAIR] == 100.0

    def test_dist_is_a_running_minimum(self, log):
        for d in (9000.0, 3000.0, 6000.0):
            log._note_conflict_state(PAIR, dist_now=d, simt=100.0, in_los=False)
        assert log.dist[PAIR] == 3000.0

    def test_intrusion_latches_true_once_seen(self, log):
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=100.0, in_los=False)
        log._note_conflict_state(PAIR, dist_now=1000.0, simt=110.0, in_los=True)
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=120.0, in_los=False)
        assert log.intrusion_occurred[PAIR] is True

    def test_tlos_records_the_first_los_tick(self, log):
        log._note_conflict_state(PAIR, dist_now=1000.0, simt=110.0, in_los=True)
        log._note_conflict_state(PAIR, dist_now=900.0, simt=120.0, in_los=True)
        assert log.tLOS[PAIR] == 110.0

    def test_conflict_without_los_records_intrusion_false(self, log):
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=100.0, in_los=False)
        assert log.intrusion_occurred[PAIR] is False
        assert PAIR not in log.tLOS


# ---------------------------------------------------------------------------
# Across events: THE BUG. Nothing may leak from event 1 into event 2.
# ---------------------------------------------------------------------------
class TestAcrossEvents:
    def test_second_event_gets_its_own_tinconf(self, log):
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=100.0, in_los=False)
        log._clear_conflict_state(PAIR)
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=500.0, in_los=False)
        assert log.tinconf[PAIR] == 500.0

    def test_intrusion_is_not_inherited(self, log):
        # event 1 goes to LoS, event 2 does not
        log._note_conflict_state(PAIR, dist_now=1000.0, simt=100.0, in_los=True)
        log._clear_conflict_state(PAIR)
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=500.0, in_los=False)
        assert log.intrusion_occurred[PAIR] is False
        assert PAIR not in log.tLOS

    def test_dist_minimum_does_not_cross_events(self, log):
        log._note_conflict_state(PAIR, dist_now=1000.0, simt=100.0, in_los=True)
        log._clear_conflict_state(PAIR)
        log._note_conflict_state(PAIR, dist_now=8000.0, simt=500.0, in_los=False)
        assert log.dist[PAIR] == 8000.0

    def test_gap_between_events_is_positive(self, log):
        """gap_to_prev_resolution_s must read as a real duration, not a
        negative number produced by a stale tinconf."""
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=100.0, in_los=False)
        log.toutconf[PAIR] = 200.0          # what update() records at event end
        log._clear_conflict_state(PAIR)
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=500.0, in_los=False)
        assert log.tinconf[PAIR] - log.toutconf[PAIR] == 300.0


# ---------------------------------------------------------------------------
# What clearing must NOT touch
# ---------------------------------------------------------------------------
class TestClearScope:
    def test_toutconf_survives(self, log):
        """The re-detection gap is measured against the PREVIOUS event's end,
        so toutconf deliberately persists across events."""
        log.toutconf[PAIR] = 200.0
        log._clear_conflict_state(PAIR)
        assert log.toutconf[PAIR] == 200.0

    def test_per_tick_severity_values_survive(self, log):
        """dcpa/dalt/tcpa/qdr are re-assigned every tick the pair is in
        confpairs, so they cannot go stale within an event and are left alone."""
        log.dcpa[PAIR], log.dalt[PAIR] = 500.0, 300.0
        log.tcpa[PAIR], log.qdr[PAIR] = -12.0, 90.0
        log._clear_conflict_state(PAIR)
        assert (log.dcpa[PAIR], log.dalt[PAIR]) == (500.0, 300.0)
        assert (log.tcpa[PAIR], log.qdr[PAIR]) == (-12.0, 90.0)

    def test_other_pairs_are_untouched(self, log):
        log._note_conflict_state(PAIR, dist_now=9000.0, simt=100.0, in_los=True)
        log._note_conflict_state(OTHER, dist_now=7000.0, simt=150.0, in_los=True)
        log._clear_conflict_state(PAIR)
        assert log.tinconf[OTHER] == 150.0
        assert log.dist[OTHER] == 7000.0
        assert log.intrusion_occurred[OTHER] is True
        assert log.tLOS[OTHER] == 150.0

    def test_clearing_an_unknown_pair_is_a_no_op(self, log):
        """update() clears on every conflict end; a pair that never reached
        the severity loop has no entries and must not raise."""
        log._clear_conflict_state(OTHER)
