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


# ===========================================================================
# BUG-02 (2026-08-24): conflicts still active at simulation end were lost.
#
# Measured on the synthetic suite: 9 of 175 runs logged ZERO conflict rows
# while the aircraft were in genuine loss of separation -- one of them for
# 4739 s. Every case was a pair still in conflict when the last aircraft was
# deleted. The failure is silent and points the wrong way: an unreleased pair
# is what a BADLY behaving resolver produces, so the worse the resolver, the
# emptier its CONFLOG, and an empty CONFLOG reads as "no conflicts".
#
# Three links, all confirmed in code:
#   1. `if traf.ntraf == 0: return` (update, ~line 521) sits BEFORE the
#      end-of-conflict logging block, so the tick that deletes the last
#      aircraft is the tick that stops conflicts ever being closed.
#   2. On the RESET path Simulation.reset() zeroes sim.simt BEFORE
#      bs.traf.reset() triggers the flush, so the flush stamps toutconf = 0.
#   3. CSVLogger.log is gated on `bs.sim.simt >= self.tlog`, and tlog is the
#      STARTLOG time -- so with simt == 0 every flushed row is discarded with
#      no error at all.
# ===========================================================================
class FakeConflog:
    """Stands in for the CSVLogger: records rows and honours the tlog gate.

    The gate is reproduced deliberately -- it is the mechanism that discarded
    the rows, so a fake without it could not catch the regression.
    """

    def __init__(self, tlog=0.0):
        self.rows = []
        self.tlog = tlog
        self.dt = 0.0
        self._open = True
        self.simt = 0.0

    def isopen(self):
        return self._open

    def log(self, *vals):
        if self._open and self.simt >= self.tlog:
            self.rows.append(vals)


@pytest.fixture
def flushlog(log):
    """A bare logger with the conflict dicts plus a fake conflog."""
    log.duration = {}
    log.init_lat1 = {}; log.init_lon1 = {}; log.init_alt1 = {}
    log.init_lat2 = {}; log.init_lon2 = {}; log.init_alt2 = {}
    log.init_hdg1 = {}; log.init_hdg2 = {}
    log.init_vs1 = {}; log.init_vs2 = {}
    log.init_vsreal1 = {}; log.init_vsreal2 = {}
    log.init_conflict_angle = {}
    log.init_mach1 = {}; log.init_mach2 = {}
    log.init_selalt1 = {}; log.init_selalt2 = {}
    log.init_tas1 = {}; log.init_tas2 = {}
    log.conflog = FakeConflog()
    log._last_simt = 0.0
    return log


class TestEndLogRunsWithoutTraffic:
    """The pair must be closed on the tick the traffic empties."""

    def test_a_pair_absent_from_confpairs_is_logged(self, flushlog):
        flushlog.duration[PAIR] = 42
        flushlog.tinconf[PAIR] = 100.0
        flushlog._last_simt = 500.0
        flushlog.conflog.simt = 500.0
        # Empty confpairs is exactly the no-traffic case.
        flushlog._log_ended_conflicts(confpairs_unique=set(), simt=500.0)
        assert len(flushlog.conflog.rows) == 1
        assert PAIR not in flushlog.duration
        assert flushlog.toutconf[PAIR] == 500.0

    def test_a_pair_still_in_conflict_is_not_logged(self, flushlog):
        flushlog.duration[PAIR] = 42
        flushlog._log_ended_conflicts(confpairs_unique={PAIR}, simt=500.0)
        assert flushlog.conflog.rows == []
        assert flushlog.duration[PAIR] == 43       # duration keeps counting

    def test_ending_a_pair_clears_its_event_state(self, flushlog):
        flushlog.duration[PAIR] = 1
        flushlog._note_conflict_state(PAIR, dist_now=1000.0, simt=100.0,
                                      in_los=True)
        flushlog._log_ended_conflicts(confpairs_unique=set(), simt=500.0)
        # Same teardown as the in-traffic path, so a re-detection starts clean.
        assert PAIR not in flushlog.tinconf
        assert PAIR not in flushlog.dist
        assert PAIR not in flushlog.intrusion_occurred

    def test_other_pairs_still_in_conflict_are_untouched(self, flushlog):
        flushlog.duration[PAIR] = 1
        flushlog.duration[OTHER] = 1
        flushlog._log_ended_conflicts(confpairs_unique={OTHER}, simt=500.0)
        assert len(flushlog.conflog.rows) == 1
        assert OTHER in flushlog.duration


class TestFlushCannotBeSilentlyDiscarded:
    """The reset-path flush must survive the periodic-logger gate."""

    def test_flush_writes_even_when_tlog_is_ahead_of_simt(self, flushlog):
        # THE REGRESSION. Simulation.reset() zeroes simt before the flush, and
        # tlog is the STARTLOG time (30 s in the paper3 batches), so the gate
        # `simt >= tlog` is False and the row vanishes without an error.
        flushlog.duration[PAIR] = 10
        flushlog._last_simt = 4800.0
        flushlog.conflog.tlog = 30.0
        flushlog.conflog.simt = 0.0                # sim.simt already zeroed
        flushlog._flush_active_conflicts()
        assert len(flushlog.conflog.rows) == 1

    def test_flush_restores_tlog(self, flushlog):
        flushlog.duration[PAIR] = 10
        flushlog.conflog.tlog = 30.0
        flushlog._flush_active_conflicts()
        assert flushlog.conflog.tlog == 30.0

    def test_flush_restores_tlog_even_if_the_write_raises(self, flushlog):
        def boom(*_):
            raise RuntimeError("disk full")
        flushlog.duration[PAIR] = 10
        flushlog.conflog.tlog = 30.0
        flushlog.conflog.log = boom
        with pytest.raises(RuntimeError):
            flushlog._flush_active_conflicts()
        assert flushlog.conflog.tlog == 30.0

    def test_flush_stamps_the_last_real_sim_time(self, flushlog):
        # Not sim.simt, which is 0.0 by the time the reset path gets here.
        flushlog.duration[PAIR] = 10
        flushlog.tinconf[PAIR] = 100.0
        flushlog._last_simt = 4800.0
        flushlog._flush_active_conflicts()
        assert flushlog.toutconf[PAIR] == 4800.0

    def test_second_flush_does_not_write_again(self, flushlog):
        flushlog.duration[PAIR] = 10
        flushlog._flush_active_conflicts()
        flushlog._flush_active_conflicts()
        assert len(flushlog.conflog.rows) == 1

    def test_flush_is_a_no_op_once_the_file_is_closed(self, flushlog):
        flushlog.duration[PAIR] = 10
        flushlog.conflog._open = False
        flushlog._flush_active_conflicts()
        assert flushlog.conflog.rows == []

    def test_flush_clears_event_state_like_the_end_log_path(self, flushlog):
        flushlog.duration[PAIR] = 10
        flushlog._note_conflict_state(PAIR, dist_now=1000.0, simt=100.0,
                                      in_los=True)
        flushlog._flush_active_conflicts()
        assert PAIR not in flushlog.tinconf
        assert PAIR not in flushlog.intrusion_occurred

    def test_a_stale_confpairs_set_would_have_left_the_bug_unfixed(self, flushlog):
        """Why update() passes an EMPTY set, not traf.cd.confpairs_unique.

        Traffic.update() returns at `if self.ntraf == 0` before reaching
        cd.update() (traffic.py:394,407), so on the no-traffic tick
        confpairs_unique still holds the pairs from the last tick that HAD
        traffic. Handing that stale set to _log_ended_conflicts closes nothing.
        """
        flushlog.duration[PAIR] = 42
        stale = {PAIR}                       # what cd still reports
        flushlog._log_ended_conflicts(stale, simt=500.0)
        assert flushlog.conflog.rows == []   # nothing closed -- the old bug
        flushlog._log_ended_conflicts(set(), simt=500.0)
        assert len(flushlog.conflog.rows) == 1
