"""
statebased_trace.py -- why did this pair leave and re-enter ``confpairs``?

Diagnostic only. ``CDMETHOD STATEBASEDREALVSTRACE``.

The question it answers. On 20250714, 53 % (vert) to 91 % (spd) of all CONFLOG
events are the same pair logged again, and in 94-98 % of those ``asas_active``
covered the WHOLE gap -- the resolver held station while the detector dropped
the pair out of ``confpairs`` and re-flagged it, median gap 1-20 s. A CONFLOG
row is written when a pair leaves ``confpairs_unique``, so the repeats are
exactly those departures. This class records, per tick and per pair, which term
of the detection test flipped.

    swconfl = swhorconf
              AND (tinconf <= toutconf)
              AND (toutconf > 0)
              AND (tinconf < dtlookahead)

with ``tinconf = max(tinver, tinhor)`` and ``toutconf = min(toutver, touthor)``
(``statebased_nat._combine_vertical``).

SAFETY PROPERTY, and the reason this is trustworthy: the returned conflict lists
come from ``super()._combine_vertical`` unchanged. The trace recomputes the
intermediates for logging only, and ASSERTS that its own reconstruction implies
exactly the pair list the parent returned. It evaluates more; it decides
identically -- the same property that made ``FTR2NAT_TRACE`` safe.

Output path from the ``CD_TRACE`` environment variable; unset, this behaves as
``StateBasedRealVS``. Intended for the 2-6 aircraft synthetic scenarios: it
writes one row per ordered pair per tick, so ``CD_TRACE_MAX_NTRAF`` (default 12)
refuses to run it on full-day traffic.
"""
import os

import numpy as np

import bluesky as bs

from bluesky.traffic.asas_nat.statebased_nat import StateBasedRealVS


class StateBasedRealVSTrace(StateBasedRealVS):
    """StateBasedRealVS plus a per-tick, per-pair record of the detection test."""

    COLUMNS = ("simt,ac1,ac2,in_confpairs,swhorconf,tinhor,touthor,"
               "tinver,toutver,tinconf,toutconf,dtlook,dist_m,dcpa_m,tcpa,"
               "dalt_m,dvs_ms,rpz_m,hpz_m,fail_term")

    def __init__(self):
        super().__init__()
        self._trace_path = os.environ.get("CD_TRACE", "")
        self._max_ntraf = int(os.environ.get("CD_TRACE_MAX_NTRAF", "12"))
        self._fh = None

    def reset(self):
        super().reset()
        self._close()

    def _close(self):
        if self._fh is not None:
            try:
                self._fh.close()
            finally:
                self._fh = None

    def _handle(self):
        if self._fh is None and self._trace_path:
            self._fh = open(self._trace_path, "a", encoding="utf-8")
            if self._fh.tell() == 0:
                self._fh.write(self.COLUMNS + "\n")
        return self._fh

    def _combine_vertical(self, ownship, intruder, hpz, dtlookahead, H,
                          vs_vec=None):
        out = super()._combine_vertical(ownship, intruder, hpz, dtlookahead, H,
                                        vs_vec=vs_vec)
        fh = self._handle()
        if fh is None or ownship.ntraf > self._max_ntraf:
            return out
        try:
            self._write(ownship, intruder, hpz, dtlookahead, H, vs_vec,
                        out)
        except Exception as exc:                      # never break a sim run
            print(f"[CD_TRACE] disabled after error: {exc}")
            self._trace_path = ""
            self._close()
        return out

    def _write(self, ownship, intruder, hpz, dtlookahead, H, vs_vec, out):
        """Reconstruct the intermediates and log them.

        Mirrors `_combine_vertical` exactly; the assertion below is what keeps
        the mirror honest.
        """
        n = ownship.ntraf
        I = H['I']
        tinhor, touthor, swhorconf = H['tinhor'], H['touthor'], H['swhorconf']
        tcpa, dcpa2, dist = H['tcpa'], H['dcpa2'], H['dist']
        rpz = H['rpz']

        dalt = (ownship.alt.reshape((1, n))
                - intruder.alt.reshape((1, n)).T + 1e9 * I)
        vs_own = ownship.vs if vs_vec is None else vs_vec
        vs_int = intruder.vs if vs_vec is None else vs_vec
        dvs = (np.asarray(vs_own).reshape(1, n)
               - np.asarray(vs_int).reshape(1, n).T)
        hpz_m = np.asarray(np.maximum(np.asmatrix(hpz), np.asmatrix(hpz).T))

        with np.errstate(divide='ignore', invalid='ignore'):
            tcrosshi = (dalt + hpz_m) / -dvs
            tcrosslo = (dalt - hpz_m) / -dvs
        tinver = np.minimum(tcrosshi, tcrosslo)
        toutver = np.maximum(tcrosshi, tcrosslo)
        tinconf = np.maximum(tinver, tinhor)
        toutconf = np.minimum(toutver, touthor)

        c1 = np.asarray(swhorconf, dtype=bool)
        c2 = tinconf <= toutconf
        c3 = toutconf > 0.0
        c4 = np.asarray(tinconf < np.asmatrix(dtlookahead).T)
        swconfl = np.array(c1 * c2 * c3 * c4 * (1.0 - I), dtype=bool)

        # The mirror must imply exactly the parent's pair list, or the trace is
        # describing a detector other than the one that ran.
        mine = {(ownship.id[i], ownship.id[j])
                for i, j in zip(*np.where(swconfl))}
        assert mine == set(out[0]), "CD trace reconstruction diverged"

        simt = float(getattr(bs.sim, "simt", 0.0))
        rows = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if not c1[i, j]:
                    fail = "swhorconf"
                elif not c2[i, j]:
                    fail = "tinconf>toutconf"
                elif not c3[i, j]:
                    fail = "toutconf<=0"
                elif not c4[i, j]:
                    fail = "tinconf>=dtlook"
                else:
                    fail = "none"
                rows.append(
                    "%.2f,%s,%s,%d,%d,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.0f,"
                    "%.1f,%.1f,%.1f,%.1f,%.4f,%.1f,%.1f,%s"
                    % (simt, ownship.id[i], ownship.id[j],
                       int(swconfl[i, j]), int(c1[i, j]),
                       tinhor[i, j], touthor[i, j],
                       tinver[i, j], toutver[i, j],
                       tinconf[i, j], toutconf[i, j],
                       float(np.asarray(dtlookahead).ravel()[i]),
                       dist[i, j], np.sqrt(max(dcpa2[i, j], 0.0)), tcpa[i, j],
                       dalt[i, j], dvs[i, j], rpz[i, j], hpz_m[i, j], fail))
        self._fh.write("\n".join(rows) + "\n")
        self._fh.flush()
