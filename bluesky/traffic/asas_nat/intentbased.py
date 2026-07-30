"""
intentbased.py — Intent-augmented conflict detection for NAT.
=============================================================

``IntentBasedNAT`` (CDMETHOD ``INTENTBASEDNAT``) extends ``StateBasedNAT`` to
catch the one conflict class the state-based method structurally misses on
fuel-optimal NAT traffic: **co-track cruise-climb vertical convergence**.

Motivation
----------
State-based CD linearly extrapolates each aircraft's *instantaneous* vertical
speed. Fuel-optimal trajectories cruise-climb at ``vs≈0`` per step, so two
co-track aircraft at slightly different optimum flight levels are projected as
staying vertically separated — until the climb has already closed the gap and
it is instantly a loss of separation (``tinconf == tLOS``, zero lead). See
``00_paperplanning/06_detection_limitation.md``.

Approach (locked with user 2026-07-30)
--------------------------------------
* **Vertical-intent only.** Horizontal geometry stays state-based; the intent
  pass only re-derives the *vertical* conflict window from planned intent.
* **UNION (intent adds only).** The full state-based result is returned
  unchanged; intent-detected pairs are *appended*. Nothing state-based finds
  is removed.
* **Intent source = route ``wpalt`` (static, set at spawn).** The fast-changing
  signals (``selalt``/instantaneous ``vs``) are deliberately NOT used, so intent
  is a stable look-ahead. The route object is read live each tick, so genuine
  route edits (waypoint-skipping in recovery, ``DIRECT``) are reflected.
* **Horizontal reused straight-line** (``StateBasedNAT._horizontal``) — accurate
  for the near-linear co-track along-track closure.

Two-layer coverage of resolution-altered trajectories: an active resolution
manoeuvre produces a large ``vs`` / off-route track, which the state-based layer
catches from the actual state; the intent layer only adds the slow *planned*
``vs≈0`` climb convergence that state-based misses. Complementary regimes.

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-07-30
"""
import numpy as np

from bluesky.tools import geo
from bluesky.tools.aero import nm
from bluesky.traffic.asas_nat.statebased_nat import StateBasedNAT


class IntentBasedNAT(StateBasedNAT):
    """State-based CD + planned cruise-climb (route-``wpalt``) vertical intent."""

    # --- tunables (class constants; no stack commands in v1) ---
    INTENT_DT = 10.0          # [s]   vertical-projection sample step
    INTENT_DALT_BAND = 2.0    # [×HPZ] only examine pairs within this vertical band now
    INTENT_MAX_NODES = 400    # safety cap on route nodes scanned per aircraft

    def detect(self, ownship, intruder, rpz, hpz, dtlookahead):
        # 1. Full state-based result (unchanged) + reusable horizontal geometry.
        H = self._horizontal(ownship, intruder, rpz, dtlookahead)
        base = self._combine_vertical(ownship, intruder, hpz, dtlookahead, H)

        n = ownship.ntraf
        if n < 2:
            return base

        (confpairs, lospairs, inconf, tcpamax,
         b_qdr, b_dist, b_dcpa, b_tcpa, b_tinconf, b_mindalt) = base

        # 2. Candidate pre-filter: horizontal breach within lookahead, currently
        #    within a vertical band, NOT already a state-based conflict.
        I = H['I']
        swhorconf = H['swhorconf']
        tinhor, touthor = H['tinhor'], H['touthor']
        hpz_m = np.asarray(np.maximum(np.asmatrix(hpz), np.asmatrix(hpz).transpose()))
        dtl_col = np.asarray(dtlookahead, dtype=float).reshape((n, 1))

        dalt_now = (ownship.alt.reshape((1, n)) - intruder.alt.reshape((1, n)).T)
        band = np.abs(dalt_now) < (self.INTENT_DALT_BAND * hpz_m)

        cand = (swhorconf & (tinhor < dtl_col) & (touthor > 0.0)
                & band & (I < 0.5))
        ii, jj = np.where(cand)
        if ii.size == 0:
            return base

        base_set = {frozenset(p) for p in confpairs}
        dtl_arr = np.asarray(dtlookahead, dtype=float).reshape(-1)
        sched_cache = {}

        def schedule(idx):
            if idx not in sched_cache:
                sched_cache[idx] = self._intent_schedule(ownship, idx)
            return sched_cache[idx]

        add = []   # list of (i, j, tinconf, mindalt)
        seen = set()
        for i, j in zip(ii.tolist(), jj.tolist()):
            key = (i, j) if i < j else (j, i)
            if key in seen:
                continue
            seen.add(key)
            if frozenset((ownship.id[i], ownship.id[j])) in base_set:
                continue
            si, sj = schedule(i), schedule(j)
            if si is None or sj is None:
                continue                              # no vertical intent → state-based only
            dtl = float(min(dtl_arr[i], dtl_arr[j]))
            res = self._intent_vertical(
                si, sj, float(ownship.gs[i]), float(ownship.gs[j]),
                float(hpz_m[i, j]), float(tinhor[i, j]), float(touthor[i, j]), dtl)
            if res is None:
                continue
            add.append((i, j, res[0], res[1]))

        if not add:
            return base

        # 3. UNION: append both directions (state-based convention) + geometry.
        confpairs = list(confpairs)
        inconf = np.array(inconf, dtype=bool)
        qdr = list(np.atleast_1d(b_qdr).astype(float))
        dist = list(np.atleast_1d(b_dist).astype(float))
        dcpa = list(np.atleast_1d(b_dcpa).astype(float))
        tcpa = list(np.atleast_1d(b_tcpa).astype(float))
        tinc = list(np.atleast_1d(b_tinconf).astype(float))
        mind = list(np.atleast_1d(b_mindalt).astype(float))

        for i, j, tinconf_ij, mindalt_ij in add:
            g_dist = float(H['dist'][i, j])
            g_dcpa = float(np.sqrt(H['dcpa2'][i, j]))
            g_tcpa = float(H['tcpa'][i, j])
            for a, b in ((i, j), (j, i)):
                confpairs.append((ownship.id[a], ownship.id[b]))
                qdr.append(float(H['qdr'][a, b]))   # bearing is direction-dependent
                dist.append(g_dist)
                dcpa.append(g_dcpa)
                tcpa.append(g_tcpa)
                tinc.append(tinconf_ij)
                mind.append(mindalt_ij)
            inconf[i] = True
            inconf[j] = True

        return (confpairs, lospairs, inconf, tcpamax,
                np.asarray(qdr), np.asarray(dist), np.asarray(dcpa),
                np.asarray(tcpa), np.asarray(tinc), np.asarray(mind))

    # ------------------------------------------------------------------ #
    # Intent helpers
    # ------------------------------------------------------------------ #
    def _intent_schedule(self, ownship, idx):
        """Planned altitude-vs-along-route-distance for aircraft ``idx``.

        Node 0 = current position/altitude; subsequent nodes = route waypoints
        with a *specified* ``wpalt`` (>=0), at their cumulative along-route
        distance from the current position. Returns ``(cumdist_m, alt_m)`` (both
        ascending-distance arrays) or ``None`` when there is no usable vertical
        intent ahead (aircraft then falls back to state-based only).
        """
        ap = getattr(ownship, 'ap', None)
        if ap is None or not hasattr(ap, 'route'):
            return None
        try:
            rte = ap.route[idx]
        except (IndexError, TypeError):
            return None
        if rte is None or getattr(rte, 'nwp', 0) == 0:
            return None

        iact = getattr(rte, 'iactwp', 0)
        if iact is None or iact < 0:
            iact = 0
        wlat = rte.wplat[iact:iact + self.INTENT_MAX_NODES]
        wlon = rte.wplon[iact:iact + self.INTENT_MAX_NODES]
        walt = rte.wpalt[iact:iact + self.INTENT_MAX_NODES]
        if len(wlat) == 0:
            return None

        cum = [0.0]
        alts = [float(ownship.alt[idx])]
        acc = 0.0
        prev_lat, prev_lon = float(ownship.lat[idx]), float(ownship.lon[idx])
        for k in range(len(wlat)):
            acc += float(geo.kwikdist(prev_lat, prev_lon,
                                      float(wlat[k]), float(wlon[k]))) * nm
            prev_lat, prev_lon = float(wlat[k]), float(wlon[k])
            wa = float(walt[k])
            if wa >= 0.0:                 # only specified altitudes are intent nodes
                cum.append(acc)
                alts.append(wa)
        if len(cum) < 2:
            return None                   # flat → no convergence to predict
        return np.asarray(cum), np.asarray(alts)

    def _intent_vertical(self, si, sj, gs_i, gs_j, hpz_ij, tinhor, touthor, dtl):
        """Sample the planned vertical gap over [0, dtl] and combine with the
        (reused) horizontal window. Returns ``(tinconf, min_dalt)`` for the first
        joint breach, or ``None`` if none within the lookahead.
        """
        if dtl <= 0.0:
            return None
        cum_i, alt_i = si
        cum_j, alt_j = sj
        ts = np.arange(0.0, dtl + 1e-6, self.INTENT_DT)
        if ts.size == 0:
            return None
        ai = np.interp(max(gs_i, 0.0) * ts, cum_i, alt_i)   # np.interp clamps at ends
        aj = np.interp(max(gs_j, 0.0) * ts, cum_j, alt_j)
        dalt_t = np.abs(ai - aj)

        vert_breach = dalt_t < hpz_ij
        horiz_in = (ts >= tinhor) & (ts <= touthor)
        both = vert_breach & horiz_in
        if not np.any(both):
            return None
        tinconf = float(ts[int(np.argmax(both))])          # first joint-breach time
        window = horiz_in if np.any(horiz_in) else np.ones_like(ts, dtype=bool)
        min_dalt = float(np.min(dalt_t[window]))
        return tinconf, min_dalt
