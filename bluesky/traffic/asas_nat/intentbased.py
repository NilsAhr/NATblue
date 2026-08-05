"""
intentbased.py — Intent-augmented conflict detection for NAT (EPP-style).
=========================================================================

``IntentBasedNAT`` (CDMETHOD ``INTENTBASEDNAT``) extends ``StateBasedNAT`` to
catch the conflict class state-based detection structurally misses on
fuel-optimal NAT traffic: **co-track cruise-climb vertical convergence** (both
aircraft climbing at ``vs≈0`` toward optima within HPZ → flagged only at LoS by
velocity-based CD, zero lead). See ``00_paperplanning/06_detection_limitation.md``.

Design (converged with user 2026-08-05)
---------------------------------------
* **UNION, state-based kept as the robust floor.** The full ``StateBasedNAT``
  result is returned unchanged; intent-detected pairs are *appended* (never
  removed) — safety cannot decrease.
* **Intent = broadcast-plan (EPP) motion projection.** For a candidate co-track
  pair, project BOTH aircraft's MOTION forward over t ∈ [0, ``T_FRAME_S``] from
  their route (``bs.traf.ap.route`` = the broadcast flight plan / EPP, unchanged
  by tactical ASAS): horizontal = move direct toward the next waypoint(s) at
  ground speed; vertical = climb at the PLANNED GRADIENT (current→next ``wpalt``)
  — instantaneous ``vs``≈0 misses it; perf-max over-predicts; the gradient is
  accurate and self-bounding (a far waypoint gives a shallow climb). Flag when
  horizontal < RPZ AND vertical < HPZ **at the same t** (co-location is temporal).
* **Time-bounded frame** (``T_FRAME_S``, not a waypoint count): a 2-h-away next
  waypoint falls outside the frame → no over-early trigger. This replaces the old
  ``tinhor < DTLOOK`` gate that fired late for same-track catch-up pairs.
* **False-positive guards:** co-track heading gate (|Δtrk| < ``INTENT_CT_ANGLE``),
  a proximity pre-filter, and **optima-within-HPZ** (if the two target FLs differ
  by ≥ HPZ the higher climbs clear → never flag).

Staleness while an aircraft is under active ASAS control is covered by the
state-based floor (actual dynamics); the intent layer reads the unchanged plan
and keeps the co-track hold engaged until horizontal clearance.

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-08-05
"""
import numpy as np

from bluesky.tools import geo
from bluesky.tools.aero import nm
from bluesky.traffic.asas_nat.statebased_nat import StateBasedNAT

_RE = 6371000.0   # earth radius [m] for the flat-earth (kwik) distance


class IntentBasedNAT(StateBasedNAT):
    """State-based CD (floor) + EPP-style route motion projection (co-track)."""

    # --- tunables (class constants) ---
    T_FRAME_S       = 1000.0   # [s]  intent look-ahead frame (~17 min; the
                               #      earliness/FP knob, matched to the ~30-40 min
                               #      cruise-climb convergence, time-bounded).
    INTENT_DT       = 15.0     # [s]  motion-projection sample step.
    INTENT_CT_ANGLE = 25.0     # [deg] co-track gate |Δtrk| (matches mvp2nat _CT_ANGLE_DEG).
    D_MAX_NM        = 130.0    # [NM] proximity pre-filter (> max closable in T_FRAME_S).
    INTENT_MAX_NODES = 400     # safety cap on route nodes scanned per aircraft.

    def detect(self, ownship, intruder, rpz, hpz, dtlookahead):
        # 1. State-based floor (unchanged) + reusable horizontal geometry.
        H = self._horizontal(ownship, intruder, rpz, dtlookahead)
        base = self._combine_vertical(ownship, intruder, hpz, dtlookahead, H)
        self._intent_pairs = set()          # frozensets flagged by intent (for logging)

        n = ownship.ntraf
        if n < 2:
            return base
        (confpairs, lospairs, inconf, tcpamax,
         b_qdr, b_dist, b_dcpa, b_tcpa, b_tinconf, b_mindalt) = base

        # 2. Candidate pre-filter: co-track + horizontally proximate + not already
        #    a state-based conflict. (Cheap boolean masks.)
        I = H['I']
        dist_now = H['dist']                                    # [m] (+1e9 on diag)
        hpz_m = np.asarray(np.maximum(np.asmatrix(hpz), np.asmatrix(hpz).transpose()))
        rpz_m = H['rpz']
        trk = np.asarray(ownship.trk, dtype=float)
        dtrk = np.abs(((trk.reshape((1, n)) - trk.reshape((n, 1)) + 180.0) % 360.0) - 180.0)
        cotrack = dtrk < self.INTENT_CT_ANGLE
        prox = dist_now < (self.D_MAX_NM * nm)
        cand = cotrack & prox & (I < 0.5)
        ii, jj = np.where(cand)
        if ii.size == 0:
            return base

        base_set = {frozenset(p) for p in confpairs}
        sched_cache = {}

        def sched(idx):
            if idx not in sched_cache:
                sched_cache[idx] = self._route_schedule(ownship, idx)
            return sched_cache[idx]

        add = []                                     # (i, j, tinconf, min_dalt, dcpa)
        seen = set()
        for i, j in zip(ii.tolist(), jj.tolist()):
            key = (i, j) if i < j else (j, i)
            if key in seen:
                continue
            seen.add(key)
            if frozenset((ownship.id[i], ownship.id[j])) in base_set:
                continue
            si, sj = sched(i), sched(j)
            if si is None or sj is None:
                continue
            # optima-within-HPZ guard: if the target FLs differ by >= HPZ the
            # higher climbs clear and they never vertically conflict → skip.
            if abs(si['target'] - sj['target']) >= float(hpz_m[i, j]):
                continue
            res = self._project_pair(si, sj, float(ownship.gs[i]), float(ownship.gs[j]),
                                     float(rpz_m[i, j]), float(hpz_m[i, j]))
            if res is None:
                continue
            add.append((i, j) + res)                 # + (tinconf, min_dalt, dcpa)

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

        for i, j, tinconf_ij, mindalt_ij, dcpa_ij in add:
            self._intent_pairs.add(frozenset((ownship.id[i], ownship.id[j])))
            g_dist = float(H['dist'][i, j])
            g_tcpa = float(H['tcpa'][i, j])
            for a, b in ((i, j), (j, i)):
                confpairs.append((ownship.id[a], ownship.id[b]))
                qdr.append(float(H['qdr'][a, b]))
                dist.append(g_dist)
                dcpa.append(dcpa_ij)                  # projected closest horizontal approach
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
    def _route_schedule(self, ownship, idx):
        """Broadcast-plan (EPP) schedule for aircraft ``idx`` as a function of
        cumulative along-route distance from its CURRENT position.

        Returns a dict of ascending-distance arrays:
          cum      [m]  position nodes (node 0 = current pos, then each waypoint)
          lat, lon [deg] at those nodes
          cum_alt  [m]  altitude nodes (node 0 = current alt, then waypoints with a
                        SPECIFIED wpalt) — the planned climb gradient is the linear
                        interpolation between them
          alt      [m]  at cum_alt
          target   [m]  the climb optimum = max planned altitude ahead
        or ``None`` when the route gives no usable plan (falls back to state-based).
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
        lat = [float(ownship.lat[idx])]
        lon = [float(ownship.lon[idx])]
        cum_alt = [0.0]
        alt = [float(ownship.alt[idx])]
        acc = 0.0
        for k in range(len(wlat)):
            acc += float(geo.kwikdist(lat[-1], lon[-1],
                                      float(wlat[k]), float(wlon[k]))) * nm
            cum.append(acc); lat.append(float(wlat[k])); lon.append(float(wlon[k]))
            wa = float(walt[k])
            if wa >= 0.0:
                cum_alt.append(acc); alt.append(wa)
        if len(cum) < 2:
            return None
        alt = np.asarray(alt)
        return {'cum': np.asarray(cum), 'lat': np.asarray(lat), 'lon': np.asarray(lon),
                'cum_alt': np.asarray(cum_alt), 'alt': alt, 'target': float(alt.max())}

    def _project_pair(self, si, sj, gs_i, gs_j, rpz_ij, hpz_ij):
        """Project both aircraft along their plans over [0, T_FRAME_S] and test
        same-time RPZ/HPZ. Returns ``(tinconf, min_dalt, dcpa)`` for the first
        joint breach, or ``None``.
        """
        ts = np.arange(0.0, self.T_FRAME_S + 1e-6, self.INTENT_DT)
        if ts.size == 0:
            return None
        di = max(gs_i, 0.0) * ts
        dj = max(gs_j, 0.0) * ts
        lat_i = np.interp(di, si['cum'], si['lat']); lon_i = np.interp(di, si['cum'], si['lon'])
        lat_j = np.interp(dj, sj['cum'], sj['lat']); lon_j = np.interp(dj, sj['cum'], sj['lon'])
        alt_i = np.interp(di, si['cum_alt'], si['alt'])
        alt_j = np.interp(dj, sj['cum_alt'], sj['alt'])
        # flat-earth (kwik) horizontal distance at each t [m]
        latm = np.radians(0.5 * (lat_i + lat_j))
        dx = _RE * np.radians(lon_j - lon_i) * np.cos(latm)
        dy = _RE * np.radians(lat_j - lat_i)
        hdist = np.hypot(dx, dy)
        vdist = np.abs(alt_i - alt_j)
        breach = (hdist < rpz_ij) & (vdist < hpz_ij)
        if not np.any(breach):
            return None
        k = int(np.argmax(breach))                   # first joint breach
        tinconf = float(ts[k])
        min_dalt = float(np.min(vdist[breach]))
        dcpa = float(np.min(hdist[breach]))
        return tinconf, min_dalt, dcpa
