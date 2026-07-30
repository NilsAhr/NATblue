"""
statebased_nat.py — State-based conflict detection for NAT.
============================================================

Full copy of ``bluesky.traffic.asas.statebased.StateBased.detect()``
so that the detection algorithm can be modified independently of the
upstream version.

Changes vs. original
--------------------
2026-07-30: ``detect()`` refactored into two reusable helpers WITHOUT any
change to its numerical output/behaviour:
  * ``_horizontal(...)``       -> horizontal CPA block (qdr/dist/tcpa/dcpa2,
                                  swhorconf, tinhor/touthor, rpz matrix, I)
  * ``_combine_vertical(...)`` -> vertical crossing + combine + conflict lists
This lets ``IntentBasedNAT`` reuse the horizontal geometry unchanged and swap
ONLY the vertical window source (see ``intentbased.py``). ``StateBasedNAT``
itself is bit-for-bit equivalent to the pre-refactor version.

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-04 (refactor 2026-07)
"""
import numpy as np
from bluesky import stack
from bluesky.tools import geo
from bluesky.tools.aero import nm
from bluesky.traffic.asas_nat.detection_nat import ConflictDetectionNAT


class StateBasedNAT(ConflictDetectionNAT):
    """State-based conflict detection — NAT variant.

    Algorithm (all-pairs O(n²) matrix):
      1. Horizontal: QDR/dist → relative velocity → tCPA → dCPA² → RPZ test
      2. Vertical:   dalt → relative VS → crossing times through ±HPZ
      3. Combined:   tin_conf = max(tin_hor, tin_ver),
                     tout_conf = min(tout_hor, tout_ver)
      4. dalt_min:   minimum vertical distance during [tin, tout]
    """

    def detect(self, ownship, intruder, rpz, hpz, dtlookahead):
        """Conflict detection between *ownship* and *intruder*."""
        H = self._horizontal(ownship, intruder, rpz, dtlookahead)
        return self._combine_vertical(ownship, intruder, hpz, dtlookahead, H)

    # ------------------------------------------------------------------ #
    # Horizontal CPA block (reused unchanged by IntentBasedNAT)
    # ------------------------------------------------------------------ #
    def _horizontal(self, ownship, intruder, rpz, dtlookahead):
        """Horizontal closest-point-of-approach geometry for all pairs.

        Returns a dict of the full O(n²) matrices needed by the vertical
        combine step (and by intent-based detection):
        ``I, qdr, dist, tcpa, dcpa2, vrel, swhorconf, tinhor, touthor, rpz``.
        ``dist`` is in metres with 1e9 added on the diagonal; ``rpz`` is the
        pairwise-max RPZ matrix.
        """
        # Identity matrix: avoid ownship–ownship pairs
        I = np.eye(ownship.ntraf)

        qdr, dist = geo.kwikqdrdist_matrix(
            np.asmatrix(ownship.lat), np.asmatrix(ownship.lon),
            np.asmatrix(intruder.lat), np.asmatrix(intruder.lon))

        qdr = np.asarray(qdr)
        dist = np.asarray(dist) * nm + 1e9 * I

        qdrrad = np.radians(qdr)
        dx = dist * np.sin(qdrrad)
        dy = dist * np.cos(qdrrad)

        owntrkrad = np.radians(ownship.trk)
        ownu = ownship.gs * np.sin(owntrkrad).reshape((1, ownship.ntraf))
        ownv = ownship.gs * np.cos(owntrkrad).reshape((1, ownship.ntraf))

        inttrkrad = np.radians(intruder.trk)
        intu = intruder.gs * np.sin(inttrkrad).reshape((1, ownship.ntraf))
        intv = intruder.gs * np.cos(inttrkrad).reshape((1, ownship.ntraf))

        du = ownu - intu.T
        dv = ownv - intv.T

        dv2 = du * du + dv * dv
        dv2 = np.where(np.abs(dv2) < 1e-6, 1e-6, dv2)
        vrel = np.sqrt(dv2)

        tcpa = -(du * dx + dv * dy) / dv2 + 1e9 * I

        dcpa2 = np.abs(dist * dist - tcpa * tcpa * dv2)

        rpz = np.asarray(np.maximum(np.asmatrix(rpz), np.asmatrix(rpz).transpose()))
        R2 = rpz * rpz
        swhorconf = dcpa2 < R2

        dxinhor = np.sqrt(np.maximum(0., R2 - dcpa2))
        dtinhor = dxinhor / vrel

        tinhor = np.where(swhorconf, tcpa - dtinhor, 1e8)
        touthor = np.where(swhorconf, tcpa + dtinhor, -1e8)

        return dict(I=I, qdr=qdr, dist=dist, tcpa=tcpa, dcpa2=dcpa2, vrel=vrel,
                    swhorconf=swhorconf, tinhor=tinhor, touthor=touthor, rpz=rpz)

    # ------------------------------------------------------------------ #
    # Vertical crossing + combine + conflict lists
    # ------------------------------------------------------------------ #
    def _combine_vertical(self, ownship, intruder, hpz, dtlookahead, H):
        """State-based vertical window + combine with horizontal window ``H``.

        Identical maths to the original monolithic ``detect()``.
        """
        I = H['I']
        tinhor, touthor = H['tinhor'], H['touthor']
        swhorconf = H['swhorconf']
        tcpa, dcpa2 = H['tcpa'], H['dcpa2']
        qdr, dist, rpz = H['qdr'], H['dist'], H['rpz']

        # === Vertical conflict ====================================================
        dalt = (ownship.alt.reshape((1, ownship.ntraf))
                - intruder.alt.reshape((1, ownship.ntraf)).T + 1e9 * I)

        dvs = (ownship.vs.reshape(1, ownship.ntraf)
               - intruder.vs.reshape(1, ownship.ntraf).T)

        hpz = np.asarray(np.maximum(np.asmatrix(hpz), np.asmatrix(hpz).transpose()))

        with np.errstate(divide='ignore', invalid='ignore'):
            tcrosshi = (dalt + hpz) / -dvs
            tcrosslo = (dalt - hpz) / -dvs

        tinver = np.minimum(tcrosshi, tcrosslo)
        toutver = np.maximum(tcrosshi, tcrosslo)

        # === Combined =============================================================
        tinconf = np.maximum(tinver, tinhor)
        toutconf = np.minimum(toutver, touthor)

        swconfl = np.array(
            swhorconf * (tinconf <= toutconf) * (toutconf > 0.0)
            * np.asarray(tinconf < np.asmatrix(dtlookahead).T) * (1.0 - I),
            dtype=bool)

        # --- Conflict lists -------------------------------------------------------
        inconf = np.any(swconfl, 1)
        tcpamax = np.max(tcpa * swconfl, 1)

        confpairs = [(ownship.id[i], ownship.id[j])
                     for i, j in zip(*np.where(swconfl))]
        swlos = (dist < rpz) * (np.abs(dalt) < hpz)
        lospairs = [(ownship.id[i], ownship.id[j])
                    for i, j in zip(*np.where(swlos))]

        # --- Minimum vertical distance during conflict window ---------------------
        with np.errstate(invalid='ignore', over='ignore'):
            dalt_entry = dalt + dvs * tinconf
            dalt_exit  = dalt + dvs * toutconf

        with np.errstate(divide='ignore', invalid='ignore'):
            t_cross = -dalt / dvs
        cross_in_window = (t_cross >= tinconf) & (t_cross <= toutconf)

        min_dalt = np.where(
            cross_in_window, 0.0,
            np.minimum(np.abs(dalt_entry), np.abs(dalt_exit)))

        return (confpairs, lospairs, inconf, tcpamax,
                qdr[swconfl], dist[swconfl], np.sqrt(dcpa2[swconfl]),
                tcpa[swconfl], tinconf[swconfl], min_dalt[swconfl])
