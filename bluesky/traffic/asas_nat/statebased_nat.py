"""
statebased_nat.py — State-based conflict detection for NAT.
============================================================

Full copy of ``bluesky.traffic.asas.statebased.StateBased.detect()``
so that the detection algorithm can be modified independently of the
upstream version.

Changes vs. original
--------------------
(none yet — exact mirror of NATBlue ``statebased.py`` as of 2026-04-08)

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-04
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
        # Identity matrix: avoid ownship–ownship pairs
        I = np.eye(ownship.ntraf)

        # === Horizontal conflict ==================================================

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
