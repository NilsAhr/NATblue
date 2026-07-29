''' NAT resume navigation: PastCPA + waypoint snap-forward.

    Identical past-CPA / horizontal-LOS / bouncing release logic as PastCPA, but
    when an aircraft transitions OUT of ASAS-active state it walks the FMS forward
    past any waypoints that lie behind it along the current route segment. This
    prevents a backward turn after a sustained vertical manoeuvre carries the
    aircraft past its "active" waypoint while the FMS still points at it.

    This is the NATBlue default resume-nav (see ResumeNavigation.select_default);
    it reproduces the behaviour previously implemented as MVP2NAT.resumenav.
'''
import numpy as np

import bluesky as bs
from bluesky.tools import geo
from bluesky.traffic.asas import PastCPA


class PastCPANAT(PastCPA):
    ''' PastCPA with NAT overflown-waypoint snap-forward recovery. '''
    def resumenav(self, conf, ownship, intruder):
        self.resopairs.update(conf.confpairs)

        delpairs = set()
        changeactive = dict()

        def anglediff(a, b):
            d = a - b
            if d > 180:
                return anglediff(a, b + 360)
            elif d < -180:
                return anglediff(a + 360, b)
            return d

        for conflict in self.resopairs:
            idx1, idx2 = bs.traf.id2idx(conflict)
            if idx1 < 0:
                delpairs.add(conflict)
                continue

            if idx2 >= 0:
                re = 6371000.
                dist = re * np.array([np.radians(intruder.lon[idx2] - ownship.lon[idx1]) *
                                      np.cos(0.5 * np.radians(intruder.lat[idx2] +
                                                              ownship.lat[idx1])),
                                      np.radians(intruder.lat[idx2] - ownship.lat[idx1])])

                vrel = np.array([intruder.gseast[idx2] - ownship.gseast[idx1],
                                 intruder.gsnorth[idx2] - ownship.gsnorth[idx1]])

                past_cpa = np.dot(dist, vrel) > 0.0

                rpz = np.max(conf.rpz[[idx1, idx2]])
                hdist = np.linalg.norm(dist)
                hor_los = hdist < rpz

                is_bouncing = \
                    abs(anglediff(ownship.trk[idx1], intruder.trk[idx2])) < 30.0 and \
                    hdist < rpz * bs.traf.cr.resofach

            if idx2 >= 0 and (not past_cpa or hor_los or is_bouncing):
                changeactive[idx1] = True
            else:
                changeactive[idx1] = changeactive.get(idx1, False)
                delpairs.add(conflict)

        for idx, active in changeactive.items():
            bs.traf.cr.active[idx] = active
            if not active:
                # NAT recovery: snap the FMS to a waypoint ahead of the aircraft,
                # walking chains of overflown waypoints (typically 0-1 skipped).
                self._advance_past_overflown(idx)

        self.resopairs -= delpairs

    @staticmethod
    def _advance_past_overflown(idx):
        ''' Find the next FMS-active waypoint, then keep advancing while it lies
            behind the aircraft along the route segment (signed along-track
            distance < 0). The first call mirrors the base findact + direct
            sequence, so behaviour is identical for aircraft that have not
            overflown their target. '''
        acrte = bs.traf.ap.route[idx]
        iwpid = acrte.findact(idx)
        if iwpid == -1:
            return
        acrte.direct(idx, acrte.wpname[iwpid])

        # Walk forward through any chain of overflown waypoints. Cap at the
        # route length so a malformed route can't loop indefinitely.
        for _ in range(acrte.nwp):
            iact = acrte.iactwp
            if iact < 0 or iact >= acrte.nwp - 1:
                return  # no next waypoint to advance to
            if iact > 0:
                qdr_leg, _ = geo.qdrdist(
                    acrte.wplat[iact - 1], acrte.wplon[iact - 1],
                    acrte.wplat[iact],     acrte.wplon[iact])
            else:
                qdr_leg = bs.traf.trk[idx]
            qdr_to_wp, _ = geo.qdrdist(
                bs.traf.lat[idx], bs.traf.lon[idx],
                acrte.wplat[iact], acrte.wplon[iact])
            delta = abs(((qdr_to_wp - qdr_leg + 180.0) % 360.0) - 180.0)
            if delta > 90.0:
                acrte.direct(idx, acrte.wpname[iact + 1])
            else:
                return
