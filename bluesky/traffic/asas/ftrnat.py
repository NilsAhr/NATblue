''' NAT Free-To-Revert: 3D revert criterion + waypoint snap-forward.

    FTR (Schaberg) tests only the HORIZONTAL closest point of approach. That is
    fine for a horizontally-resolving ASAS, but NATBlue's resolvers are
    vertical-first: MVP2NAT and every MVP2NAT_HYBRID3_* sweep variant run with
    swresovert=True / swresohoriz=False. Two aircraft on near-parallel NAT
    tracks resolved by climbing one still converge horizontally to well inside
    RPZ, so the plain FTR criterion is effectively never satisfied and ASAS
    would hold the resolution indefinitely.

    FTRNAT therefore replaces the horizontal-only test with the full 3D
    question that conflict DETECTION already asks: would the reverted state be
    predicted in conflict at all? A pair is free to revert when the horizontal
    and vertical violation windows do not overlap inside DTLOOK -- i.e. exactly
    the negation of bluesky.traffic.asas.statebased's conflict criterion.

    It also inherits PastCPANAT's recovery: on release the FMS is walked past
    any waypoints overflown during a sustained vertical manoeuvre, instead of
    the stock findact + direct which can leave the aircraft turning backwards.

    Selected with `RESNAV FTRNAT`, which MUST come after `RESO <method>` in a
    scenario -- RESO calls ResumeNavigation.select_default() and would
    otherwise overwrite the choice with PastCPANAT.
'''
import numpy as np

import bluesky as bs
from bluesky.tools.aero import ft
from bluesky.traffic.asas.ftr import FTR
from bluesky.traffic.asas.pastcpanat import PastCPANAT


class FTRNAT(FTR):
    ''' Free-To-Revert with a 3D release criterion and NAT waypoint recovery. '''

    # Reuse PastCPANAT's overflown-waypoint recovery verbatim rather than
    # re-deriving it, so the two NAT resume-nav methods cannot drift apart.
    _advance_past_overflown = staticmethod(PastCPANAT._advance_past_overflown)

    def resumenav(self, conf, ownship, intruder):
        '''
            Decide for each aircraft in the conflict list whether the ASAS
            should be followed or not, based on whether reverting to the
            autopilot state -- horizontally AND vertically -- would leave the
            pair out of conflict.
        '''
        # Record the assumed intent on the tick a pair becomes active, before
        # the pair is merged into resopairs, so it is the velocity the intruder
        # was flying when this conflict was first detected. Unlike FTR this
        # stores the vertical rate too, since the criterion needs it.
        for pair in conf.confpairs:
            if pair not in self.resopairs:
                idx2 = bs.traf.id2idx(pair[1])
                if idx2 >= 0:
                    self.assumed[pair] = (intruder.gseast[idx2],
                                          intruder.gsnorth[idx2],
                                          intruder.vs[idx2])

        # Add new conflicts to resopairs
        self.resopairs.update(conf.confpairs)

        # Conflict pairs to be deleted
        delpairs = set()
        changeactive = dict()

        # Look at all conflicts, also the ones that are solved but for which
        # reverting would not yet be safe
        for conflict in self.resopairs:
            idx1, idx2 = bs.traf.id2idx(conflict)
            # If the ownship aircraft is deleted remove its conflict from the list
            if idx1 < 0:
                delpairs.add(conflict)
                continue

            if idx2 >= 0:
                # Distance vector using flat earth approximation. Every relative
                # quantity below is intruder-minus-ownship, matching statebased.
                re = 6371000.
                dist = re * np.array([np.radians(intruder.lon[idx2] - ownship.lon[idx1]) *
                                      np.cos(0.5 * np.radians(intruder.lat[idx2] +
                                                              ownship.lat[idx1])),
                                      np.radians(intruder.lat[idx2] - ownship.lat[idx1])])

                rpz = np.max(conf.rpz[[idx1, idx2]])
                hpz = np.max(conf.hpz[[idx1, idx2]])
                dtlook = np.max(conf.dtlookahead[[idx1, idx2]])

                alt_own = ownship.alt[idx1]
                apalt_own = ownship.ap.alt[idx1]
                alt_intr = intruder.alt[idx2]
                dalt = alt_intr - alt_own

                # The ownship's own side of both criteria: what the autopilot
                # would fly the moment ASAS lets go, horizontally and vertically.
                vown = self.desired_spd_trk(ownship, idx1)
                vsown = self.desired_vs(ownship, idx1)

                # Criterion 1: the intruder holds its current, observed state
                # while the ownship reverts and settles at its commanded level.
                vintr = np.array([intruder.gseast[idx2], intruder.gsnorth[idx2]])
                ftr = not self.revert_conflict(
                    dist, vintr - vown, rpz, hpz, dtlook,
                    dalt, intruder.vs[idx2] - vsown,
                    alt_intr - apalt_own, intruder.vs[idx2])

                # Criterion 2: the intruder reverts to its desired state too
                if ftr and self.intent != 'OFF':
                    vrevert = vsrevert = None
                    if self.intent == 'ASSUMED':
                        assumed = self.assumed.get(conflict)
                        if assumed is not None:
                            vrevert = np.array(assumed[:2])
                            vsrevert = assumed[2]
                            # ASSUMED gives a velocity, not a target level, so
                            # the intruder keeps that rate once ownship settles.
                            dalt_settled = alt_intr - apalt_own
                            dvs_settled = vsrevert
                    else:
                        vrevert = self.desired_spd_trk(intruder, idx2)
                        vsrevert = self.desired_vs(intruder, idx2)
                        # DECLARED: both settle at their commanded levels.
                        dalt_settled = intruder.ap.alt[idx2] - apalt_own
                        dvs_settled = 0.0
                    if vrevert is not None:
                        ftr = not self.revert_conflict(
                            dist, vrevert - vown, rpz, hpz, dtlook,
                            dalt, vsrevert - vsown,
                            dalt_settled, dvs_settled)

            # Keep resolving for ownship while reverting would not clear the
            # intruder.
            if idx2 >= 0 and not ftr:
                # Enable ASAS for this aircraft
                changeactive[idx1] = True
            else:
                # Switch ASAS off for ownship if there are no other conflicts
                # that this aircraft is involved in.
                changeactive[idx1] = changeactive.get(idx1, False)
                # If reverting is safe, remove the pair from the resopairs list
                delpairs.add(conflict)

        for idx, active in changeactive.items():
            # Loop a second time: this is to avoid that ASAS resolution is
            # turned off for an aircraft that is involved simultaneously in
            # multiple conflicts, where the first, but not all conflicts are
            # resolved.
            bs.traf.cr.active[idx] = active
            if not active:
                # NAT recovery: snap the FMS to a waypoint ahead of the
                # aircraft, walking chains of overflown waypoints.
                self._advance_past_overflown(idx)

        # Remove pairs that are free to revert or have deleted aircraft
        self.resopairs -= delpairs
        for conflict in delpairs:
            self.assumed.pop(conflict, None)

    @classmethod
    def revert_conflict(cls, dist, vrel, rpz, hpz, dtlook,
                        dalt_now, dvs_now, dalt_settled, dvs_settled):
        '''
            Whether reverting leaves the pair in conflict, given that the
            reverted VERTICAL trajectory is a ramp followed by a hold.

            conflict_predicted extrapolates a constant vertical rate forever.
            That is exact while an aircraft is still capturing its commanded
            level, and badly wrong afterwards: an aircraft 11 ft below cruise
            reports the autopilot default of 1500 fpm, which extrapolates to
            1000 ft of separation in 40 s when in reality it levels off almost
            immediately. Believing that phantom escape releases every co-level
            pair on the spot.

            Rather than integrate the piecewise profile, bracket it with the two
            linear cases that are each exact over half of it and hold the
            resolution if either is in conflict:

              - "now"     : current levels closing at the current relative rate
                            -- exact until the aircraft capture their levels.
              - "settled" : where the pair ends up once the ownship has
                            captured its commanded level -- exact from then on.

            The result is conservative where the two overlap, which is the right
            direction of error for a release criterion.
        '''
        return (cls.conflict_predicted(dist, vrel, dalt_now, dvs_now,
                                       rpz, hpz, dtlook)
                or cls.conflict_predicted(dist, vrel, dalt_settled, dvs_settled,
                                          rpz, hpz, dtlook))

    @staticmethod
    def desired_vs(traf, idx):
        '''
            Signed vertical speed [m/s] the autopilot of aircraft idx would fly
            right now: climb or descend toward its commanded altitude at its
            default rate, or hold level if it is already there.

            Autopilot.vs is NOT a statement of intent -- it falls back to vsdef
            (1500 fpm, autopilot.py:101) whenever no vertical speed is selected,
            so a level cruising aircraft still reports 7.62 m/s. Taking that at
            face value fabricates a vertical escape for every co-level pair and
            releases the resolution almost immediately.

            This mirrors Traffic.update_airspeed (traffic.py:463-465):
                target_vs = swaltsel * sign(delta_alt) * |vs|
            using the same altitude-capture dead band, floored at the 10 ft the
            older fixed-band version used so the result does not depend on the
            current time step.
        '''
        delta_alt = float(traf.ap.alt[idx]) - float(traf.alt[idx])
        vsmag = abs(float(traf.ap.vs[idx]))
        # bs.sim is absent when this is exercised outside a running sim.
        simdt = abs(float(getattr(bs.sim, 'simdt', 0.0) or 0.0))
        band = max(10.0 * ft,
                   1.05 * max(simdt * vsmag, simdt * abs(float(traf.vs[idx]))))
        if abs(delta_alt) <= band:
            return 0.0
        return np.sign(delta_alt) * vsmag

    @staticmethod
    def conflict_predicted(dist, vrel, dalt, dvs, rpz, hpz, dtlook):
        '''
            Whether a pair in the given relative state is predicted in conflict.

            Single-pair mirror of the matrix algebra in
            bluesky/traffic/asas/statebased.py (StateBased.detect), so that
            "free to revert" is the exact negation of "would be detected as a
            conflict". All arguments are intruder-minus-ownship:

                dist  [m]   horizontal separation vector (east, north)
                vrel  [m/s] horizontal relative velocity (east, north)
                dalt  [m]   vertical separation
                dvs   [m/s] relative vertical speed
        '''
        # Horizontal ---------------------------------------------------------
        dv2 = float(np.dot(vrel, vrel))
        if abs(dv2) < 1e-6:
            dv2 = 1e-6  # limit lower absolute value, as StateBased does
        vabs = np.sqrt(dv2)
        dist2 = float(np.dot(dist, dist))

        tcpa = -float(np.dot(dist, vrel)) / dv2
        dcpa2 = abs(dist2 - tcpa * tcpa * dv2)

        R2 = rpz * rpz
        if not dcpa2 < R2:
            # The pair never comes within RPZ: no conflict, whatever it does
            # vertically.
            return False

        dtinhor = np.sqrt(max(0.0, R2 - dcpa2)) / vabs
        tinhor = tcpa - dtinhor
        touthor = tcpa + dtinhor

        # Vertical -----------------------------------------------------------
        # StateBased divides by -dvs and relies on the resulting +/-inf; do it
        # explicitly here so the degenerate case does not depend on the sign of
        # a floating-point zero.
        if abs(dvs) < 1e-6:
            if abs(dalt) >= hpz:
                return False        # never vertically overlapping
            tinver, toutver = -np.inf, np.inf   # overlapping for all time
        else:
            tcrosshi = (dalt + hpz) / -dvs
            tcrosslo = (dalt - hpz) / -dvs
            tinver = min(tcrosshi, tcrosslo)
            toutver = max(tcrosshi, tcrosslo)

        # Combine ------------------------------------------------------------
        tinconf = max(tinver, tinhor)
        toutconf = min(toutver, touthor)

        return bool(tinconf <= toutconf and toutconf > 0.0 and tinconf < dtlook)
