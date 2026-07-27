''' Conflict resolution based on the Modified Voltage Potential algorithm.
    mvp2nat: starts as a clean copy of stock MVP. Step 1 (vs_stop) was removed —
    the near-ceiling climb-into-PZ case it targeted does not appear in the NAT
    scenario, and the trigger fired for normal cruise pairs and blocked MVP's
    legitimate vertical separation commands. '''
from collections import deque

import numpy as np
import bluesky as bs
from bluesky import stack
from bluesky.tools import geo
from bluesky.tools.aero import ft, fpm
from bluesky.traffic.asas import ConflictResolution


# Resolver safety clamps for degenerate geometries (in-trail, near-tangent).
# Without these, divisions by near-zero geometric factors emit numerical
# garbage that BADA perf.limits cannot always catch in a useful way.
_ERRATUM_FLOOR = np.sin(np.radians(15.0))  # ~0.259; caps (alpha-beta) at 75 deg
_MAX_DH_M = 2000.0 * ft                    # cap |asasalt - ownship.alt|

# ── Co-altitude vertical symmetry breaker (ported from mvp_nat.MVPNAT) ────────
# When two aircraft are at ~the same flight level, the relative vertical speed
# carries no usable sign, so the plain MVP formula hands BOTH aircraft the same
# dv3 sign (both descend) and never separates them. Below _VREL_VERT_FLOOR we
# instead assign climb/descend deterministically by aircraft index.
_VREL_VERT_FLOOR     = 0.1          # m/s (~20 fpm): below this, treat pair as co-altitude
_SYM_BREAK_TSOLV_MAX = 60.0         # s: cap on vertical solve-time in the degenerate case
_CEIL_MARGIN_M       = 500.0 * ft   # keep the aircraft assigned to climb this far below its ceiling

# ── Heading-fallback configuration ───────────────────────────────────────────
# H2 redesign (2026-07-27, Stage-C iter 1): the fallback commands the MVP
# HORIZONTAL resolution vector directly (RPZ-sized), gated only on co-altitude
# + still-approaching geometry — no tunable time/offset/class constants. All
# gated by self.swfallbackhdg; when False the fallback path is dead code.


class MVP2NAT(ConflictResolution):
    ''' Conflict resolution using the Modified Voltage Potential Method.
        Functionally identical to stock MVP; renamed to allow side-by-side
        comparison and to host NAT-specific safety clamps (tasactive in
        vertical-only mode, erratum floor, dh cap).

        Default resolution domain is VERTICAL (stock MVP defaults to
        horizontal). RMETHH2NAT / RMETHV2NAT stack commands silently
        fail under BlueSky's replaceable-Entity machinery: when RESO
        MVP2NAT is called from a scenario, no MVP2NAT instance yet
        exists. FuncObjects bind to the class, the command sets class
        attributes that are then shadowed by the instance __init__ once
        the first conflict materialises. Confirmed empirically on
        2026-05-12 — all four domain settings produced byte-identical
        results. For domain comparison studies use the dedicated
        derived classes below (MVP2NAT_VERT / MVP2NAT_SPD / MVP2NAT_HDG
        / MVP2NAT_BOTH) and select via RESO MVP2NAT_<DOMAIN>. '''
    def __init__(self):
        super().__init__()
        # NAT scenario default: vertical-only.
        self.swresohoriz = False
        self.swresospd = False
        self.swresohdg = False
        self.swresovert = True

        # Opt-in heading-fallback on top of vertical resolution. When False
        # the pure-vertical code path is byte-equivalent to legacy mvp2nat
        # (aside from the resumenav waypoint-snap fix). See MVP2NAT_HYBRID.
        self.swfallbackhdg = False

        # Heading-fallback state (H2). Rebuilt every resolve() tick:
        # _fallback_state[cs] = (commanded_track_deg, |dv_horiz|^2) for each
        # ownship being actively turned THIS tick. Keyed by callsign so aircraft
        # index compaction on traf.delete() cannot corrupt entries.
        self._fallback_state = {}

    def reset(self):
        super().reset()
        self._fallback_state.clear()

    def setprio(self, flag=None, priocode=''):
        '''Set the prio switch and the type of prio '''
        if flag is None:
            return True, "PRIORULES [ON/OFF] [PRIOCODE]" + \
                            "\nAvailable priority codes: " + \
                            "\n     FF1:  Free Flight Primary (No Prio) " + \
                            "\n     FF2:  Free Flight Secondary (Cruising has priority)" + \
                            "\n     FF3:  Free Flight Tertiary (Climbing/descending has priority)" + \
                            "\n     LAY1: Layers Primary (Cruising has priority + horizontal resolutions)" + \
                            "\n     LAY2: Layers Secondary (Climbing/descending has priority + horizontal resolutions)" + \
                            "\nPriority is currently " + ("ON" if self.swprio else "OFF") + \
                            "\nPriority code is currently: " + \
                str(self.priocode)
        options = ["FF1", "FF2", "FF3", "LAY1", "LAY2"]
        if priocode not in options:
            return False, "Priority code Not Understood. Available Options: " + str(options)
        return super().setprio(flag, priocode)

    @stack.command(name="RMETHH2NAT")
    def setresometh(self, value:'txt'=''):
        """ Processes the RMETHH command. Sets swresovert = False"""
        # Acceptable arguments for this command
        options = ["BOTH", "SPD", "HDG", "NONE", "ON", "OFF", "OF"]
        if not value:
            return True, "RMETHH2NAT [ON / BOTH / OFF / NONE / SPD / HDG]" + \
                         "\nHorizontal resolution limitation is currently " + ("ON" if self.swresohoriz else "OFF") + \
                         "\nSpeed resolution limitation is currently " + ("ON" if self.swresospd else "OFF") + \
                         "\nHeading resolution limitation is currently " + \
                ("ON" if self.swresohdg else "OFF")
        if value not in options:
            return False, "RMETH Not Understood" + "\nRMETHH2NAT [ON / BOTH / OFF / NONE / SPD / HDG]"
        else:
            if value == "ON" or value == "BOTH":
                self.swresohoriz = True
                self.swresospd = True
                self.swresohdg = True
                self.swresovert = False
            elif value == "OFF" or value == "OF" or value == "NONE":
                # Do NOT swtich off self.swresovert if value == OFF
                self.swresohoriz = False
                self.swresospd = False
                self.swresohdg = False
            elif value == "SPD":
                self.swresohoriz = True
                self.swresospd = True
                self.swresohdg = False
                self.swresovert = False
            elif value == "HDG":
                self.swresohoriz = True
                self.swresospd = False
                self.swresohdg = True
                self.swresovert = False

    @stack.command(name='RMETHV2NAT')
    def setresometv(self, value:'txt'=''):
        """ Processes the RMETHV command. Sets swresohoriz = False."""
        # Acceptable arguments for this command
        options = ["NONE", "ON", "OFF", "OF", "V/S"]
        if not value:
            return True, "RMETHV2NAT [ON / V/S / OFF / NONE]" + \
                         "\nVertical resolution limitation is currently " + \
                ("ON" if self.swresovert else "OFF")
        if value not in options:
            return False, f"RMETHV2NAT '{value}' Not Understood\nRMETHV2NAT [ON / V/S / OFF / NONE]"

        if value == "ON" or value == "V/S":
            self.swresovert = True
            self.swresohoriz = False
            self.swresospd = False
            self.swresohdg = False
        elif value == "OFF" or value == "OF" or value == "NONE":
            # Do NOT swtich off self.swresohoriz if value == OFF
            self.swresovert = False

    @property
    def tasactive(self):
        ''' In vertical-only mode the autopilot retains Mach-hold control.
            Without this override, tasactive inherits self.active (True for any
            conflicting aircraft), causing aporasas to command a fixed TAS during
            the climb — which raises Mach as altitude increases and produces MMO
            violations. Returning False here keeps the speed channel with the AP. '''
        if self.swresovert and not self.swresohoriz:
            return np.zeros(bs.traf.ntraf, dtype=bool)
        return self.active

    @property
    def hdgactive(self):
        ''' When swfallbackhdg is OFF this returns the base behaviour
            (self.active). When ON, heading is ASAS-controlled ONLY while
            the 60 s fallback hold is in progress for that aircraft — outside
            the hold window the FMS owns the heading channel, so the aircraft
            navigates to its next waypoint as soon as the offset expires. '''
        if not self.swfallbackhdg:
            return super().hdgactive
        arr = np.zeros(bs.traf.ntraf, dtype=bool)
        for cs in self._fallback_state:
            idx = bs.traf.id2idx(cs)
            if idx >= 0:
                arr[idx] = True
        return arr

    @property
    def fallback_hdg_active(self):
        ''' Per-aircraft boolean: True while a heading fallback turn is being
            commanded this tick. Always False when swfallbackhdg is OFF.
            Surfaced for loggerff_asas. '''
        arr = np.zeros(bs.traf.ntraf, dtype=bool)
        if not self.swfallbackhdg:
            return arr
        for cs in self._fallback_state:
            idx = bs.traf.id2idx(cs)
            if idx >= 0:
                arr[idx] = True
        return arr

    def applyprio(self, dv_mvp, dv1, dv2, vs1, vs2):
        ''' Apply the desired priority setting to the resolution '''

        # Primary Free Flight prio rules (no priority)
        if self.priocode == 'FF1':
            # since cooperative, the vertical resolution component can be halved, and then dv_mvp can be added
            dv_mvp[2] = dv_mvp[2] / 2.0
            dv1 = dv1 - dv_mvp
            dv2 = dv2 + dv_mvp

        # Secondary Free Flight (Cruising aircraft has priority, combined resolutions)
        if self.priocode == 'FF2':
            # since cooperative, the vertical resolution component can be halved, and then dv_mvp can be added
            dv_mvp[2] = dv_mvp[2]/2.0
            # If aircraft 1 is cruising, and aircraft 2 is climbing/descending -> aircraft 2 solves conflict
            if abs(vs1) < 0.1 and abs(vs2) > 0.1:
                dv2 = dv2 + dv_mvp
            # If aircraft 2 is cruising, and aircraft 1 is climbing -> aircraft 1 solves conflict
            elif abs(vs2) < 0.1 and abs(vs1) > 0.1:
                dv1 = dv1 - dv_mvp
            else:  # both are climbing/descending/cruising -> both aircraft solves the conflict
                dv1 = dv1 - dv_mvp
                dv2 = dv2 + dv_mvp

        # Tertiary Free Flight (Climbing/descending aircraft have priority and crusing solves with horizontal resolutions)
        elif self.priocode == 'FF3':
            # If aircraft 1 is cruising, and aircraft 2 is climbing/descending -> aircraft 1 solves conflict horizontally
            if abs(vs1) < 0.1 and abs(vs2) > 0.1:
                dv_mvp[2] = 0.0
                dv1 = dv1 - dv_mvp
            # If aircraft 2 is cruising, and aircraft 1 is climbing -> aircraft 2 solves conflict horizontally
            elif abs(vs2) < 0.1 and abs(vs1) > 0.1:
                dv_mvp[2] = 0.0
                dv2 = dv2 + dv_mvp
            else:  # both are climbing/descending/cruising -> both aircraft solves the conflict, combined
                dv_mvp[2] = dv_mvp[2]/2.0
                dv1 = dv1 - dv_mvp
                dv2 = dv2 + dv_mvp

        # Primary Layers (Cruising aircraft has priority and clmibing/descending solves. All conflicts solved horizontally)
        elif self.priocode == 'LAY1':
            dv_mvp[2] = 0.0
            # If aircraft 1 is cruising, and aircraft 2 is climbing/descending -> aircraft 2 solves conflict horizontally
            if abs(vs1) < 0.1 and abs(vs2) > 0.1:
                dv2 = dv2 + dv_mvp
            # If aircraft 2 is cruising, and aircraft 1 is climbing -> aircraft 1 solves conflict horizontally
            elif abs(vs2) < 0.1 and abs(vs1) > 0.1:
                dv1 = dv1 - dv_mvp
            else:  # both are climbing/descending/cruising -> both aircraft solves the conflict horizontally
                dv1 = dv1 - dv_mvp
                dv2 = dv2 + dv_mvp

        # Secondary Layers (Climbing/descending aircraft has priority and cruising solves. All conflicts solved horizontally)
        elif self.priocode == 'LAY2':
            dv_mvp[2] = 0.0
            # If aircraft 1 is cruising, and aircraft 2 is climbing/descending -> aircraft 1 solves conflict horizontally
            if abs(vs1) < 0.1 and abs(vs2) > 0.1:
                dv1 = dv1 - dv_mvp
            # If aircraft 2 is cruising, and aircraft 1 is climbing -> aircraft 2 solves conflict horizontally
            elif abs(vs2) < 0.1 and abs(vs1) > 0.1:
                dv2 = dv2 + dv_mvp
            else:  # both are climbing/descending/cruising -> both aircraft solves the conflic horizontally
                dv1 = dv1 - dv_mvp
                dv2 = dv2 + dv_mvp

        return dv1, dv2


    def resolve(self, conf, ownship, intruder):
        ''' Resolve all current conflicts '''
        # Initialize an array to store the resolution velocity vector for all A/C
        dv = np.zeros((ownship.ntraf, 3))

        # Initialize an array to store time needed to resolve vertically
        timesolveV = np.ones(ownship.ntraf) * 1e9

        fb_enabled = self.swfallbackhdg and self.swresovert
        # Heading fallback (H2): rebuild the per-tick map of ownship callsigns
        # being actively turned. An entry exists only while that aircraft is
        # co-altitude with and still approaching a conflict, so the MVP
        # horizontal turn is held through the approach and released past CPA.
        if fb_enabled:
            self._fallback_state = {}

        # Call MVP function to resolve conflicts-----------------------------------
        for ((ac1, ac2), qdr, dist, tcpa, tLOS) in zip(conf.confpairs, conf.qdr, conf.dist, conf.tcpa, conf.tLOS):
            idx1 = ownship.id.index(ac1)
            idx2 = intruder.id.index(ac2)

            # If A/C indexes are found, then apply MVP on this conflict pair
            # Because ADSB is ON, this is done for each aircraft separately
            if idx1 >-1 and idx2 > -1:
                dv_mvp, tsolV, dabsH = self.MVP(ownship, intruder, conf, qdr, dist, tcpa, tLOS, idx1, idx2)
                if tsolV < timesolveV[idx1]:
                    timesolveV[idx1] = tsolV

                # Use priority rules if activated
                if self.swprio:
                    dv[idx1], _ = self.applyprio(dv_mvp, dv[idx1], dv[idx2], ownship.vs[idx1], intruder.vs[idx2])
                else:
                    # since cooperative, the vertical resolution component can be halved, and then dv_mvp can be added
                    dv_mvp[2] = 0.5 * dv_mvp[2]
                    dv[idx1] = dv[idx1] - dv_mvp

                # Check the noreso aircraft. Nobody avoids noreso aircraft.
                # But noreso aircraft will avoid other aircraft
                if self.noresoac[idx2]:
                    dv[idx1] = dv[idx1] + dv_mvp

                # Check the resooff aircraft. These aircraft will not do resolutions.
                if self.resooffac[idx1]:
                    dv[idx1] = 0.0

                # Heading fallback (only enabled in vert + swfallbackhdg mode):
                # command the discarded MVP horizontal resolution vector for
                # co-altitude pairs vertical has not yet separated.
                if fb_enabled:
                    self._eval_horizontal_fallback(
                        ownship, intruder, ac1, ac2, idx1, idx2, dv_mvp, conf)

        # Determine new speed and limit resolution direction for all aicraft-------

        # Resolution vector for all aircraft, cartesian coordinates
        dv = np.transpose(dv)

        # The old speed vector, cartesian coordinates
        v = np.array([ownship.gseast, ownship.gsnorth, ownship.vs])

        # The new speed vector, cartesian coordinates
        newv = v + dv

        # Limit resolution direction if required-----------------------------------

        # [mvp2nat step 2] Generalised performance limits check.
        # The original implementation capped ground-speed against a
        # CAS-/TAS-based vmin/vmax window — a unit mismatch that, under
        # NAT tailwinds, capped GS against TAS limits and produced
        # spurious accelerations/decelerations. Here we instead route
        # the magnitude through ``bs.traf.perf.limits`` (BADA-/OpenAP-
        # aware) and treat it as TAS. The trade-off is that wind is no
        # longer subtracted on the way out, so the returned speed is
        # technically TAS — incompatible with downstream modules that
        # expect GS, but consistent with how perf.limits operates.
        # Multiple GS↔TAS conversions caused singularities at NAT-aligned
        # angles; doing it once here avoids them.
        if self.swresohoriz: # horizontal resolutions
            if self.swresospd and not self.swresohdg: # SPD only
                newtrack = ownship.trk
                newtas   = np.sqrt(newv[0,:]**2 + newv[1,:]**2)
                newvs    = ownship.vs
            elif self.swresohdg and not self.swresospd: # HDG only
                newtrack = (np.arctan2(newv[0,:],newv[1,:])*180/np.pi) % 360
                newtas   = ownship.tas
                newvs    = ownship.vs
            else: # SPD + HDG
                newtrack = (np.arctan2(newv[0,:],newv[1,:])*180/np.pi) %360
                newtas   = np.sqrt(newv[0,:]**2 + newv[1,:]**2)
                newvs    = ownship.vs
        elif self.swresovert: # vertical resolutions
            newtrack = ownship.trk
            newtas   = ownship.tas
            newvs    = newv[2,:]
        else: # horizontal + vertical
            newtrack = (np.arctan2(newv[0,:],newv[1,:])*180/np.pi) %360
            newtas   = np.sqrt(newv[0,:]**2 + newv[1,:]**2)
            newvs    = newv[2,:]

        # Determine ASAS module commands for all aircraft--------------------------

        # Single, BADA-/OpenAP-aware cap applied in TAS domain.
        # Returns: allowed TAS, capped VS, (third value unused).
        allowed_tas, vscapped, _ = ownship.perf.limits(
            newtas, newvs, ownship.alt, bs.traf.ax)

        # Calculate if Autopilot selected altitude should be followed. This avoids ASAS from
        # climbing or descending longer than it needs to if the autopilot leveloff
        # altitude also resolves the conflict. Because asasalttemp is calculated using
        # the time to resolve, it may result in climbing or descending more than the selected
        # altitude.
        # Cap |dh| at 5000 ft. In-trail same-altitude geometries can yield
        # tiny vrel[2] -> tsolV collapsed to tLOS, and tLOS * vsmax can push
        # the altitude target far beyond any sensible resolution. The cap
        # lets the resolver give up gracefully on geometries it cannot handle.
        dh = np.clip(vscapped * timesolveV, -_MAX_DH_M, _MAX_DH_M)
        asasalttemp = dh + ownship.alt
        signdvs = np.sign(vscapped - ownship.ap.vs * np.sign(ownship.selalt - ownship.alt))
        signalt = np.sign(asasalttemp - ownship.selalt)
        alt = np.where(np.logical_or(signdvs == 0, signdvs == signalt), asasalttemp, ownship.selalt)

        # To compute asas alt, timesolveV is used. timesolveV is a really big value (1e9)
        # when there is no conflict. Therefore asas alt is only updated when its
        # value is less than the look-ahead time, because for those aircraft are in conflict
        altCondition = np.logical_and(timesolveV<conf.dtlookahead, np.abs(dv[2,:])>0.0)
        alt[altCondition] = asasalttemp[altCondition]

        # If resolutions are limited in the horizontal direction, then asasalt should
        # be equal to auto pilot alt (aalt). This is to prevent a new asasalt being computed
        # using the auto pilot vertical speed (ownship.avs) using the code in line 106 (asasalttemp) when only
        # horizontal resolutions are allowed.
        alt = alt * (1 - self.swresohoriz) + ownship.selalt * self.swresohoriz

        # Heading-fallback override: aircraft turned this tick get their
        # commanded track replaced by the MVP horizontal resolution track. The
        # vertical channel is left untouched, so MVP keeps building vertical
        # separation in parallel.
        if fb_enabled and self._fallback_state:
            newtrack = np.array(newtrack, dtype=float, copy=True)
            for cs, (trk, _mag) in self._fallback_state.items():
                idx = bs.traf.id2idx(cs)
                if idx >= 0:
                    newtrack[idx] = trk

        return newtrack, allowed_tas, vscapped, alt

    def MVP(self, ownship, intruder, conf, qdr, dist, tcpa, tLOS, idx1, idx2):
        """Modified Voltage Potential (MVP) resolution method"""
        # Preliminary calculations-------------------------------------------------
        # Determine largest RPZ and HPZ of the conflict pair, use lookahead of ownship
        rpz_m = np.max(conf.rpz[[idx1, idx2]] * self.resofach)
        hpz_m = np.max(conf.hpz[[idx1, idx2]] * self.resofacv)
        dtlook = conf.dtlookahead[idx1]
        # Convert qdr from degrees to radians
        qdr = np.radians(qdr)

        # Relative position vector between id1 and id2
        drel = np.array([np.sin(qdr) * dist, \
                        np.cos(qdr) * dist, \
                        intruder.alt[idx2] - ownship.alt[idx1]])

        # Write velocities as vectors and find relative velocity vector
        v1 = np.array([ownship.gseast[idx1], ownship.gsnorth[idx1], ownship.vs[idx1]])
        v2 = np.array([intruder.gseast[idx2], intruder.gsnorth[idx2], intruder.vs[idx2]])
        vrel = v2 - v1


        # Horizontal resolution----------------------------------------------------

        # Find horizontal distance at the tcpa (min horizontal distance)
        dcpa  = drel + vrel*tcpa
        dabsH = np.sqrt(dcpa[0] * dcpa[0] + dcpa[1] * dcpa[1])

        # Compute horizontal intrusion
        iH = rpz_m - dabsH

        # Exception handlers for head-on conflicts
        # This is done to prevent division by zero in the next step
        if dabsH <= 10.:
            dabsH = 10.
            dcpa[0] = drel[1] / dist * dabsH
            dcpa[1] = -drel[0] / dist * dabsH

        # If intruder is outside the ownship PZ, then apply extra factor
        # to make sure that resolution does not graze IPZ
        if rpz_m < dist and dabsH < dist:
            # Compute the resolution velocity vector in horizontal direction.
            # abs(tcpa) because it bcomes negative during intrusion.
            erratum = np.cos(np.arcsin(rpz_m / dist)-np.arcsin(dabsH / dist))
            # Floor erratum at sin(15 deg). Near-tangent geometries
            # (rpz_m ~ dist with dabsH ~ 0, e.g. in-trail) drive
            # (alpha-beta) toward pi/2 and erratum toward 0, which makes
            # rpz_m/erratum blow up. The floor caps the angle at 75 deg
            # so the resolver backs out instead of emitting garbage.
            erratum = max(erratum, _ERRATUM_FLOOR)
            dv1 = ((rpz_m / erratum - dabsH) * dcpa[0]) / (abs(tcpa) * dabsH)
            dv2 = ((rpz_m / erratum - dabsH) * dcpa[1]) / (abs(tcpa) * dabsH)
        else:
            dv1 = (iH * dcpa[0]) / (abs(tcpa) * dabsH)
            dv2 = (iH * dcpa[1]) / (abs(tcpa) * dabsH)

        # Vertical resolution------------------------------------------------------
        # Two regimes, split on the relative vertical speed:
        #   |vrel[2]| >  _VREL_VERT_FLOOR : normal MVP — the direction comes from
        #                                   the sign of the closing vertical rate.
        #   |vrel[2]| <= _VREL_VERT_FLOOR : co-altitude degenerate case. The sign
        #                                   term -vrel[2]/|vrel[2]| is undefined,
        #                                   so plain MVP hands BOTH aircraft the
        #                                   same dv3 sign and they never separate.
        #                                   Use a deterministic symmetry breaker.
        if abs(vrel[2]) > _VREL_VERT_FLOOR:
            # --- Normal case: meaningful relative vertical speed ---
            iV    = hpz_m
            tsolV = abs(drel[2] / vrel[2])
            # If solving vertically would take longer than the look-ahead,
            # collapse to the horizontal loss-of-separation time instead.
            if tsolV > dtlook:
                tsolV = tLOS
            if tsolV <= 0.0:
                # Already at LOS vertically — the dh clamp downstream catches it.
                dv3 = 0.0
            else:
                dv3 = (iV / tsolV) * (-vrel[2] / abs(vrel[2]))
        else:
            # --- Co-altitude symmetry breaker (ported from mvp_nat.MVPNAT) ---
            # Assign climb/descend deterministically by aircraft index so the two
            # ordered evaluations of the pair, (A,B) and (B,A), produce OPPOSITE
            # signs. With the downstream accumulation dv[idx1] -= dv_mvp, the
            # effective ownship VS is -dv3, so:
            #   idx1 < idx2 -> dv3 > 0 -> ownship (lower index) DESCENDS
            #   idx1 > idx2 -> dv3 < 0 -> ownship (higher index) CLIMBS
            # tsolV is capped at _SYM_BREAK_TSOLV_MAX (not tLOS) so ~1000 ft of
            # separation can actually build before CPA; a tLOS-based reference is
            # far too weak (~100 fpm per aircraft).
            iV    = hpz_m
            tsolV = min(max(abs(tcpa), 1.0), _SYM_BREAK_TSOLV_MAX)
            magnitude = iV / tsolV
            dv3 = magnitude if idx1 < idx2 else -magnitude

            # Ceiling guard: the aircraft assigned to CLIMB (the higher index)
            # must not be driven above its service ceiling. If it is already
            # within _CEIL_MARGIN_M of the ceiling, flip the pair so the climber
            # descends and its partner climbs instead. Evaluated identically in
            # both orderings, so the flip stays consistent across (A,B)/(B,A).
            if idx1 < idx2:
                climber, climber_idx = intruder, idx2
            else:
                climber, climber_idx = ownship, idx1
            ceil = self._service_ceiling(climber, climber_idx)
            if climber.alt[climber_idx] > ceil - _CEIL_MARGIN_M:
                dv3 = -dv3

        # It is necessary to cap dv3 to prevent that a vertical conflict
        # is solved in 1 timestep, leading to a vertical separation that is too
        # high (high vs assumed in traf). If vertical dynamics are included to
        # aircraft  model in traffic.py, the below three lines should be deleted.
        #    mindv3 = -400*fpm# ~ 2.016 [m/s]
        #    maxdv3 = 400*fpm
        #    dv3 = np.maximum(mindv3,np.minimum(maxdv3,dv3))


        # Combine resolutions------------------------------------------------------

        # combine the dv components
        dv = np.array([dv1, dv2, dv3])

        # dabsH is returned so resolve() can feed it into the heading-fallback
        # dcpa-trend tracker. Stock MVP discards it after the local branch.
        return dv, tsolV, dabsH

    @staticmethod
    def _service_ceiling(ac, idx):
        ''' Service ceiling [m] for aircraft index idx, used by the co-altitude
            ceiling guard. Returns the larger of the certified max operating
            altitude (hmo) and the OPF max altitude (hmax): hmax alone is the
            MTOW-weight-limited ceiling and can sit below cruise, so taking the
            max of the two is the safe reference. Falls back to hmax if hmo is
            unavailable (e.g. a non-BADA performance model). '''
        hmax = ac.perf.hmax[idx]
        hmo = getattr(ac.perf, 'hmo', None)
        if hmo is not None:
            try:
                return max(float(hmax), float(hmo[idx]))
            except (TypeError, IndexError):
                pass
        return float(hmax)


    # ─────────────────────────────────────────────────────────────────────────
    # Recovery — snap past waypoints that were overflown during the maneuver
    # ─────────────────────────────────────────────────────────────────────────
    def resumenav(self, conf, ownship, intruder):
        ''' Mirror of the base implementation, with one addition: when an
            aircraft transitions out of ASAS-active state, we walk the FMS
            forward past any waypoints that lie behind the aircraft along
            the current route segment. This prevents a backward turn after
            a sustained vertical maneuver carries the aircraft past its
            "active" waypoint while the FMS still pointed at it.

            Constraint motivation: in conftest.scn every waypoint carries
            altitude and speed constraints (granular cruise climb), so a
            "do not skip constrained waypoints" rule would disable this fix
            entirely. We deliberately do skip — the next waypoint's
            constraints become active on direct() and cruise-climb
            continuity is preserved by the FMS. '''
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
                    hdist < rpz * self.resofach

            if idx2 >= 0 and (not past_cpa or hor_los or is_bouncing):
                changeactive[idx1] = True
            else:
                changeactive[idx1] = changeactive.get(idx1, False)
                delpairs.add(conflict)

        for idx, active in changeactive.items():
            self.active[idx] = active
            if not active:
                # Snap the FMS to a waypoint that is ahead of the aircraft
                # along the current route segment. _advance_past_overflown
                # walks chains of overflown waypoints; in the typical NAT case
                # zero or one waypoints are skipped.
                self._advance_past_overflown(idx)

        self.resopairs -= delpairs

    def _advance_past_overflown(self, idx):
        ''' Find the next FMS-active waypoint, then keep advancing while it
            lies behind the aircraft along the route segment (signed
            along-track distance < 0). The first call mirrors the base
            class's findact + direct sequence, so behaviour is identical
            for aircraft that have not overflown their target. '''
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
            # Bearing of the segment leading TO the active waypoint.
            # If no previous waypoint exists, fall back to aircraft track.
            if iact > 0:
                qdr_leg, _ = geo.qdrdist(
                    acrte.wplat[iact - 1], acrte.wplon[iact - 1],
                    acrte.wplat[iact],     acrte.wplon[iact])
            else:
                qdr_leg = bs.traf.trk[idx]
            # Bearing from aircraft to the active waypoint.
            qdr_to_wp, _ = geo.qdrdist(
                bs.traf.lat[idx], bs.traf.lon[idx],
                acrte.wplat[iact], acrte.wplon[iact])
            # Signed along-track sign: if |delta| > 90 deg, the waypoint is
            # behind us projected onto the segment direction.
            delta = abs(((qdr_to_wp - qdr_leg + 180.0) % 360.0) - 180.0)
            if delta > 90.0:
                acrte.direct(idx, acrte.wpname[iact + 1])
            else:
                return

    # ─────────────────────────────────────────────────────────────────────────
    # Heading fallback (H2, Stage-C iter 1). Dead code when swfallbackhdg=False.
    # ─────────────────────────────────────────────────────────────────────────
    def _eval_horizontal_fallback(self, ownship, intruder, ac1, ac2,
                                   idx1, idx2, dv_mvp, conf):
        ''' Command the ownship the MVP HORIZONTAL resolution vector
            (dv_mvp[:2]) — the RPZ-sized lateral avoidance MVP computes but
            discards in vertical-only mode — whenever the pair is (a) still
            co-altitude (|Δalt| < HPZ, so vertical has not yet separated) and
            (b) still approaching (not past CPA). No class gate, no time gate;
            re-evaluated every tick, so the turn is held through the whole
            approach and released once past CPA, and re-fires while closing.
            Runs in parallel with the vertical channel (only newtrack is
            overridden). Replaces the fire-once, 90 s-gated, class-gated ±5°
            offset (ticket 2026-07-16 backlog #2: that turn was <1 NM vs 5 NM
            RPZ, fired too late, and never triggered on head-on). '''
        # Respect noreso / resooff exactly as the vertical channel does.
        if self.resooffac[idx1] or self.noresoac[idx2]:
            return
        # (a) co-altitude gate: only step in while vertical has NOT built HPZ.
        hpz_m = np.max(conf.hpz[[idx1, idx2]])
        if abs(intruder.alt[idx2] - ownship.alt[idx1]) >= hpz_m:
            return
        # (b) still-approaching gate: horizontal relative-position · relative-
        # velocity > 0 means the pair is past CPA (separating) -> release.
        re = 6371000.0
        dx = re * np.radians(intruder.lon[idx2] - ownship.lon[idx1]) * \
            np.cos(0.5 * np.radians(intruder.lat[idx2] + ownship.lat[idx1]))
        dy = re * np.radians(intruder.lat[idx2] - ownship.lat[idx1])
        dvx = intruder.gseast[idx2] - ownship.gseast[idx1]
        dvy = intruder.gsnorth[idx2] - ownship.gsnorth[idx1]
        if (dx * dvx + dy * dvy) > 0.0:
            return
        # Command the MVP horizontal resolution. Same sign convention as the
        # vertical accumulation dv[idx1] -= dv_mvp in resolve(): ownship target
        # horizontal velocity = current - dv_mvp[:2].
        tgt_e = ownship.gseast[idx1] - dv_mvp[0]
        tgt_n = ownship.gsnorth[idx1] - dv_mvp[1]
        if tgt_e == 0.0 and tgt_n == 0.0:
            return
        new_trk = np.degrees(np.arctan2(tgt_e, tgt_n)) % 360.0
        # With multiple simultaneous conflicts, keep the largest lateral
        # resolution (biggest |dv_horiz|) for this ownship.
        mag = dv_mvp[0] * dv_mvp[0] + dv_mvp[1] * dv_mvp[1]
        prev = self._fallback_state.get(ac1)
        if prev is None or mag > prev[1]:
            self._fallback_state[ac1] = (new_trk, mag)


# ─────────────────────────────────────────────────────────────────────────────
# Domain-specific derived classes for resolver comparison studies.
# Each only differs from MVP2NAT in __init__ defaults. They exist because
# RMETHH2NAT / RMETHV2NAT silently fail (see MVP2NAT docstring), so domain
# selection must be done at class level. Activate via RESO MVP2NAT_<DOMAIN>.
# ─────────────────────────────────────────────────────────────────────────────
class MVP2NAT_VERT(MVP2NAT):
    ''' Vertical-only resolution (same as MVP2NAT default). '''
    def __init__(self):
        super().__init__()
        self.swresohoriz = False
        self.swresospd   = False
        self.swresohdg   = False
        self.swresovert  = True


class MVP2NAT_SPD(MVP2NAT):
    ''' Horizontal speed-only resolution. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz = True
        self.swresospd   = True
        self.swresohdg   = False
        self.swresovert  = False


class MVP2NAT_HDG(MVP2NAT):
    ''' Horizontal heading-only resolution. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz = True
        self.swresospd   = False
        self.swresohdg   = True
        self.swresovert  = False


class MVP2NAT_BOTH(MVP2NAT):
    ''' Horizontal resolution using both speed and heading. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz = True
        self.swresospd   = True
        self.swresohdg   = True
        self.swresovert  = False


class MVP2NAT_HYBRID(MVP2NAT):
    ''' Vertical-primary resolution + opt-in ±5° heading-fallback after a
        90 s ASAS-active period if dcpa is still decreasing, |dalt| < HPZ,
        and the geometry is class A (vertical convergence) or class D
        (near in-trail). Targets the residual LoS clusters that pure
        vertical leaves behind. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz   = False
        self.swresospd     = False
        self.swresohdg     = False
        self.swresovert    = True
        self.swfallbackhdg = True
