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

# ── Co-track speed reduction (Stage-C iter 2) ────────────────────────────────
# The residual worst-day LoS are adjacent-FL co-track pairs (angle<2 deg, vs=0,
# ~1000 ft = HPZ apart): overtakes (Mach differs) + same-speed parallel locks.
# No vertical/heading command separates them (already at the HPZ boundary;
# sin-floor horizontal). Slowing ONE aircraft opens ALONG-TRACK separation,
# clearing both. MMO-safe (reducing Mach). Gated by self.swspeedcotrack.
_CT_ANGLE_DEG   = 25.0              # max |trk1-trk2| to count as co-track/parallel
_CT_VS_FLOOR    = 100.0 * fpm       # both aircraft |vs| below this = level [m/s]
_CT_DALT_MAX    = 2000.0 * ft       # within 2 FL (catches the 1000 ft adjacent-FL) [m]
_CT_SPD_FACTOR  = 0.90             # slow the chosen aircraft to this fraction of its TAS
_CT_VSREL_MAX   = 500.0 * fpm      # max |vs1-vs2| to count as co-track (not a fast crossing) [m/s]

# ── Co-track FL-hold (iter 3) ────────────────────────────────────────────────
# The TRAILING aircraft of a sustained co-track pair is held at a fixed FL below
# the leader AND slowed, until the pair is horizontally clear (> RPZ x release
# factor). Stable target (no per-tick decrement) => no runaway dive.
_HOLD_HPZ_MARGIN   = 1.10           # hold the trailing aircraft this x HPZ below the leader
_HOLD_VS_MAX       = 1500.0 * fpm   # cap the VS used to reach the hold target [m/s]
_HOLD_RELEASE_FACT = 1.10           # release when horizontal dist > this x RPZ


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

        # Co-track speed reduction (iter 2). When True, populate _speed_state
        # each tick with {cs: target_TAS} for co-track pairs to slow; applied to
        # newtas + tasactive. Off by default (base/domain/HYBRID classes).
        self.swspeedcotrack = False
        self._speed_state = {}

        # Co-track FL-hold + speed (iter 3). When swaltholdcotrack is True, a
        # co-track (sustained same-track cruise-climb competition) pair puts its
        # TRAILING aircraft into a PERSISTENT hold: fly a fixed lower FL AND slow
        # down, held until the pair is horizontally clear (> RPZ). Persistent (NOT
        # rebuilt per tick) so it survives the conflict clearing vertically -> no
        # rebound. _hold[trailing_cs] = {'partner','alt','tas'}.
        self.swaltholdcotrack = False
        self._hold = {}

    def reset(self):
        super().reset()
        self._fallback_state.clear()
        self._speed_state.clear()
        self._hold.clear()

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
            violations. Returning False here keeps the speed channel with the AP,
            EXCEPT for aircraft the co-track speed reduction is actively slowing
            this tick (a REDUCTION, so MMO-safe). '''
        if self.swresovert and not self.swresohoriz:
            arr = np.zeros(bs.traf.ntraf, dtype=bool)
            if self.swspeedcotrack:
                for cs in self._speed_state:
                    idx = bs.traf.id2idx(cs)
                    if idx >= 0:
                        arr[idx] = True
            if self.swaltholdcotrack:
                for cs in self._hold:            # trailing aircraft is slowed
                    idx = bs.traf.id2idx(cs)
                    if idx >= 0:
                        arr[idx] = True
            return arr
        return self.active

    @property
    def altactive(self):
        ''' Base returns self.active (alt controlled while in conflict). For the
            co-track FL-hold we ALSO force it True for held (trailing) aircraft
            even after the conflict clears vertically, so the AP cannot climb them
            back into the conflict (rebound) before they are horizontally clear. '''
        base = self.active
        if not self.swaltholdcotrack or not self._hold:
            return base
        arr = np.array(base, dtype=bool, copy=True)
        for cs in self._hold:
            idx = bs.traf.id2idx(cs)
            if idx >= 0:
                arr[idx] = True
        return arr

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
        if self.swspeedcotrack:
            self._speed_state = {}

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

                # Co-track speed reduction (iter 2): slow one aircraft of a
                # co-track co-altitude pair to open along-track separation.
                if self.swspeedcotrack:
                    self._eval_speed_cotrack(
                        ownship, intruder, ac1, ac2, idx1, idx2, conf)

                # Co-track FL-hold + speed (iter 3): put the TRAILING aircraft of a
                # sustained co-track pair into a persistent hold (fixed lower FL +
                # slow) so both eventually reach optimal FL without LoS/runaway.
                if self.swaltholdcotrack:
                    self._eval_althold_cotrack(
                        ownship, intruder, ac1, ac2, idx1, idx2, conf)

        # Maintain the persistent co-track hold: refresh targets, release pairs
        # that are now horizontally clear (> RPZ). Runs for ALL held aircraft,
        # including any whose conflict has dropped out of confpairs (rebound guard).
        if self.swaltholdcotrack:
            self._maintain_hold(ownship, conf)

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

        # Co-track speed override: replace the commanded TAS of aircraft being
        # slowed this tick with the reduced target (tasactive is True for them).
        if self.swspeedcotrack and self._speed_state:
            newtas = np.array(newtas, dtype=float, copy=True)
            for cs, ttas in self._speed_state.items():
                idx = bs.traf.id2idx(cs)
                if idx >= 0:
                    newtas[idx] = ttas

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

        # Co-track FL-hold override (iter 3): held (trailing) aircraft fly a FIXED
        # target altitude (bounded VS toward it, so no runaway) and a reduced TAS,
        # held until horizontally clear. altactive/tasactive are forced True for them.
        if self.swaltholdcotrack and self._hold:
            alt = np.array(alt, dtype=float, copy=True)
            allowed_tas = np.array(allowed_tas, dtype=float, copy=True)
            vscapped = np.array(vscapped, dtype=float, copy=True)
            for cs, h in self._hold.items():
                idx = bs.traf.id2idx(cs)
                if idx < 0:
                    continue
                alt[idx] = h['alt']
                allowed_tas[idx] = h['tas']
                # bounded VS toward the fixed target (reach over ~60 s, capped)
                vscapped[idx] = np.clip((h['alt'] - ownship.alt[idx]) / 60.0,
                                        -_HOLD_VS_MAX, _HOLD_VS_MAX)

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


    # Recovery / resume-navigation moved out of the resolver (2026-07-29): now a
    # pluggable ResumeNavigation Entity. NAT waypoint snap-forward lives in
    # bluesky/traffic/asas/pastcpanat.py (PastCPANAT), auto-selected via RESO.

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
        # Don't ALSO turn an aircraft already under co-track FL-hold+speed (the
        # turn interferes with the along-track speed separation — iter-2 finding).
        if self.swaltholdcotrack and (ac1 in self._hold or ac2 in self._hold):
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

    # ─────────────────────────────────────────────────────────────────────────
    # Co-track speed reduction (iter 2). Dead code when swspeedcotrack=False.
    # ─────────────────────────────────────────────────────────────────────────
    def _eval_speed_cotrack(self, ownship, intruder, ac1, ac2, idx1, idx2, conf):
        ''' For a co-track (near-parallel, level, within ~2 FL) pair that is
            still approaching, slow ONE aircraft (the faster; lower callsign if
            ~equal) to _CT_SPD_FACTOR of its TAS. This opens ALONG-TRACK
            separation, which is the only thing that clears adjacent-FL co-track
            conflicts — vertical is stuck at the HPZ boundary and horizontal is
            sin-floor-degenerate at ~0 deg. A REDUCTION, so it cannot breach MMO.
            Re-evaluated each tick (held while approaching, released past CPA). '''
        if self.resooffac[idx1] or self.noresoac[idx2]:
            return
        # level: both aircraft near-zero vertical speed
        if abs(ownship.vs[idx1]) >= _CT_VS_FLOOR or abs(intruder.vs[idx2]) >= _CT_VS_FLOOR:
            return
        # within ~2 FL (co-altitude/adjacent-FL, incl. the 1000 ft residual)
        if abs(intruder.alt[idx2] - ownship.alt[idx1]) >= _CT_DALT_MAX:
            return
        # co-track: relative track angle small (parallel / in-trail)
        rel_trk = abs(((ownship.trk[idx1] - intruder.trk[idx2] + 180.0) % 360.0) - 180.0)
        if rel_trk >= _CT_ANGLE_DEG:
            return
        # still approaching (hold until past CPA)
        re = 6371000.0
        dx = re * np.radians(intruder.lon[idx2] - ownship.lon[idx1]) * \
            np.cos(0.5 * np.radians(intruder.lat[idx2] + ownship.lat[idx1]))
        dy = re * np.radians(intruder.lat[idx2] - ownship.lat[idx1])
        dvx = intruder.gseast[idx2] - ownship.gseast[idx1]
        dvy = intruder.gsnorth[idx2] - ownship.gsnorth[idx1]
        if (dx * dvx + dy * dvy) > 0.0:
            return
        # Choose which aircraft to slow: the faster; lower callsign if ~equal.
        # Deterministic across both (A,B)/(B,A) orderings of the pair.
        t1, t2 = float(ownship.tas[idx1]), float(intruder.tas[idx2])
        if abs(t1 - t2) < 1.0:                      # ~equal speed
            slow_cs, slow_tas = (ac1, t1) if str(ac1) < str(ac2) else (ac2, t2)
        else:
            slow_cs, slow_tas = (ac1, t1) if t1 > t2 else (ac2, t2)
        self._speed_state[slow_cs] = _CT_SPD_FACTOR * slow_tas

    # ─────────────────────────────────────────────────────────────────────────
    # Co-track FL-hold + speed (iter 3). Dead code when swaltholdcotrack=False.
    # ─────────────────────────────────────────────────────────────────────────
    def _eval_althold_cotrack(self, ownship, intruder, ac1, ac2, idx1, idx2, conf):
        ''' Detect a sustained co-track pair and put its TRAILING aircraft into a
            PERSISTENT hold (fixed lower FL + reduced TAS). Idempotent: does
            nothing if either aircraft is already held. Trailing = the one behind
            along-track (the other is ahead of it in its own velocity direction). '''
        if self.resooffac[idx1] or self.noresoac[idx2]:
            return
        if ac1 in self._hold or ac2 in self._hold:
            return                                  # already handled
        # NOT a fast vertical crossing (that is the tactical regime) — but DO
        # allow a slowly cruise-climbing trailing aircraft (the residual regime
        # IS a cruise-climb competition, so requiring both level missed it).
        if abs(intruder.vs[idx2] - ownship.vs[idx1]) >= _CT_VSREL_MAX:
            return
        # co-altitude / adjacent-FL and near-parallel
        if abs(intruder.alt[idx2] - ownship.alt[idx1]) >= _CT_DALT_MAX:
            return
        rel_trk = abs(((ownship.trk[idx1] - intruder.trk[idx2] + 180.0) % 360.0) - 180.0)
        if rel_trk >= _CT_ANGLE_DEG:
            return
        # along-track ordering: is the intruder ahead of the ownship?
        re = 6371000.0
        dx = re * np.radians(intruder.lon[idx2] - ownship.lon[idx1]) * \
            np.cos(0.5 * np.radians(intruder.lat[idx2] + ownship.lat[idx1]))
        dy = re * np.radians(intruder.lat[idx2] - ownship.lat[idx1])
        ahead = (dx * ownship.gseast[idx1] + dy * ownship.gsnorth[idx1]) > 0.0
        # trailing aircraft = the one with the other ahead of it (ownship and
        # intruder are the same ADSB traffic array, so index either via ownship).
        if ahead:                                   # intruder ahead -> ownship trails
            trail, trail_i, lead_i = ac1, idx1, idx2
        else:
            trail, trail_i, lead_i = ac2, idx2, idx1
        hpz_m = np.max(conf.hpz[[idx1, idx2]])
        target_alt = min(ownship.alt[trail_i],
                         ownship.alt[lead_i] - _HOLD_HPZ_MARGIN * hpz_m)
        self._hold[trail] = {'partner': ac2 if trail == ac1 else ac1,
                             'alt': float(target_alt),
                             'tas': float(_CT_SPD_FACTOR * ownship.tas[trail_i])}

    def _maintain_hold(self, ownship, conf):
        ''' Refresh each persistent hold and release it once the pair is
            horizontally clear (> RPZ x release factor). Runs every tick for ALL
            held aircraft — including pairs no longer in confpairs — so the AP
            cannot rebound the trailing aircraft into the conflict. '''
        re = 6371000.0
        for cs in list(self._hold.keys()):
            i = bs.traf.id2idx(cs)
            j = bs.traf.id2idx(self._hold[cs]['partner'])
            if i < 0 or j < 0:                       # an aircraft left the sim
                del self._hold[cs]
                continue
            dx = re * np.radians(bs.traf.lon[j] - bs.traf.lon[i]) * \
                np.cos(0.5 * np.radians(bs.traf.lat[j] + bs.traf.lat[i]))
            dy = re * np.radians(bs.traf.lat[j] - bs.traf.lat[i])
            hdist = np.hypot(dx, dy)
            rpz_m = np.max(conf.rpz[[i, j]])
            # Release only when the pair is PAST CPA (separating along-track) AND
            # horizontally clear. During the approach they can be far apart yet
            # still closing — releasing then lets the trailing ac rebound (bug in
            # iter 3: released on hdist>RPZ alone before the hold did anything).
            dvx = bs.traf.gseast[j] - bs.traf.gseast[i]
            dvy = bs.traf.gsnorth[j] - bs.traf.gsnorth[i]
            separating = (dx * dvx + dy * dvy) > 0.0
            if separating and hdist > _HOLD_RELEASE_FACT * rpz_m:
                del self._hold[cs]
                continue
            # refresh the target below the (possibly climbed) leader
            hpz_m = np.max(conf.hpz[[i, j]])
            self._hold[cs]['alt'] = float(min(self._hold[cs]['alt'],
                                              bs.traf.alt[j] - _HOLD_HPZ_MARGIN * hpz_m))


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
    ''' Vertical-primary resolution + MVP-horizontal-vector heading fallback for
        co-altitude, still-approaching pairs (iter 1 redesign). Targets the
        residual LoS clusters that pure vertical leaves behind. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz   = False
        self.swresospd     = False
        self.swresohdg     = False
        self.swresovert    = True
        self.swfallbackhdg = True


class MVP2NAT_HYBRID2(MVP2NAT):
    ''' Approach A (Stage-C iter 2): H2 hybrid (vertical + MVP-horizontal-vector
        fallback on all co-altitude approaching pairs) PLUS a co-track speed
        reduction. Layers the along-track-separation fix for the adjacent-FL
        co-track residual on top of the working vertical+horizontal channels. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz    = False
        self.swresospd      = False
        self.swresohdg      = False
        self.swresovert     = True
        self.swfallbackhdg  = True     # broad MVP-horizontal fallback (H2)
        self.swspeedcotrack = True     # + co-track speed reduction


class MVP2NAT_CLUSTER(MVP2NAT):
    ''' Approach B (Stage-C iter 2): conflict-cluster resolver. Applies the
        domain that physically fits each conflict class instead of a
        vertical-primary + secondary trigger: vertical is the default best
        domain for head-on / crossing / vertical-crossing (Stage-B: IPR>=0.98),
        and the co-track class (near-parallel, level, adjacent-FL) — which
        vertical cannot open — is routed to a SPEED reduction that builds
        along-track separation. NO broad horizontal fallback (heading was never
        uniquely best for any class). = vertical + targeted co-track speed. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz    = False
        self.swresospd      = False
        self.swresohdg      = False
        self.swresovert     = True
        self.swfallbackhdg  = False    # no broad horizontal fallback
        self.swspeedcotrack = True     # co-track class -> speed


class MVP2NAT_HYBRID3(MVP2NAT):
    ''' Stage-C iter 3: two-regime resolver. Tactical conflicts (crossing /
        head-on / vertical-crossing) use the H2 vertical + MVP-horizontal
        fallback. Sustained CO-TRACK cruise-climb competitions (the residual: two
        same-track aircraft both VNAV-climbing to the same FL) are handled by a
        PERSISTENT FL-hold + speed on the TRAILING aircraft: hold it a fixed FL
        below the leader (stable target -> no runaway dive) and slow it to open
        along-track separation, with authority over VNAV (forced altactive), held
        until horizontally clear (> RPZ) so it cannot rebound. The trailing
        aircraft then climbs + re-accelerates to its optimum behind the leader. '''
    def __init__(self):
        super().__init__()
        self.swresohoriz    = False
        self.swresospd      = False
        self.swresohdg      = False
        self.swresovert     = True
        self.swfallbackhdg  = True     # tactical: MVP-horizontal fallback (H2)
        self.swspeedcotrack = False    # co-track handled by FL-hold below, not the iter-2 speed
        self.swaltholdcotrack = True   # co-track: persistent FL-hold + speed on trailing ac
