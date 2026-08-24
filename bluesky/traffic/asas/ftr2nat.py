''' FTR2NAT — a release criterion that lets go early WITHOUT the conflict
    re-forming.

    Measured motivation (synthetic suite, 2026-08-21/24). Over 17 scenarios,
    MVP2NAT_VERT produced 0 s of loss of separation under PastCPANAT and 113 s
    under FTRNAT. That is not a safety/efficiency trade: it is two methods
    wrong in opposite directions.

    * `PastCPANAT` holds while `not past_cpa or hor_los or is_bouncing`. On
      co-track or still-converging NAT geometry that is permanently true, so
      the aircraft is frozen and nothing can re-form. Safety by immobility --
      and the reason it holds aircraft for 2637 s and never releases at all in
      three of the lifecycle scenarios.
    * `FTRNAT` releases correctly by the instantaneous separation test (every
      measured release had >= RPZ or >= HPZ) and then the conflict comes back.
      One aircraft was released six times, re-engaging after 19, 35, 14, 43 and
      30 s. Each cycle re-commands the aircraft: on L04 the ownship's altitude
      range was 1175 ft under FTRNAT against 259 ft under PastCPANAT, and that
      movement is what drove it into a third aircraft sitting 1041 ft below.

    So the missing property is not a dimension -- FTRNAT is already 3D, and its
    `conflict_predicted` is a faithful mirror of StateBased.detect. The missing
    property is that a release should be *stable*: let go as soon as the
    geometry is genuinely clear, and not one tick earlier.

    Two classes, so a run can attribute an improvement to the right cause:

    `FTR2NAT_FIX`  FTRNAT plus the two defects the audit found, nothing else.
    `FTR2NAT`      the above plus margin, dwell, and the horizontal bracket.

    FTRNAT and PastCPANAT are deliberately left untouched, so all four methods
    run in one batch at a single SHA.
'''
import numpy as np

import bluesky as bs
from bluesky.traffic.asas.ftrnat import FTRNAT


class FTR2NAT_FIX(FTRNAT):
    ''' FTRNAT with the two audited defects fixed, and nothing else.

        Exists so that "FTRNAT was buggy" and "the release needs a margin" are
        separable in the results table rather than confounded in one class.

        DEFECT 1 -- the ASSUMED-intent path fabricates a vertical escape.
        `FTRINTENT` defaults to ASSUMED, and FTRNAT records the intruder's raw
        instantaneous `vs` on the tick the pair is first detected. If the
        intruder was already manoeuvring under ASAS, that is a RESOLUTION rate,
        not intent -- typically the autopilot default of 1500 fpm. FTRNAT then
        sets `dvs_settled` to it and extrapolates it forever, which over a
        300 s lookahead invents ~7500 ft of vertical separation and silently
        turns criterion 2 into a no-op. This is precisely the bug class that
        `desired_vs` and `revert_conflict` were written to eliminate, let back
        in through the ASSUMED path. Fixed by recording the intruder's
        AUTOPILOT target (`desired_vs`, which reads ap.alt/ap.vs and is
        therefore unaffected by ASAS) and by settling at zero rate.

        DEFECT 2 -- the predicate is directionally asymmetric. For (A,B) the
        relative velocity is `v_B_actual - v_A_desired`; for (B,A) it is
        `v_A_actual - v_B_desired`. Those are not negatives of each other, so a
        pair can release in one direction and hold in the other. Measured: on
        L02/L03/L04 FTRNAT left 1 of 2 aircraft engaged where PastCPANAT left
        2 of 2. Fixed by deciding per UNORDERED pair -- release only when both
        directions clear.
    '''

    def _record_assumed(self, pair, intruder, idx2):
        ''' What the intruder would fly if IT reverted -- its own intent.

            FTRNAT stores `intruder.vs[idx2]`, the rate it happens to be flying
            right now, which under an active resolver is the resolver's
            command. `desired_vs` reads the autopilot instead.
        '''
        self.assumed[pair] = (intruder.gseast[idx2],
                              intruder.gsnorth[idx2],
                              self.desired_vs(intruder, idx2))

    def _intruder_settled_alt(self, intruder, idx2, alt_intr):
        ''' The level the intruder ends up at. Its CURRENT one, here. '''
        return alt_intr

    def _assumed_settled(self, alt_intr, apalt_own, vsrevert):
        ''' Where the pair ends up once BOTH have settled, under ASSUMED.

            FTRNAT returns `(alt_intr - apalt_own, vsrevert)`, i.e. it keeps
            the intruder climbing at the assumed rate for ever. Nothing settles
            for ever; returning a zero rate makes the "settled" branch of
            revert_conflict mean what its name says, and is the conservative
            direction for a release criterion.
        '''
        return alt_intr - apalt_own, 0.0

    # ------------------------------------------------------------------
    def _pair_clears(self, conf, ownship, intruder, idx1, idx2, conflict):
        ''' The one-directional FTRNAT criterion, factored out unchanged.

            Returns True when reverting the OWNSHIP is predicted to leave the
            pair out of conflict. Subclasses override the zone margins and the
            bracketing; the state wiring stays here so there is one place where
            "what does the ownship revert to" is decided.
        '''
        re = 6371000.
        dist = re * np.array([np.radians(intruder.lon[idx2] - ownship.lon[idx1]) *
                              np.cos(0.5 * np.radians(intruder.lat[idx2] +
                                                      ownship.lat[idx1])),
                              np.radians(intruder.lat[idx2] - ownship.lat[idx1])])

        rpz, hpz = self._release_zones(conf, idx1, idx2)
        dtlook = np.max(conf.dtlookahead[[idx1, idx2]])

        alt_own = ownship.alt[idx1]
        apalt_own = ownship.ap.alt[idx1]
        alt_intr = intruder.alt[idx2]
        dalt = alt_intr - alt_own

        vown = self.desired_spd_trk(ownship, idx1)
        vsown = self.desired_vs(ownship, idx1)
        vnow = np.array([intruder.gseast[idx2], intruder.gsnorth[idx2]])
        # The CURRENT relative velocity -- both sides actual. This is what the
        # aircraft is flying while it is still turning back onto its route, and
        # it is what detection sees during that transient.
        vrel_now = vnow - np.array([ownship.gseast[idx1], ownship.gsnorth[idx1]])

        # How long the ownship would still be climbing/descending after it
        # reverts, before it captures its commanded level. Past that moment the
        # "now" rate is fiction -- see _ramp_horizon.
        t_cap = self._ramp_horizon(alt_own, apalt_own, vsown, dtlook)

        # Criterion 1: the intruder holds its current state while the ownship
        # reverts and settles at its commanded level.
        #
        # ...except that a rate the RESOLVER is commanding is not a statement
        # about what the intruder will do. This is the same defect as the
        # ASSUMED-intent one, one criterion further up: FTRNAT feeds
        # `intruder.vs` in raw, so when both aircraft are parked off their
        # route levels by a vertical resolver, criterion 1 models the ownship
        # climbing back while the intruder stays put -- a convergence that
        # never happens, because the intruder climbs back too.
        #
        # Measured on R01 under MVP2NAT_VERT: both aircraft frozen ~1600 and
        # ~2000 ft below their routes, 2100 ft apart, criterion 2 clear and
        # criterion 1 holding, so the pair was never released in 4002 s.
        vs_intr = (self.desired_vs(intruder, idx2)
                   if self._asas_owns(idx2) else intruder.vs[idx2])
        # Where the intruder ENDS UP -- again, not where the resolver is
        # currently holding it (see _intruder_settled_alt).
        alt_intr_end = self._intruder_settled_alt(intruder, idx2, alt_intr)
        clear = not self._revert_conflict(
            dist, vnow - vown, vrel_now, rpz, hpz, dtlook,
            dalt, vs_intr - vsown,
            alt_intr_end - apalt_own, vs_intr, t_cap)

        # Criterion 2: the intruder reverts to its desired state too.
        if clear and self.intent != 'OFF':
            vrevert = vsrevert = None
            if self.intent == 'ASSUMED' and self._asas_owns(idx2):
                # The ASSUMED record is taken on the tick the pair is first
                # detected and never refreshed. That is fine for an intruder
                # flying its own plan, but useless for one the resolver has
                # since taken over and parked off its route level: the record
                # says "level" while the aircraft now intends 1500 fpm back up.
                # Measured on R01: criteria 1 cleared, criterion 2 held on a
                # stale record, and the pair was never released in 4002 s.
                # Where ASAS owns the intruder, read its live intent instead.
                vrevert = self.desired_spd_trk(intruder, idx2)
                vsrevert = self.desired_vs(intruder, idx2)
                dalt_settled = alt_intr_end - apalt_own
                dvs_settled = 0.0
            elif self.intent == 'ASSUMED':
                assumed = self.assumed.get(conflict)
                if assumed is None:
                    # FTRNAT drops criterion 2 for the pair's whole lifetime if
                    # the record is missing. Record it now instead.
                    self._record_assumed(conflict, intruder, idx2)
                    assumed = self.assumed.get(conflict)
                if assumed is not None:
                    vrevert = np.array(assumed[:2])
                    vsrevert = assumed[2]
                    dalt_settled, dvs_settled = self._assumed_settled(
                        alt_intr_end, apalt_own, vsrevert)
            else:
                vrevert = self.desired_spd_trk(intruder, idx2)
                vsrevert = self.desired_vs(intruder, idx2)
                dalt_settled = intruder.ap.alt[idx2] - apalt_own
                dvs_settled = 0.0
            if vrevert is not None:
                clear = not self._revert_conflict(
                    dist, vrevert - vown, vrel_now, rpz, hpz, dtlook,
                    dalt, vsrevert - vsown,
                    dalt_settled, dvs_settled, t_cap)
        return clear

    def _release_zones(self, conf, idx1, idx2):
        ''' The zone sizes the release is judged against. Plain here. '''
        return (np.max(conf.rpz[[idx1, idx2]]),
                np.max(conf.hpz[[idx1, idx2]]))

    @staticmethod
    def _ramp_horizon(alt, apalt, vs_desired, dtlook):
        ''' How long the reverted vertical ramp actually lasts, in seconds.

            revert_conflict brackets the ramp-then-hold profile by ORing a
            "now" case (current levels closing at the current rate) with a
            "settled" case. The now case extrapolates that rate FOREVER, which
            is what the bracket is for -- but it means any transient
            convergence during level capture holds the pair indefinitely.
            Bounding the now case by the capture time is the missing half of
            the bracket: a predicted conflict that only happens AFTER the
            aircraft has levelled off is not a conflict.
        '''
        if abs(vs_desired) < 1e-6:
            return dtlook
        return min(float(dtlook), abs(apalt - alt) / abs(vs_desired))

    def _revert_conflict(self, dist, vrel_revert, vrel_now, rpz, hpz, dtlook,
                         dalt_now, dvs_now, dalt_settled, dvs_settled,
                         dtlook_now=None):
        ''' Hold-or-release for one direction. `vrel_now` and `dtlook_now` are
            unused here and are taken up by FTR2NAT. '''
        return FTRNAT.revert_conflict(dist, vrel_revert, rpz, hpz, dtlook,
                                      dalt_now, dvs_now,
                                      dalt_settled, dvs_settled)

    def _release_ok(self, key, clear, dt):
        ''' Turn a per-tick predicate into a release decision. Immediate here;
            FTR2NAT requires it to hold for a dwell time. '''
        return clear

    # ------------------------------------------------------------------
    def resumenav(self, conf, ownship, intruder):
        ''' As FTRNAT, but the release decision is taken per UNORDERED pair.

            Both directions are evaluated and a pair is released only when both
            clear, so `cr.active` can never be dropped for one aircraft of a
            pair while the other is still resolving.
        '''
        for pair in conf.confpairs:
            if pair not in self.resopairs:
                idx2 = bs.traf.id2idx(pair[1])
                if idx2 >= 0:
                    self._record_assumed(pair, intruder, idx2)

        self.resopairs.update(conf.confpairs)

        delpairs = set()
        changeactive = dict()
        dt = self._resnav_dt()

        # First pass: the per-direction predicate, collected by unordered pair.
        clear_by_pair = {}
        seen = {}
        for conflict in self.resopairs:
            idx1, idx2 = bs.traf.id2idx(conflict)
            if idx1 < 0:
                delpairs.add(conflict)
                continue
            key = frozenset(conflict)
            seen.setdefault(key, []).append((conflict, idx1, idx2))
            if idx2 < 0:
                # The intruder is gone: nothing to hold against.
                clear_by_pair[key] = clear_by_pair.get(key, True) and True
                continue
            one = self._pair_clears(conf, ownship, intruder, idx1, idx2, conflict)
            clear_by_pair[key] = clear_by_pair.get(key, True) and one

        # Second pass: apply the pair-level decision to both directions.
        for key, entries in seen.items():
            clear = clear_by_pair.get(key, True)
            release = self._release_ok(key, clear, dt)
            for conflict, idx1, idx2 in entries:
                if idx2 >= 0 and not release:
                    changeactive[idx1] = True
                else:
                    changeactive[idx1] = changeactive.get(idx1, False)
                    delpairs.add(conflict)

        for idx, active in changeactive.items():
            bs.traf.cr.active[idx] = active
            if not active:
                self._advance_past_overflown(idx)

        self.resopairs -= delpairs
        for conflict in delpairs:
            self.assumed.pop(conflict, None)
            self._forget(frozenset(conflict))

    def _intruder_settled_alt(self, intruder, idx2, alt_intr):
        ''' Where the intruder ENDS UP once everyone has been released.

            FTRNAT uses its CURRENT altitude, but under a vertical resolver
            that is a held value, not a destination. Measured on R01: the
            intruder was frozen at 36069 ft while its own route commanded
            38081 ft, so the settled gap read 538 ft against the ownship's
            35531 ft -- inside HPZ, so the pair was held for ever. The true
            settled gap was 2550 ft, comfortably clear.

            Only substituted while ASAS actually owns the intruder; otherwise
            its current altitude IS where it is going.
        '''
        if self._asas_owns(idx2):
            try:
                return float(intruder.ap.alt[idx2])
            except Exception:
                return alt_intr
        return alt_intr

    @staticmethod
    def _asas_owns(idx):
        ''' Is the resolver currently commanding this aircraft?

            If so its instantaneous state is a RESOLUTION, not intent, and must
            not be extrapolated as a prediction of future behaviour.
        '''
        try:
            act = bs.traf.cr.active
            return bool(0 <= idx < len(act) and act[idx])
        except Exception:
            return False

    @staticmethod
    def _resnav_dt():
        ''' Seconds of simulated time between two resumenav calls.

            NOT bs.sim.simdt. resumenav is driven by traffic.py's `asastimer`
            (traffic.py:82,407), so it runs once per ASAS interval -- 1.0 s by
            default -- while simdt is the integration step, 0.05 s here. Using
            simdt made a 15 s dwell into a 300 s one, a factor of 20.
        '''
        dt = getattr(bs.settings, 'asas_dt', None)
        if not dt:
            dt = getattr(bs.sim, 'simdt', 1.0)
        return float(dt) or 1.0

    def _forget(self, key):
        ''' Drop any per-pair release state. Nothing to drop here. '''
        return


class FTR2NAT(FTR2NAT_FIX):
    ''' FTR2NAT_FIX plus the three things that make a release STABLE.

        1. MARGIN, per plane. FTRNAT judges the release on plain `conf.rpz` /
           `conf.hpz` -- the bare detection zone -- while the resolver is
           pushing toward `rpz*resofach` / `hpz*resofacv`. It therefore lets go
           at exactly the boundary detection re-triggers on, and any wobble
           brings the conflict straight back. Releasing on the same margin the
           resolver aims at removes that.

        2. DWELL. `revert_conflict` is a pure instantaneous predicate with no
           hysteresis, re-evaluated every tick. A predicate flapping at its own
           boundary produces exactly the observed chatter, so the clearance has
           to persist for `_rel_dwell_s` before the aircraft is handed back.

        3. HORIZONTAL BRACKET. `revert_conflict` ORs a "now" and a "settled"
           vertical case -- but both branches use the SAME horizontal terms, so
           the bracket is vertical-only. The horizontal side assumes the
           ownship attains its autopilot track and speed instantaneously, when
           in fact it ramps at the turn rate. During that transient the ACTUAL
           state is what detection sees, which is mechanism 3 behind the
           re-engagements. Holding when EITHER the current or the reverted
           horizontal velocity conflicts is the honest reading of "the release
           must check both dimensions".
    '''

    def __init__(self):
        super().__init__()
        # Default to the resolver's own margins, so release and resolution
        # agree on where "clear" is. bs.traf.cr may not exist at construction.
        self._rel_fach = None      # None => follow cr.resofach at run time
        self._rel_facv = None      # None => follow cr.resofacv at run time
        self._rel_dwell_s = 15.0
        self._dwell = {}

    def reset(self):
        super().reset()
        self._dwell.clear()
        self._rel_fach = None
        self._rel_facv = None
        self._rel_dwell_s = 15.0

    def _margins(self):
        fach, facv = self._rel_fach, self._rel_facv
        if fach is None or facv is None:
            cr = getattr(bs.traf, 'cr', None)
            fach = getattr(cr, 'resofach', 1.0) if fach is None else fach
            facv = getattr(cr, 'resofacv', 1.0) if facv is None else facv
        return float(fach), float(facv)

    def _release_zones(self, conf, idx1, idx2):
        fach, facv = self._margins()
        return (np.max(conf.rpz[[idx1, idx2]]) * fach,
                np.max(conf.hpz[[idx1, idx2]]) * facv)

    def _revert_conflict(self, dist, vrel_revert, vrel_now, rpz, hpz, dtlook,
                         dalt_now, dvs_now, dalt_settled, dvs_settled,
                         dtlook_now=None):
        """Bracket the reversion in BOTH planes, and bound the ramp in time.

        Two brackets, ORed as FTRNAT's is -- hold if any branch conflicts.

        HORIZONTAL: the reverted velocity is exact once the aircraft has
        finished turning; the current velocity is exact until it starts. FTRNAT
        only ever uses the reverted one, so it assumes the turn is
        instantaneous.

        VERTICAL, and this is the part that decides whether the vertical domain
        can release at all: the "now" case is bounded by the CAPTURE time. The
        resolver parks an aircraft off its route level, so desired_vs is
        +/-1500 fpm back toward it, and an unbounded now-case reads that as a
        permanent closure and holds the pair for ever. A predicted conflict
        that only happens after the aircraft has levelled off is not a
        conflict. The settled case keeps the full lookahead, so a genuine
        convergence is still caught.
        """
        tl_now = dtlook if dtlook_now is None else min(dtlook, dtlook_now)
        for vrel in (vrel_revert, vrel_now):
            if FTRNAT.conflict_predicted(dist, vrel, dalt_now, dvs_now,
                                         rpz, hpz, tl_now):
                return True
            if FTRNAT.conflict_predicted(dist, vrel, dalt_settled, dvs_settled,
                                         rpz, hpz, dtlook):
                return True
        return False

    def _release_ok(self, key, clear, dt):
        ''' Release only after the clearance has held for the dwell time. '''
        if not clear:
            self._dwell.pop(key, None)
            return False
        held = self._dwell.get(key, 0.0) + dt
        self._dwell[key] = held
        return held >= self._rel_dwell_s

    def _forget(self, key):
        self._dwell.pop(key, None)


class FTR2NAT_NOMARGIN(FTR2NAT):
    ''' Ablation: dwell + horizontal bracket, but no zone margin. '''
    def __init__(self):
        super().__init__()
        self._rel_fach = 1.0
        self._rel_facv = 1.0


class FTR2NAT_NODWELL(FTR2NAT):
    ''' Ablation: margin + horizontal bracket, but no dwell. '''
    def __init__(self):
        super().__init__()
        self._rel_dwell_s = 0.0
