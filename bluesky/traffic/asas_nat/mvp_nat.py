"""
mvp_nat.py — Modified Voltage Potential resolution for NAT en-route traffic.
=============================================================================

Inherits from the stock ``MVP`` so it appears in the standard ``RESO``
command (``RESO MVPNAT``).  All changes are NAT-specific; the upstream
MVP remains untouched.

Design
------
Default resolution dimension : **vertical** (original MVP defaults to
horizontal).  A control-logic layer decides *per conflict pair* which
dimension to use:

  1. **Vertical first** — preferred for NAT (fuel-efficient, small
     deviation from track).
  2. **Ceiling guard** — if the ASAS altitude would exceed the aircraft
     service ceiling (``perf.hmax``), suppress the vertical component
     for that aircraft and fall back to horizontal.
  3. **MMO guard** — if horizontal speed resolution would push an
     aircraft past MMO, suppress the speed component and use heading
     (or vertical).
  4. **Symmetry breaker** — when both aircraft are at the same FL
     (``vrel_vert ≈ 0``), deterministically assign one aircraft UP and
     the other DOWN, instead of both getting the same sign.
  5. **VS-stop** — if two aircraft are already within horizontal PZ and
     a vertical manoeuvre (climb/descent at waypoint) would breach the
     vertical PZ next timestep, command VS = 0 to prevent intrusion.
  6. **Anti-bouncing for parallel tracks** — if two aircraft are on
     near-parallel headings (Δtrack < 30°) and in sustained conflict,
     keep vertical resolution active until well past CPA instead of
     toggling on/off.

Recovery algorithm (Schaberg 2-criteria, NAT-adapted)
-----------------------------------------------------
Replaces the stock 1-criterion past-CPA check.  Two criteria must
both pass before ASAS hands control back to the autopilot:

  **Criterion 1**: Dcpa(Vo_desired, Vi_current_reso) > RPZ × factor
  **Criterion 2**: Dcpa(Vo_desired, Vi_initial)      > RPZ × factor

Additional NAT guards:
  - min hold time (default 30 s) to prevent premature recovery
  - max hold time (default 0 = off) as safety valve for same-dir
    traffic where Dcpa may never exceed threshold

Parameters (configurable via stack commands):
  - ``nat_ceiling_margin_ft``: altitude buffer below perf.hmax [ft]
  - ``nat_min_hold_t``: minimum ASAS hold time [s]
  - ``nat_max_hold_t``: maximum ASAS hold time [s]  (0 = unlimited)
  - ``nat_dcpa_margin``: Dcpa threshold multiplier on RPZ

Stack commands registered
-------------------------
  RESO MVPNAT          — activate this resolver
  RMETHH_NAT / RMETHV_NAT  — override dimension switches (testing)
  NATCEILMARGIN [ft]   — set ceiling margin
  NATDOMAIN [domain]   — set resolution domain strategy:
                          VERT        vertical only (default)
                          HORIZ_BOTH  horizontal speed + heading
                          HORIZ_SPD   horizontal speed only
                          HORIZ_HDG   horizontal heading only
                          VERT_FIRST  try vertical, fallback horizontal
                          HORIZ_FIRST try horizontal, fallback vertical
  NATRECOVERY [setting] [value] — configure 2-criteria recovery:
                          MINHOLD s   min hold time (default 30)
                          MAXHOLD s   max hold time (default 0=off)
                          DCPAFAC x   Dcpa factor on RPZ (default 1.0)

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-04
"""
import numpy as np

import bluesky as bs
from bluesky import stack
from bluesky.tools.aero import ft, nm, fpm, vtas2mach
from bluesky.traffic.asas.mvp import MVP

# --- DEBUG ---
_MVPNAT_VERSION = '2026-04-10b'  # bump on every edit
print(f'[MVPNAT] loaded version {_MVPNAT_VERSION}')
_dbg_count = 0   # limit debug output to first N calls
_DBG_MAX   = 40  # max debug lines

# --- Domain strategy constants ---
_VALID_DOMAINS = {
    'VERT',          # vertical only (default)
    'HORIZ_BOTH',    # horizontal speed + heading
    'HORIZ_SPD',     # horizontal speed only
    'HORIZ_HDG',     # horizontal heading only
    'VERT_FIRST',    # try vertical, fallback horizontal if guard triggers
    'HORIZ_FIRST',   # try horizontal, fallback vertical if guard triggers
}

_dbg_recov_count = 0  # separate debug counter for recovery


def _compute_dcpa_h(dist_h, vrel_h):
    """Horizontal distance at CPA given 2-D position & relative velocity.

    Parameters
    ----------
    dist_h : ndarray (2,)   east/north distance vector [m]
    vrel_h : ndarray (2,)   east/north relative velocity [m/s]

    Returns
    -------
    float   distance at CPA [m].  If tcpa < 0 (diverging), return
            current distance (aircraft will never be closer).
    """
    vrel_sq = float(np.dot(vrel_h, vrel_h))
    if vrel_sq < 1e-6:
        # Nearly zero relative velocity → distance is constant.
        return float(np.linalg.norm(dist_h))
    tcpa = -float(np.dot(dist_h, vrel_h)) / vrel_sq
    if tcpa < 0:
        # CPA is in the past → aircraft are diverging.
        return float(np.linalg.norm(dist_h))
    dcpa_vec = dist_h + vrel_h * tcpa
    return float(np.linalg.norm(dcpa_vec))


class MVPNAT(MVP):
    """Modified Voltage Potential — NAT variant.

    Key differences from stock MVP
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - Default dimension is **vertical** (swresovert=True).
    - Per-pair control logic: ceiling guard, MMO guard, symmetry
      breaker, VS-stop, anti-bouncing.
    - Altitude computation clamps to service ceiling and floor.
    """

    def __init__(self):
        super().__init__()
        # ------ Override defaults: vertical first ------
        self.swresohoriz = False
        self.swresospd   = False
        self.swresohdg   = False
        self.swresovert  = True

        # ------ NAT-specific tunables ------
        self.nat_ceiling_margin_ft = 500.0  # buffer below hmax [ft]
        self.nat_domain = 'VERT'            # resolution domain strategy

        # ------ Schaberg 2-criteria recovery ------
        self._preconf_state = {}   # acid → {gseast, gsnorth, vs, alt}
        self._reso_start_t  = {}   # (ac1, ac2) → simt at first detection
        self.nat_min_hold_t  = 30.0  # [s] min ASAS hold before recovery check
        self.nat_max_hold_t  = 0.0   # [s] max ASAS hold, force recovery (0 = off)
        self.nat_dcpa_margin = 1.0   # Dcpa threshold = RPZ × this factor

    def reset(self):
        super().reset()
        self._preconf_state.clear()
        self._reso_start_t.clear()

    # ==================================================================
    #  Stack commands for NAT tunables
    # ==================================================================
    @stack.command(name='NATCEILMARGIN')
    def set_ceil_margin(self, value: float = -1.0):
        """Set altitude margin below service ceiling [ft]."""
        if value < 0:
            return True, (f'NATCEILMARGIN [ft]\n'
                          f'Current: {self.nat_ceiling_margin_ft:.0f} ft')
        self.nat_ceiling_margin_ft = value
        return True, f'NAT ceiling margin set to {value:.0f} ft'

    @stack.command(name='NATDOMAIN')
    def set_nat_domain(self, domain: 'txt' = ''):
        """Set resolution domain strategy.

        Usage: NATDOMAIN [domain]

        Domains
        -------
        VERT        : Vertical only (default). Preferred for NAT.
        HORIZ_BOTH  : Horizontal speed + heading combined.
        HORIZ_SPD   : Horizontal speed only.
        HORIZ_HDG   : Horizontal heading only.
        VERT_FIRST  : Try vertical; fallback horizontal if ceiling guard
                      blocks both aircraft.
        HORIZ_FIRST : Try horizontal; fallback vertical if MMO guard
                      blocks the speed change.

        Without argument: print current domain.
        """
        if not domain:
            return True, (f'NATDOMAIN [domain]\n'
                          f'Current: {self.nat_domain}\n'
                          f'Options: {", ".join(sorted(_VALID_DOMAINS))}')
        domain = domain.upper().strip()
        if domain not in _VALID_DOMAINS:
            return False, (f'Unknown domain "{domain}". '
                           f'Valid: {", ".join(sorted(_VALID_DOMAINS))}')

        self.nat_domain = domain

        # Update the base-class dimension switches to match
        if domain == 'VERT':
            self.swresovert = True
            self.swresohoriz = False
            self.swresospd = False
            self.swresohdg = False
        elif domain == 'HORIZ_BOTH':
            self.swresovert = False
            self.swresohoriz = True
            self.swresospd = True
            self.swresohdg = True
        elif domain == 'HORIZ_SPD':
            self.swresovert = False
            self.swresohoriz = True
            self.swresospd = True
            self.swresohdg = False
        elif domain == 'HORIZ_HDG':
            self.swresovert = False
            self.swresohoriz = True
            self.swresospd = False
            self.swresohdg = True
        elif domain == 'VERT_FIRST':
            self.swresovert = True
            self.swresohoriz = False
            self.swresospd = True
            self.swresohdg = True
        elif domain == 'HORIZ_FIRST':
            self.swresovert = False
            self.swresohoriz = True
            self.swresospd = True
            self.swresohdg = True

        print(f'[MVPNAT] Domain set to {self.nat_domain}  '
              f'(swresovert={self.swresovert}, swresohoriz={self.swresohoriz}, '
              f'swresospd={self.swresospd}, swresohdg={self.swresohdg})')
        return True, f'NAT domain set to {self.nat_domain}'

    @stack.command(name='NATRECOVERY')
    def set_recovery(self, setting: 'txt' = '', value: float = -1.0):
        """Configure Schaberg 2-criteria recovery parameters.

        Usage: NATRECOVERY [setting] [value]

        Settings
        --------
        MINHOLD  [s]   Minimum hold time before recovery check (default 30).
        MAXHOLD  [s]   Maximum hold time, force recovery. 0 = off (default 0).
        DCPAFAC  [x]   Dcpa threshold as multiplier of RPZ (default 1.0).

        Without arguments: print current values.
        """
        if not setting:
            return True, (
                f'NATRECOVERY [MINHOLD|MAXHOLD|DCPAFAC] [value]\n'
                f'  MINHOLD  = {self.nat_min_hold_t:.0f} s\n'
                f'  MAXHOLD  = {self.nat_max_hold_t:.0f} s  (0 = unlimited)\n'
                f'  DCPAFAC  = {self.nat_dcpa_margin:.2f}\n'
                f'  Pre-conflict states stored: {len(self._preconf_state)}\n'
                f'  Active reso pairs tracked : {len(self._reso_start_t)}')

        setting = setting.upper().strip()
        if value < 0:
            return False, f'NATRECOVERY {setting} [value]  — value must be >= 0'

        if setting == 'MINHOLD':
            self.nat_min_hold_t = value
            return True, f'NAT recovery min hold set to {value:.0f} s'
        elif setting == 'MAXHOLD':
            self.nat_max_hold_t = value
            return True, f'NAT recovery max hold set to {value:.0f} s (0 = unlimited)'
        elif setting == 'DCPAFAC':
            self.nat_dcpa_margin = max(value, 0.5)
            return True, f'NAT recovery Dcpa factor set to {self.nat_dcpa_margin:.2f}'
        else:
            return False, (f'Unknown setting "{setting}". '
                           f'Valid: MINHOLD, MAXHOLD, DCPAFAC')

    # Re-register dimension switches under _NAT names so they don't
    # collide with the original MVP commands when both are loaded.
    @stack.command(name="RMETHH_NAT")
    def setresometh_nat(self, value: 'txt' = ''):
        """Horizontal resolution method for MVPNAT."""
        return super().setresometh(value)

    @stack.command(name='RMETHV_NAT')
    def setresometv_nat(self, value: 'txt' = ''):
        """Vertical resolution method for MVPNAT."""
        return super().setresometv(value)

    # ==================================================================
    #  resolve() — main entry: loops all conflict pairs
    # ==================================================================
    def resolve(self, conf, ownship, intruder):
        """Resolve all current conflicts with NAT control logic."""
        ntraf = ownship.ntraf
        dv = np.zeros((ntraf, 3))
        timesolveV = np.ones(ntraf) * 1e9

        # Pre-compute per-aircraft limits once
        # Use hmo (max operating altitude, ~FL431 for B77W) instead of
        # hmax (max alt at MTOW, ~FL343 for B77W!).  hmax is the
        # weight-limited ceiling which can be BELOW cruise altitude.
        hmax_raw = ownship.perf.hmax
        hmo = getattr(ownship.perf, 'hmo', hmax_raw)  # BADA has hmo
        # Use the higher of hmo and hmax (some models only set hmax)
        hmax = np.maximum(hmo, hmax_raw)
        ceil_margin_m = self.nat_ceiling_margin_ft * ft  # [m]
        mmo = getattr(ownship.perf, 'mmo',
                      np.full(ntraf, 0.87))             # MMO array

        # ----------------------------------------------------------
        #  Per-pair loop
        # ----------------------------------------------------------
        for (ac1, ac2), qdr, dist, tcpa, tLOS in zip(
                conf.confpairs, conf.qdr, conf.dist,
                conf.tcpa, conf.tLOS):

            idx1 = ownship.id.index(ac1)
            idx2 = intruder.id.index(ac2)
            if idx1 < 0 or idx2 < 0:
                continue

            # --- Raw MVP vector (all 3 components) ---
            dv_mvp, tsolV = self.MVP_NAT(
                ownship, intruder, conf, qdr, dist, tcpa, tLOS,
                idx1, idx2)

            # ===========================================================
            #  CONTROL LOGIC — per-pair dimension selection
            # ===========================================================
            # --- Initial dimension from domain strategy ---
            if self.nat_domain == 'VERT':
                use_vert  = True
                use_horiz = False
            elif self.nat_domain in ('HORIZ_BOTH', 'HORIZ_SPD', 'HORIZ_HDG'):
                use_vert  = False
                use_horiz = True
            elif self.nat_domain == 'VERT_FIRST':
                use_vert  = True
                use_horiz = False
            elif self.nat_domain == 'HORIZ_FIRST':
                use_vert  = False
                use_horiz = True
            else:
                # Fallback to vertical
                use_vert  = True
                use_horiz = False

            alt1 = ownship.alt[idx1]
            alt2 = intruder.alt[idx2]

            # --- 1) Ceiling guard ---
            # After cooperative halving, ownship gets  dv_eff = -dv_mvp/2.
            # Positive dv_eff[2] → ownship climbs.
            # dv_mvp[2] < 0  →  ownship climbs  (because -(-) = +)
            # dv_mvp[2] > 0  →  ownship descends (because -(+) = -)
            #
            # Intruder (from *its* loop iteration where it is ownship)
            # gets  dv_eff = -(-dv_mvp)/2 = +dv_mvp/2.
            # So intruder climbs when dv_mvp[2] > 0.
            near_ceil_1 = alt1 > (hmax[idx1] - ceil_margin_m)
            near_ceil_2 = alt2 > (hmax[idx2] - ceil_margin_m)

            ownship_would_climb = dv_mvp[2] < 0  # ownship gets -dv_mvp
            intruder_would_climb = dv_mvp[2] > 0  # intruder gets +dv_mvp

            if (near_ceil_1 and ownship_would_climb) or \
               (near_ceil_2 and intruder_would_climb):
                # At least one aircraft would breach ceiling
                if not (near_ceil_1 and near_ceil_2):
                    # Only one is near ceiling -> flip who climbs/descends
                    dv_mvp[2] = -dv_mvp[2]
                else:
                    # Both near ceiling
                    if use_vert and not use_horiz:
                        if self.nat_domain == 'VERT_FIRST':
                            # Fallback: switch to horizontal
                            use_vert  = False
                            use_horiz = True
                        else:
                            # Force both descend (VERT-only, no fallback)
                            dv_mvp[2] = abs(dv_mvp[2])
                    else:
                        # Already horizontal or combined → keep as is
                        dv_mvp[2] = abs(dv_mvp[2])

            # --- 2) MMO guard (only relevant if horizontal is used) ---
            if use_horiz:
                gs1 = max(ownship.gs[idx1], 1.0)
                newgs1 = np.sqrt(
                    (ownship.gseast[idx1] - dv_mvp[0] * 0.5)**2
                    + (ownship.gsnorth[idx1] - dv_mvp[1] * 0.5)**2)
                mach_ratio = newgs1 / gs1
                mach_new = vtas2mach(ownship.tas[idx1], alt1) * mach_ratio
                if mach_new > mmo[idx1]:
                    if self.nat_domain == 'HORIZ_SPD':
                        # Speed-only domain but MMO violated:
                        # nothing to do here — perf.limits will cap
                        pass
                    elif self.nat_domain == 'HORIZ_HDG':
                        # Heading-only: MMO irrelevant (no speed change)
                        pass
                    elif self.nat_domain == 'HORIZ_FIRST':
                        # Fallback: switch to vertical
                        if not (near_ceil_1 and near_ceil_2):
                            use_vert  = True
                            use_horiz = False
                        # else: stay horizontal, perf.limits will cap MMO
                    else:
                        # HORIZ_BOTH or fallback from VERT_FIRST
                        # Suppress speed, keep heading (or fall back to vert)
                        if not (near_ceil_1 and near_ceil_2):
                            use_vert = True
                            use_horiz = False
                        # else: keep horizontal but perf.limits will cap MMO

            # --- 3) Apply dimension mask ---
            if use_vert and not use_horiz:
                dv_mvp[0] = 0.0
                dv_mvp[1] = 0.0
            elif use_horiz and not use_vert:
                dv_mvp[2] = 0.0
                # Sub-mode: suppress speed or heading component
                if self.nat_domain == 'HORIZ_SPD':
                    # Keep speed (east/north magnitude), suppress heading
                    # by projecting dv onto the current velocity direction
                    ge = ownship.gseast[idx1]
                    gn = ownship.gsnorth[idx1]
                    gs_sq = ge**2 + gn**2
                    if gs_sq > 1.0:
                        # Project dv_h onto velocity unit vector → speed only
                        proj = (dv_mvp[0]*ge + dv_mvp[1]*gn) / gs_sq
                        dv_mvp[0] = proj * ge
                        dv_mvp[1] = proj * gn
                elif self.nat_domain == 'HORIZ_HDG':
                    # Keep heading (perpendicular component), suppress speed
                    ge = ownship.gseast[idx1]
                    gn = ownship.gsnorth[idx1]
                    gs_sq = ge**2 + gn**2
                    if gs_sq > 1.0:
                        # Remove the parallel component → heading only
                        proj = (dv_mvp[0]*ge + dv_mvp[1]*gn) / gs_sq
                        dv_mvp[0] -= proj * ge
                        dv_mvp[1] -= proj * gn
            # else: both True → combined (rare fallback)

            # --- 4) VS-stop: prevent vertical PZ breach from AP ---
            rpz_m = np.max(conf.rpz[[idx1, idx2]])
            hpz_m = np.max(conf.hpz[[idx1, idx2]])
            if dist < rpz_m:
                # Within horizontal PZ — check if AP-commanded VS
                # would breach vertical PZ next timestep
                dalt_now = abs(alt1 - alt2)
                # Use autopilot VS to predict next-step altitude
                selalt1 = ownship.selalt[idx1]
                selalt2 = intruder.selalt[idx2]
                vs1_ap = abs(ownship.ap.vs[idx1]) * np.sign(selalt1 - alt1)
                vs2_ap = abs(intruder.ap.vs[idx2]) * np.sign(selalt2 - alt2)
                dalt_next = abs((alt1 + vs1_ap * bs.sim.simdt)
                                - (alt2 + vs2_ap * bs.sim.simdt))
                if dalt_now >= hpz_m and dalt_next < hpz_m:
                    # Imminent vertical PZ breach → zero VS component
                    dv_mvp[2] = 0.0

            # --- 5) Anti-bouncing for near-parallel tracks ---
            # If headings are nearly parallel (|Δtrk| < 30°) and conflict
            # is long-lived (tLOS is small → already deep in conflict),
            # keep vertical resolution strong and suppress horizontal
            # to avoid lateral oscillation.
            dtrk = abs(ownship.trk[idx1] - intruder.trk[idx2])
            if dtrk > 180:
                dtrk = 360 - dtrk
            if dtrk < 30.0 and dist < rpz_m * self.resofach:
                # Near-parallel and close → force vertical, suppress horiz
                dv_mvp[0] = 0.0
                dv_mvp[1] = 0.0

            # ===========================================================
            #  Accumulate dv (cooperative halving, same as stock MVP)
            # ===========================================================
            if tsolV < timesolveV[idx1]:
                timesolveV[idx1] = tsolV

            # --- DEBUG: log pre-halving values ---
            global _dbg_count
            if _dbg_count < _DBG_MAX:
                _dbg_count += 1
                print(f'[MVPNAT] t={bs.sim.simt:.1f} pair=({ac1},{ac2}) '
                      f'idx=({idx1},{idx2}) '
                      f'dv_mvp=[{dv_mvp[0]:.4f},{dv_mvp[1]:.4f},{dv_mvp[2]:.4f}] '
                      f'tsolV={tsolV:.2f} '
                      f'use_v={use_vert} use_h={use_horiz} '
                      f'domain={self.nat_domain} '
                      f'alt={alt1/ft:.0f}/{alt2/ft:.0f}ft '
                      f'hmax={hmax[idx1]/ft:.0f}/{hmax[idx2]/ft:.0f}ft '
                      f'near_c={near_ceil_1}/{near_ceil_2}')

            if self.swprio:
                dv[idx1], _ = self.applyprio(
                    dv_mvp, dv[idx1], dv[idx2],
                    ownship.vs[idx1], intruder.vs[idx2])
            else:
                dv_mvp[2] = 0.5 * dv_mvp[2]
                dv[idx1] = dv[idx1] - dv_mvp

            # NORESO / RESOOFF
            if self.noresoac[idx2]:
                dv[idx1] = dv[idx1] + dv_mvp
            if self.resooffac[idx1]:
                dv[idx1] = 0.0

        # ==============================================================
        #  Convert dv to track / gs / vs / alt
        # ==============================================================
        dv = np.transpose(dv)
        v = np.array([ownship.gseast, ownship.gsnorth, ownship.vs])
        newv = v + dv

        # --- Apply global dimension switches --------------------------
        # MVPNAT returns TAS directly (not GS), avoiding the
        # singularities from multiple GS↔TAS conversions at NAT
        # flight angles.  For horizontal fallback the GS-domain
        # value is kept (accepted trade-off for rare fallback path).
        #
        # For *_FIRST domains, different aircraft may use different
        # dimensions (decided per-pair in the loop above), so we
        # must compute all three output channels from dv.
        if self.nat_domain in ('VERT_FIRST', 'HORIZ_FIRST'):
            # Mixed-dimension mode: extract all channels from dv.
            # Aircraft with dv_h zeroed keep their AP track/TAS.
            # Aircraft with dv_v zeroed keep their AP VS.
            has_dv_h = (np.abs(dv[0, :]) + np.abs(dv[1, :])) > 1e-8
            has_dv_v = np.abs(dv[2, :]) > 1e-8

            newtrack = np.where(
                has_dv_h,
                (np.arctan2(newv[0, :], newv[1, :]) * 180 / np.pi) % 360,
                ownship.trk)
            newtas = np.where(
                has_dv_h,
                np.sqrt(newv[0, :]**2 + newv[1, :]**2),
                ownship.tas)
            newvs = np.where(has_dv_v, newv[2, :], ownship.vs)
        elif self.swresohoriz:
            if self.swresospd and not self.swresohdg:
                newtrack = ownship.trk
                newtas = np.sqrt(newv[0, :]**2 + newv[1, :]**2)
                newvs = ownship.vs
            elif self.swresohdg and not self.swresospd:
                newtrack = (np.arctan2(newv[0, :], newv[1, :])
                            * 180 / np.pi) % 360
                newtas = ownship.gs
                newvs = ownship.vs
            else:
                newtrack = (np.arctan2(newv[0, :], newv[1, :])
                            * 180 / np.pi) % 360
                newtas = np.sqrt(newv[0, :]**2 + newv[1, :]**2)
                newvs = ownship.vs
        elif self.swresovert:
            newtrack = ownship.trk
            newtas = ownship.tas  # actual TAS, not GS
            newvs = newv[2, :]
        else:
            newtrack = (np.arctan2(newv[0, :], newv[1, :])
                        * 180 / np.pi) % 360
            newtas = np.sqrt(newv[0, :]**2 + newv[1, :]**2)
            newvs = newv[2, :]

        # --- Cap to performance limits --------------------------------
        # Use perf.limits() for proper TAS-domain limiting instead of
        # the old manual min/max capping that compared ground speed
        # against CAS-based vmin/vmax (unit mismatch).  perf.limits()
        # correctly handles MMO, vmin, vmo, and VS limits and returns
        # TAS directly.
        allowed_tas, vscapped, _ = ownship.perf.limits(
            newtas, newvs, ownship.alt, bs.traf.ax)

        # --- Altitude computation -------------------------------------
        asasalttemp = vscapped * timesolveV + ownship.alt

        # Ceiling clamp: never exceed service ceiling
        asasalttemp = np.minimum(asasalttemp, hmax - ceil_margin_m)
        # Floor clamp: never go below 0
        asasalttemp = np.maximum(asasalttemp, 0.0)

        signdvs = np.sign(
            vscapped
            - ownship.ap.vs * np.sign(ownship.selalt - ownship.alt))
        signalt = np.sign(asasalttemp - ownship.selalt)
        alt = np.where(
            np.logical_or(signdvs == 0, signdvs == signalt),
            asasalttemp, ownship.selalt)

        altCondition = np.logical_and(
            timesolveV < conf.dtlookahead, np.abs(dv[2, :]) > 0.0)
        alt[altCondition] = asasalttemp[altCondition]

        # Horizontal-only: keep AP altitude
        if self.nat_domain in ('VERT_FIRST', 'HORIZ_FIRST'):
            # Per-aircraft: keep AP alt for aircraft that used horizontal
            horiz_only = has_dv_h & ~has_dv_v
            alt = np.where(horiz_only, ownship.selalt, alt)
        else:
            alt = alt * (1 - self.swresohoriz) + ownship.selalt * self.swresohoriz

        # Final clamps
        alt = np.minimum(alt, hmax - ceil_margin_m)
        alt = np.maximum(alt, 0.0)

        # --- DEBUG: final output per aircraft in conflict ---
        if _dbg_count <= _DBG_MAX:
            for i in range(ntraf):
                if timesolveV[i] < 1e8:
                    print(f'[MVPNAT-OUT] t={bs.sim.simt:.1f} ac={ownship.id[i]} '
                          f'dv2={dv[2,i]:.4f} newvs={newvs[i]/fpm:.1f}fpm '
                          f'vscapped={vscapped[i]/fpm:.1f}fpm '
                          f'tsolV={timesolveV[i]:.2f} '
                          f'asasalt={asasalttemp[i]/ft:.0f}ft '
                          f'alt_out={alt[i]/ft:.0f}ft '
                          f'curr_alt={ownship.alt[i]/ft:.0f}ft '
                          f'selalt={ownship.selalt[i]/ft:.0f}ft '
                          f'swresovert={self.swresovert}')

        return newtrack, allowed_tas, vscapped, alt

    # ==================================================================
    #  resumenav() — Schaberg 2-criteria recovery (NAT adaptation)
    # ==================================================================
    def resumenav(self, conf, ownship, intruder):
        """Decide per aircraft whether ASAS should remain active.

        Replaces the stock 1-criterion past-CPA check with the Schaberg
        2-criteria algorithm (Schaberg, 2021) adapted for NAT:

        **Criterion 1** — *intruder keeps current resolution velocity*:
            Vrel = Vi_current - Vo_desired
            |Dcpa(Vrel)| > RPZ × dcpa_margin

        **Criterion 2** — *intruder reverts to pre-conflict velocity*:
            Vrel = Vi_initial - Vo_desired
            |Dcpa(Vrel)| > RPZ × dcpa_margin

        Recovery is allowed only when **both** criteria pass.

        Additional NAT guards:
          - **min hold time** : ASAS stays active for at least
            ``nat_min_hold_t`` seconds after first detection.
          - **max hold time** : if > 0, forces recovery after this
            many seconds (safety valve for same-direction traffic
            where Dcpa may never exceed threshold).
          - **hor_los gate** : ASAS stays active while inside RPZ,
            regardless of criteria.
        """
        global _dbg_recov_count

        # ----------------------------------------------------------
        #  1. Capture pre-conflict state for NEW conflict pairs
        # ----------------------------------------------------------
        for conflict in conf.confpairs:
            if conflict not in self.resopairs:
                ac1, ac2 = conflict
                idx1, idx2 = bs.traf.id2idx(conflict)
                # Store pre-conflict velocity — first time only
                if ac1 not in self._preconf_state and idx1 >= 0:
                    self._preconf_state[ac1] = dict(
                        gseast=float(ownship.gseast[idx1]),
                        gsnorth=float(ownship.gsnorth[idx1]),
                        vs=float(ownship.vs[idx1]),
                        alt=float(ownship.alt[idx1]))
                if ac2 not in self._preconf_state and idx2 >= 0:
                    self._preconf_state[ac2] = dict(
                        gseast=float(intruder.gseast[idx2]),
                        gsnorth=float(intruder.gsnorth[idx2]),
                        vs=float(intruder.vs[idx2]),
                        alt=float(intruder.alt[idx2]))
                # Timestamp for hold-time logic
                self._reso_start_t[conflict] = bs.sim.simt

        # ----------------------------------------------------------
        #  2. Track all active + newly-detected pairs
        # ----------------------------------------------------------
        self.resopairs.update(conf.confpairs)

        delpairs = set()
        changeactive = dict()

        for conflict in self.resopairs:
            ac1, ac2 = conflict
            idx1, idx2 = bs.traf.id2idx(conflict)

            # Ownship deleted → remove pair
            if idx1 < 0:
                delpairs.add(conflict)
                self._reso_start_t.pop(conflict, None)
                continue

            # Intruder deleted → safe to recover
            if idx2 < 0:
                changeactive[idx1] = changeactive.get(idx1, False)
                delpairs.add(conflict)
                self._reso_start_t.pop(conflict, None)
                continue

            # ---- Hold time gates ------------------------------------
            t_start = self._reso_start_t.get(conflict, bs.sim.simt)
            dt_active = bs.sim.simt - t_start

            # Min hold: too early to check criteria
            if dt_active < self.nat_min_hold_t:
                changeactive[idx1] = True
                continue

            # Max hold: force recovery (safety valve)
            if self.nat_max_hold_t > 0 and dt_active > self.nat_max_hold_t:
                changeactive[idx1] = changeactive.get(idx1, False)
                delpairs.add(conflict)
                self._reso_start_t.pop(conflict, None)
                if _dbg_recov_count < _DBG_MAX:
                    _dbg_recov_count += 1
                    print(f'[MVPNAT-RECOV] t={bs.sim.simt:.1f} '
                          f'MAXHOLD pair=({ac1},{ac2}) '
                          f'dt={dt_active:.0f}s')
                continue

            # ---- Current horizontal geometry ------------------------
            re = 6371000.0
            dist = re * np.array([
                np.radians(intruder.lon[idx2] - ownship.lon[idx1])
                * np.cos(0.5 * np.radians(
                    intruder.lat[idx2] + ownship.lat[idx1])),
                np.radians(intruder.lat[idx2] - ownship.lat[idx1])])
            hdist = float(np.linalg.norm(dist))
            rpz = float(np.max(conf.rpz[[idx1, idx2]]))
            dcpa_thr = rpz * self.nat_dcpa_margin

            # Immediate safety: still inside horizontal PZ → stay active
            if hdist < rpz:
                changeactive[idx1] = True
                continue

            # ---- Retrieve pre-conflict velocities -------------------
            pre1 = self._preconf_state.get(ac1)
            pre2 = self._preconf_state.get(ac2)

            if pre1 is None or pre2 is None:
                # Fallback: stock past-CPA check
                vrel = np.array([
                    intruder.gseast[idx2] - ownship.gseast[idx1],
                    intruder.gsnorth[idx2] - ownship.gsnorth[idx1]])
                if np.dot(dist, vrel) > 0.0:
                    changeactive[idx1] = changeactive.get(idx1, False)
                    delpairs.add(conflict)
                    self._reso_start_t.pop(conflict, None)
                else:
                    changeactive[idx1] = True
                continue

            # ---- Criterion 1: intruder keeps resolution velocity ----
            vo_des = np.array([pre1['gseast'], pre1['gsnorth']])
            vi_cur = np.array([intruder.gseast[idx2],
                               intruder.gsnorth[idx2]])
            vrel1  = vi_cur - vo_des
            dcpa1  = _compute_dcpa_h(dist, vrel1)

            # ---- Criterion 2: intruder reverts to initial velocity --
            vi_init = np.array([pre2['gseast'], pre2['gsnorth']])
            vrel2   = vi_init - vo_des
            dcpa2   = _compute_dcpa_h(dist, vrel2)

            # ---- Combined decision ----------------------------------
            safe = (dcpa1 > dcpa_thr) and (dcpa2 > dcpa_thr)

            if safe:
                changeactive[idx1] = changeactive.get(idx1, False)
                delpairs.add(conflict)
                self._reso_start_t.pop(conflict, None)
            else:
                changeactive[idx1] = True

            # Debug
            if _dbg_recov_count < _DBG_MAX:
                _dbg_recov_count += 1
                tag = 'RECOVERED' if safe else 'HOLD'
                print(f'[MVPNAT-RECOV] t={bs.sim.simt:.1f} {tag} '
                      f'pair=({ac1},{ac2}) '
                      f'dcpa1={dcpa1/nm:.2f}NM '
                      f'dcpa2={dcpa2/nm:.2f}NM '
                      f'thr={dcpa_thr/nm:.2f}NM '
                      f'hdist={hdist/nm:.2f}NM '
                      f'dt={dt_active:.0f}s')

        # ----------------------------------------------------------
        #  3. Apply active flags & clean up
        # ----------------------------------------------------------
        for idx, active in changeactive.items():
            self.active[idx] = active
            if not active:
                # Clean up pre-conflict state for this aircraft
                acid = bs.traf.id[idx]
                self._preconf_state.pop(acid, None)
                # Waypoint recovery: go to next active waypoint
                iwpid = bs.traf.ap.route[idx].findact(idx)
                if iwpid != -1:
                    bs.traf.ap.route[idx].direct(
                        idx, bs.traf.ap.route[idx].wpname[iwpid])

        self.resopairs -= delpairs

    # ==================================================================
    #  MVP_NAT — single conflict pair with symmetry breaker
    # ==================================================================
    def MVP_NAT(self, ownship, intruder, conf, qdr, dist, tcpa, tLOS,
                idx1, idx2):
        """MVP for one conflict pair with NAT symmetry breaker.

        When ``|vrel_vert| < VREL_THRESHOLD`` (both at ~same FL), the
        original MVP gives both aircraft the *same* dv3 sign → both
        descend.  This variant deterministically assigns:

          - Lower-index aircraft → descend (dv3 > 0, because
            ``dv[idx] -= dv_mvp`` → effective VS = -dv3 < 0)
          - Higher-index aircraft → climb (gets +dv_mvp from its own
            loop iteration where it is "ownship")
        """
        rpz_m  = np.max(conf.rpz[[idx1, idx2]] * self.resofach)
        hpz_m  = np.max(conf.hpz[[idx1, idx2]] * self.resofacv)
        dtlook = conf.dtlookahead[idx1]

        qdr_rad = np.radians(qdr)

        drel = np.array([np.sin(qdr_rad) * dist,
                         np.cos(qdr_rad) * dist,
                         intruder.alt[idx2] - ownship.alt[idx1]])

        v1 = np.array([ownship.gseast[idx1], ownship.gsnorth[idx1],
                        ownship.vs[idx1]])
        v2 = np.array([intruder.gseast[idx2], intruder.gsnorth[idx2],
                        intruder.vs[idx2]])
        vrel = v2 - v1

        # === Horizontal resolution ====================================
        dcpa = drel + vrel * tcpa
        dabsH = np.sqrt(dcpa[0]**2 + dcpa[1]**2)
        iH = rpz_m - dabsH

        if dabsH <= 10.:
            dabsH = 10.
            dcpa[0] = drel[1] / dist * dabsH
            dcpa[1] = -drel[0] / dist * dabsH

        if rpz_m < dist and dabsH < dist:
            erratum = np.cos(np.arcsin(rpz_m / dist)
                             - np.arcsin(dabsH / dist))
            dv1 = ((rpz_m / erratum - dabsH) * dcpa[0]) / (abs(tcpa) * dabsH)
            dv2 = ((rpz_m / erratum - dabsH) * dcpa[1]) / (abs(tcpa) * dabsH)
        else:
            dv1 = (iH * dcpa[0]) / (abs(tcpa) * dabsH)
            dv2 = (iH * dcpa[1]) / (abs(tcpa) * dabsH)

        # === Vertical resolution (with symmetry breaker) ===============
        VREL_THRESHOLD = 0.1  # m/s (~20 fpm) — below this, treat as "same FL"

        if abs(vrel[2]) > VREL_THRESHOLD:
            # --- Normal case: meaningful relative vertical speed ---
            iV = hpz_m
            tsolV = abs(drel[2] / vrel[2])
            if tsolV > dtlook:
                tsolV = tLOS
                iV = hpz_m
            dv3 = (iV / tsolV) * (-vrel[2] / abs(vrel[2]))
            if _dbg_count < _DBG_MAX:
                print(f'[MVP_NAT] NORMAL path: vrel[2]={vrel[2]:.6f} drel[2]={drel[2]:.3f} '
                      f'tsolV={tsolV:.2f} dv3={dv3:.4f}({dv3/fpm:.1f}fpm) idx=({idx1},{idx2})')
        else:
            # --- SYMMETRY BREAKER ---
            # Both at ~same FL (|vrel_vert| < threshold).
            #
            # Assign climb/descend deterministically by index:
            #   Lower-index  --> dv3 > 0 --> effective = -dv3 (descend)
            #   Higher-index --> dv3 < 0 --> effective = +|dv3| (climb)
            #
            # Use aggressive time reference: min(|tcpa|, 60 s).
            # The old tLOS-based reference (~280 s) produced only
            # ~112 fpm per aircraft — far too weak to build 1000 ft
            # separation before CPA.  Capping at 60 s gives ~1050 fpm
            # per aircraft (close to vsmax after capping).
            iV = hpz_m                                   # full zone height
            tsolV = min(max(abs(tcpa), 1.0), 60.0)       # aggressive: <= 60 s

            magnitude = iV / tsolV
            dv3 = magnitude if idx1 < idx2 else -magnitude
            if _dbg_count < _DBG_MAX:
                print(f'[MVP_NAT] SYM-BREAK: vrel[2]={vrel[2]:.6f} drel[2]={drel[2]:.3f} '
                      f'tcpa={tcpa:.1f} tsolV={tsolV:.2f} mag={magnitude:.4f}({magnitude/fpm:.1f}fpm) '
                      f'dv3={dv3:.4f} idx=({idx1},{idx2}) idx1<idx2={idx1<idx2}')

        return np.array([dv1, dv2, dv3]), tsolV
