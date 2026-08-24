""" BlueSky loggerff_asas — extended flight-state & conflict logger
    with ASAS diagnostic columns for NAT CD&R testing.
    
    Created 2026-04
    @author: ahre_ni, Nils Ahrenhold
"""
import atexit
import numpy as np
from bluesky import core, stack, traf, sim, navdb
from bluesky.core import Entity
from bluesky.tools import datalog, geo
from bluesky.tools.aero import ft, kts, nm, fpm, vtas2cas, vtas2mach
from bluesky.tools.position import txt2pos
import bluesky as bs

# Window over which the realized vertical rate (vs_real column) is measured.
# Matches StateBasedRealVS.VS_SMOOTH_S so the logged column equals the rate
# CDMETHOD STATEBASEDREALVS feeds into its vertical conflict projection.
VS_SMOOTH_S = 20.0

# ======================================================================
# FLIGHT-STATE HEADER  (loggerff columns + ASAS diagnostic columns)
# ======================================================================
flstheader = (
    # --- loggerff core columns (same order) ---
    'simt,'
    'flightid,'
    'ac_type,'
    'flighttime,'
    'currentmass,'
    'distanceflown,'
    'totalfuel,'
    'latitude,'
    'longitude,'
    'spawntime,'
    'actualdistance2D,'
    'actualdistance3D,'
    'workdone,'
    'positivefuelflow,'
    'rawfuelflow,'
    'thrust,'
    'altitude,'
    'tas,'
    'gs,'
    'cas,'
    'mach,'
    'vs,'
    'vs_real,'              # [fpm] realized d(alt)/dt over VS_SMOOTH_S window
    'heading,'
    'track,'
    # --- ASAS diagnostic columns ---
    'asas_active,'          # bool: is ASAS controlling this aircraft?
    'asas_fallback_hdg,'    # bool: heading-fallback offset currently held (mvp2nat hybrid mode)
    'asas_trk,'             # [deg] ASAS-commanded track
    'asas_tas,'             # [kts] ASAS-commanded TAS
    'asas_vs,'              # [fpm] ASAS-commanded VS
    'asas_alt,'             # [ft]  ASAS-commanded altitude
    'ap_trk,'               # [deg] autopilot desired track (planned route)
    'ap_tas,'               # [kts] autopilot desired TAS
    'ap_vs,'                # [fpm] autopilot desired VS
    'ap_alt,'               # [ft]  autopilot desired altitude
    'pilot_trk,'            # [deg] final selected track (after aporasas)
    'pilot_tas,'            # [kts] final selected TAS
    'pilot_vs,'             # [fpm] final selected VS
    'pilot_alt,'            # [ft]  final selected altitude
    'selalt,'               # [ft]  FMS selected altitude
    'n_conflicts,'          # [-]   number of conflicts this aircraft is in
    # --- performance margin columns ---
    'perf_vmin,'            # [kts] minimum speed (from perf model)
    'perf_vmax,'            # [kts] maximum speed (from perf model)
    'perf_vsmin,'           # [fpm] minimum VS (from perf model)
    'perf_vsmax,'           # [fpm] maximum VS (from perf model)
    # --- mode flags ---
    'swlnav,'
    'swvnav,'
    'swvnavspd,'
    # --- wind ---
    'gsnorth,'
    'gseast,'
    'windnorth,'
    'windeast,'
    'headwind,'
    'crosswind,'
    # --- atmosphere ---
    'temp,'
    'rho,'
    'flight_phase,'
    'drag\n'
)

# ======================================================================
# CONFLICT LOG HEADER  (loggerff columns + ASAS diagnostic columns)
# ======================================================================
confheader = (
    # --- loggerff core columns ---
    'simt[s],'
    'ac1,'
    'ac2,'
    'latitude_ac1[deg],'
    'longitude_ac1[deg],'
    'altitude_ac1[ft],'
    'latitude_ac2[deg],'
    'longitude_ac2[deg],'
    'altitude_ac2[ft],'
    'heading_ac1[deg],'
    'heading_ac2[deg],'
    'vs_ac1[fpm],'
    'vs_ac2[fpm],'
    'vs_real_ac1[fpm],'    # [fpm] windowed d(alt)/dt at conflict start
    'vs_real_ac2[fpm],'
    'dcpa[nm],'
    'tcpa[sec],'
    'tLOS[sec],'
    'qdr[deg],'
    'dist[nm],'
    'dalt_min[ft],'
    'tinconf[sec],'
    'toutconf[sec],'
    'duration[s],'
    'intrusion,'
    # --- ASAS diagnostic columns ---
    'conflict_angle[deg],'  # relative heading angle at conflict start
    'asas_active_ac1,'      # bool: was ASAS active for ac1 at end?
    'asas_active_ac2,'      # bool: was ASAS active for ac2 at end?
    'mach_ac1,'             # [-] Mach of ac1 at conflict start
    'mach_ac2,'             # [-] Mach of ac2 at conflict start
    'selalt_ac1[ft],'       # [ft] FMS alt target ac1 at conflict start
    'selalt_ac2[ft],'       # [ft] FMS alt target ac2 at conflict start
    'tas_ac1[kts],'         # [kts] TAS of ac1 at conflict start
    'tas_ac2[kts],'         # [kts] TAS of ac2 at conflict start
    'los_severity\n'        # min_dist / rpz  (0 = centre, 1 = boundary)
)


# ======================================================================
# Plugin entry point
# ======================================================================
def init_plugin():
    """Plugin initialisation function."""
    global loggerff_asas
    loggerff_asas = LoggerffAsas()

    config = {
        'plugin_name': 'LOGGERFF_ASAS',
        'plugin_type': 'sim',
    }

    stackfunctions = {
        'STARTLOG': [
            'STARTLOG',
            '',
            loggerff_asas.start_log,
            'Start the extended ASAS flight-state and conflict logger'
        ],
        'SETMASS': [
            'SETMASS acid,mass',
            '[acid,float]',
            loggerff_asas.setmass,
            'Set aircraft mass [kg] after creation'
        ]
    }

    return config, stackfunctions


# ======================================================================
# Logger entity
# ======================================================================
class LoggerffAsas(Entity):
    """Extended logger with ASAS diagnostic columns for NAT CD&R testing."""

    def __init__(self):
        super().__init__()

        # Destination proximity for auto-delete [nm]
        self.d = 10

        # --- Conflict tracking dicts (same as loggerff) ---
        self.duration = {}
        self.init_lat1 = {}
        self.init_lon1 = {}
        self.init_alt1 = {}
        self.init_hdg1 = {}
        self.init_vs1  = {}
        self.init_vsreal1 = {}
        self.init_lat2 = {}
        self.init_lon2 = {}
        self.init_alt2 = {}
        self.init_hdg2 = {}
        self.init_vs2  = {}
        self.init_vsreal2 = {}
        self.intrusion_occurred = {}

        # Severity & timing
        self.dcpa     = {}
        self.dalt     = {}
        self.tLOS     = {}
        self.dist     = {}
        self.qdr      = {}
        self.tcpa     = {}
        self.tinconf  = {}
        self.toutconf = {}

        # --- NEW: ASAS diagnostic dicts per conflict pair ---
        self.init_conflict_angle = {}   # relative heading angle at conflict start
        self.init_mach1  = {}           # Mach at conflict start
        self.init_mach2  = {}
        self.init_selalt1 = {}          # FMS selected alt at conflict start
        self.init_selalt2 = {}
        self.init_tas1   = {}           # TAS at conflict start
        self.init_tas2   = {}

        # --- Realized vertical rate (vs_real column) ---
        # traf.vs is the integrator's instantaneous VS state, and it reads 0
        # whenever traffic.py takes the swaltsel "snap" branch (alt is assigned
        # from aporasas.alt instead of integrated from vs).  On shallow VNAV
        # cruise climbs that is most ticks, so traf.vs is not the derivative of
        # the logged altitude.  vs_real is that derivative, computed the same
        # way as StateBasedRealVS._real_vs (asas_nat/statebased_nat.py) so the
        # column is directly comparable to what CDMETHOD STATEBASEDREALVS uses.
        self._vs_alt_ref = {}           # acid -> (t_ref [s], alt_ref [m])
        self._vs_real    = {}           # acid -> last computed real vs [m/s]

        self.sim_name = stack.get_scenname()

        # Log files
        self.flst    = datalog.crelog('FLSTLOG_ASAS', None, flstheader)
        self.conflog = datalog.crelog('CONFLOG_ASAS', None, confheader)

        # Last sim time seen by update(). The RESET path flushes AFTER
        # Simulation.reset() has zeroed sim.simt (simulation.py:206 runs before
        # bs.traf.reset() at :216, which is what triggers our reset()), so
        # sim.simt is useless as an end-of-conflict stamp there. This is.
        self._last_simt = 0.0

        # Hook into conflog.reset() so active conflicts are flushed
        # before the CSV file is closed.  This covers both RESET and
        # QUIT paths (datalog.reset() calls conflog.reset()).
        #
        # Installed ONCE per logger, not once per LoggerffAsas. crelog returns
        # the EXISTING logger when the name is already registered
        # (datalog.py:45), so a second __init__ would otherwise wrap its own
        # wrapper, nesting closures that each hold a stale `self`. The wrapper
        # delegates to whichever instance currently owns the logger instead.
        conflog = self.conflog
        conflog._flush_owner = self
        if not getattr(conflog, '_flush_hooked', False):
            _original_conflog_reset = conflog.reset

            def _conflog_reset_with_flush():
                owner = getattr(conflog, '_flush_owner', None)
                if owner is not None:
                    owner._flush_active_conflicts()
                _original_conflog_reset()

            conflog.reset = _conflog_reset_with_flush
            conflog._flush_hooked = True

            # Safety net for abnormal shutdowns. Also once per logger: an
            # atexit registration per instance would replay the same flush.
            atexit.register(
                lambda: getattr(conflog, '_flush_owner', None)
                and conflog._flush_owner._flush_active_conflicts())

        with self.settrafarrays():
            self.distance2D       = np.array([])
            self.distance3D       = np.array([])
            self.create_time      = np.array([])
            self.total_fuel       = np.array([])
            self.last_update_time = np.array([])

    # ------------------------------------------------------------------
    def _asas_flags(self, ac1, ac2):
        """ASAS-active flag per aircraft, or -1 if it has been deleted."""
        try:
            ntraf = traf.ntraf
            i1, i2 = traf.id2idx(ac1), traf.id2idx(ac2)
            return (int(traf.cr.active[i1]) if 0 <= i1 < ntraf else -1,
                    int(traf.cr.active[i2]) if 0 <= i2 < ntraf else -1)
        except Exception:
            return -1, -1

    def _write_conflict_row(self, cpf, ac1, ac2, toutconf, duration):
        """Emit one CONFLOG row for a finished conflict event.

        Shared by the end-of-conflict path and by the end-of-simulation flush.
        The two used to carry 34 duplicated arguments each, which is how they
        drifted apart in the first place.

        Every field comes from the per-event dicts, never from live traffic, so
        a row is complete even when both aircraft have already been deleted.
        """
        try:
            rpz_val = bs.traf.cd.rpz_def if hasattr(bs.traf.cd, 'rpz_def') else (5.0 * nm)
        except Exception:
            rpz_val = 5.0 * nm
        los_sev = self.dist.get(cpf, rpz_val) / rpz_val if rpz_val > 0 else -1
        asas1, asas2 = self._asas_flags(ac1, ac2)

        self.conflog.log(
            ac1, ac2,
            self.init_lat1.get(cpf, np.nan),
            self.init_lon1.get(cpf, np.nan),
            self.init_alt1.get(cpf, np.nan) / ft,
            self.init_lat2.get(cpf, np.nan),
            self.init_lon2.get(cpf, np.nan),
            self.init_alt2.get(cpf, np.nan) / ft,
            self.init_hdg1.get(cpf, np.nan),
            self.init_hdg2.get(cpf, np.nan),
            self.init_vs1.get(cpf, np.nan),
            self.init_vs2.get(cpf, np.nan),
            self.init_vsreal1.get(cpf, np.nan),
            self.init_vsreal2.get(cpf, np.nan),
            self.dcpa.get(cpf, np.nan) / nm,
            self.tcpa.get(cpf, np.nan),
            self.tLOS.get(cpf, -1),
            self.qdr.get(cpf, np.nan),
            self.dist.get(cpf, np.nan) / nm,
            self.dalt.get(cpf, np.nan) / ft,
            self.tinconf.get(cpf, np.nan),
            toutconf,
            duration,
            int(self.intrusion_occurred.get(cpf, False)),
            # --- ASAS diagnostic columns ---
            self.init_conflict_angle.get(cpf, np.nan),
            asas1,
            asas2,
            self.init_mach1.get(cpf, np.nan),
            self.init_mach2.get(cpf, np.nan),
            self.init_selalt1.get(cpf, np.nan),
            self.init_selalt2.get(cpf, np.nan),
            self.init_tas1.get(cpf, np.nan),
            self.init_tas2.get(cpf, np.nan),
            los_sev,
        )

    def _close_conflict(self, cpf, toutconf):
        """Log a finished event and drop the state it accumulated."""
        # self.duration is keyed by frozenset, so this pair order is arbitrary.
        # Harmless -- the postprocessing dedup key is the SORTED pair -- and
        # sorting here would change the output of every existing run.
        ac1, ac2 = tuple(cpf)
        self._write_conflict_row(cpf, ac1, ac2, toutconf, self.duration[cpf])
        del self.duration[cpf]
        self.toutconf[cpf] = toutconf
        # The event is over: drop its accumulated state so a re-detection of
        # this pair starts clean (2026-08-20).
        self._clear_conflict_state(cpf)

    def _log_ended_conflicts(self, confpairs_unique, simt):
        """Advance every tracked pair: tick its duration, or close it.

        Split out of update() for BUG-02 (2026-08-24). It used to be inline
        AFTER the `if traf.ntraf == 0: return` guard, so the tick that deleted
        the last aircraft was also the tick that stopped conflicts ever being
        closed -- and any pair still open was orphaned in self.duration.
        """
        for cpf in list(self.duration):
            if cpf in confpairs_unique:
                self.duration[cpf] += 1
            else:
                self._close_conflict(cpf, simt)

    # ------------------------------------------------------------------
    def _flush_active_conflicts(self):
        """Write all still-active conflicts to conflog before reset/quit.

        The net for the paths update() cannot cover: HOLD, QUIT, and a hard
        shutdown. With BUG-02 fixed, the ordinary end-of-simulation case is
        closed by _log_ended_conflicts() at the real sim time instead, so this
        should now write nothing on a normal run.

        Two traps this has to work around, both outside the plugin:

        * `sim.simt` is already 0.0 here on the RESET path -- Simulation.reset()
          zeroes it (simulation.py:206) before bs.traf.reset() (:216) triggers
          our reset(). So the end stamp comes from self._last_simt.
        * CSVLogger.log is gated on `bs.sim.simt >= self.tlog` (datalog.py:170),
          and tlog is the STARTLOG time -- 30 s in the paper3 batches. With simt
          zeroed the gate is False, there is no else branch and no warning, and
          every flushed row is discarded in silence. That gate is meant for
          PERIODIC loggers; CONFLOG is an event logger (dt=None -> 0.0). We drop
          tlog around the write and restore it, rather than patching core
          BlueSky shared with upstream.

        Safe to call repeatedly -- self.duration is emptied, so a second call
        (e.g. atexit after reset) writes nothing.
        """
        if not self.duration:
            return  # nothing to flush

        # Guard: conflog file must still be open
        if not (hasattr(self, 'conflog') and self.conflog.isopen()):
            return

        toutconf = getattr(self, '_last_simt', 0.0)
        saved_tlog = getattr(self.conflog, 'tlog', None)
        n_flushed = 0
        try:
            if saved_tlog is not None:
                self.conflog.tlog = float('-inf')
            for cpf in list(self.duration):
                self._close_conflict(cpf, toutconf)
                n_flushed += 1
        finally:
            if saved_tlog is not None:
                self.conflog.tlog = saved_tlog

        print(f"LOGGERFF_ASAS: Flushed {n_flushed} active "
              f"conflict(s) to conflog at simt={toutconf:.1f}s")

    # ------------------------------------------------------------------
    def reset(self):
        # Flush conflicts that are still active before clearing everything
        self._flush_active_conflicts()

        super().reset()
        self.duration.clear()
        self.init_lat1.clear(); self.init_lon1.clear(); self.init_alt1.clear()
        self.init_hdg1.clear(); self.init_vs1.clear(); self.init_vsreal1.clear()
        self.init_lat2.clear(); self.init_lon2.clear(); self.init_alt2.clear()
        self.init_hdg2.clear(); self.init_vs2.clear(); self.init_vsreal2.clear()
        self.intrusion_occurred.clear()
        self.dcpa.clear(); self.dalt.clear(); self.tLOS.clear()
        self.dist.clear(); self.qdr.clear(); self.tcpa.clear()
        self.tinconf.clear(); self.toutconf.clear()
        self.init_conflict_angle.clear()
        self.init_mach1.clear(); self.init_mach2.clear()
        self.init_selalt1.clear(); self.init_selalt2.clear()
        self.init_tas1.clear(); self.init_tas2.clear()
        self._vs_alt_ref.clear(); self._vs_real.clear()
        self.sim_name = None

    # ------------------------------------------------------------------
    def create(self, n=1):
        super().create(n)
        self.create_time[-n:]      = sim.simt
        self.total_fuel[-n:]       = 0.0
        self.last_update_time[-n:] = sim.simt

    # ------------------------------------------------------------------
    def _real_vs(self, n):
        """Realized vertical rate [m/s] for the first *n* aircraft.

        (alt(t) - alt(t-dt)) / dt over a ~VS_SMOOTH_S window, held between
        updates.  Keyed by callsign so it survives index compaction when
        traf.delete() removes an aircraft.  Mirrors
        StateBasedRealVS._real_vs in asas_nat/statebased_nat.py.
        """
        t   = float(sim.simt)
        alt = np.asarray(traf.alt[:n], dtype=float)
        out = np.zeros(n, dtype=float)
        for k in range(n):
            acid = traf.id[k]
            ref  = self._vs_alt_ref.get(acid)
            if ref is None:
                self._vs_alt_ref[acid] = (t, alt[k])
                out[k] = self._vs_real.get(acid, 0.0)
                continue
            t0, a0 = ref
            dt_ref = t - t0
            if dt_ref >= VS_SMOOTH_S:
                v = (alt[k] - a0) / dt_ref
                self._vs_real[acid]    = v
                self._vs_alt_ref[acid] = (t, alt[k])
                out[k] = v
            else:
                out[k] = self._vs_real.get(acid, 0.0)
        return out

    # ------------------------------------------------------------------
    # Per-EVENT conflict state
    #
    # A conflict EVENT is one contiguous span of confpairs membership for a
    # pair. The same pair can produce several events in a run, so the state an
    # event accumulates must be scoped to that event -- see
    # _clear_conflict_state for what happens when it is not.
    # ------------------------------------------------------------------
    def _note_conflict_state(self, cpf, dist_now, simt, in_los):
        """Accumulate one tick of the current event for pair `cpf`.

        tinconf and tLOS are write-once (the FIRST tick, and the first tick in
        LoS); dist is the event MINIMUM; intrusion_occurred latches True.
        """
        if cpf not in self.dist:
            self.dist[cpf] = dist_now
        else:
            self.dist[cpf] = min(self.dist[cpf], dist_now)
        if cpf not in self.tinconf:
            self.tinconf[cpf] = simt
        # Intrusion tracking
        if in_los:
            self.intrusion_occurred[cpf] = True
            if cpf not in self.tLOS:
                self.tLOS[cpf] = simt
        elif cpf not in self.intrusion_occurred:
            self.intrusion_occurred[cpf] = False

    def _clear_conflict_state(self, cpf):
        """Discard the per-event state of a conflict that has just ended.

        Every value _note_conflict_state writes is guarded write-once or
        accumulated, so without this a pair re-entering confpairs logs its
        SECOND row carrying the FIRST event's tinconf, an inherited intrusion
        flag, and a dist minimised across both events. Downstream that corrupts
        n_conflicts (the postprocessing dedup key is (pair, tinconf), so
        re-detections collapse into one row), n_intrusions, and
        gap_to_prev_resolution_s (which goes negative for exactly the
        re-detections it measures). Fixed 2026-08-20; logs written before that
        date still carry the defect.

        Deliberately NOT cleared:
          toutconf            persists across events on purpose -- the
                              re-detection gap is tinconf_k - toutconf_{k-1}.
          dcpa/dalt/tcpa/qdr  re-assigned every tick the pair is in confpairs,
                              so they cannot go stale within an event.
          init_*              already per-event: their capture is guarded by
                              `up not in self.duration`, and duration IS
                              deleted at the end of every event.
        """
        self.tinconf.pop(cpf, None)
        self.dist.pop(cpf, None)
        self.tLOS.pop(cpf, None)
        self.intrusion_occurred.pop(cpf, None)

    # ------------------------------------------------------------------
    @core.timed_function(name='LOGGERFF_ASAS', dt=1.0)
    def update(self, dt):
        """
        1. Aircraft deletion (destination proximity + zombies)
        2. Early-exit check
        3. Distance & fuel integration
        4. Conflict processing (extended)
        5. Flight-state logging (extended)
        """
        # The last sim time we actually saw. The reset-path flush runs after
        # Simulation.reset() has zeroed sim.simt, so this is the only honest
        # end-of-conflict stamp available there.
        self._last_simt = sim.simt

        # ==============================================================
        # 1. AIRCRAFT DELETION
        # ==============================================================
        if traf.ntraf > 0:
            lat_ac = np.array(traf.lat)
            lon_ac = np.array(traf.lon)
            dest_names = np.array(traf.ap.dest)

            lat_dest   = np.zeros(traf.ntraf)
            lon_dest   = np.zeros(traf.ntraf)
            valid_dest = np.ones(traf.ntraf, dtype=bool)

            for idx, apname in enumerate(dest_names):
                apidx = bs.navdb.getaptidx(apname)
                if apidx < 0:
                    if bs.traf.ap.route[idx].nwp > 0:
                        lat_dest[idx] = bs.traf.ap.route[idx].wplat[-1]
                        lon_dest[idx] = bs.traf.ap.route[idx].wplon[-1]
                    else:
                        lat_dest[idx] = traf.lat[idx]
                        lon_dest[idx] = traf.lon[idx]
                    success, posobj = txt2pos(apname, lat_dest[idx], lon_dest[idx])
                    if success:
                        lat_dest[idx] = posobj.lat
                        lon_dest[idx] = posobj.lon
                    else:
                        valid_dest[idx] = False
                else:
                    lat_dest[idx] = bs.navdb.aptlat[apidx]
                    lon_dest[idx] = bs.navdb.aptlon[apidx]

            _, distances = geo.qdrdist(
                lat_ac[valid_dest], lon_ac[valid_dest],
                lat_dest[valid_dest], lon_dest[valid_dest])
            to_delete_valid = np.where(distances < self.d)[0]
            to_delete = np.nonzero(valid_dest)[0][to_delete_valid]

            for idx in sorted(to_delete, reverse=True):
                cs = traf.id[idx]
                traf.delete(idx)
                print(f"LOGGERFF_ASAS - {self.sim_name}: {cs} landed "
                      f"at {sim.simt}; active: {traf.ntraf}")

            # Zombie cleanup (>24 h or invalid lat)
            if traf.ntraf > 0:
                flight_ages = sim.simt - self.create_time[:traf.ntraf]
                max_s = 24.0 * 3600.0
                zombie = (flight_ages > max_s) | (np.abs(traf.lat[:traf.ntraf]) > 90.0)
                for idx in sorted(np.where(zombie)[0], reverse=True):
                    cs = traf.id[idx]
                    traf.delete(idx)
                    print(f"LOGGERFF_ASAS - {self.sim_name}: ZOMBIE {cs} deleted; "
                          f"active: {traf.ntraf}")

        # ==============================================================
        # 2. EARLY-EXIT
        # ==============================================================
        if sim.simt >= 48 * 3600 and traf.ntraf == 0:
            name = self.sim_name or "UNNAMED"
            print(f"END of simulation: {name} at {sim.simt}s = {sim.simt/3600:.1f}h")
            stack.stack('reset')
            return

        if traf.ntraf == 0:
            # BUG-02. No traffic means no confpairs, so every pair still in
            # self.duration has ended by definition -- close them HERE, at the
            # real sim time. This guard used to sit in front of the
            # end-of-conflict block, so the tick that deleted the last aircraft
            # was also the tick that stopped conflicts ever being closed, and
            # an open pair was orphaned until the reset-path flush (which then
            # stamped toutconf = 0 and was silently discarded by the periodic
            # logger gate). Measured: 9 of 175 synthetic runs logged ZERO
            # conflict rows against up to 4739 s of real loss of separation.
            #
            # An EMPTY set, deliberately, not traf.cd.confpairs_unique:
            # Traffic.update() returns at `if self.ntraf == 0` before reaching
            # cd.update() (traffic.py:394,407), so confpairs_unique is stale
            # here -- it still holds the pairs from the last tick that had
            # traffic. Passing it would leave every pair looking "still in
            # conflict" and the bug unfixed.
            self._log_ended_conflicts(set(), sim.simt)
            return

        n = traf.ntraf

        # ==============================================================
        # 3. DISTANCE & FUEL INTEGRATION
        # ==============================================================
        resultantspd = np.sqrt(traf.gs[:n]**2 + traf.vs[:n]**2)
        self.distance2D[:n] += dt * traf.gs[:n]
        self.distance3D[:n] += dt * resultantspd

        try:
            raw_fuelflow = traf.perf.fuelflow[:n]
            pos_fuelflow = np.maximum(0, raw_fuelflow)
            self.total_fuel[:n] += pos_fuelflow * dt
            self.last_update_time[:n] = sim.simt
        except (AttributeError, IndexError):
            raw_fuelflow = np.full(n, np.nan)
            pos_fuelflow = np.zeros(n)

        # Thrust
        thrust = np.zeros(n)
        if hasattr(traf.perf, 'thrust_effective') and len(traf.perf.thrust_effective) >= n:
            thrust = traf.perf.thrust_effective[:n]
        elif hasattr(traf.perf, 'thrust') and len(traf.perf.thrust) >= n:
            thrust = traf.perf.thrust[:n]

        # Current mass
        current_mass = np.zeros(n)
        if hasattr(traf.perf, 'mass') and len(traf.perf.mass) >= n:
            current_mass = traf.perf.mass[:n]

        # Realized vertical rate. Computed BEFORE conflict processing so a
        # conflict starting on this tick records the current window value.
        vs_real = self._real_vs(n)

        # ==============================================================
        # 4. CONFLICT PROCESSING  (with ASAS diagnostic extras)
        # ==============================================================
        # Build index map: frozenset(pair) → position in confpairs list
        idxdict = {frozenset(v): i for i, v in enumerate(traf.cd.confpairs)}

        # --- Count per-aircraft conflicts ---
        ac_nconf = np.zeros(n, dtype=int)
        for pair in traf.cd.confpairs:
            i1 = traf.id2idx(pair[0])
            if 0 <= i1 < n:
                ac_nconf[i1] += 1

        # --- Update severity values for ongoing conflicts ---
        for cpf in traf.cd.confpairs_unique:
            if cpf in idxdict:
                i = idxdict[cpf]
                self.dcpa[cpf] = np.asarray(traf.cd.dcpa)[i]
                self.dalt[cpf] = np.asarray(traf.cd.dalt)[i]
                self.tcpa[cpf] = np.asarray(traf.cd.tcpa)[i]
                self.qdr[cpf]  = np.asarray(traf.cd.qdr)[i]
                # dcpa/dalt/tcpa/qdr above are per-TICK (last value wins);
                # everything below is per-EVENT and is discarded when the
                # event ends -- see _note_conflict_state/_clear_conflict_state.
                self._note_conflict_state(
                    cpf,
                    dist_now=np.asarray(traf.cd.dist)[i],
                    simt=sim.simt,
                    in_los=cpf in traf.cd.lospairs_unique,
                )

        # --- New conflict pair detection & initial state capture ---
        for pair in traf.cd.confpairs:
            up = frozenset(pair)
            if up not in self.duration:
                self.duration[up] = 1
                ac1, ac2 = tuple(pair)
                i1, i2 = traf.id2idx(ac1), traf.id2idx(ac2)
                if i1 < 0 or i2 < 0:
                    continue

                # Standard initial-state capture
                self.init_lat1[up] = traf.lat[i1]
                self.init_lon1[up] = traf.lon[i1]
                self.init_alt1[up] = traf.alt[i1]
                self.init_hdg1[up] = traf.hdg[i1]
                self.init_vs1[up]  = traf.vs[i1] / fpm
                self.init_vsreal1[up] = self._vs_real.get(traf.id[i1], 0.0) / fpm
                self.init_lat2[up] = traf.lat[i2]
                self.init_lon2[up] = traf.lon[i2]
                self.init_alt2[up] = traf.alt[i2]
                self.init_hdg2[up] = traf.hdg[i2]
                self.init_vs2[up]  = traf.vs[i2] / fpm
                self.init_vsreal2[up] = self._vs_real.get(traf.id[i2], 0.0) / fpm

                # NEW: conflict geometry at start
                dhdg = (traf.hdg[i1] - traf.hdg[i2] + 180) % 360 - 180
                self.init_conflict_angle[up] = abs(dhdg)
                self.init_mach1[up]   = traf.M[i1]
                self.init_mach2[up]   = traf.M[i2]
                self.init_selalt1[up] = traf.selalt[i1] / ft
                self.init_selalt2[up] = traf.selalt[i2] / ft
                self.init_tas1[up]    = traf.tas[i1] / kts
                self.init_tas2[up]    = traf.tas[i2] / kts

        # --- Duration update & end-of-conflict logging ---
        self._log_ended_conflicts(traf.cd.confpairs_unique, sim.simt)

        # ==============================================================
        # 5. FLIGHT-STATE LOGGING  (extended)
        # ==============================================================
        if n == 0:
            return

        # --- Pre-compute derived quantities ---
        trk_rad  = np.radians(traf.trk[:n])
        wn       = traf.windnorth[:n]
        we       = traf.windeast[:n]
        headwind  = -(wn * np.cos(trk_rad) + we * np.sin(trk_rad))
        crosswind = -wn * np.sin(trk_rad) + we * np.cos(trk_rad)

        # Performance limits (safe access)
        perf_vmin  = traf.perf.vmin[:n] / kts  if hasattr(traf.perf, 'vmin')  and len(traf.perf.vmin)  >= n else np.full(n, np.nan)
        perf_vmax  = traf.perf.vmax[:n] / kts  if hasattr(traf.perf, 'vmax')  and len(traf.perf.vmax)  >= n else np.full(n, np.nan)
        perf_vsmin = traf.perf.vsmin[:n] / fpm if hasattr(traf.perf, 'vsmin') and len(traf.perf.vsmin) >= n else np.full(n, np.nan)
        perf_vsmax = traf.perf.vsmax[:n] / fpm if hasattr(traf.perf, 'vsmax') and len(traf.perf.vsmax) >= n else np.full(n, np.nan)

        # Mode flags
        swlnav    = traf.swlnav[:n].astype(int)   if len(traf.swlnav)    >= n else np.zeros(n, dtype=int)
        swvnav    = traf.swvnav[:n].astype(int)    if len(traf.swvnav)    >= n else np.zeros(n, dtype=int)
        swvnavspd = traf.swvnavspd[:n].astype(int) if len(traf.swvnavspd) >= n else np.zeros(n, dtype=int)

        # Atmosphere & performance
        Temp         = traf.Temp[:n]
        rho          = traf.rho[:n]
        flight_phase = traf.perf.phase[:n].astype(int) if hasattr(traf.perf, 'phase') and len(traf.perf.phase) >= n else np.zeros(n, dtype=int)
        drag         = traf.perf.D[:n] if hasattr(traf.perf, 'D') and len(traf.perf.D) >= n else np.zeros(n)

        self.flst.log(
            # --- core columns ---
            traf.id[:n],                                    # flightid
            traf.type[:n],                                  # ac_type
            sim.simt - self.create_time[:n],                # flighttime [s]
            current_mass,                                   # currentmass [kg]
            traf.distflown[:n] / nm,                        # distanceflown [nm]
            self.total_fuel[:n],                            # totalfuel [kg]
            traf.lat[:n],                                   # latitude [deg]
            traf.lon[:n],                                   # longitude [deg]
            self.create_time[:n],                           # spawntime [s]
            self.distance2D[:n],                            # actualdistance2D [m]
            self.distance3D[:n],                            # actualdistance3D [m]
            traf.work[:n] * 1e-6,                           # workdone [MJ]
            pos_fuelflow,                                   # positivefuelflow [kg/s]
            raw_fuelflow,                                   # rawfuelflow [kg/s]
            thrust,                                         # thrust [N]
            traf.alt[:n] / ft,                              # altitude [ft]
            traf.tas[:n] / kts,                             # tas [kts]
            traf.gs[:n] / kts,                              # gs [kts]
            traf.cas[:n] / kts,                             # cas [kts]
            traf.M[:n],                                     # mach [-]
            traf.vs[:n] / fpm,                              # vs [fpm]
            vs_real / fpm,                                  # vs_real [fpm]
            traf.hdg[:n],                                   # heading [deg]
            traf.trk[:n],                                   # track [deg]
            # --- ASAS diagnostic columns ---
            traf.cr.active[:n].astype(int),                 # asas_active [0/1]
            (getattr(traf.cr, 'fallback_hdg_active',
                     np.zeros(n, dtype=bool))[:n]).astype(int),  # asas_fallback_hdg [0/1]
            traf.cr.trk[:n],                                # asas_trk [deg]
            traf.cr.tas[:n] / kts,                          # asas_tas [kts]
            traf.cr.vs[:n] / fpm,                           # asas_vs [fpm]
            traf.cr.alt[:n] / ft,                           # asas_alt [ft]
            traf.ap.trk[:n],                                # ap_trk [deg]
            traf.ap.tas[:n] / kts,                          # ap_tas [kts]
            traf.ap.vs[:n] / fpm,                           # ap_vs [fpm]
            traf.ap.alt[:n] / ft,                           # ap_alt [ft]
            traf.aporasas.trk[:n],                          # pilot_trk [deg]
            traf.aporasas.tas[:n] / kts,                    # pilot_tas [kts]
            traf.aporasas.vs[:n] / fpm,                     # pilot_vs [fpm]
            traf.aporasas.alt[:n] / ft,                     # pilot_alt [ft]
            traf.selalt[:n] / ft,                           # selalt [ft]
            ac_nconf,                                       # n_conflicts [-]
            # --- performance margin ---
            perf_vmin,                                      # perf_vmin [kts]
            perf_vmax,                                      # perf_vmax [kts]
            perf_vsmin,                                     # perf_vsmin [fpm]
            perf_vsmax,                                     # perf_vsmax [fpm]
            # --- mode flags ---
            swlnav,                                         # swlnav [0/1]
            swvnav,                                         # swvnav [0/1]
            swvnavspd,                                      # swvnavspd [0/1]
            # --- wind ---
            traf.gsnorth[:n] / kts,                         # gsnorth [kts]
            traf.gseast[:n] / kts,                          # gseast [kts]
            traf.windnorth[:n] / kts,                       # windnorth [kts]
            traf.windeast[:n] / kts,                        # windeast [kts]
            headwind / kts,                                 # headwind [kts]
            crosswind / kts,                                # crosswind [kts]
            # --- atmosphere ---
            Temp,                                           # temp [K]
            rho,                                            # rho [kg/m³]
            flight_phase,                                   # flight_phase [1-6]
            drag,                                           # drag [N]
        )

    # ------------------------------------------------------------------
    def start_log(self):
        """Start both FLST and CONF logs."""
        print(f"LOGGERFF_ASAS: FLST and CONF logger started: {self.sim_name}")
        self.flst.start()
        self.conflog.start()
        return True, (f'LOGGERFF_ASAS: FLST and CONF logger is ON '
                       f'for simulation: {self.sim_name}.')

    # ------------------------------------------------------------------
    def setmass(self, idx, mass_kg):
        """Set aircraft mass [kg] in the BADA performance model."""
        if mass_kg <= 0:
            return False, f'SETMASS: Mass must be positive, got {mass_kg} kg.'
        acid = traf.id[idx]
        if not hasattr(traf.perf, 'mass') or len(traf.perf.mass) == 0:
            return False, (f'SETMASS {acid}: No performance model with mass '
                           f'data active.')
        if hasattr(traf.perf, 'mmin') and hasattr(traf.perf, 'mmax'):
            mmin = traf.perf.mmin[idx]
            mmax = traf.perf.mmax[idx]
            if mass_kg < mmin:
                return False, (f'SETMASS {acid}: Mass {mass_kg:.1f} kg below '
                               f'OEW ({mmin:.1f} kg).')
            if mass_kg > mmax:
                return False, (f'SETMASS {acid}: Mass {mass_kg:.1f} kg exceeds '
                               f'MTOW ({mmax:.1f} kg).')
        old_mass = traf.perf.mass[idx]
        traf.perf.mass[idx] = mass_kg
        print(f'SETMASS {acid}: Mass set from {old_mass:.1f} to {mass_kg:.1f} kg')
        return True, (f'SETMASS {acid}: Mass set to {mass_kg:.1f} kg '
                       f'(was {old_mass:.1f} kg).')
