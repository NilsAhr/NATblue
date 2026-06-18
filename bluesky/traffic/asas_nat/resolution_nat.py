"""
resolution_nat.py — NAT-aware conflict resolution base class.
==============================================================

Full copy of ``bluesky.traffic.asas.resolution.ConflictResolution``
so that resumenav() and the active-flag logic can be modified
independently for NAT operations.

Changes vs. original
--------------------
(none yet — exact mirror as of 2026-04-08)

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-04
"""
import numpy as np

import bluesky as bs
from bluesky.core import Entity
from bluesky.stack import command
from bluesky.tools.aero import nm, ft


bs.settings.set_variable_defaults(asas_marh=1.01, asas_marv=1.01)


class ConflictResolutionNAT(Entity, replaceable=True):
    """NAT-specific conflict resolution base class.

    Mirrors ``ConflictResolution`` so we can modify ``resumenav()``
    and the priority/active logic without touching the upstream code.
    """

    def __init__(self):
        super().__init__()
        self.swprio = False
        self.priocode = ''
        self.resopairs = set()

        self.resofach = bs.settings.asas_marh
        self.resofacv = bs.settings.asas_marv

        self.resodhrelative = True
        self.resorrelative  = True

        with self.settrafarrays():
            self.resooffac = np.array([], dtype=bool)
            self.noresoac  = np.array([], dtype=bool)
            self.active = np.array([], dtype=bool)
            self.trk = np.array([])
            self.tas = np.array([])
            self.alt = np.array([])
            self.vs  = np.array([])

    def reset(self):
        super().reset()
        self.swprio = False
        self.priocode = ''
        self.resopairs.clear()
        self.resofach = bs.settings.asas_marh
        self.resofacv = bs.settings.asas_marv
        self.resodhrelative = True
        self.resorrelative  = True

    # ----- Channel properties (can be overloaded) -----
    @property
    def hdgactive(self):
        return self.active

    @property
    def vsactive(self):
        return self.active

    @property
    def altactive(self):
        return self.active

    @property
    def tasactive(self):
        return self.active

    def resolve(self, conf, ownship, intruder):
        """Resolve all current conflicts (override in subclass)."""
        return ownship.ap.trk, ownship.ap.tas, ownship.ap.vs, ownship.ap.alt

    def update(self, conf, ownship, intruder):
        """Perform one CR update step."""
        if ConflictResolutionNAT.selected() is not ConflictResolutionNAT:
            if conf.confpairs:
                self.trk, self.tas, self.vs, self.alt = \
                    self.resolve(conf, ownship, intruder)
            self.resumenav(conf, ownship, intruder)

    # ------------------------------------------------------------------
    #  resumenav — decide when to hand control back to autopilot
    # ------------------------------------------------------------------
    def resumenav(self, conf, ownship, intruder):
        """Decide per aircraft whether ASAS or autopilot should be followed."""
        self.resopairs.update(conf.confpairs)

        delpairs = set()
        changeactive = dict()

        def anglediff(a, b):
            d = a - b
            if d > 180:
                return anglediff(a, b + 360)
            elif d < -180:
                return anglediff(a + 360, b)
            else:
                return d

        for conflict in self.resopairs:
            idx1, idx2 = bs.traf.id2idx(conflict)
            if idx1 < 0:
                delpairs.add(conflict)
                continue

            if idx2 >= 0:
                re = 6371000.
                dist = re * np.array([
                    np.radians(intruder.lon[idx2] - ownship.lon[idx1])
                    * np.cos(0.5 * np.radians(intruder.lat[idx2]
                                               + ownship.lat[idx1])),
                    np.radians(intruder.lat[idx2] - ownship.lat[idx1])])

                vrel = np.array([
                    intruder.gseast[idx2] - ownship.gseast[idx1],
                    intruder.gsnorth[idx2] - ownship.gsnorth[idx1]])

                past_cpa = np.dot(dist, vrel) > 0.0

                rpz = np.max(conf.rpz[[idx1, idx2]])
                hdist = np.linalg.norm(dist)
                hor_los = hdist < rpz

                is_bouncing = (
                    abs(anglediff(ownship.trk[idx1], intruder.trk[idx2])) < 30.0
                    and hdist < rpz * self.resofach)

            if idx2 >= 0 and (not past_cpa or hor_los or is_bouncing):
                changeactive[idx1] = True
            else:
                changeactive[idx1] = changeactive.get(idx1, False)
                delpairs.add(conflict)

        for idx, active in changeactive.items():
            self.active[idx] = active
            if not active:
                iwpid = bs.traf.ap.route[idx].findact(idx)
                if iwpid != -1:
                    bs.traf.ap.route[idx].direct(
                        idx, bs.traf.ap.route[idx].wpname[iwpid])

        self.resopairs -= delpairs

    # ------------------------------------------------------------------
    #  Stack commands (same names + _NAT suffix to avoid clashes)
    # ------------------------------------------------------------------
    @command(name='PRIORULES_NAT')
    def setprio(self, flag: bool = None, priocode=''):
        """Define priority rules for NAT conflict resolution."""
        if flag is None:
            if self.__class__ is ConflictResolutionNAT:
                return False, 'No NAT conflict resolution enabled.'
            return False, (f'Resolution algorithm {self.__class__.__name__} '
                           "hasn't implemented priority.")
        self.swprio = flag
        self.priocode = priocode
        return True

    @command(name='NORESO_NAT')
    def setnoreso(self, *idx: 'acid'):
        """Aircraft that nobody will avoid (NAT variant)."""
        if not idx:
            return True, ('NORESO_NAT [ACID, ...]\nCurrent list: '
                          + ', '.join(np.array(bs.traf.id)[self.noresoac]))
        idx = list(idx)
        self.noresoac[idx] = np.logical_not(self.noresoac[idx])
        return True

    @command(name='RESOOFF_NAT')
    def setresooff(self, *idx: 'acid'):
        """Aircraft that will not resolve (NAT variant)."""
        if not idx:
            return True, ('RESOOFF_NAT [ACID, ...]\nCurrent list: '
                          + ', '.join(np.array(bs.traf.id)[self.resooffac]))
        idx = list(idx)
        self.resooffac[idx] = np.logical_not(self.resooffac[idx])
        return True

    @command(name='RFACH_NAT', aliases=('RESOFACH_NAT',))
    def setresofach(self, factor: float = None):
        """Set horizontal resolution factor (NAT)."""
        if factor is None:
            return True, f'RFACH_NAT [FACTOR]\nCurrent: {self.resofach}'
        self.resofach = factor
        self.resorrelative = True
        return True, f'Horizontal resolution factor set to {self.resofach}'

    @command(name='RFACV_NAT', aliases=('RESOFACV_NAT',))
    def setresofacv(self, factor: float = None):
        """Set vertical resolution factor (NAT)."""
        if factor is None:
            return True, f'RFACV_NAT [FACTOR]\nCurrent: {self.resofacv}'
        self.resofacv = factor
        self.resodhrelative = True
        return True, f'Vertical resolution factor set to {self.resofacv}'

    @staticmethod
    @command(name='RESO_NAT')
    def setmethod(name: 'txt' = ''):
        """Select a NAT conflict resolution method."""
        methods = ConflictResolutionNAT.derived()
        names = ['OFF' if n == 'CONFLICTRESOLUTIONNAT' else n
                 for n in methods]
        if not name:
            curname = ('OFF'
                       if ConflictResolutionNAT.selected() is ConflictResolutionNAT
                       else ConflictResolutionNAT.selected().__name__)
            return True, (f'Current NAT CR method: {curname}\n'
                          f'Available: {", ".join(names)}')
        if name == 'OFF':
            ConflictResolutionNAT.select()
            return True, 'NAT Conflict Resolution turned off.'
        method = methods.get(name, None)
        if method is None:
            return False, (f'{name} not found.\n'
                           f'Available: {", ".join(names)}')
        method.select()
        return True, f'Selected {method.__name__} as NAT CR method.'
