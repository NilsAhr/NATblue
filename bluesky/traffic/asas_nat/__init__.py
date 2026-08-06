"""
asas_nat — ASAS (CD&R) package optimised for North Atlantic en-route traffic.
=============================================================================

This package extends the standard ``bluesky.traffic.asas`` classes with
NAT-specific improvements.  The classes inherit from the upstream bases
so they register in the same ``replaceable`` Entity lineage and become
available via the standard ``RESO`` / ``CDMETHOD`` stack commands:

    RESO MVPNAT          — NAT-optimised Modified Voltage Potential
    CDMETHOD STATEBASEDNAT — NAT state-based detection (mirror, future use)

Modules
-------
mvp_nat
    ``MVPNAT`` — vertical-first MVP with ceiling guard, MMO guard,
    symmetry breaker, VS-stop, and anti-bouncing logic.
statebased_nat
    ``StateBasedNAT`` — mirror of upstream StateBased (for future
    NAT-specific detection tweaks).

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-04
"""
from .mvp_nat import MVPNAT
from .mvp2nat import MVP2NAT
from .statebased_nat import StateBasedNAT, StateBasedRealVS
from .intentbased import IntentBasedNAT
