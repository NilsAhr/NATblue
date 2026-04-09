"""
detection_nat.py — NAT-aware Conflict Detection base class.
============================================================

Currently a thin wrapper around the standard ``ConflictDetection``.
This exists so that NAT-specific detection tweaks (e.g. extended
lookahead, asymmetric PZ, RNP-based PZ scaling) can be added here
without touching the upstream class.

@author : Nils Ahrenhold (ahre_ni)
@date   : 2026-04
"""
import numpy as np

import bluesky as bs
from bluesky.tools.aero import ft, nm
from bluesky.core import Entity
from bluesky.stack import command

# Inherit from the *original* base so BlueSky's Entity replacement
# mechanism treats this as a separate lineage.
from bluesky.traffic.asas.detection import ConflictDetection


class ConflictDetectionNAT(ConflictDetection):
    """NAT-specific conflict detection base.

    Inherits all behaviour from ``ConflictDetection``.
    Override methods here for NAT-specific tweaks.
    """

    def __init__(self):
        super().__init__()
        # --- NAT-specific additions go below ---
