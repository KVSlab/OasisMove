"""Top-level package for OasisMove."""
from importlib.metadata import metadata

from .run_oasis import main

# Solvers
from .solvers.NSfracStep import IPCS_ABCN_Move

# Common
from .common import io
from .common import utilities

# Problems
from .problems import __init__ as problems_init
from .problems.NSfracStep import MovingCommon
from .problems.NSfracStep import MovingCylinder
from .problems.NSfracStep import MovingTaylorGreen3D
from .problems.NSfracStep import MovingVortex
from .problems.NSfracStep import MovingWall
from .problems.NSfracStep import Probe


meta = metadata("oasismove")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]

__all__ = [
    "IPCS_ABCN_Move",
    "problems_init",
    "MovingWall",
    "MovingVortex",
    "MovingCommon",
    "MovingTaylorGreen3D",
    "MovingCylinder",
    "Probe",
    "io",
    "utilities"
]
