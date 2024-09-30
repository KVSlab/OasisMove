"""Top-level package for OasisMove."""

from importlib.metadata import metadata

# Common
from .common import io, utilities

# Problems
from .problems import __init__ as problems_init
from .problems.NSfracStep import (
    MovingCommon,
    MovingCylinder,
    MovingTaylorGreen3D,
    MovingVortex,
    MovingWall,
)
from .run_oasis import main

# Solvers
from .solvers.NSfracStep import IPCS_ABCN_Move

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
    "io",
    "utilities",
]
