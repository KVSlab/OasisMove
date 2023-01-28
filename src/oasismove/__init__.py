"""Top-level package for OasisMove."""
from importlib.metadata import metadata

from .run_oasis import main
from .solvers.NSfracStep import IPCS_ABCN_Move

meta = metadata("oasismove")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]

__all__ = [
    "IPCS_ABCN_Move",
]
