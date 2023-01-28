"""Top-level package for OasisMove."""
from importlib.metadata import metadata

from .run_oasis import main


meta = metadata("oasismove")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]
