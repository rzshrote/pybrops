"""
Module containing items related to population genetics.
"""

__all__ = [
    "gmap",
    "gmat",
    "bvmat",
    "cmat",
]

# order dependent

# genetic map components
from pybrops.popgen import gmap

# matrix components
from pybrops.popgen import gmat
from pybrops.popgen import bvmat

# other matrix components
# from pybrops.popgen import ldmat
from pybrops.popgen import cmat
