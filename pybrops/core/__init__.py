"""
Module containing core, generalized utilities for ``PyBrOpS``.

These utilities are:

- ``error`` - Provides error checking routines.
- ``io`` - Provides input/output interfaces.
- ``mat`` - Provides general matrix interfaces and implementations.
- ``random`` - Provides random number generation.
- ``util`` - Provides general, miscellaneous utilities.
"""

__all__ = [
    "error",
    "util",
    "random",
    "io",
    "mat",
]

# order dependent import

# all error functions must go first!
from pybrops.core import error

# utility functions
from pybrops.core import util

# random number generator interface
from pybrops.core import random

# input/output interfaces: must be imported before ``mat``
from pybrops.core import io

# base matrix interfaces
from pybrops.core import mat
