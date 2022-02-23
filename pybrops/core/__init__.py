"""
Module containing core, generalized utilities for ``PyBrOpS``.

These utilities are:

- ``df`` - Provides general dataframe support.
- ``error`` - Provides error checking routines.
- ``io`` - Provides input/output interfaces.
- ``mat`` - Provides general matrix interfaces and implementations.
- ``random`` - Provides random number generation.
- ``util`` - Provides general, miscellaneous utilities.
"""

# order dependent import

# all error functions must go first!
from . import error

# utility functions
from . import util

# random number generator interface
from . import random

# input/output interfaces
from . import io

# base matrix interfaces
from . import mat
