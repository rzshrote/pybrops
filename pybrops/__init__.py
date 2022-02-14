"""
PyBrOpS: a software package for numerical optimization and simulation of
breeding programs.
"""

# order dependent import of submodules

# core before all other submodules
from . import core

# algorithms
from . import algo

# order dependent submodules
from . import popgen
from . import model
from . import breed
