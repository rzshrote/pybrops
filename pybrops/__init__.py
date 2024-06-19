"""
PyBrOpS: a software package for numerical optimization and simulation of
breeding programs.
"""

__all__ = [
    "core",
    "opt",
    "popgen",
    "model",
    "breed",
    "test",
]

# define version
__version__ = '1.0.4'

# order dependent import of submodules

# core before all other submodules
from pybrops import core

# algorithms
from pybrops import opt

# order dependent submodules
from pybrops import popgen
from pybrops import model
from pybrops import breed

# testing
from pybrops import test