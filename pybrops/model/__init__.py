"""
Module containing breeding models. This is not restricted to simply genomic
prediction models, but can include other types of models like variance
estimations.
"""

__all__ = [
    "gmod",
    "vmat"
]

# imports are order dependent!!!!

# submodules
from pybrops.model import gmod
# from pybrops.model import emat
from pybrops.model import vmat
