"""
Module containing Progeny Mean Estimated Breeding Value (pmEBV) matrix classes.
"""

__all__ = [
    "ProgenyMeanEstimatedBreedingValueMatrix",
    "DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix",
]

# imports are order dependent!

# abstract interfaces
from pybrops.model.pmebvmat import ProgenyMeanEstimatedBreedingValueMatrix

# implmentations
from pybrops.model.pmebvmat import DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix
