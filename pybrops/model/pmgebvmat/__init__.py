"""
Module containing Progeny Mean Genomic Estimated Breeding Value (pmGEBV) matrix classes.
"""

__all__ = [
    "ProgenyMeanGenomicEstimatedBreedingValueMatrix",
    "DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix",
]

# order dependent imports!

# abstract interfaces
from pybrops.model.pmgebvmat import ProgenyMeanGenomicEstimatedBreedingValueMatrix

# implementations
from pybrops.model.pmgebvmat import DenseTwoWayProgenyMeanGenomicEstimatedBreedingValueMatrix
