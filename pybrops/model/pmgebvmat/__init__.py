"""
Module containing Progeny Mean Genomic Estimated Breeding Value (pmGEBV) matrix classes.
"""

__all__ = [
    "ParentalMeanGenomicEstimatedBreedingValueMatrix",
    "DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix",
]

# order dependent imports!

# abstract interfaces
from pybrops.model.pmgebvmat import ParentalMeanGenomicEstimatedBreedingValueMatrix

# implementations
from pybrops.model.pmgebvmat import DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix
