"""
Module containing Weighted Genomic Estimated Breeding Value (wGEBV) matrix 
classes.

These are considered models because they include a classmethod which allows 
them to be constructed from Genomic Models.
"""

__all__ = [
    "WeightedGenomicEstimatedBreedingValueMatrix",
    "DenseWeightedGenomicEstimatedBreedingValueMatrix",
]

# imports are order dependent!

# abstract interfaces
from pybrops.model.wgebvmat import WeightedGenomicEstimatedBreedingValueMatrix

# implementations
from pybrops.model.wgebvmat import DenseWeightedGenomicEstimatedBreedingValueMatrix
