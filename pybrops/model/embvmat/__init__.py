"""
Module containing Expected Maximum Breeding Value (EMBV) matrix classes.

These are considered models because they include a classmethod which allows 
them to be constructed from Genomic Models.
"""

__all__ = [
    "ExpectedMaximumBreedingValueMatrix",
    "DenseExpectedMaximumBreedingValueMatrix",
]

# imports are order dependent!

# abstract interfaces
from pybrops.model.embvmat import ExpectedMaximumBreedingValueMatrix

# implementations
from pybrops.model.embvmat import DenseExpectedMaximumBreedingValueMatrix
