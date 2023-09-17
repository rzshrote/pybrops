"""
Module containing genotyping protocols.
"""

__all__ = [
    "GenotypingProtocol",
    "DenseUnphasedGenotyping",
]

# order dependent imports

# abstract classes
from pybrops.breed.prot.gt import GenotypingProtocol

# concrete classes
from pybrops.breed.prot.gt import DenseUnphasedGenotyping
