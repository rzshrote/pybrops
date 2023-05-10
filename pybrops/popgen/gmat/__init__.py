"""
Module containing genotype matrix containers
"""

__all__ = [
    "GenotypeMatrix",
    "HaplotypeMatrix",
    "PhasedGenotypeMatrix",
    "PhasedHaplotypeMatrix",
    "DenseGenotypeMatrix",
    "DensePhasedGenotypeMatrix"
]

# imports are order dependent!!!

# abstract classes
from pybrops.popgen.gmat import GenotypeMatrix
from pybrops.popgen.gmat import HaplotypeMatrix
from pybrops.popgen.gmat import PhasedGenotypeMatrix
from pybrops.popgen.gmat import PhasedHaplotypeMatrix

# implemented dense classes
from pybrops.popgen.gmat import DenseGenotypeMatrix
from pybrops.popgen.gmat import DensePhasedGenotypeMatrix
