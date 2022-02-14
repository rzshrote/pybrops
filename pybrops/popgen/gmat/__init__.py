"""
Module containing genotype matrix containers
"""
# imports are order dependent!!!

# abstract classes
from . import GenotypeMatrix
from . import HaplotypeMatrix
from . import PhasedGenotypeMatrix
from . import PhasedHaplotypeMatrix

# implemented dense classes
from . import DenseGenotypeMatrix
from . import DensePhasedGenotypeMatrix
