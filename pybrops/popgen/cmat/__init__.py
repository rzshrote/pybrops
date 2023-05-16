"""
Module for calculating coancestry and kinship relationship matrices.
"""

__all__ = [
    "CoancestryMatrix",
    "DenseCoancestryMatrix",
    "DenseGeneralizedWeightedCoancestryMatrix",
    "DenseMolecularCoancestryMatrix",
    "DenseVanRadenCoancestryMatrix",
    "DenseYangCoancestryMatrix"
]

# abstract interfaces
from pybrops.popgen.cmat import CoancestryMatrix

# semi-abstract classes
from pybrops.popgen.cmat import DenseCoancestryMatrix

# concrete classes
from pybrops.popgen.cmat import DenseGeneralizedWeightedCoancestryMatrix
from pybrops.popgen.cmat import DenseMolecularCoancestryMatrix
from pybrops.popgen.cmat import DenseVanRadenCoancestryMatrix
from pybrops.popgen.cmat import DenseYangCoancestryMatrix
