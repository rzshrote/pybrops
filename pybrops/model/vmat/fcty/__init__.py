"""
Module containing Genetic and Genic Variance Matrix Factory interfaces and implementations
"""

__all__ = [
    "GeneticVarianceMatrixFactory",
    "AdditiveGeneticVarianceMatrixFactory",
    "GenicVarianceMatrixFactory",
    "AdditiveGenicVarianceMatrixFactory",
    "DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory",
    "DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory",
    "DenseFourWayDHAdditiveGeneticVarianceMatrixFactory",
    "DenseDihybridDHAdditiveGeneticVarianceMatrixFactory",
    "DenseTwoWayDHAdditiveGenicVarianceMatrixFactory"
]

##### abstract interfaces #####

# genetic variance matrix factories
from pybrops.model.vmat.fcty import GeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty import AdditiveGeneticVarianceMatrixFactory
# genic variance matrix factories
from pybrops.model.vmat.fcty import GenicVarianceMatrixFactory
from pybrops.model.vmat.fcty import AdditiveGenicVarianceMatrixFactory

##### concrete classes #####

# genetic variance matrix factories
from pybrops.model.vmat.fcty import DenseTwoWayDHAdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty import DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty import DenseFourWayDHAdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty import DenseDihybridDHAdditiveGeneticVarianceMatrixFactory

# genic variance matrix factories
from pybrops.model.vmat.fcty import DenseTwoWayDHAdditiveGenicVarianceMatrixFactory