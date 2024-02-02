"""
Experimental module containing variance matrix items.
"""

__all__ = [
    "util",
    "GeneticVarianceMatrix",
    "AdditiveGeneticVarianceMatrix",
    "DenseGeneticVarianceMatrix",
    "DenseAdditiveGeneticVarianceMatrix",
    "DenseTwoWayDHAdditiveGeneticVarianceMatrix",
    "DenseThreeWayDHAdditiveGeneticVarianceMatrix",
    "DenseFourWayDHAdditiveGeneticVarianceMatrix",
    "DenseDihybridDHAdditiveGeneticVarianceMatrix",
    "GenicVarianceMatrix",
    "AdditiveGenicVarianceMatrix",
    "DenseGenicVarianceMatrix",
    "DenseAdditiveGenicVarianceMatrix",
    "DenseTwoWayDHAdditiveGenicVarianceMatrix",
]

# order dependent imports

# level 0 utilities
from pybrops.model.vmat import util

# GeneticVarianceMatrix family
from pybrops.model.vmat import GeneticVarianceMatrix
from pybrops.model.vmat import AdditiveGeneticVarianceMatrix
from pybrops.model.vmat import DenseGeneticVarianceMatrix
from pybrops.model.vmat import DenseAdditiveGeneticVarianceMatrix
from pybrops.model.vmat import DenseTwoWayDHAdditiveGeneticVarianceMatrix
from pybrops.model.vmat import DenseThreeWayDHAdditiveGeneticVarianceMatrix
from pybrops.model.vmat import DenseFourWayDHAdditiveGeneticVarianceMatrix
from pybrops.model.vmat import DenseDihybridDHAdditiveGeneticVarianceMatrix

# GenicVarianceMatrix family
from pybrops.model.vmat import GenicVarianceMatrix
from pybrops.model.vmat import AdditiveGenicVarianceMatrix
from pybrops.model.vmat import DenseGenicVarianceMatrix
from pybrops.model.vmat import DenseAdditiveGenicVarianceMatrix
from pybrops.model.vmat import DenseTwoWayDHAdditiveGenicVarianceMatrix
