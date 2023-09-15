"""
Module providing matrix interfaces and basic matrix functionality.
"""

__all__ = [
    "util",
    "Matrix",
    "PrunableMatrix",
    "SquareMatrix",
    "MutableMatrix",
    "PhasedMatrix",
    "SortableMatrix",
    "TraitMatrix",
    "GroupableMatrix",
    "TaxaMatrix",
    "VariantMatrix",
    "SquareTaxaMatrix",
    "TaxaTraitMatrix",
    "TaxaVariantMatrix",
    "PhasedTaxaVariantMatrix",
    "DenseMatrix",
    "DenseSquareMatrix",
    "DenseMutableMatrix",
    "DensePhasedMatrix",
    "DenseTraitMatrix",
    "DenseTaxaMatrix",
    "DenseVariantMatrix",
    "DenseSquareTaxaMatrix",
    "DenseTaxaTraitMatrix",
    "DenseTaxaVariantMatrix",
    "DensePhasedTaxaVariantMatrix",
]

# order dependent!

# Utilities
from pybrops.core.mat import util

# level 0 interface
from pybrops.core.mat import Matrix           # order 0
from pybrops.core.mat import PrunableMatrix   # order 1
from pybrops.core.mat import SquareMatrix     # order 1

# level 1 interface
from pybrops.core.mat import MutableMatrix    # order 0
from pybrops.core.mat import PhasedMatrix     # order 1

# level 2 interface
from pybrops.core.mat import SortableMatrix   # order 0
from pybrops.core.mat import TraitMatrix      # order 1

# level 3 interface
from pybrops.core.mat import GroupableMatrix  # order 0
from pybrops.core.mat import TaxaMatrix       # order 1
from pybrops.core.mat import VariantMatrix    # order 1

# level 4 interface
from pybrops.core.mat import SquareTaxaMatrix
from pybrops.core.mat import TaxaTraitMatrix      # order 0
from pybrops.core.mat import TaxaVariantMatrix    # order 0

# level 5 interface
from pybrops.core.mat import PhasedTaxaVariantMatrix  # order 0

################################################################################
################################################################################
################################################################################

# level 0 implementation
from pybrops.core.mat import DenseMatrix          # order 0
from pybrops.core.mat import DenseSquareMatrix    # order 1

# level 1 implementation
from pybrops.core.mat import DenseMutableMatrix   # order 0
from pybrops.core.mat import DensePhasedMatrix    # order 1

# level 2 implementation
from pybrops.core.mat import DenseTraitMatrix

# level 3 implementation
from pybrops.core.mat import DenseTaxaMatrix      # order 0
from pybrops.core.mat import DenseVariantMatrix   # order 0

# level 4 implementation
from pybrops.core.mat import DenseSquareTaxaMatrix
from pybrops.core.mat import DenseTaxaTraitMatrix         # order 0
from pybrops.core.mat import DenseTaxaVariantMatrix       # order 0
from pybrops.core.mat import DensePhasedTaxaVariantMatrix # order 1
