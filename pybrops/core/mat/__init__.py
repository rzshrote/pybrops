"""
Module providing base matrix functionality.
"""

# order dependent!

# Utilities
from . import util

# level 0 interface
from . import Matrix           # order 0
from . import PrunableMatrix   # order 1
from . import SquareMatrix     # order 1

# level 1 interface
from . import MutableMatrix    # order 0
from . import PhasedMatrix     # order 1

# level 2 interface
from . import SortableMatrix   # order 0
from . import TraitMatrix      # order 1

# level 3 interface
from . import GroupableMatrix  # order 0
from . import TaxaMatrix       # order 1
from . import VariantMatrix    # order 1

# level 4 interface
from . import SquareTaxaMatrix
from . import TaxaTraitMatrix      # order 0
from . import TaxaVariantMatrix    # order 0

# level 5 interface
from . import PhasedTaxaVariantMatrix  # order 0

################################################################################
################################################################################
################################################################################

# level 0 implementation
from . import DenseMatrix          # order 0
from . import DenseSquareMatrix    # order 1

# level 1 implementation
from . import DenseMutableMatrix   # order 0
from . import DensePhasedMatrix    # order 1

# level 2 implementation
from . import DenseTraitMatrix

# level 3 implementation
from . import DenseTaxaMatrix      # order 0
from . import DenseVariantMatrix   # order 0

# level 4 implementation
from . import DenseSquareTaxaMatrix
from . import DenseTaxaTraitMatrix         # order 0
from . import DenseTaxaVariantMatrix       # order 0
from . import DensePhasedTaxaVariantMatrix # order 1
