# order dependent!

# Utilities
from .mat_util import *

# level 0 interface
from .Matrix import *           # order 0
from .PrunableMatrix import *   # order 1

# level 1 interface
from .MutableMatrix import *    # order 0
from .PhasedMatrix import *     # order 1

# level 2 interface
from .SortableMatrix import *   # order 0
from .TraitMatrix import *      # order 1

# level 3 interface
from .GroupableMatrix import *  # order 0
from .TaxaMatrix import *       # order 1
from .VariantMatrix import *    # order 1

# level 4 interface
from .TaxaTraitMatrix import *      # order 0
from .TaxaVariantMatrix import *    # order 0

# level 5 interface
from .PhasedTaxaVariantMatrix import *  # order 0

################################################################################
################################################################################
################################################################################

# level 0 implementation
from .DenseMatrix import *  # order 0

# level 1 implementation
from .DenseMutableMatrix import *   # order 0
from .DensePhasedMatrix import *    # order 1

# level 2 implementation
from .DenseTraitMatrix import *

# level 3 implementation
from .DenseTaxaMatrix import *      # order 0
from .DenseVariantMatrix import *   # order 0

# level 4 implementation
from .DenseTaxaTraitMatrix import *         # order 0
from .DenseTaxaVariantMatrix import *       # order 0
from .DensePhasedTaxaVariantMatrix import * # order 1
