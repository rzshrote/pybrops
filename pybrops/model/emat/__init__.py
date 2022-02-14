"""
Experimental module. Still under development.
"""
# imports are order dependent!!!

# abstract classes
# level 1 interface
from . import EffectMatrix

# level 2 interface
from . import GenotypeEffectMatrix
from . import HaplotypeEffectMatrix

# level 3 interface
from . import GenotypeEffectVariantMatrix
from . import HaplotypeEffectVariantMatrix

# concrete classes
# level 1 implementation
from . import PhasedGenotypeEffectMatrix
from . import PhasedHaplotypeEffectMatrix

# level 2 implementation
from . import PhasedGenotypeEffectVariantMatrix
from . import PhasedHaplotypeEffectVariantMatrix
