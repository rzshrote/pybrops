"""
Module containing selection configuration definitions.
"""

__all__ = [
    "SampledSelectionConfigurationMixin",
    "MateSelectionConfiguration",
    "SelectionConfiguration",
    "SimpleSelectionConfiguration",
    "BinarySelectionConfiguration",
    "IntegerSelectionConfiguration",
    "RealSelectionConfiguration",
    "SubsetSelectionConfiguration",
    "SimpleMateSelectionConfiguration",
    "BinaryMateSelectionConfiguration",
    "IntegerMateSelectionConfiguration",
    "RealMateSelectionConfiguration",
    "SubsetMateSelectionConfiguration"
]

# order dependent imports
# mixin semi-abstract classes
from pybrops.breed.prot.sel.cfg import SampledSelectionConfigurationMixin

# basal semi-abstract class
from pybrops.breed.prot.sel.cfg import SelectionConfiguration
from pybrops.breed.prot.sel.cfg import MateSelectionConfiguration

# concrete classes
from pybrops.breed.prot.sel.cfg import SimpleSelectionConfiguration
from pybrops.breed.prot.sel.cfg import SimpleMateSelectionConfiguration

# Sampled selection configurations
from pybrops.breed.prot.sel.cfg import BinarySelectionConfiguration
from pybrops.breed.prot.sel.cfg import IntegerSelectionConfiguration
from pybrops.breed.prot.sel.cfg import RealSelectionConfiguration
from pybrops.breed.prot.sel.cfg import SubsetSelectionConfiguration

from pybrops.breed.prot.sel.cfg import BinaryMateSelectionConfiguration
from pybrops.breed.prot.sel.cfg import IntegerMateSelectionConfiguration
from pybrops.breed.prot.sel.cfg import RealMateSelectionConfiguration
from pybrops.breed.prot.sel.cfg import SubsetMateSelectionConfiguration
