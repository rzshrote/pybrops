"""
Module containing selection configuration definitions.
"""

__all__ = [
    "SampledSelectionConfigurationMixin",
    "CrossMapSelectionConfigurationMixin",
    "SelectionConfiguration",
    "SimpleSelectionConfiguration",
    "BinarySelectionConfiguration",
    "IntegerSelectionConfiguration",
    "RealSelectionConfiguration",
    "SubsetSelectionConfiguration"
]

# order dependent imports
# mixin semi-abstract classes
from pybrops.breed.prot.sel.cfg import SampledSelectionConfigurationMixin
from pybrops.breed.prot.sel.cfg import CrossMapSelectionConfigurationMixin


# basal semi-abstract class
from pybrops.breed.prot.sel.cfg import SelectionConfiguration
from pybrops.breed.prot.sel.cfg import SimpleSelectionConfiguration

# Sampled selection configurations
from pybrops.breed.prot.sel.cfg import BinarySelectionConfiguration
from pybrops.breed.prot.sel.cfg import IntegerSelectionConfiguration
from pybrops.breed.prot.sel.cfg import RealSelectionConfiguration
from pybrops.breed.prot.sel.cfg import SubsetSelectionConfiguration
