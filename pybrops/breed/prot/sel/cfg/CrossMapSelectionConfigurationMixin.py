"""
Mixin class to provide functionality for selection configurations that require cross maps.
"""

from abc import ABCMeta, abstractmethod

import numpy

from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_ndim


class CrossMapSelectionConfigurationMixin(metaclass=ABCMeta):
    """
    A mixin class to provide functionality for selection configurations which 
    require random samples to be drawn.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################
    @property
    def xconfig_xmap(self) -> numpy.ndarray:
        """Cross map corresponding to the decision space."""
        return self._xconfig_xmap
    @xconfig_xmap.setter
    def xconfig_xmap(self, value: numpy.ndarray) -> None:
        """Set xconfig_xmap."""
        check_is_ndarray(value, "xconfig_xmap")
        check_ndarray_ndim(value, "xconfig_xmap", 2)
        self._xconfig_xmap = value

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
