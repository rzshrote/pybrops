"""
Mixin class to provide functionality for selection configurations that require cross maps.
"""

from abc import ABCMeta, abstractmethod

import numpy
from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len, check_ndarray_ndim


class MateSelectionConfiguration(SelectionConfiguration,metaclass=ABCMeta):
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
        check_ndarray_dtype_is_integer(value, "xconfig_xmap")
        check_ndarray_ndim(value, "xconfig_xmap", 2)
        check_ndarray_axis_len(value, "xconfig_xmap", 1, self.nparent)
        self._xconfig_xmap = value

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_MateSelectionConfiguration(v: object, vname: str) -> None:
    """
    Check if object is of type MateSelectionConfiguration, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, MateSelectionConfiguration):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,MateSelectionConfiguration.__name__,type(v).__name__))
