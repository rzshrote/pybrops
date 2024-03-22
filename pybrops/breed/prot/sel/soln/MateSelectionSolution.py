"""
Module containing mate selection solution class interface.
"""

from abc import ABCMeta

import numpy

from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_ndim


class MateSelectionSolution(SelectionSolution,metaclass=ABCMeta):
    """
    Semi-abstract interface 
    """
    ########################## Special Object Methods ##########################
    # do not implement __init__()

    ############################ Object Properties #############################
    @property
    def decn_space_xmap(self) -> numpy.ndarray:
        """decn_space_xmap."""
        return self._decn_space_xmap
    @decn_space_xmap.setter
    def decn_space_xmap(self, value: numpy.ndarray) -> None:
        """Set decn_space_xmap."""
        check_is_ndarray(value, "decn_space_xmap")
        check_ndarray_ndim(value, "decn_space_xmap", 2)
        check_ndarray_dtype_is_integer(value, "decn_space_xmap")
        self._decn_space_xmap = value



################################## Utilities ###################################
def check_is_MateSelectionSolution(v: object, vname: str) -> None:
    """
    Check if object is of type MateSelectionSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, MateSelectionSolution):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,MateSelectionSolution.__name__,type(v).__name__))
