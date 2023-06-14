"""
Module defining interfaces for mate selection problems.
"""

from abc import ABCMeta

import numpy
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_ndim


class MateSelectionProblem(SelectionProblem,metaclass=ABCMeta):
    """
    docstring for MateSelectionProblem.
    """

    ########################## Special Object Methods ##########################
    # do not define __init__() since this is a semi-abstract/mixin class

    ############################ Object Properties #############################
    @property
    def decn_space_xmap(self) -> numpy.ndarray:
        """Decision space cross map."""
        return self._decn_space_xmap
    @decn_space_xmap.setter
    def decn_space_xmap(self, value: numpy.ndarray) -> None:
        """Set decision space cross map."""
        check_is_ndarray(value, "decn_space_xmap")
        check_ndarray_ndim(value, "decn_space_xmap", 2)
        check_ndarray_dtype_is_integer(value, "decn_space_xmap")
        self._decn_space_xmap = value
