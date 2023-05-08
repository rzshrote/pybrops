"""
Partial implementation of the SetSolution interface.
"""

# list of all public imports in the module
__all__ = [
    "IntegerSolution",
    "check_is_IntegerSolution"
]

# imports
from numbers import Integral
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_shape_eq
from pybrops.opt.soln.Solution import Solution

class IntegerSolution(Solution):
    """
    Class for optimization problem solutions with integer decision variables.
    """

    ########################## Special Object Methods ##########################
    
    # implementation of abstract method
    def __init__(
            self,
            ndecn: Integral,
            decn_space: numpy.ndarray,
            decn_space_lower: numpy.ndarray,
            decn_space_upper: numpy.ndarray,
            nobj: Integral,
            obj_wt: numpy.ndarray,
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            nsoln: Integral,
            soln_decn: numpy.ndarray,
            soln_obj: numpy.ndarray,
            soln_ineqcv: numpy.ndarray,
            soln_eqcv: numpy.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for IntegerSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(IntegerSolution, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            nsoln = nsoln,
            soln_decn = soln_decn,
            soln_obj = soln_obj,
            soln_ineqcv = soln_ineqcv,
            soln_eqcv = soln_eqcv,
            **kwargs
        )

    ############################ Object Properties #############################
    
    # override decn_space setter properties
    @Solution.decn_space.setter
    def decn_space(self, value: Union[numpy.ndarray,None]) -> None:
        """Set decision space boundaries."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_shape_eq(value, "decn_space", (2,self.ndecn))
            check_ndarray_dtype_is_integer(value, "decn_space")
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space' must be of type numpy.ndarray or None")
        self._decn_space = value




################################################################################
################################## Utilities ###################################
################################################################################
def check_is_IntegerSolution(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerSolution):
        raise TypeError("'{0}' must be of type IntegerSolution.".format(vname))
