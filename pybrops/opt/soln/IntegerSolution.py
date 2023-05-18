"""
Partial implementation of the SetSolution interface.
"""

# list of all public imports in the module
__all__ = [
    "IntegerSolution",
    "check_is_IntegerSolution"
]

# imports
from numbers import Integral, Real
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_shape_eq
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
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Integral,None],
            decn_space_upper: Union[numpy.ndarray,Integral,None],
            nobj: Integral,
            obj_wt: Union[numpy.ndarray,Real,None],
            nineqcv: Union[Integral,None],
            ineqcv_wt: Union[numpy.ndarray,Real,None],
            neqcv: Union[Integral,None],
            eqcv_wt: Union[numpy.ndarray,Real,None],
            nsoln: Integral,
            soln_decn: numpy.ndarray,
            soln_obj: numpy.ndarray,
            soln_ineqcv: Union[numpy.ndarray,None],
            soln_eqcv: Union[numpy.ndarray,None],
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
            check_ndarray_dtype_is_integer(value, "decn_space")
            check_ndarray_shape_eq(value, "decn_space", (2,self.ndecn))
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space' must be of type numpy.ndarray or None")
        self._decn_space = value

    # override decn_space setter properties
    @Solution.decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,Integral,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_integer(value, "decn_space_lower")
            check_ndarray_len_eq(value, "decn_space_lower", self.ndecn)
        elif isinstance(value, Integral):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_lower' must be of type numpy.ndarray, an Integral type, or None")
        self._decn_space_lower = value

    # override decn_space setter properties
    @Solution.decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,Integral,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_integer(value, "decn_space_upper")
            check_ndarray_len_eq(value, "decn_space_upper", self.ndecn)
        elif isinstance(value, Integral):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_upper' must be of type numpy.ndarray, an Integral type, or None")
        self._decn_space_upper = value




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
