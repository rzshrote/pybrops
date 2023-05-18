"""
Partial implementation of the SetSolution interface.
"""

# list of all public imports in the module
__all__ = [
    "BinarySolution",
    "check_is_BinarySolution"
]

# imports
from numbers import Integral, Real
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_bool_or_integer
from pybrops.core.error.error_value_numpy import check_ndarray_is_binary, check_ndarray_len_eq, check_ndarray_shape_eq
from pybrops.opt.soln.Solution import Solution

class BinarySolution(Solution):
    """
    Class for optimization problem solutions with binary decision variables.
    """

    ########################## Special Object Methods ##########################
    
    # implementation of abstract method
    def __init__(
            self,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Integral,bool,None],
            decn_space_upper: Union[numpy.ndarray,Integral,bool,None],
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
        Constructor for BinarySolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(BinarySolution, self).__init__(
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
            check_ndarray_dtype_is_bool_or_integer(value, "decn_space")
            check_ndarray_shape_eq(value, "decn_space", (2,self.ndecn))
            check_ndarray_is_binary(value, "decn_space")
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space' must be of type numpy.ndarray or None")
        self._decn_space = value

    # override decn_space setter properties
    @Solution.decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,Integral,bool,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_bool_or_integer(value, "decn_space_lower")
            check_ndarray_len_eq(value, "decn_space_lower", self.ndecn)
            check_ndarray_is_binary(value, "decn_space")
        elif isinstance(value, Integral):
            value = numpy.repeat(value, self.ndecn)
        elif isinstance(value, bool):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_lower' must be of type numpy.ndarray, an Integral type, bool, or None")
        self._decn_space_lower = value

    # override decn_space setter properties
    @Solution.decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,Integral,bool,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_bool_or_integer(value, "decn_space_upper")
            check_ndarray_len_eq(value, "decn_space_upper", self.ndecn)
            check_ndarray_is_binary(value, "decn_space")
        elif isinstance(value, Integral):
            value = numpy.repeat(value, self.ndecn)
        elif isinstance(value, bool):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_upper' must be of type numpy.ndarray, an Integral type, bool, or None")
        self._decn_space_upper = value





################################################################################
################################## Utilities ###################################
################################################################################
def check_is_BinarySolution(v: object, vname: str) -> None:
    """
    Check if object is of type BinarySolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BinarySolution):
        raise TypeError("'{0}' must be of type BinarySolution.".format(vname))
