"""
Partial implementation of the SetSolution interface.
"""

# list of all public imports in the module
__all__ = [
    "SubsetSolution",
    "check_is_SubsetSolution"
]

# imports
from numbers import Integral, Real
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_len_gteq, check_ndarray_ndim
from pybrops.opt.soln.Solution import Solution

class SubsetSolution(Solution):
    """
    Class for optimization problem solutions with nominal decision variables.
    """

    ########################## Special Object Methods ##########################
    
    # implementation of abstract method
    def __init__(
            self,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,object,None],
            decn_space_upper: Union[numpy.ndarray,object,None],
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
        Constructor for SubsetSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # order dependent assignments
        self.ndecn = ndecn
        self.decn_space = decn_space
        self.decn_space_lower = decn_space_lower
        self.decn_space_upper = decn_space_upper
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt
        self.nsoln = nsoln
        self.soln_decn = soln_decn
        self.soln_obj = soln_obj
        self.soln_ineqcv = soln_ineqcv
        self.soln_eqcv = soln_eqcv

    ############################ Object Properties #############################
    
    # override decn_space setter properties
    @Solution.decn_space.setter
    def decn_space(self, value: Union[numpy.ndarray,None]) -> None:
        """Set decision space boundaries."""
        check_is_ndarray(value, "decn_space")
        check_ndarray_ndim(value, "decn_space", 1)
        check_ndarray_len_gteq(value, "decn_space", self.ndecn)
        self._decn_space = value

    # override decn_space setter properties
    @Solution.decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,object,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "decn_space_lower", self.ndecn)
        elif isinstance(value, object):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_lower' must be of type numpy.ndarray, object, or None")
        self._decn_space_lower = value

    # override decn_space setter properties
    @Solution.decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,object,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "decn_space_upper", self.ndecn)
        elif isinstance(value, object):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_upper' must be of type numpy.ndarray, object, or None")
        self._decn_space_upper = value



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SubsetSolution(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetSolution):
        raise TypeError("'{0}' must be of type SubsetSolution.".format(vname))
