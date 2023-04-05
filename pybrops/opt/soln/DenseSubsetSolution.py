"""
Partial implementation of the SetSolution interface.
"""

# list of all public imports in the module
__all__ = [
    "DenseSubsetSolution",
    "check_is_DenseSubsetSolution"
]

# imports
from numbers import Integral
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_is_1d, check_ndarray_len_gteq
from pybrops.opt.soln.DenseSolution import DenseSolution
from pybrops.opt.soln.SubsetSolution import SubsetSolution

class DenseSubsetSolution(DenseSolution,SubsetSolution):
    """
    Partially implemented class for optimization problems with nominal decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            ndecn: Integral,
            decn_space: numpy.ndarray,
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
        Constructor for DenseSubsetSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseSubsetSolution, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
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

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    # override decn_space setter properties
    @DenseSolution.decn_space.setter
    def decn_space(self, value: numpy.ndarray) -> None:
        """Set decision space boundaries."""
        check_is_ndarray(value, "decn_space")
        check_ndarray_is_1d(value, "decn_space")
        check_ndarray_len_gteq(value, "decn_space", self.ndecn)
        self._decn_space = value



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseSubsetSolution(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSubsetSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSubsetSolution):
        raise TypeError("'{0}' must be of type DenseSubsetSolution.".format(vname))
