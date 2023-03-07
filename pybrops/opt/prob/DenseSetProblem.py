"""
Partial implementation of the SetProblem interface.
"""

# list of public objects in this module
__all__ = [
    "DenseSetProblem",
    "check_is_DenseSetProblem"
]

# imports
from numbers import Integral
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_is_1d, check_ndarray_len_gteq
from pybrops.opt.prob.DenseProblem import DenseProblem
from pybrops.opt.prob.SetProblem import SetProblem

# inheritance ordering is important for method resolution order
class DenseSetProblem(DenseProblem,SetProblem):
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
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSetProblem.
        
        Parameters
        ----------
        ndecn : Integral
            Number of decision variables.
        decn_space : numpy.ndarray
            A 1d array containing the set of available elements
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseSetProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            nobj = nobj,
            obj_wt = obj_wt,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    # override decn_space setter properties
    @DenseProblem.decn_space.setter
    def decn_space(self, value: numpy.ndarray) -> None:
        """Set decision space boundaries."""
        check_is_ndarray(value, "decn_space")
        check_ndarray_is_1d(value, "decn_space")
        check_ndarray_len_gteq(value, "decn_space", self.ndecn)
        self._decn_space = value




################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseSetProblem(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSetProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSetProblem):
        raise TypeError("'{0}' must be of type DenseSetProblem.".format(vname))
