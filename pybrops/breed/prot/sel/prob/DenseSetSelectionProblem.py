"""
Module partially implementing the SetSelectionProblem interface.
"""

# list of public objects in this module
__all__ = [
    "DenseSetSelectionProblem",
    "check_is_DenseSetSelectionProblem"
]

# imports
from numbers import Integral
from typing import Callable, Tuple

import numpy
from pybrops.breed.prot.sel.prob.DenseSelectionProblem import DenseSelectionProblem
from pybrops.breed.prot.sel.prob.SetSelectionProblem import SetSelectionProblem
from pybrops.opt.prob.DenseSetProblem import DenseSetProblem

# inheritance ordering is important here to avoid circular dependency/method resolution issues
class DenseSetSelectionProblem(DenseSetProblem,DenseSelectionProblem,SetSelectionProblem):
    """
    docstring for DenseSetSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            ndecn: Integral,
            decn_space: numpy.ndarray,
            encode_trans: Callable[[numpy.ndarray,dict],Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]],
            encode_trans_kwargs: dict,
            nobj: Integral,
            obj_wt: numpy.ndarray,
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSetSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # do not use super constructor since there are two dense constructors
        # order dependent
        self.ndecn = ndecn
        self.decn_space = decn_space
        self.encode_trans = encode_trans
        self.encode_trans_kwargs = encode_trans_kwargs
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    # leave encodefn abstract
    # leave evalfn abstract



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseSetSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSetSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSetSelectionProblem):
        raise TypeError("'{0}' must be of type DenseSetSelectionProblem.".format(vname))
