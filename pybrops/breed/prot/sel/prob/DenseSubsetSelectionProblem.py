"""
Module partially implementing the SetSelectionProblem interface.
"""

# list of public objects in this module
__all__ = [
    "DenseSubsetSelectionProblem",
    "check_is_DenseSubsetSelectionProblem"
]

# imports
import numpy
from numbers import Integral, Number
from typing import Callable, Iterable, Optional, Sequence, Tuple, Union
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation

from pybrops.breed.prot.sel.prob.DenseSelectionProblem import DenseSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.opt.prob.DenseSubsetProblem import DenseSubsetProblem

# inheritance ordering is important here to avoid circular dependency/method resolution issues
class DenseSubsetSelectionProblem(DenseSubsetProblem,DenseSelectionProblem,SubsetSelectionProblem):
    """
    docstring for DenseSubsetSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            encode_trans: Callable[[numpy.ndarray,dict],Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]],
            encode_trans_kwargs: dict,
            nobj: Integral,
            obj_wt: numpy.ndarray,
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            vtype: Optional[type] = None,
            vars: Optional[Sequence] = None,
            elementwise: bool = True,
            elementwise_func: type = ElementwiseEvaluationFunction,
            elementwise_runner: Callable = LoopedElementwiseEvaluation(),
            replace_nan_values_by: Optional[Number] = None,
            exclude_from_serialization: Optional[Iterable] = None,
            callback: Optional[Callable] = None,
            strict: bool = True,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSubsetSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseSubsetSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            encode_trans = encode_trans,
            encode_trans_kwargs = encode_trans_kwargs,
            nobj = nobj,
            obj_wt = obj_wt,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            vtype = vtype,
            vars = vars,
            elementwise = elementwise,
            elementwise_func = elementwise_func,
            elementwise_runner = elementwise_runner,
            replace_nan_values_by = replace_nan_values_by,
            exclude_from_serialization = exclude_from_serialization,
            callback = callback,
            strict = strict,
            **kwargs
        )

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    # leave encodefn abstract
    # evalfn defined by DenseSelectionProblem
    # _evaluate defined by DenseSelectionProblem



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseSubsetSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSubsetSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSubsetSelectionProblem):
        raise TypeError("'{0}' must be of type DenseSubsetSelectionProblem.".format(vname))
