"""
Module partially implementing the SelectionProblem class
"""

# list of public objects in this module
__all__ = [
    "DenseSelectionProblem",
    "check_is_DenseSelectionProblem"
]

# imports
from numbers import Integral, Number
from typing import Callable, Iterable, Optional, Sequence, Tuple, Union
import numpy
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.core.error.error_type_python import check_is_Callable, check_is_dict
from pybrops.opt.prob.DenseProblem import DenseProblem
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation

# inheritance order is important for method resolution order
class DenseSelectionProblem(DenseProblem,SelectionProblem):
    """
    docstring for DenseSelectionProblem.
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
            elementwise: bool = False,
            elementwise_func: type = ElementwiseEvaluationFunction,
            elementwise_runner: Callable = LoopedElementwiseEvaluation(),
            replace_nan_values_by: Optional[Number] = None,
            exclude_from_serialization: Optional[Iterable] = None,
            callback: Optional[Callable] = None,
            strict: bool = True,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseSelectionProblem, self).__init__(
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
        # make assignments
        self.encode_trans = encode_trans
        self.encode_trans_kwargs = encode_trans_kwargs

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def encode_trans(self) -> Callable[[numpy.ndarray,dict],Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]]:
        """Function which transforms outputs from ``encodefn`` to a tuple ``(obj,ineqcv,eqcv)``."""
        return self._encd_trans
    @encode_trans.setter
    def encode_trans(self, value: Callable[[numpy.ndarray,dict],Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]]) -> None:
        """Set ``encodefn`` output transformation function."""
        check_is_Callable(value, "encd_trans")
        # TODO: check signature
        self._encd_trans = value
    @encode_trans.deleter
    def encode_trans(self) -> None:
        """Delete ``encodefn`` output transformation function."""
        del self._encd_trans
    
    @property
    def encode_trans_kwargs(self) -> dict:
        """``encodefn`` output transformation function keyword arguments."""
        return self._encd_trans_kwargs
    @encode_trans_kwargs.setter
    def encode_trans_kwargs(self, value: dict) -> None:
        """Set ``encodefn`` output transformation function keyword arguments."""
        check_is_dict(value, "encd_trans_kwargs")
        self._encd_trans_kwargs = value
    @encode_trans_kwargs.deleter
    def encode_trans_kwargs(self) -> None:
        """Delete ``encodefn`` output transformation function keyword arguments."""
        del self._encd_trans_kwargs

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    # leave encodefn abstract
    # leave evalfn abstract
    # leave _evaluate abstract



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSelectionProblem):
        raise TypeError("'{0}' must be of type DenseSelectionProblem.".format(vname))
