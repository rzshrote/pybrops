"""
Partial implementation of the IntegerProblem interface.
"""

# list of public objects in this module
__all__ = [
    "IntegerProblem",
    "check_is_IntegerProblem"
]

# imports
from numbers import Integral, Real
from typing import Callable, Iterable, Optional, Sequence, Union
import numpy
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_integer
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_shape_eq
from pybrops.opt.prob.Problem import Problem

# inheritance ordering is important for method resolution order
class IntegerProblem(Problem):
    """
    Partially implemented class for optimization problems with nominal decision 
    variables where the goal is to select an optimal subset.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            vtype: Optional[type] = None,
            vars: Optional[Sequence] = None,
            elementwise: bool = False,
            elementwise_func: type = ElementwiseEvaluationFunction,
            elementwise_runner: Callable = LoopedElementwiseEvaluation(),
            replace_nan_values_by: Optional[Real] = None,
            exclude_from_serialization: Optional[Iterable] = None,
            callback: Optional[Callable] = None,
            strict: bool = True,
            **kwargs: dict
        ) -> None:
        """
        Constructor for IntegerProblem.
        
        Parameters
        ----------
        ndecn : Integral
            Number of decision variables.
        decn_space : numpy.ndarray
            A 1d array containing the set of available elements
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(IntegerProblem, self).__init__(
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

    ############################ Object Properties #############################
    # override decn_space setter properties
    @Problem.decn_space.setter
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
    @Problem.decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_integer(value, "decn_space_lower")
            check_ndarray_len_eq(value, "decn_space_lower", self.ndecn)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_lower' must be of type numpy.ndarray, Real, or None")
        self._decn_space_lower = value
        self._xl = value

    # override decn_space setter properties
    @Problem.decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_integer(value, "decn_space_upper")
            check_ndarray_len_eq(value, "decn_space_upper", self.ndecn)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_upper' must be of type numpy.ndarray, Real, or None")
        self._decn_space_upper = value
        self._xu = value



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_IntegerProblem(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerProblem):
        raise TypeError("'{0}' must be of type IntegerProblem.".format(vname))
