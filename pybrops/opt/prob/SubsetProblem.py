"""
Partial implementation of the SubsetProblem interface.
"""

# list of public objects in this module
__all__ = [
    "SubsetProblem",
    "check_is_SubsetProblem"
]

# imports
from numbers import Integral, Real
from typing import Callable, Iterable, Optional, Sequence, Union
import numpy
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_len_gteq, check_ndarray_ndim
from pybrops.opt.prob.Problem import Problem

# inheritance ordering is important for method resolution order
class SubsetProblem(Problem):
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
        Constructor for SubsetProblem.
        
        Parameters
        ----------
        ndecn : Integral
            Number of decision variables.
        decn_space : numpy.ndarray
            A 1d array containing the set of available elements
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # order dependent assignments (for PyBrOpS interface)
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

        # call PyMOO constructor to set things its way (for PyMOO interface)
        super(Problem, self).__init__(
            n_var = ndecn,
            n_obj = nobj,
            n_ieq_constr = nineqcv,
            n_eq_constr = neqcv,
            xl = decn_space_lower,
            xu = decn_space_upper,
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
    def decn_space(self, value: numpy.ndarray) -> None:
        """Set decision space boundaries."""
        check_is_ndarray(value, "decn_space")
        check_ndarray_ndim(value, "decn_space", 1)
        check_ndarray_len_gteq(value, "decn_space", self.ndecn)
        self._decn_space = value



################################## Utilities ###################################
def check_is_SubsetProblem(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetProblem):
        raise TypeError("'{0}' must be of type SubsetProblem.".format(vname))
