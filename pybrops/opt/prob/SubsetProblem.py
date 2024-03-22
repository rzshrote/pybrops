"""
Partial implementation of the SubsetProblem interface.
"""

# list of public objects in this module
__all__ = [
    "SubsetProblem",
    "check_is_SubsetProblem",
]

# imports
from numbers import Integral
from numbers import Real
from typing import Callable
from typing import Iterable
from typing import Optional
from typing import Sequence
from typing import Union
import numpy
from pymoo.core.problem import ElementwiseEvaluationFunction
from pymoo.core.problem import LoopedElementwiseEvaluation
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_len_gteq
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.opt.prob.Problem import Problem

# inheritance ordering is important for method resolution order
class SubsetProblem(
        Problem,
    ):
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
            elementwise: bool = True,
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
        
        decn_space_lower : numpy.ndarray
            A 1d array representing the lower bound of the decision variables.
        
        decn_space_upper : numpy.ndarray
            A 1d array representing the upper bound of the decision variables.
        
        nobj : Integral
            Number of objectives for the problem.
        
        obj_wt : numpy.ndarray, Real, None
            Objective function weights. Weights from this vector are applied 
            to objective function values via the Hadamard product. If values 
            are ``1.0`` or ``-1.0``, this can be used to specify minimizing 
            and maximizing objectives, respectively.

            If ``obj_wt`` is ``numpy.ndarray``, then the array must be of shape 
            ``(nobj,)``.

            If ``obj_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``.

            If ``obj_wt`` is ``None``, then the value ``1.0`` is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``. This assumes that all 
            objectives are to be minimized.
        
        nineqcv : Integral, None
            Number of inequality constraint violation functions. This is 
            equivalent to the vector length returned by the ``ineqcv_trans`` 
            function. Must be ``Integral`` greater than or equal to zero.

            If ``nineqcv`` is ``None``, then set to zero.

        ineqcv_wt : numpy.ndarray, None
            Inequality constraint violation function weights. Weights from this 
            vector are applied to inequality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``ineqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(nineqcv,)``. This assumes that 
            all constraints are to be minimized.

        neqcv : Integral, None
            Number of equality constraint violations. This is equivalent to the 
            vector length returned by the ``eqcv_trans`` function. Must be 
            ``Integral`` greater than or equal to zero.
        
            If ``neqcv`` is ``None``, then set to zero.

        eqcv_wt : numpy.ndarray, None
            Equality constraint violation function weights. Weights from this 
            vector are applied to equality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``eqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(neqcv,)``. This assumes that 
            all constraints are to be minimized.

        vtype : type, None
            Used by PyMOO interface to construct a PyMOO Problem object.

        vars : Sequence, None
            Used by PyMOO interface to construct a PyMOO Problem object.

        elementwise : bool
            Used by PyMOO interface to construct a PyMOO Problem object.
        
        elementwise_func : type
            Used by PyMOO interface to construct a PyMOO Problem object.
        
        elementwise_runner : Callable
            Used by PyMOO interface to construct a PyMOO Problem object.
        
        replace_nan_values_by : Real, None
            Used by PyMOO interface to construct a PyMOO Problem object.
        
        exclude_from_serialization : Iterable, None
            Used by PyMOO interface to construct a PyMOO Problem object.
        
        callback : Callable, None
            Used by PyMOO interface to construct a PyMOO Problem object.
        
        strict : bool
            Used by PyMOO interface to construct a PyMOO Problem object.

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
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'.".format(vname,SubsetProblem.__name__,type(v).__name__))

def check_SubsetProblem_is_single_objective(v: SubsetProblem, vname: str) -> None:
    """
    Check if a SubsetProblem is single objective in nature, otherwise raise TypeError.

    Parameters
    ----------
    v : SubsetProblem
        A SubsetProblem for which to check the number of objectives.
    vname : str
        Name of variable to print in TypeError message.
    """
    if v.nobj != 1:
        raise TypeError("{0} '{1}' must be single objective in nature but received {1}.nobj == {2}".format(SubsetProblem.__name__,vname,v.nobj))

def check_SubsetProblem_is_multi_objective(v: SubsetProblem, vname: str) -> None:
    """
    Check if a SubsetProblem is multi objective in nature, otherwise raise TypeError.

    Parameters
    ----------
    v : SubsetProblem
        A SubsetProblem for which to check the number of objectives.
    vname : str
        Name of variable to print in TypeError message.
    """
    if v.nobj <= 1:
        raise TypeError("{0} '{1}' must be multi objective in nature but received {1}.nobj == {2}".format(SubsetProblem.__name__,vname,v.nobj))
