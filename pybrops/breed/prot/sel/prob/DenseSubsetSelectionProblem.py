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
from typing import Callable, Iterable, Optional, Sequence, Union
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
            nobj: Integral,
            obj_wt: numpy.ndarray,
            obj_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            eqcv_trans_kwargs: Optional[dict],
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
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Number, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Number is provided, construct an array of shape ``(ndecn,)`` containing the Number.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Number, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Number is provided, construct an array of shape ``(ndecn,)`` containing the Number.
            If None, do not set a upper limit for the decision variables.
        nobj: Integral
            Number of objectives.
        obj_wt: numpy.ndarray
            Objective function weights.
        obj_trans: Callable, None
            A transformation function transforming a latent space vector to an objective space vector.
            The transformation function must be of the form: ``obj_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the identity transformation function: copy the latent space vector to the objective space vector.
        obj_trans_kwargs: dict, None
            Keyword arguments for the latent space to objective space transformation function.
            If None, an empty dictionary is used.
        nineqcv: Integral,
            Number of inequality constraints.
        ineqcv_wt: numpy.ndarray,
            Inequality constraint violation weights.
        ineqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an inequality constraint violation vector.
            The transformation function must be of the form: ``ineqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        ineqcv_trans_kwargs: Optional[dict],
            Keyword arguments for the latent space to inequality constraint violation space transformation function.
            If None, an empty dictionary is used.
        neqcv: Integral
            Number of equality constraints.
        eqcv_wt: numpy.ndarray
            Equality constraint violation weights.
        eqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an equality constraint violation vector.
            The transformation function must be of the form: ``eqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        eqcv_trans_kwargs: dict, None
            Keyword arguments for the latent space to equality constraint violation space transformation function.
            If None, an empty dictionary is used.
        vtype: type, None
            The variable type. So far, just used as a type hint. See PyMOO documentation.
        vars: Sequence, None
            Variables provided in their explicit form. See PyMOO documentation.
        elementwise: bool
            Whether the evaluation function should be run elementwise. See PyMOO documentation.
        elementwise_func: type
            A class that creates the function that evaluates a single individual. See PyMOO documentation.
        elementwise_runner: Callable
            A function that runs the function that evaluates a single individual. See PyMOO documentation.
        replace_nan_values_by: Number, None
            Value for which to replace NaN values. See PyMOO documentation.
        exclude_from_serialization: Iterable, None
            Attributes which are excluded from being serialized. See PyMOO documentation.
        callback: Callable, None
            A callback function to be called after every evaluation. See PyMOO documentation.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance. See PyMOO documentation.
        """
        # call the DenseSubsetProblem constructor
        super(DenseSubsetSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans, # somehow required??? TODO: understand why
            obj_trans_kwargs = obj_trans_kwargs, # somehow required??? TODO: understand why
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans, # somehow required??? TODO: understand why
            ineqcv_trans_kwargs = ineqcv_trans_kwargs, # somehow required??? TODO: understand why
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans, # somehow required??? TODO: understand why
            eqcv_trans_kwargs = eqcv_trans_kwargs, # somehow required??? TODO: understand why
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
        self.obj_trans = obj_trans
        self.obj_trans_kwargs = obj_trans_kwargs
        self.ineqcv_trans = ineqcv_trans
        self.ineqcv_trans_kwargs = ineqcv_trans_kwargs
        self.eqcv_trans = eqcv_trans
        self.eqcv_trans_kwargs = eqcv_trans_kwargs

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    # leave nlatent property abstract

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    # leave latentfn abstract
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
