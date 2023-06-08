"""
Module partially implementing the SelectionProblem class.
"""

# list of public objects in this module
__all__ = [
    "SelectionProblem",
    "check_is_SelectionProblem"
]

# imports
from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Iterable, Optional, Sequence, Tuple, Union
import numpy
from pybrops.breed.prot.sel.prob.trans import trans_empty, trans_identity
from pybrops.core.error.error_type_python import check_is_Callable, check_is_dict
from pybrops.opt.prob.Problem import Problem
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation

# inheritance order is important for method resolution order
class SelectionProblem(Problem,metaclass=ABCMeta):
    """
    Basal interface for all genotype selection problem specifications.

    All selection optimization problems have the form:

    .. math::

        \\min_{\\mathbf{x}} \\mathbf{w_{obj} \\odot T_{obj}(L(x))}

    Such that:

    .. math::

        \\mathbf{w_{ineqcv} \\odot T_{ineqcv}(L(x)) \\leq 0}

        \\mathbf{w_{eqcv} \\odot T_{eqcv}(L(x)) = 0}
    
    Where:

        - :math:`\\mathbf{x}` is a selection decision vector.
        - :math:`L(\\cdot)` is a latent vector encoding function transforming 
            the decision vector into a latent space. The ``latentfn`` function 
            defines this function in this interface.
        - :math:`w_{obj}` is an objective function weight vector.
        - :math:`T_{obj}(\\cdot)` is a function transforming a latent space 
            vector to an objective space vector.
        - :math:`w_{ineqcv}` is an inequality constraint violation function 
            weight vector.
        - :math:`T_{ineqcv}(\\cdot)` is a function transforming a latent space 
            vector to an inequality constraint violation vector.
        - :math:`w_{eqcv}` is an equality constraint violation function 
            weight vector.
        - :math:`T_{eqcv}(\\cdot)` is a function transforming a latent space 
            vector to an equality constraint violation vector.
    """

    ########################## Special Object Methods ##########################
    @abstractmethod
    def __init__(
            self,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
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
        Constructor for SelectionProblem.
        
        Parameters
        ----------
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
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
        replace_nan_values_by: Real, None
            Value for which to replace NaN values. See PyMOO documentation.
        exclude_from_serialization: Iterable, None
            Attributes which are excluded from being serialized. See PyMOO documentation.
        callback: Callable, None
            A callback function to be called after every evaluation. See PyMOO documentation.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance. See PyMOO documentation.
        """
        super(SelectionProblem, self).__init__(
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
    @property
    @abstractmethod
    def nlatent(self) -> Integral:
        """Size of the latent vector for the optimization problem."""
        raise NotImplementedError("property is abstract")

    @property
    def obj_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to objective function values."""
        return self._obj_trans
    @obj_trans.setter
    def obj_trans(self, value: Union[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to objective space transformation function."""
        if value is None:
            value = trans_identity
        check_is_Callable(value, "obj_trans")
        self._obj_trans = value
    
    @property
    def obj_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to objective space transformation function."""
        return self._obj_trans_kwargs
    @obj_trans_kwargs.setter
    def obj_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to objective space transformation function."""
        if value is None:
            value = {}
        check_is_dict(value, "obj_trans_kwargs")
        self._obj_trans_kwargs = value
    
    @property
    def ineqcv_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to inequality constraint violation values."""
        return self._ineqcv_trans
    @ineqcv_trans.setter
    def ineqcv_trans(self, value: Union[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to inequality constraint violation transformation function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "ineqcv_trans")
        self._ineqcv_trans = value
    
    @property
    def ineqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to inequality constraint violation transformation function."""
        return self._ineqcv_trans_kwargs
    @ineqcv_trans_kwargs.setter
    def ineqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to inequality constraint violation transformation function."""
        if value is None:
            value = {}
        check_is_dict(value, "ineqcv_trans_kwargs")
        self._ineqcv_trans_kwargs = value
    
    @property
    def eqcv_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to equality constraint violation values."""
        return self._eqcv_trans
    @eqcv_trans.setter
    def eqcv_trans(self, value: Union[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to equality constraint violation transformation function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "eqcv_trans")
        self._eqcv_trans = value 
    
    @property
    def eqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to equality constraint violation transformation function."""
        return self._eqcv_trans_kwargs
    @eqcv_trans_kwargs.setter
    def eqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to equality constraint violation transformation function."""
        if value is None:
            value = {}
        check_is_dict(value, "eqcv_trans_kwargs")
        self._eqcv_trans_kwargs = value

    ############################## Object Methods ##############################
    # leave latentfn abstract
    @abstractmethod
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Encode a candidate solution for the given Problem into an ``l`` 
        dimensional latent evaluation space.
        
        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(l,)`` containing latent evaluation values.
        """
        raise NotImplementedError("method is abstract")

    # provide generic implementation of this function; users may override if needed
    def evalfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> Tuple:
        """
        Evaluate a candidate solution for the given Problem.
        
        This calculates three vectors which are to be minimized:

        .. math::

            \\mathbf{v_{obj}} = \\mathbf{w_{obj} \\odot T_{obj}(L(x))} \\
            \\mathbf{v_{ineqcv}} = \\mathbf{w_{ineqcv} \\odot T_{ineqcv}(L(x))} \\
            \\mathbf{v_{eqcv}} = \\mathbf{w_{eqcv} \\odot T_{eqcv}(L(x))}

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : tuple
            A tuple ``(obj, ineqcv, eqcv)``.
            
            Where:
            
            - ``obj`` is a numpy.ndarray of shape ``(nobj,)`` that contains 
                objective function evaluations.
                This is equivalent to :math:`\\mathbf{v_{obj}}`
            - ``ineqcv`` is a numpy.ndarray of shape ``(nineqcv,)`` that contains 
                inequality constraint violation values.
                This is equivalent to :math:`\\mathbf{v_{ineqcv}}`
            - ``eqcv`` is a numpy.ndarray of shape ``(neqcv,)`` that contains 
                equality constraint violation values.
                This is equivalent to :math:`\\mathbf{v_{eqcv}}`
        """
        # get latent space evaluation
        latent = self.latentfn(x, *args, **kwargs)
        # transform latent values into evaluated values
        obj = self.obj_wt * self.obj_trans(x, latent, **self.obj_trans_kwargs)
        ineqcv = self.ineqcv_wt * self.ineqcv_trans(x, latent, **self.ineqcv_trans_kwargs)
        eqcv = self.eqcv_wt * self.eqcv_trans(x, latent, **self.eqcv_trans_kwargs)
        # return results
        return obj, ineqcv, eqcv

    # provide generic implementation of this function; users may override if needed
    def _evaluate(
            self, 
            x: numpy.ndarray, 
            out: dict, 
            *args: tuple, 
            **kwargs: dict
        ) -> None:
        """
        Evaluate a set of candidate solutions for the given Problem.
        This function is used to interface with PyMOO.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(nsoln,ndecn)``.
            Where ``nsoln`` is the number of candidates solutions and ``ndecn``
            is the number of decision variables.
        out : dict
            Dictionary to which to output function evaluations.
        args : tuple
            Additional arguments.
        kwargs : dict
            Additional keyword arguments.
        """
        if x.ndim == 1:
            # get evaluations
            vals = self.evalfn(x, *args, **kwargs)
            # create temporary dictionary
            tmp = {key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0}
            # update output dictionary
            out.update(tmp)
        else:
            # create lists for accumulating variables
            objs = []
            ineqcvs = []
            eqcvs = []
            # for each row in x
            for v in x:
                # get evaluations
                obj, ineqcv, eqcv = self.evalfn(v, *args, **kwargs)
                # append values to lists
                objs.append(obj)
                ineqcvs.append(ineqcv)
                eqcvs.append(eqcv)
            # stack outputs
            objs = numpy.stack(objs)
            ineqcvs = numpy.stack(ineqcvs)
            eqcvs = numpy.stack(eqcvs)
            # create temporary dictionary
            tmp = {key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0}
            # update output dictionary
            out.update(tmp)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type SelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelectionProblem):
        raise TypeError("'{0}' must be of type SelectionProblem.".format(vname))
