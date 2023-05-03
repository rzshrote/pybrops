"""
Module defining optimization problems.
"""

# list of public objects in this module
__all__ = [
    "ProblemType",
    "check_is_ProblemType"
]

# imports
from abc import abstractmethod
from numbers import Integral, Real
from typing import Callable, Iterable, Sequence, Tuple, Union
import numpy
import pymoo.core.problem

class ProblemType(pymoo.core.problem.Problem):
    """
    Base class for all optimization problems. This basal abstract class extends
    the PyMOO Problem class.

    The general formulation for an optimization problem should be:

    .. math::

        \\min_{\\mathbf{x}} \\mathbf{w_F \\odot F(x)}

    Such that:

    .. math::

        \\mathbf{w_G \\odot G(x) \\leq 0}

        \\mathbf{w_H \\odot H(x) = 0}
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for ProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(ProblemType, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ########################################################
    # Variables inherited from pymoo.core.problem.Problem ##
    ########################################################
    @property
    def n_var(self) -> Integral:
        """Number of decision variables."""
        return self._n_var
    @n_var.setter
    def n_var(self, value: Integral) -> None:
        """Set number of decision variables."""
        self._n_var = value

    @property
    def n_obj(self) -> Integral:
        """Number of objectives."""
        return self._n_obj
    @n_obj.setter
    def n_obj(self, value: Integral) -> None:
        """Set number of objectives."""
        self._n_obj = value
    
    @property
    def n_ieq_constr(self) -> Integral:
        """Number of inequality constraints."""
        return self._n_ieq_constr
    @n_ieq_constr.setter
    def n_ieq_constr(self, value: Integral) -> None:
        """Set number of inequality constraints."""
        self._n_ieq_constr = value
    
    @property
    def n_eq_constr(self) -> Integral:
        """Number of equality constraints."""
        return self._n_eq_constr
    @n_eq_constr.setter
    def n_eq_constr(self, value: Integral) -> None:
        """Set number of equality constraints."""
        self._n_eq_constr = value
    
    @property
    def xl(self) -> Union[numpy.ndarray,None]:
        """Lower boundary of the decision space."""
        return self._xl
    @xl.setter
    def xl(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set lower boundary of the decision space."""
        self._xl = value
    
    @property
    def xu(self) -> Union[numpy.ndarray,None]:
        """Upper boundary of the decision space."""
        return self._xu
    @xu.setter
    def xu(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set upper boundary of the decision space."""
        self._xu = value

    @property
    def vtype(self) -> Union[type,None]:
        """The variable type. So far, just used as a type hint."""
        return self._vtype
    @vtype.setter
    def vtype(self, value: Union[type,None]) -> None:
        """Set the variable type."""
        self._vtype = value
    
    @property
    def vars(self) -> Union[Sequence,None]:
        """Variables provided in their explicit form."""
        return self._vars
    @vars.setter
    def vars(self, value: Union[Sequence,None]) -> None:
        """Set variables provided in their explicit form."""
        self._vars = value
    
    @property
    def elementwise(self) -> bool:
        """Whether the evaluation function should be run elementwise."""
        return self._elementwise
    @elementwise.setter
    def elementwise(self, value: bool) -> None:
        """Set whether the evaluation function should be run elementwise."""
        self._elementwise = value
    
    @property
    def elementwise_func(self) -> type:
        """A class that creates the function that evaluates a single individual."""
        return self._elementwise_func
    @elementwise_func.setter
    def elementwise_func(self, value: type) -> None:
        """Set the class that creates the function that evaluates a single individual."""
        self._elementwise_func = value
    
    @property
    def elementwise_runner(self) -> Callable:
        """A function that runs the function that evaluates a single individual."""
        return self._elementwise_runner
    @elementwise_runner.setter
    def elementwise_runner(self, value: Callable) -> None:
        """Set the function that runs the function that evaluates a single individual."""
        self._elementwise_runner = value
    
    @property
    def replace_nan_values_by(self) -> Union[Real,None]:
        """Value for which to replace NaN values."""
        return self._replace_nan_values_by
    @replace_nan_values_by.setter
    def replace_nan_values_by(self, value: Union[Real,None]) -> None:
        """Set value for which to replace NaN values."""
        self._replace_nan_values_by = value
    
    @property
    def exclude_from_serialization(self) -> Union[Iterable,None]:
        """Attributes which are excluded from being serialized."""
        return self._exclude_from_serialization
    @exclude_from_serialization.setter
    def exclude_from_serialization(self, value: Union[Iterable,None]) -> None:
        """Set attributes which are excluded from being serialized."""
        self._exclude_from_serialization = value
    
    @property
    def callback(self) -> Union[Callable,None]:
        """A callback function to be called after every evaluation."""
        return self._callback
    @callback.setter
    def callback(self, value: Union[Callable,None]) -> None:
        """Set a callback function to be called after every evaluation."""
        self._callback = value
    
    @property
    def strict(self) -> bool:
        """Whether the shapes are checked strictly."""
        return self._strict
    @strict.setter
    def strict(self, value: bool) -> None:
        """Set whether the shapes are checked strictly."""
        self._strict = value
    
    @property
    def data(self) -> dict:
        """Type of the variable to be evaluated."""
        return self._data
    @data.setter
    def data(self, value: dict) -> None:
        """Set type of the variable to be evaluated."""
        self._data = value

    ########################################################
    ##### Properties unique to this ProblemType class ######
    ########################################################
    @property
    def ndecn(self) -> Integral:
        """Number of decision variables."""
        raise NotImplementedError("property is abstract")
    @ndecn.setter
    def ndecn(self, value: Integral) -> None:
        """Set number of decision variables."""
        raise NotImplementedError("property is abstract")
    
    @property
    def decn_space(self) -> numpy.ndarray:
        """Decision space boundaries."""
        raise NotImplementedError("property is abstract")
    @decn_space.setter
    def decn_space(self, value: numpy.ndarray) -> None:
        """Set decision space boundaries."""
        raise NotImplementedError("property is abstract")

    @property
    def decn_space_lower(self) -> Union[numpy.ndarray,None]:
        """Lower boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set lower boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    
    @property
    def decn_space_upper(self) -> Union[numpy.ndarray,None]:
        """Upper boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set upper boundary of the decision space."""
        raise NotImplementedError("property is abstract")

    @property
    def nobj(self) -> Integral:
        """Number of objectives."""
        raise NotImplementedError("property is abstract")
    @nobj.setter
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        raise NotImplementedError("property is abstract")
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        raise NotImplementedError("property is abstract")
    @obj_wt.setter
    def obj_wt(self, value: numpy.ndarray) -> None:
        """Set objective function weights."""
        raise NotImplementedError("property is abstract")

    @property
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")
    @nineqcv.setter
    def nineqcv(self, value: Integral) -> None:
        """Set number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: numpy.ndarray) -> None:
        """Set inequality constraint violation function weights."""
        raise NotImplementedError("property is abstract")

    @property
    def neqcv(self) -> Integral:
        """Number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    @neqcv.setter
    def neqcv(self, value: Integral) -> None:
        """Set number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @eqcv_wt.setter
    def eqcv_wt(self, value: numpy.ndarray) -> None:
        """Set equality constraint violation function weights."""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @abstractmethod
    def evalfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]:
        """
        Evaluate a candidate solution for the given ProblemType.

        This calculates three vectors which are to be minimized:

        .. math::

            \\mathbf{v_{obj}} = \\mathbf{w_{obj} \\odot F_{obj}(x)} \\
            \\mathbf{v_{ineqcv}} = \\mathbf{w_{ineqcv} \\odot G_{ineqcv}(x)} \\
            \\mathbf{v_{eqcv}} = \\mathbf{w_{eqcv} \\odot H_{eqcv}(x)}
        
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
            - ``ineqcv`` is a numpy.ndarray of shape ``(nineqcv,)`` that contains 
                inequality constraint violation values.
            - ``eqcv`` is a numpy.ndarray of shape ``(neqcv,)`` that contains 
                equality constraint violation values.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def _evaluate(
            self, 
            x: numpy.ndarray, 
            out: dict, 
            *args: tuple, 
            **kwargs: dict
        ) -> None:
        """
        Evaluate a set of candidate solutions for the given ProblemType.

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
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_ProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type ProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ProblemType):
        raise TypeError("'{0}' must be of type ProblemType.".format(vname))
