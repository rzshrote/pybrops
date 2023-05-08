"""
Partial implementation of the Problem interface.
"""

# list of public objects in this module
__all__ = [
    "Problem",
    "check_is_Problem"
]

# imports
from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Container, Iterable, Optional, Tuple, Union
import numpy
from pybrops.core.error.error_type_python import check_is_Integral, check_is_type
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_ndim, check_ndarray_shape_eq
from pybrops.core.error.error_value_python import check_is_gteq
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation
import pymoo.core.problem

class Problem(pymoo.core.problem.Problem,metaclass=ABCMeta):
    """
    A semi-abstract base class for representing all optimization problems.
    This basal semi-abstract class extends the PyMOO Problem class.

    The general formulation for an optimization problem should be:

    .. math::

        \\min_{\\mathbf{x}} \\mathbf{w_F \\odot F(x)}

    Such that:

    .. math::

        \\mathbf{w_G \\odot G(x) \\leq 0}

        \\mathbf{w_H \\odot H(x) = 0}

    A user must implement the following abstract methods in derivatives:
        1) ``__init__``
        2) ``evalfn``
        3) ``_evaluate``

    Notes:
        1) It is possible to call the constructor of this semi-abstract from a
           derived class.
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
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            vtype: Optional[type] = None,
            vars: Optional[Container] = None,
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
        Constructor for Problem.
        
        Parameters
        ----------
        n_var : int
            Number of variables.
        n_obj : int
            Number of objectives.
        n_ieq_constr : int
            Number of Inequality Constraints
        n_eq_constr : int
            Number of Equality Constraints
        xl : np.array, float, int
            Lower bounds for the variables. if integer all lower bounds are equal.
        xu : np.array, float, int
            Upper bounds for the variable. if integer all upper bounds are equal.
        vtype : type
            The variable type. So far, just used as a type hint.
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

        # call PyMOO constructor to set things their way (for PyMOO interface)
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

    # Variables inherited from pymoo.core.problem.Problem ##
    @property
    def n_var(self) -> Integral:
        """Number of decision variables."""
        return self._n_var
    @n_var.setter
    def n_var(self, value: Integral) -> None:
        """Set number of decision variables."""
        check_is_Integral(value, "n_var")
        self._n_var = value
        self._ndecn = value # set ndecn too; used for easy separation from PyMOO
    
    @property
    def n_obj(self) -> Integral:
        """Number of objectives."""
        return self._n_obj
    @n_obj.setter
    def n_obj(self, value: Integral) -> None:
        """Set number of objectives."""
        check_is_Integral(value, "n_obj")
        check_is_gteq(value, "n_obj", 1)
        self._n_obj = value
        self._nobj = value # set nobj too; used for easy separation from PyMOO
    
    @property
    def n_ieq_constr(self) -> Integral:
        """Number of inequality constraints."""
        return self._n_ieq_constr
    @n_ieq_constr.setter
    def n_ieq_constr(self, value: Union[Integral,None]) -> None:
        """Set number of inequality constraints."""
        if value is None:
            value = 0
        check_is_Integral(value, "n_ieq_constr")
        check_is_gteq(value, "n_ieq_constr", 0)
        self._n_ieq_constr = value
        self._nineqcv = value # used for easy separation from PyMOO
    
    @property
    def n_eq_constr(self) -> Integral:
        """n_eq_constr."""
        return self._n_eq_constr
    @n_eq_constr.setter
    def n_eq_constr(self, value: Union[Integral,None]) -> None:
        """Set n_eq_constr."""
        if value is None:
            value = 0
        check_is_Integral(value, "n_eq_constr")
        check_is_gteq(value, "n_eq_constr", 0)
        self._n_eq_constr = value
        self._neqcv = value # used for easy separation from PyMOO
    
    @property
    def xl(self) -> Union[numpy.ndarray,None]:
        """Lower boundary of the decision space."""
        return self._xl
    @xl.setter
    def xl(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "xl", self.n_var)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.n_var)
        elif value is None:
            pass
        else:
            raise TypeError("'xl' must be of type numpy.ndarray, Real, or None")
        self._xl = value
        self._decn_space_lower = value
    
    @property
    def xu(self) -> Union[numpy.ndarray,None]:
        """Upper boundary of the decision space."""
        return self._xu
    @xu.setter
    def xu(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "xl", self.n_var)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.n_var)
        elif value is None:
            pass
        else:
            raise TypeError("'xu' must be of type numpy.ndarray, Real, or None")
        self._xu = value
        self._decn_space_upper = value

    @property
    def vtype(self) -> Union[type,None]:
        """The variable type. So far, just used as a type hint."""
        return self._vtype
    @vtype.setter
    def vtype(self, value: Union[type,None]) -> None:
        """Set the variable type."""
        if (not isinstance(value, type)) and (value is not None):
            raise TypeError("'vtype' must be of type type or None")
        self._vtype = value
    
    @property
    def vars(self) -> Union[Container,None]:
        """Variables provided in their explicit form."""
        return self._vars
    @vars.setter
    def vars(self, value: Union[Container,None]) -> None:
        """Set variables provided in their explicit form."""
        if (not isinstance(value, Container)) and (value is not None):
            raise TypeError("'vars' must be a Container or None")
        self._vars = value
    
    @property
    def elementwise(self) -> bool:
        """Whether the evaluation function should be run elementwise."""
        return self._elementwise
    @elementwise.setter
    def elementwise(self, value: bool) -> None:
        """Set whether the evaluation function should be run elementwise."""
        if not isinstance(value, bool):
            raise TypeError("'elementwise' must be of type bool")
        self._elementwise = value
    
    @property
    def elementwise_func(self) -> type:
        """A class that creates the function that evaluates a single individual."""
        return self._elementwise_func
    @elementwise_func.setter
    def elementwise_func(self, value: type) -> None:
        """Set the class that creates the function that evaluates a single individual."""
        if not isinstance(value, type):
            raise TypeError("'elementwise_func' must be a type")
        check_is_type(value, "elementwise_func")
        self._elementwise_func = value
    
    @property
    def elementwise_runner(self) -> Callable:
        """A function that runs the function that evaluates a single individual."""
        return self._elementwise_runner
    @elementwise_runner.setter
    def elementwise_runner(self, value: Callable) -> None:
        """Set the function that runs the function that evaluates a single individual."""
        if not isinstance(value, Callable):
            raise TypeError("'elementwise_runner' must be a Callable type")
        self._elementwise_runner = value
    
    @property
    def replace_nan_values_by(self) -> Union[Real,None]:
        """replace_nan_values_by."""
        return self._replace_nan_values_by
    @replace_nan_values_by.setter
    def replace_nan_values_by(self, value: Union[Real,None]) -> None:
        """Set replace_nan_values_by."""
        if (not isinstance(value, Real)) and (value is not None):
            raise TypeError("'replace_nan_values_by' must be a Real type")
        self._replace_nan_values_by = value
    
    @property
    def exclude_from_serialization(self) -> Union[Iterable,None]:
        """attributes which are excluded from being serialized."""
        return self._exclude_from_serialization
    @exclude_from_serialization.setter
    def exclude_from_serialization(self, value: Union[Iterable,None]) -> None:
        """Set attributes which are excluded from being serialized."""
        if (not isinstance(value, Iterable)) and (value is not None):
            raise TypeError("'exclude_from_serialization' must be an Iterable type or None")
        self._exclude_from_serialization = value
    
    @property
    def callback(self) -> Union[Callable,None]:
        """A callback function to be called after every evaluation."""
        return self._callback
    @callback.setter
    def callback(self, value: Union[Callable,None]) -> None:
        """Set a callback function to be called after every evaluation."""
        if (not isinstance(value, Callable)) and (value is not None):
            raise TypeError("'callback' must be a Callable type or None")
        self._callback = value
    
    @property
    def strict(self) -> bool:
        """Whether the shapes are checked strictly."""
        return self._strict
    @strict.setter
    def strict(self, value: bool) -> None:
        """Set whether the shapes are checked strictly."""
        if not isinstance(value, bool):
            raise TypeError("'strict' must be of type bool")
        self._strict = value
    
    @property
    def data(self) -> dict:
        """Type of the variable to be evaluated."""
        return self._data
    @data.setter
    def data(self, value: dict) -> None:
        """Set type of the variable to be evaluated."""
        if not isinstance(value, dict):
            raise TypeError("'data' must be of type dict")
        self._data = value

    ####### Properties unique to this Problem class ########

    ############## Decision space properties ###############
    @property
    def ndecn(self) -> Integral:
        """Number of decision variables."""
        return self._ndecn
    @ndecn.setter
    def ndecn(self, value: Integral) -> None:
        """Set number of decision variables."""
        check_is_Integral(value, "ndecn")
        #check_is_gteq(value, "ndecn", 1) # PyMOO allows for -1
        self._ndecn = value
        self._n_var = value # used for easy separation from PyMOO
    
    @property
    def decn_space(self) -> Union[numpy.ndarray,None]:
        """Decision space boundaries."""
        return self._decn_space
    @decn_space.setter
    def decn_space(self, value: Union[numpy.ndarray,None]) -> None:
        """Set decision space boundaries."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_shape_eq(value, "decn_space", (2,self.ndecn))
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space' must be of type numpy.ndarray or None")
        self._decn_space = value

    @property
    def decn_space_lower(self) -> Union[numpy.ndarray,None]:
        """Lower boundary of the decision space."""
        return self._decn_space_lower
    @decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "xl", self.ndecn)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_lower' must be of type numpy.ndarray, Real, or None")
        self._decn_space_lower = value
        self._xl = value
    
    @property
    def decn_space_upper(self) -> Union[numpy.ndarray,None]:
        """Upper boundary of the decision space."""
        return self._decn_space_upper
    @decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "xl", self.ndecn)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_upper' must be of type numpy.ndarray, Real, or None")
        self._decn_space_upper = value
        self._xu = value

    ############## Objective space properties ##############
    @property
    def nobj(self) -> Integral:
        """Number of objectives."""
        return self._nobj
    @nobj.setter
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        check_is_Integral(value, "nobj")
        check_is_gteq(value, "nobj", 1)     # cannot have 0 objectives
        self._nobj = value
        self._n_obj = value # for easy separation from PyMOO
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        return self._obj_wt
    @obj_wt.setter
    def obj_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set objective function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "obj_wt", 1)
            check_ndarray_len_eq(value, "obj_wt", self.nobj)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.nobj)
        elif value is None:
            value = numpy.repeat(1.0, self.nobj)
        else:
            raise TypeError("'obj_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._obj_wt = value

    ######## Inequality constraint space properties ########
    @property
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violation functions."""
        return self._nineqcv
    @nineqcv.setter
    def nineqcv(self, value: Union[Integral,None]) -> None:
        """Set number of inequality constraint violation functions."""
        if value is None:
            value = 0
        check_is_Integral(value, "nineqcv")
        check_is_gteq(value, "nineqcv", 0)  # possible to have 0 inequality constraints
        self._nineqcv = value
        self._n_ieq_constr = value # for easy separation from PyMOO

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        return self._ineqcv_wt
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set inequality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "ineqcv_wt", 1)
            check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.nineqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.nineqcv)
        else:
            raise TypeError("'ineqcv_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._ineqcv_wt = value

    ######### Equality constraint space properties #########
    @property
    def neqcv(self) -> Integral:
        """Number of equality constraint violations."""
        return self._neqcv
    @neqcv.setter
    def neqcv(self, value: Union[Integral,None]) -> None:
        """Set number of equality constraint violations."""
        if value is None:
            value = 0
        check_is_Integral(value, "neqcv")
        check_is_gteq(value, "neqcv", 0)    # possible to have 0 equality constraints
        self._neqcv = value
        self._n_eq_constr = value # for easy separation from PyMOO
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        return self._eqcv_wt
    @eqcv_wt.setter
    def eqcv_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set equality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "eqcv_wt", 1)
            check_ndarray_len_eq(value, "eqcv_wt", self.neqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.neqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.neqcv)
        else:
            raise TypeError("'eqcv_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._eqcv_wt = value

    ############################## Object Methods ##############################

    ### method required by PyBrOpS interface ###
    @abstractmethod
    def evalfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]:
        """
        Evaluate a candidate solution for the given Problem.

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

    ### method required by PyMOO interface ###
    @abstractmethod
    def _evaluate(
            self, 
            x: numpy.ndarray, 
            out: dict, 
            *args: tuple, 
            **kwargs: dict
        ) -> None:
        """
        Evaluate a set of candidate solutions for the given Problem.

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



################################## Utilities ###################################
def check_is_Problem(v: object, vname: str) -> None:
    """
    Check if object is of type Problem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, Problem):
        raise TypeError("'{0}' must be of type Problem.".format(vname))
