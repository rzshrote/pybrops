"""
Partial implementation of the Problem interface.
"""

# list of public objects in this module
__all__ = [
    "Problem",
    "check_is_Problem"
]

# imports
from numbers import Integral, Real
from typing import Callable, Container, Iterable, Optional, Union
import numpy
from pybrops.core.error.error_type_python import check_is_Integral, check_is_type
from pybrops.core.error.error_value_numpy import check_ndarray_is_1d, check_ndarray_len_eq, check_ndarray_shape_eq
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.opt.prob.ProblemType import ProblemType
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation

class Problem(ProblemType):
    """
    docstring for Problem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
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
        check_is_Integral(value, "n_var")
        self._n_var = value
        self._ndecn = value # set ndecn too; used for easy separation from PyMOO
    @n_var.deleter
    def n_var(self) -> None:
        """Delete number of decision variables."""
        del self._n_var
    
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
    @n_obj.deleter
    def n_obj(self) -> None:
        """Delete number of objectives."""
        del self._n_obj
    
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
    @n_ieq_constr.deleter
    def n_ieq_constr(self) -> None:
        """Delete number of inequality constraints."""
        del self._n_ieq_constr
    
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
    @n_eq_constr.deleter
    def n_eq_constr(self) -> None:
        """Delete n_eq_constr."""
        del self._n_eq_constr
    
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
    @xl.deleter
    def xl(self) -> None:
        """Delete lower boundary of the decision space."""
        del self._xl
    
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
    @xu.deleter
    def xu(self) -> None:
        """Delete upper boundary of the decision space."""
        del self._xu

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
    @vtype.deleter
    def vtype(self) -> None:
        """Delete the variable type."""
        del self._vtype
    
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
    @vars.deleter
    def vars(self) -> None:
        """Delete variables provided in their explicit form."""
        del self._vars
    
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
    @elementwise.deleter
    def elementwise(self) -> None:
        """Delete whether the evaluation function should be run elementwise."""
        del self._elementwise
    
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
    @elementwise_func.deleter
    def elementwise_func(self) -> None:
        """Delete the class that creates the function that evaluates a single individual."""
        del self._elementwise_func
    
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
    @elementwise_runner.deleter
    def elementwise_runner(self) -> None:
        """Delete the function that runs the function that evaluates a single individual."""
        del self._elementwise_runner
    
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
    @replace_nan_values_by.deleter
    def replace_nan_values_by(self) -> None:
        """Delete replace_nan_values_by."""
        del self._replace_nan_values_by
    
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
    @exclude_from_serialization.deleter
    def exclude_from_serialization(self) -> None:
        """Delete attributes which are excluded from being serialized."""
        del self._exclude_from_serialization
    
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
    @callback.deleter
    def callback(self) -> None:
        """Delete a callback function to be called after every evaluation."""
        del self._callback
    
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
    @strict.deleter
    def strict(self) -> None:
        """Delete whether the shapes are checked strictly."""
        del self._strict
    
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
    @data.deleter
    def data(self) -> None:
        """Delete type of the variable to be evaluated."""
        del self._data

    ########################################################
    ####### Properties unique to this Problem class ########
    ########################################################
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
    @ndecn.deleter
    def ndecn(self) -> None:
        """Delete number of decision variables."""
        del self._ndecn
    
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
            raise TypeError("'decn_space' must be of type numpy.ndarrray or None")
        self._decn_space = value
    @decn_space.deleter
    def decn_space(self) -> None:
        """Delete decision space boundaries."""
        del self._decn_space

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
    @decn_space_lower.deleter
    def decn_space_lower(self) -> None:
        """Delete lower boundary of the decision space."""
        del self._decn_space_lower
    
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
    @decn_space_upper.deleter
    def decn_space_upper(self) -> None:
        """Delete upper boundary of the decision space."""
        del self._decn_space_upper

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
    @nobj.deleter
    def nobj(self) -> None:
        """Delete number of objectives."""
        del self._nobj
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        return self._obj_wt
    @obj_wt.setter
    def obj_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set objective function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "obj_wt")
            check_ndarray_len_eq(value, "obj_wt", self.nobj)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.nobj)
        elif value is None:
            value = numpy.repeat(1.0, self.nobj)
        else:
            raise TypeError("'obj_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._obj_wt = value
    @obj_wt.deleter
    def obj_wt(self) -> None:
        """Delete objective function weights."""
        del self._obj_wt

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
    @nineqcv.deleter
    def nineqcv(self) -> None:
        """Delete number of inequality constraint violation functions."""
        del self._nineqcv

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        return self._ineqcv_wt
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set inequality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "ineqcv_wt")
            check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.nineqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.nineqcv)
        else:
            raise TypeError("'ineqcv_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._ineqcv_wt = value
    @ineqcv_wt.deleter
    def ineqcv_wt(self) -> None:
        """Delete inequality constraint violation function weights."""
        del self._ineqcv_wt

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
    @neqcv.deleter
    def neqcv(self) -> None:
        """Delete number of equality constraint violations."""
        del self._neqcv
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        return self._eqcv_wt
    @eqcv_wt.setter
    def eqcv_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set equality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "eqcv_wt")
            check_ndarray_len_eq(value, "eqcv_wt", self.neqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.neqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.neqcv)
        else:
            raise TypeError("'eqcv_wt' must be of type numpy.ndarray, a numeric type, or None")
        self._eqcv_wt = value
    @eqcv_wt.deleter
    def eqcv_wt(self) -> None:
        """Delete equality constraint violation function weights."""
        del self._eqcv_wt

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
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
