"""
Implementation of the Solution interface.
"""

# list of all public objects in this module
__all__ = [
    "Solution",
    "check_is_Solution",
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_len_eq, check_ndarray_ndim, check_ndarray_shape_eq
from pybrops.core.error.error_value_python import check_is_gteq


class Solution(metaclass=ABCMeta):
    """
    A semi-abstract class defining a Solution interface and implementing several
    essential properties for all Solution classes.

    A user must implement the following abstract methods in derivatives:
        1) ``__init__``

    Notes:
        1) It is possible to call the constructor of this semi-abstract from a
           derived class.
    """

    ########################## Special Object Methods ##########################
    # do not implement __init__() since this is an interface/mixin class

    ############################ Object Properties #############################

    ############## Decision space properties ###############
    @property
    def ndecn(self) -> Integral:
        """Number of decision variables."""
        return self._ndecn
    @ndecn.setter
    def ndecn(self, value: Integral) -> None:
        """Set number of decision variables."""
        check_is_Integral(value, "ndecn")
        check_is_gteq(value, "ndecn", 1)
        self._ndecn = value
    
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
            raise TypeError("'obj_wt' must be of type numpy.ndarray or a Real type")
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

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        return self._ineqcv_wt
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set inequality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "ineqcv_wt", 1)
            check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.nineqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.nineqcv)
        else:
            raise TypeError("'ineqcv_wt' must be of type numpy.ndarray or a Real type")
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
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        return self._eqcv_wt
    @eqcv_wt.setter
    def eqcv_wt(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set equality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "eqcv_wt", 1)
            check_ndarray_len_eq(value, "eqcv_wt", self.neqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.neqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.neqcv)
        else:
            raise TypeError("'eqcv_wt' must be of type numpy.ndarray or a Real type")
        self._eqcv_wt = value

    ################# Solution properties ##################
    @property
    def nsoln(self) -> Integral:
        """Number of solutions to the problem."""
        return self._nsoln
    @nsoln.setter
    def nsoln(self, value: Integral) -> None:
        """Set number of solutions to the problem."""
        check_is_Integral(value, "nsoln")
        check_is_gteq(value, "nsoln", 0)    # possible to have no solutions
        self._nsoln = value

    @property
    def soln_decn(self) -> numpy.ndarray:
        """Matrix of solution vectors in the decision space."""
        return self._soln
    @soln_decn.setter
    def soln_decn(self, value: numpy.ndarray) -> None:
        """Set matrix of solution vectors in the decision space."""
        check_is_ndarray(value, "soln")
        check_ndarray_ndim(value, "soln", 2)
        check_ndarray_shape_eq(value, "soln", (self.nsoln,self.ndecn))
        self._soln = value
    
    @property
    def soln_obj(self) -> numpy.ndarray:
        """Solution objective function values."""
        return self._soln_obj
    @soln_obj.setter
    def soln_obj(self, value: numpy.ndarray) -> None:
        """Set solution objective function values."""
        check_is_ndarray(value, "soln_obj")
        check_ndarray_ndim(value, "soln_obj", 2)
        check_ndarray_shape_eq(value, "soln", (self.nsoln,self.nobj))
        self._soln_obj = value
    
    @property
    def soln_ineqcv(self) -> numpy.ndarray:
        """Solution inequality constraint violation function values."""
        return self._soln_ineqcv
    @soln_ineqcv.setter
    def soln_ineqcv(self, value: Union[numpy.ndarray,None]) -> None:
        """Set solution inequality constraint violation function values."""
        if value is None:
            value = numpy.zeros((self.nsoln,self.nineqcv), dtype = self.soln_obj.dtype)
        check_is_ndarray(value, "soln_ineqcv")
        check_ndarray_ndim(value, "soln_ineqcv", 2)
        check_ndarray_shape_eq(value, "soln_ineqcv", (self.nsoln,self.nineqcv))
        self._soln_ineqcv = value
    
    @property
    def soln_eqcv(self) -> numpy.ndarray:
        """Solution equality constraint violation function values."""
        return self._soln_eqcv
    @soln_eqcv.setter
    def soln_eqcv(self, value: Union[numpy.ndarray,None]) -> None:
        """Set solution equality constraint violation function values."""
        if value is None:
            value = numpy.zeros((self.nsoln,self.nineqcv), dtype = self.soln_obj.dtype)
        check_is_ndarray(value, "soln_eqcv")
        check_ndarray_ndim(value, "soln_eqcv", 2)
        check_ndarray_shape_eq(value, "soln_eqcv", (self.nsoln,self.neqcv))
        self._soln_eqcv = value



################################## Utilities ###################################
def check_is_Solution(v: object, vname: str) -> None:
    """
    Check if object is of type Solution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, Solution):
        raise TypeError("'{0}' must be of type Solution.".format(vname))
