"""
Partial implementation of the Problem interface.
"""

# list of public objects in this module
__all__ = [
    "DenseProblem",
    "check_is_DenseProblem"
]

# imports
from numbers import Integral
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_is_1d, check_ndarray_len_eq, check_ndarray_len_gteq
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.opt.prob.Problem import Problem

class DenseProblem(Problem):
    """
    docstring for DenseProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            ndecn: Integral,
            decn_space: numpy.ndarray,
            nobj: Integral,
            obj_wt: numpy.ndarray,
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseProblem, self).__init__(**kwargs)
        # order dependent assignments
        self.ndecn = ndecn
        self.decn_space = decn_space
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
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
    @ndecn.deleter
    def ndecn(self) -> None:
        """Delete number of decision variables."""
        del self._ndecn
    
    @property
    def decn_space(self) -> numpy.ndarray:
        """Decision space boundaries."""
        return self._decn_space
    @decn_space.setter
    def decn_space(self, value: numpy.ndarray) -> None:
        """Set decision space boundaries."""
        check_is_ndarray(value, "decn_space")
        self._decn_space = value
    @decn_space.deleter
    def decn_space(self) -> None:
        """Delete decision space boundaries."""
        del self._decn_space

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
    @nobj.deleter
    def nobj(self) -> None:
        """Delete number of objectives."""
        del self._nobj
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        return self._obj_wt
    @obj_wt.setter
    def obj_wt(self, value: numpy.ndarray) -> None:
        """Set objective function weights."""
        check_is_ndarray(value, "obj_wt")
        check_ndarray_is_1d(value, "obj_wt")
        check_ndarray_len_eq(value, "obj_wt", self.nobj)
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
    def nineqcv(self, value: Integral) -> None:
        """Set number of inequality constraint violation functions."""
        check_is_Integral(value, "nineqcv")
        check_is_gteq(value, "nineqcv", 0)  # possible to have 0 inequality constraints
        self._nineqcv = value
    @nineqcv.deleter
    def nineqcv(self) -> None:
        """Delete number of inequality constraint violation functions."""
        del self._nineqcv

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        return self._ineqcv_wt
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: numpy.ndarray) -> None:
        """Set inequality constraint violation function weights."""
        check_is_ndarray(value, "ineqcv_wt")
        check_ndarray_is_1d(value, "ineqcv_wt")
        check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
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
    def neqcv(self, value: Integral) -> None:
        """Set number of equality constraint violations."""
        check_is_Integral(value, "neqcv")
        check_is_gteq(value, "neqcv", 0)    # possible to have 0 equality constraints
        self._neqcv = value
    @neqcv.deleter
    def neqcv(self) -> None:
        """Delete number of equality constraint violations."""
        del self._neqcv
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        return self._eqcv_wt
    @eqcv_wt.setter
    def eqcv_wt(self, value: numpy.ndarray) -> None:
        """Set equality constraint violation function weights."""
        check_is_ndarray(value, "eqcv_wt")
        check_ndarray_is_1d(value, "eqcv_wt")
        check_is_gteq(value, "eqcv_wt", self.neqcv)
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
def check_is_DenseProblem(v: object, vname: str) -> None:
    """
    Check if object is of type DenseProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseProblem):
        raise TypeError("'{0}' must be of type DenseProblem.".format(vname))
