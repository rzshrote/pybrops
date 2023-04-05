"""
Implementation of the Solution interface.
"""

# list of all public objects in this module
__all__ = [
    "DenseSolution",
    "check_is_DenseSolution"
]

from numbers import Integral, Number
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_is_1d, check_ndarray_is_2d, check_ndarray_len_eq, check_ndarray_shape_eq
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.opt.soln.Solution import Solution


class DenseSolution(Solution):
    """
    docstring for DenseSolution.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            ndecn: Integral,
            decn_space: numpy.ndarray,
            decn_space_lower: numpy.ndarray,
            decn_space_upper: numpy.ndarray,
            nobj: Integral,
            obj_wt: numpy.ndarray,
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            nsoln: Integral,
            soln_decn: numpy.ndarray,
            soln_obj: numpy.ndarray,
            soln_ineqcv: numpy.ndarray,
            soln_eqcv: numpy.ndarray,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseSolution, self).__init__(**kwargs)
        # order dependent assignments
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
        self.nsoln = nsoln
        self.soln_decn = soln_decn
        self.soln_obj = soln_obj
        self.soln_ineqcv = soln_ineqcv
        self.soln_eqcv = soln_eqcv

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
    def decn_space_lower(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set lower boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "xl", self.ndecn)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_lower' must be of type numpy.ndarray, Number, or None")
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
    def decn_space_upper(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set upper boundary of the decision space."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_len_eq(value, "xl", self.ndecn)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.ndecn)
        elif value is None:
            pass
        else:
            raise TypeError("'decn_space_upper' must be of type numpy.ndarray, Number, or None")
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
    @nobj.deleter
    def nobj(self) -> None:
        """Delete number of objectives."""
        del self._nobj
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        return self._obj_wt
    @obj_wt.setter
    def obj_wt(self, value: Union[numpy.ndarray,Number]) -> None:
        """Set objective function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "obj_wt")
            check_ndarray_len_eq(value, "obj_wt", self.nobj)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.nobj)
        else:
            raise TypeError("'obj_wt' must be of type numpy.ndarray or a numeric type")
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
    def ineqcv_wt(self, value: Union[numpy.ndarray,Number]) -> None:
        """Set inequality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "ineqcv_wt")
            check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.nineqcv)
        else:
            raise TypeError("'ineqcv_wt' must be of type numpy.ndarray or a numeric type")
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
    def eqcv_wt(self, value: Union[numpy.ndarray,Number]) -> None:
        """Set equality constraint violation function weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_is_1d(value, "eqcv_wt")
            check_ndarray_len_eq(value, "eqcv_wt", self.neqcv)
        elif isinstance(value, Number):
            value = numpy.repeat(value, self.neqcv)
        else:
            raise TypeError("'eqcv_wt' must be of type numpy.ndarray or a numeric type")
        self._eqcv_wt = value
    @eqcv_wt.deleter
    def eqcv_wt(self) -> None:
        """Delete equality constraint violation function weights."""
        del self._eqcv_wt

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
    @nsoln.deleter
    def nsoln(self) -> None:
        """Delete number of solutions to the problem."""
        del self._nsoln

    @property
    def soln_decn(self) -> numpy.ndarray:
        """Matrix of solution vectors in the decision space."""
        return self._soln
    @soln_decn.setter
    def soln_decn(self, value: numpy.ndarray) -> None:
        """Set matrix of solution vectors in the decision space."""
        check_is_ndarray(value, "soln")
        check_ndarray_is_2d(value, "soln")
        check_ndarray_shape_eq(value, "soln", (self.nsoln,self.ndecn))
        self._soln = value
    @soln_decn.deleter
    def soln_decn(self) -> None:
        """Delete matrix of solution vectors in the decision space."""
        del self._soln
    
    @property
    def soln_obj(self) -> numpy.ndarray:
        """Solution objective function values."""
        return self._soln_obj
    @soln_obj.setter
    def soln_obj(self, value: numpy.ndarray) -> None:
        """Set solution objective function values."""
        check_is_ndarray(value, "soln_obj")
        check_ndarray_is_2d(value, "soln_obj")
        check_ndarray_shape_eq(value, "soln", (self.nsoln,self.nobj))
        self._soln_obj = value
    @soln_obj.deleter
    def soln_obj(self) -> None:
        """Delete solution objective function values."""
        del self._soln_obj
    
    @property
    def soln_ineqcv(self) -> numpy.ndarray:
        """Solution inequality constraint violation function values."""
        return self._soln_ineqcv
    @soln_ineqcv.setter
    def soln_ineqcv(self, value: numpy.ndarray) -> None:
        """Set solution inequality constraint violation function values."""
        check_is_ndarray(value, "soln_ineqcv")
        check_ndarray_is_2d(value, "soln_ineqcv")
        check_ndarray_shape_eq(value, "soln_ineqcv", (self.nsoln,self.nineqcv))
        self._soln_ineqcv = value
    @soln_ineqcv.deleter
    def soln_ineqcv(self) -> None:
        """Delete solution inequality constraint violation function values."""
        del self._soln_ineqcv
    
    @property
    def soln_eqcv(self) -> numpy.ndarray:
        """Solution equality constraint violation function values."""
        return self._soln_eqcv
    @soln_eqcv.setter
    def soln_eqcv(self, value: numpy.ndarray) -> None:
        """Set solution equality constraint violation function values."""
        check_is_ndarray(value, "soln_eqcv")
        check_ndarray_is_2d(value, "soln_eqcv")
        check_ndarray_shape_eq(value, "soln_eqcv", (self.nsoln,self.neqcv))
        self._soln_eqcv = value
    @soln_eqcv.deleter
    def soln_eqcv(self) -> None:
        """Delete solution equality constraint violation function values."""
        del self._soln_eqcv

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
def check_is_DenseSolution(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSolution):
        raise TypeError("'{0}' must be of type DenseSolution.".format(vname))
