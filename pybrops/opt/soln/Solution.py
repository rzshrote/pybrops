"""
Module defining optimization solutions.
"""

# list of all public imports in the module
__all__ = [
    "Solution",
    "check_is_Solution"
]

# imports
from numbers import Integral, Number
from typing import Union
import numpy

class Solution:
    """
    Base class for all optimization problem solutions.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for Solution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(Solution, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def ndecn(self) -> Integral:
        """Number of decision variables."""
        raise NotImplementedError("property is abstract")
    @ndecn.setter
    def ndecn(self, value: Integral) -> None:
        """Set number of decision variables."""
        raise NotImplementedError("property is abstract")
    @ndecn.deleter
    def ndecn(self) -> None:
        """Delete number of decision variables."""
        raise NotImplementedError("property is abstract")
    
    @property
    def decn_space(self) -> numpy.ndarray:
        """Decision space boundaries."""
        raise NotImplementedError("property is abstract")
    @decn_space.setter
    def decn_space(self, value: numpy.ndarray) -> None:
        """Set decision space boundaries."""
        raise NotImplementedError("property is abstract")
    @decn_space.deleter
    def decn_space(self) -> None:
        """Delete decision space boundaries."""
        raise NotImplementedError("property is abstract")

    @property
    def decn_space_lower(self) -> Union[numpy.ndarray,None]:
        """Lower boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_lower.setter
    def decn_space_lower(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set lower boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_lower.deleter
    def decn_space_lower(self) -> None:
        """Delete lower boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    
    @property
    def decn_space_upper(self) -> Union[numpy.ndarray,None]:
        """Upper boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_upper.setter
    def decn_space_upper(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set upper boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_upper.deleter
    def decn_space_upper(self) -> None:
        """Delete upper boundary of the decision space."""
        raise NotImplementedError("property is abstract")

    @property
    def nobj(self) -> Integral:
        """Number of objectives."""
        raise NotImplementedError("property is abstract")
    @nobj.setter
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        raise NotImplementedError("property is abstract")
    @nobj.deleter
    def nobj(self) -> None:
        """Delete number of objectives."""
        raise NotImplementedError("property is abstract")
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        raise NotImplementedError("property is abstract")
    @obj_wt.setter
    def obj_wt(self, value: numpy.ndarray) -> None:
        """Set objective function weights."""
        raise NotImplementedError("property is abstract")
    @obj_wt.deleter
    def obj_wt(self) -> None:
        """Delete objective function weights."""
        raise NotImplementedError("property is abstract")

    @property
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")
    @nineqcv.setter
    def nineqcv(self, value: Integral) -> None:
        """Set number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")
    @nineqcv.deleter
    def nineqcv(self) -> None:
        """Delete number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: numpy.ndarray) -> None:
        """Set inequality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @ineqcv_wt.deleter
    def ineqcv_wt(self) -> None:
        """Delete inequality constraint violation function weights."""
        raise NotImplementedError("property is abstract")

    @property
    def neqcv(self) -> Integral:
        """Number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    @neqcv.setter
    def neqcv(self, value: Integral) -> None:
        """Set number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    @neqcv.deleter
    def neqcv(self) -> None:
        """Delete number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @eqcv_wt.setter
    def eqcv_wt(self, value: numpy.ndarray) -> None:
        """Set equality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @eqcv_wt.deleter
    def eqcv_wt(self) -> None:
        """Delete equality constraint violation function weights."""
        raise NotImplementedError("property is abstract")

    @property
    def nsoln(self) -> Integral:
        """Number of solutions to the problem."""
        raise NotImplementedError("property is abstract")
    @nsoln.setter
    def nsoln(self, value: Integral) -> None:
        """Set number of solutions to the problem."""
        raise NotImplementedError("property is abstract")
    @nsoln.deleter
    def nsoln(self) -> None:
        """Delete number of solutions to the problem."""
        raise NotImplementedError("property is abstract")
    
    @property
    def soln_decn(self) -> numpy.ndarray:
        """Matrix of solution vectors in the decision space."""
        raise NotImplementedError("property is abstract")
    @soln_decn.setter
    def soln_decn(self, value: numpy.ndarray) -> None:
        """Set matrix of solution vectors in the decision space."""
        raise NotImplementedError("property is abstract")
    @soln_decn.deleter
    def soln_decn(self) -> None:
        """Delete matrix of solution vectors in the decision space."""
        raise NotImplementedError("property is abstract")
    
    @property
    def soln_obj(self) -> numpy.ndarray:
        """Solution objective function values."""
        raise NotImplementedError("property is abstract")
    @soln_obj.setter
    def soln_obj(self, value: numpy.ndarray) -> None:
        """Set solution objective function values."""
        raise NotImplementedError("property is abstract")
    @soln_obj.deleter
    def soln_obj(self) -> None:
        """Delete solution objective function values."""
        raise NotImplementedError("property is abstract")
    
    @property
    def soln_ineqcv(self) -> numpy.ndarray:
        """Solution inequality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    @soln_ineqcv.setter
    def soln_ineqcv(self, value: numpy.ndarray) -> None:
        """Set solution inequality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    @soln_ineqcv.deleter
    def soln_ineqcv(self) -> None:
        """Delete solution inequality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    
    @property
    def soln_eqcv(self) -> numpy.ndarray:
        """Solution equality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    @soln_eqcv.setter
    def soln_eqcv(self, value: numpy.ndarray) -> None:
        """Set solution equality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    @soln_eqcv.deleter
    def soln_eqcv(self) -> None:
        """Delete solution equality constraint violation function values."""
        raise NotImplementedError("property is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
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
