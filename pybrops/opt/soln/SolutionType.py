"""
Module defining optimization solutions.
"""

# list of all public imports in the module
__all__ = [
    "SolutionType",
    "check_is_SolutionType"
]

# imports
from abc import ABCMeta, abstractmethod
from numbers import Integral, Number
from typing import Union
import numpy

class SolutionType(metaclass=ABCMeta):
    """
    Abstract Base Class for all optimization problem solutions.
    """
    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    @abstractmethod
    def ndecn(self) -> Integral:
        """Number of decision variables."""
        raise NotImplementedError("property is abstract")
    @ndecn.setter
    @abstractmethod
    def ndecn(self, value: Integral) -> None:
        """Set number of decision variables."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def decn_space(self) -> numpy.ndarray:
        """Decision space boundaries."""
        raise NotImplementedError("property is abstract")
    @decn_space.setter
    @abstractmethod
    def decn_space(self, value: numpy.ndarray) -> None:
        """Set decision space boundaries."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def decn_space_lower(self) -> Union[numpy.ndarray,None]:
        """Lower boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_lower.setter
    @abstractmethod
    def decn_space_lower(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set lower boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def decn_space_upper(self) -> Union[numpy.ndarray,None]:
        """Upper boundary of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_space_upper.setter
    @abstractmethod
    def decn_space_upper(self, value: Union[numpy.ndarray,Number,None]) -> None:
        """Set upper boundary of the decision space."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nobj(self) -> Integral:
        """Number of objectives."""
        raise NotImplementedError("property is abstract")
    @nobj.setter
    @abstractmethod
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        raise NotImplementedError("property is abstract")
    @obj_wt.setter
    @abstractmethod
    def obj_wt(self, value: numpy.ndarray) -> None:
        """Set objective function weights."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")
    @nineqcv.setter
    @abstractmethod
    def nineqcv(self, value: Integral) -> None:
        """Set number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @ineqcv_wt.setter
    @abstractmethod
    def ineqcv_wt(self, value: numpy.ndarray) -> None:
        """Set inequality constraint violation function weights."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def neqcv(self) -> Integral:
        """Number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    @neqcv.setter
    @abstractmethod
    def neqcv(self, value: Integral) -> None:
        """Set number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        raise NotImplementedError("property is abstract")
    @eqcv_wt.setter
    @abstractmethod
    def eqcv_wt(self, value: numpy.ndarray) -> None:
        """Set equality constraint violation function weights."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def nsoln(self) -> Integral:
        """Number of solutions to the problem."""
        raise NotImplementedError("property is abstract")
    @nsoln.setter
    @abstractmethod
    def nsoln(self, value: Integral) -> None:
        """Set number of solutions to the problem."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def soln_decn(self) -> numpy.ndarray:
        """Matrix of solution vectors in the decision space."""
        raise NotImplementedError("property is abstract")
    @soln_decn.setter
    @abstractmethod
    def soln_decn(self, value: numpy.ndarray) -> None:
        """Set matrix of solution vectors in the decision space."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def soln_obj(self) -> numpy.ndarray:
        """Solution objective function values."""
        raise NotImplementedError("property is abstract")
    @soln_obj.setter
    @abstractmethod
    def soln_obj(self, value: numpy.ndarray) -> None:
        """Set solution objective function values."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def soln_ineqcv(self) -> numpy.ndarray:
        """Solution inequality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    @soln_ineqcv.setter
    @abstractmethod
    def soln_ineqcv(self, value: numpy.ndarray) -> None:
        """Set solution inequality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def soln_eqcv(self) -> numpy.ndarray:
        """Solution equality constraint violation function values."""
        raise NotImplementedError("property is abstract")
    @soln_eqcv.setter
    @abstractmethod
    def soln_eqcv(self, value: numpy.ndarray) -> None:
        """Set solution equality constraint violation function values."""
        raise NotImplementedError("property is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SolutionType(v: object, vname: str) -> None:
    """
    Check if object is of type ABCSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SolutionType):
        raise TypeError("'{0}' must be of type ABCSolution.".format(vname))
