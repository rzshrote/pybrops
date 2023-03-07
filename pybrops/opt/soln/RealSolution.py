"""
Module for defining optimization problem solutions with real decision variables.
"""

# list of all public imports in the module
__all__ = [
    "RealSolution",
    "check_is_RealSolution"
]

# imports
import numpy
from pybrops.opt.soln.Solution import Solution

class RealSolution(Solution):
    """
    Base class for all optimization problem solutions with real decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(RealSolution, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def decn_lower(self) -> numpy.ndarray:
        """Lower bound of decision space."""
        raise NotImplementedError("property is abstract")
    @decn_lower.setter
    def decn_lower(self, value: numpy.ndarray) -> None:
        """Set lower bound of decision space."""
        raise NotImplementedError("property is abstract")
    @decn_lower.deleter
    def decn_lower(self) -> None:
        """Delete lower bound of decision space."""
        raise NotImplementedError("property is abstract")
    
    @property
    def decn_upper(self) -> numpy.ndarray:
        """Upper bound of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_upper.setter
    def decn_upper(self, value: numpy.ndarray) -> None:
        """Set upper bound of the decision space."""
        raise NotImplementedError("property is abstract")
    @decn_upper.deleter
    def decn_upper(self) -> None:
        """Delete upper bound of the decision space."""
        raise NotImplementedError("property is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_RealSolution(v: object, vname: str) -> None:
    """
    Check if object is of type RealSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealSolution):
        raise TypeError("'{0}' must be of type RealSolution.".format(vname))
