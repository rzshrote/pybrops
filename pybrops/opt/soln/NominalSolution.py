"""
Module for defining optimization problem solutions with nominal decision variables.
"""

from pybrops.opt.soln.Solution import Solution

class NominalSolution(Solution):
    """
    Base class for all optimization problem solutions with nominal decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for NominalSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(NominalSolution, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_NominalSolution(v: object, vname: str) -> None:
    """
    Check if object is of type NominalSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, NominalSolution):
        raise TypeError("'{0}' must be of type NominalSolution.".format(vname))
