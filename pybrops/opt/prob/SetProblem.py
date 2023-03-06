"""
Module for defining optimization problems with nominal decision variables.
"""

from pybrops.opt.prob.Problem import Problem

class SetProblem(Problem):
    """
    Base class for all optimization problems with nominal decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SetProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SetProblem, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SetProblem(v: object, vname: str) -> None:
    """
    Check if object is of type SetProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SetProblem):
        raise TypeError("'{0}' must be of type SetProblem.".format(vname))
